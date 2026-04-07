#!/usr/bin/env python3
"""
01i Trajectory Analysis - CLI (Hit Validation)
================================================
Extracts per-frame data from MD trajectories using cpptraj and ProLIF.

Input:  05_results/{campaign}/01g_mmpbsa_decomp/{name}/
Output: 05_results/{campaign}/01i_trajectory_analysis/{name}/

Usage:
    # All molecules:
    python 02_scripts/01i_trajectory_analysis.py \
        --config 03_configs/01i_trajectory_analysis.yaml \
        --campaigns 04_data/campaigns/druglikeness_1/campaign_config.yaml

    # Specific molecules:
    python 02_scripts/01i_trajectory_analysis.py \
        --config 03_configs/01i_trajectory_analysis.yaml \
        --campaigns 04_data/campaigns/druglikeness_1/campaign_config.yaml \
        --molecules PubChem-146072367 PubChem-59068465

    # Skip ProLIF (cpptraj only):
    python 02_scripts/01i_trajectory_analysis.py \
        --config 03_configs/01i_trajectory_analysis.yaml \
        --campaigns 04_data/campaigns/druglikeness_1/campaign_config.yaml \
        --no-prolif

Project: hit_validation
Module: 01i
Version: 1.0
"""

import argparse
import logging
import sys
import yaml
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-8s | %(message)s")
sys.path.insert(0, str(Path(__file__).parent.parent / "01_src"))

from hit_validation.m01_docking.trajectory_analysis import (
    run_trajectory_analysis,
    run_trajectory_analysis_batch,
)

logger = logging.getLogger(__name__)


def load_yaml(path):
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def setup_log_file(log_path: Path, log_level: str = "INFO"):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(str(log_path), encoding="utf-8")
    fh.setLevel(getattr(logging, log_level.upper(), logging.INFO))
    fh.setFormatter(logging.Formatter("%(asctime)s | %(levelname)-8s | %(message)s"))
    logging.getLogger().addHandler(fh)


def main():
    parser = argparse.ArgumentParser(
        description="01i Trajectory Analysis -- cpptraj + ProLIF",
    )
    parser.add_argument("--config", "-c", type=str, required=True, help="Module YAML")
    parser.add_argument("--campaigns", type=str, required=True, help="Campaign YAML")
    parser.add_argument("--output", "-o", type=str, default=None)
    parser.add_argument("--molecules", nargs="*", default=None,
                        help="Specific molecules (default: all)")
    parser.add_argument("--no-prolif", action="store_true",
                        help="Skip ProLIF analysis (cpptraj only)")
    parser.add_argument("--reparse-only", action="store_true",
                        help="Skip cpptraj, re-parse existing .dat files only")
    parser.add_argument("--residue-mapping", type=str, default=None,
                        help="Path to residue_mapping.csv (default: auto-detect from 04b)")
    parser.add_argument("--log-level", type=str, default=None,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    args = parser.parse_args()

    # =========================================================================
    # LOAD CONFIGS
    # =========================================================================
    cc = load_yaml(args.campaigns)
    campaign_dir = Path(args.campaigns).parent
    campaign_id = cc.get("campaign_id", campaign_dir.name)
    mc = load_yaml(args.config)
    params = mc.get("parameters", {})

    # =========================================================================
    # RESOLVE PATHS
    # =========================================================================
    results_base = Path("05_results") / campaign_id
    output_subdir = mc.get("outputs", {}).get("subdir", "01i_trajectory_analysis")
    output_dir = Path(args.output) if args.output else results_base / output_subdir

    mmpbsa_dir = results_base / "01g_mmpbsa_decomp"

    if args.residue_mapping:
        residue_mapping_csv = args.residue_mapping
    else:
        footprint_dir = results_base / "04_dock6_analysis" / "04b_footprint_analysis"
        residue_mapping_csv = str(footprint_dir / "residue_mapping.csv")
    if not Path(residue_mapping_csv).exists():
        logger.warning(f"Residue mapping not found: {residue_mapping_csv}")
        residue_mapping_csv = None

    # =========================================================================
    # MERGE PARAMETERS
    # =========================================================================
    lig_resname = params.get("lig_resname", "UNL")
    n_protein_res = params.get("n_protein_res", 706)
    cpptraj_timeout = params.get("cpptraj_timeout", 600)
    cpptraj_solvated_timeout = params.get("cpptraj_solvated_timeout", 3600)
    do_prolif = not args.no_prolif and params.get("run_prolif", True)

    log_level = args.log_level or params.get("log_level", "INFO")
    logging.getLogger().setLevel(getattr(logging, log_level.upper()))

    # =========================================================================
    # SETUP LOGGING
    # =========================================================================
    setup_log_file(output_dir / "01i_trajectory_analysis.log", log_level)

    logger.info("=" * 60)
    logger.info("  HIT_VALIDATION - Module 01i: Trajectory Analysis")
    logger.info("=" * 60)
    logger.info(f"Campaign:   {campaign_id}")
    logger.info(f"Molecules:  {args.molecules or 'ALL'}")
    logger.info(f"ProLIF:     {'yes' if do_prolif else 'no'}")
    logger.info(f"Output:     {output_dir}")

    # =========================================================================
    # EXECUTE
    # =========================================================================
    if args.molecules and len(args.molecules) == 1:
        mol_name = args.molecules[0]
        mol_mmpbsa = mmpbsa_dir / mol_name
        mol_output = output_dir / mol_name

        result = run_trajectory_analysis(
            mol_name=mol_name,
            mmpbsa_dir=str(mol_mmpbsa),
            output_dir=str(mol_output),
            campaign_config=cc,
            residue_mapping_csv=residue_mapping_csv,
            lig_resname=lig_resname,
            n_protein_res=n_protein_res,
            cpptraj_timeout=cpptraj_timeout,
            cpptraj_solvated_timeout=cpptraj_solvated_timeout,
            run_prolif=do_prolif,
            reparse_only=args.reparse_only,
        )
    else:
        result = run_trajectory_analysis_batch(
            campaign_id=campaign_id,
            mmpbsa_base_dir=str(mmpbsa_dir),
            output_base_dir=str(output_dir),
            campaign_config=cc,
            residue_mapping_csv=residue_mapping_csv,
            molecules=args.molecules,
            lig_resname=lig_resname,
            n_protein_res=n_protein_res,
            run_prolif=do_prolif,
            cpptraj_timeout=cpptraj_timeout,
            cpptraj_solvated_timeout=cpptraj_solvated_timeout,
            reparse_only=args.reparse_only,
        )

    if not result.get("success"):
        logger.error(f"Error: {result.get('error')}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())