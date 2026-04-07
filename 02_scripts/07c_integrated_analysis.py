#!/usr/bin/env python3
"""
07c Integrated Residue Analysis - CLI (Hit Validation)
=======================================================
Per-residue, per-molecule report integrating all pipeline evidence.

Input:  04b footprint, 01h MMPBSA, 03a PLIP (optional), 01i MD (optional)
Output: 05_results/{campaign}/07c_integrated_analysis/

Usage:
    python 02_scripts/07c_integrated_analysis.py \
        --config 03_configs/07c_integrated_analysis.yaml \
        --campaigns 04_data/campaigns/druglikeness_1/campaign_config.yaml

Project: hit_validation
Module: 07c
Version: 1.0
"""

import argparse
import logging
import sys
import yaml
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-8s | %(message)s")
sys.path.insert(0, str(Path(__file__).parent.parent / "01_src"))

from hit_validation.m07_decision_report.integrated_analysis import run_integrated_analysis

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
        description="07c Integrated Residue Analysis",
    )
    parser.add_argument("--config", "-c", type=str, required=True, help="Module YAML")
    parser.add_argument("--campaigns", type=str, required=True, help="Campaign YAML")
    parser.add_argument("--output", "-o", type=str, default=None)
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
    output_subdir = mc.get("outputs", {}).get("subdir", "07c_integrated_analysis")
    output_dir = Path(args.output) if args.output else results_base / output_subdir

    # Required inputs
    footprint_csv = results_base / "04_dock6_analysis" / "04b_footprint_analysis" / "footprint_per_molecule.csv"
    mmpbsa_global_csv = results_base / "01h_mmpbsa_analysis" / "consolidated_mmpbsa.csv"

    if not footprint_csv.exists():
        logger.error(f"Required: {footprint_csv}")
        return 1
    if not mmpbsa_global_csv.exists():
        logger.error(f"Required: {mmpbsa_global_csv}")
        return 1

    # Optional inputs (directories)
    mmpbsa_decomp_dir = results_base / "01h_mmpbsa_analysis"
    plip_dir = results_base / "03a_plip_analysis"
    md_dir = results_base / "01i_trajectory_analysis"

    mmpbsa_decomp_dir = str(mmpbsa_decomp_dir) if mmpbsa_decomp_dir.exists() else None
    plip_dir = str(plip_dir) if plip_dir.exists() else None
    md_dir = str(md_dir) if md_dir.exists() else None

    # Residue mapping for H-bond/water bridge matching (PDB <-> sequential)
    residue_mapping_csv = results_base / "04_dock6_analysis" / "04b_footprint_analysis" / "residue_mapping.csv"
    residue_mapping_csv = str(residue_mapping_csv) if residue_mapping_csv.exists() else None

    # =========================================================================
    # MERGE PARAMETERS
    # =========================================================================
    log_level = args.log_level or params.get("log_level", "INFO")
    logging.getLogger().setLevel(getattr(logging, log_level.upper()))

    # =========================================================================
    # SETUP LOGGING
    # =========================================================================
    setup_log_file(output_dir / "07c_integrated_analysis.log", log_level)

    logger.info("=" * 60)
    logger.info("  HIT_VALIDATION - Module 07c: Integrated Analysis")
    logger.info("=" * 60)
    logger.info(f"Campaign:       {campaign_id}")
    logger.info(f"Footprint:      {footprint_csv}")
    logger.info(f"MMPBSA global:  {mmpbsa_global_csv}")
    logger.info(f"MMPBSA decomp:  {mmpbsa_decomp_dir or 'not found'}")
    logger.info(f"PLIP:           {plip_dir or 'not found'}")
    logger.info(f"MD (01i):       {md_dir or 'not found'}")
    logger.info(f"Output:         {output_dir}")

    # =========================================================================
    # EXECUTE
    # =========================================================================
    result = run_integrated_analysis(
        campaign_id=campaign_id,
        footprint_csv=str(footprint_csv),
        mmpbsa_global_csv=str(mmpbsa_global_csv),
        output_dir=str(output_dir),
        campaign_config=cc,
        mmpbsa_decomp_dir=mmpbsa_decomp_dir,
        plip_dir=plip_dir,
        md_dir=md_dir,
        residue_mapping_csv=residue_mapping_csv,
        energy_cutoff=params.get("energy_cutoff", 0.5),
        reject_threshold=params.get("mmpbsa_reject_threshold", 0.0),
        weak_threshold=params.get("mmpbsa_weak_threshold", -15.0),
        unstable_threshold=params.get("rmsd_unstable_threshold", 4.0),
        mobile_threshold=params.get("rmsd_mobile_threshold", 2.0),
    )

    if not result.get("success"):
        logger.error(f"Error: {result.get('error')}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())