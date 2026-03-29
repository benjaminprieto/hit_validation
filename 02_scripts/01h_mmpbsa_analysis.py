#!/usr/bin/env python3
"""
01h MMPBSA Analysis - CLI (Hit Validation)
=============================================
Batch analysis of MMPBSA per-residue decomposition output from 01g.

Parses FINAL_DECOMP_MMPBSA.dat for each molecule, maps to PDB numbering,
compares with footprint, and produces consolidated comparison tables.

Input:  05_results/{campaign}/01g_mmpbsa_decomp/{name}/mmpbsa/FINAL_DECOMP_MMPBSA.dat
Output: 05_results/{campaign}/01h_mmpbsa_analysis/

Usage:
    python 02_scripts/01h_mmpbsa_analysis.py \\
        --config 03_configs/01h_mmpbsa_analysis.yaml \\
        --campaigns 04_data/campaigns/NAMIKI_top20_pH63/campaign_config.yaml

    # With explicit footprint directory:
    python 02_scripts/01h_mmpbsa_analysis.py \\
        --config 03_configs/01h_mmpbsa_analysis.yaml \\
        --campaigns 04_data/campaigns/NAMIKI_top20_pH63/campaign_config.yaml \\
        --footprint-dir 05_results/.../04b_footprint_analysis/

Project: hit_validation
Module: 01h
Version: 1.0 (2026-03-28)
"""

import argparse
import logging
import sys
import yaml
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-8s | %(message)s")
sys.path.insert(0, str(Path(__file__).parent.parent / "01_src"))

from hit_validation.m01_docking.mmpbsa_analysis import run_mmpbsa_batch_analysis

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
        description="01h MMPBSA Analysis — batch parse + footprint comparison",
    )
    parser.add_argument("--config", "-c", type=str, required=True, help="Module YAML")
    parser.add_argument("--campaigns", type=str, required=True, help="Campaign YAML")
    parser.add_argument("--output", "-o", type=str, default=None)
    parser.add_argument("--decomp-dir", type=str, default=None,
                        help="Override 01g output directory")
    parser.add_argument("--footprint-dir", type=str, default=None,
                        help="Path to 04b footprint analysis output")
    parser.add_argument("--no-compare-footprint", action="store_true",
                        help="Skip footprint comparison")
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
    output_subdir = mc.get("outputs", {}).get("subdir", "01h_mmpbsa_analysis")
    output_dir = Path(args.output) if args.output else results_base / output_subdir

    # 01g output directory
    decomp_base = Path(args.decomp_dir) if args.decomp_dir else results_base / "01g_mmpbsa_decomp"

    if not decomp_base.exists():
        logger.error(f"MMPBSA results directory not found: {decomp_base}")
        logger.error("Run module 01g first.")
        return 1

    # Receptor PDB (from 00b)
    rec_pdb = results_base / "00b_receptor_preparation" / "receptor_protonated.pdb"
    if not rec_pdb.exists():
        rec_pdb = results_base / "00b_receptor_preparation" / "rec_noH.pdb"
    if not rec_pdb.exists():
        logger.error("Receptor PDB not found in 00b_receptor_preparation/")
        return 1

    # Footprint directory (auto-resolve from 04b)
    compare_fp = params.get("compare_footprint", True) and not args.no_compare_footprint
    footprint_dir = args.footprint_dir

    if not footprint_dir and compare_fp:
        fp_candidates = [
            results_base / "04b_footprint_analysis",
            results_base / "04_dock6_analysis" / "04b_footprint_analysis",
        ]
        for fp in fp_candidates:
            if fp.exists():
                footprint_dir = str(fp)
                break

    # Reference context from campaign config
    reference_context = cc.get("reference_context")

    # Zones from campaign config
    zones = {}
    raw_zones = cc.get("zones", {})
    for zone_id, zdef in raw_zones.items():
        zones[zone_id] = {
            "residues": set(zdef.get("residues", [])),
            "label": zdef.get("label", zone_id),
            "description": zdef.get("description", ""),
        }

    # =========================================================================
    # SETUP
    # =========================================================================
    log_level = args.log_level or params.get("log_level", "INFO")
    logging.getLogger().setLevel(getattr(logging, log_level.upper()))
    setup_log_file(output_dir / "01h_mmpbsa_analysis.log", log_level)

    logger.info("=" * 60)
    logger.info("  HIT_VALIDATION - Module 01h: MMPBSA Batch Analysis")
    logger.info("=" * 60)
    logger.info(f"Campaign:   {campaign_id}")
    logger.info(f"Decomp dir: {decomp_base}")
    logger.info(f"Receptor:   {rec_pdb.name}")
    logger.info(f"Footprint:  {footprint_dir or 'not found'}")
    logger.info(f"Output:     {output_dir}")

    # =========================================================================
    # EXECUTE
    # =========================================================================
    result = run_mmpbsa_batch_analysis(
        mmpbsa_results_dir=str(decomp_base),
        receptor_pdb=str(rec_pdb),
        output_dir=str(output_dir),
        reference_context=reference_context,
        compare_footprint=compare_fp,
        footprint_dir=footprint_dir,
        campaign_id=campaign_id,
        zones=zones if zones else None,
    )

    if not result.get("success"):
        logger.error(f"Error: {result.get('error')}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
