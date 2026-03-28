#!/usr/bin/env python3
"""
03a PLIP Interaction Analysis - CLI (Hit Validation)
=======================================================
Analyze protein-ligand interactions for docked poses using PLIP.

Unlike reference_docking (crystal complex), hit_validation runs PLIP
on docked poses (receptor + best scored pose from 01c).

Usage:
    python 02_scripts/03a_plip_interaction_analysis.py \
        --campaigns 04_data/campaigns/NAMIKI_top20_pH63/campaign_config.yaml \
        --config 03_configs/03a_plip_interaction_analysis.yaml

Project: hit_validation
Module: 03a
Version: 2.0 (2026-03-27)
"""

import argparse
import logging
import sys
import yaml
from pathlib import Path

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s | %(levelname)-8s | %(message)s")

sys.path.insert(0, str(Path(__file__).parent.parent / "01_src"))

from hit_validation.m03_interaction_analysis.plip_interaction_analysis import (
    run_plip_batch_analysis,
)

logger = logging.getLogger(__name__)


def load_yaml(path):
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def setup_log_file(log_path: Path, log_level: str = "INFO"):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(str(log_path), mode="w", encoding="utf-8")
    fh.setLevel(getattr(logging, log_level.upper(), logging.INFO))
    fh.setFormatter(logging.Formatter("%(asctime)s | %(levelname)-8s | %(message)s"))
    logging.getLogger().addHandler(fh)


def main():
    parser = argparse.ArgumentParser(
        description="03a PLIP Interaction Analysis — docked pose interactions",
    )
    parser.add_argument("--campaigns", type=str, required=True)
    parser.add_argument("--config", "-c", type=str, default=None)
    parser.add_argument("--output", "-o", type=str, default=None)
    parser.add_argument("--log-level", type=str, default=None,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    args = parser.parse_args()

    # --- Config defaults ---
    pose_source = "01c_dock6_run"
    pose_selection = "best_score"
    output_name = "interactions"
    log_level = "INFO"

    if args.config:
        mc = load_yaml(args.config)
        params = mc.get("parameters", {})
        pose_source = params.get("pose_source", pose_source)
        pose_selection = params.get("pose_selection", pose_selection)
        output_name = params.get("output_name", output_name)
        log_level = params.get("log_level", log_level)

    if args.log_level:
        log_level = args.log_level
    logging.getLogger().setLevel(getattr(logging, log_level))

    # --- Resolve paths from campaign config ---
    cc = load_yaml(args.campaigns)
    campaign_dir = Path(args.campaigns).parent
    campaign_id = cc.get("campaign_id", campaign_dir.name)
    results_base = Path("05_results") / campaign_id

    # Docking output directory (from 01c)
    docking_dir = str(results_base / pose_source)

    # Receptor: protonated from 00b
    receptor_path = None
    candidates = [
        results_base / "00b_receptor_preparation" / "receptor_protonated.pdb",
        results_base / "00b_receptor_preparation" / "receptor_clean.pdb",
    ]
    for c in candidates:
        if c.exists():
            receptor_path = str(c)
            break

    if receptor_path is None:
        # Fallback to campaign receptor PDB
        rec_cfg = cc.get("receptor", {})
        rec_pdb = rec_cfg.get("pdb")
        if rec_pdb:
            candidate = campaign_dir / rec_pdb
            if candidate.exists():
                receptor_path = str(candidate)

    # Output directory
    output_subdir = "03a_plip_analysis"
    if args.config:
        mc = load_yaml(args.config)
        output_subdir = mc.get("outputs", {}).get("subdir", output_subdir)
    output_dir = args.output or str(results_base / output_subdir)

    # --- Validate ---
    if receptor_path is None:
        logger.error("Receptor PDB not found. Run 00b first.")
        sys.exit(1)
    if not Path(docking_dir).exists():
        logger.error(f"Docking dir not found: {docking_dir}")
        logger.error("Run 01c first.")
        sys.exit(1)

    # Log file
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    setup_log_file(out_path / "03a_plip_analysis.log", log_level)

    logger.info("=" * 60)
    logger.info("  HIT_VALIDATION - Module 03a: PLIP Batch Analysis")
    logger.info("=" * 60)
    logger.info(f"  Campaign:    {campaign_id}")
    logger.info(f"  Docking dir: {docking_dir}")
    logger.info(f"  Receptor:    {receptor_path}")
    logger.info(f"  Output:      {output_dir}")

    # --- Run ---
    result = run_plip_batch_analysis(
        docking_dir=docking_dir,
        output_dir=output_dir,
        receptor_pdb=receptor_path,
        pose_source=pose_source,
        pose_selection=pose_selection,
        output_name=output_name,
    )

    if result.get("success"):
        logger.info(f"\nPLIP analysis complete: {result['n_ok']}/{result['n_total']}")
        sys.exit(0)
    else:
        logger.error(f"Failed: {result.get('error')}")
        sys.exit(1)


if __name__ == "__main__":
    main()