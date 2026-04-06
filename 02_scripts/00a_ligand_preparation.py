#!/usr/bin/env python3
"""
00a Ligand Preparation - CLI (Hit Validation)
===============================================
Prepares screening hit ligands for DOCK6 validation docking.
Uses antechamber with AM1-BCC charges + Sybyl atom types.

Usage:
    python 02_scripts/00a_ligand_preparation.py \
        --config 03_configs/00a_ligand_preparation.yaml \
        --campaigns 04_data/campaigns/NAMIKI_top20_pH63/campaign_config.yaml

Project: hit_validation
Module: 00a
Version: 1.0 (2026-03-27)
"""

import argparse
import logging
import sys
import yaml
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-8s | %(message)s")
sys.path.insert(0, str(Path(__file__).parent.parent / "01_src"))

from hit_validation.m00_preparation.ligand_preparation import run_ligand_preparation

logger = logging.getLogger(__name__)


def load_yaml(path):
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def main():
    parser = argparse.ArgumentParser(description="00a Ligand Preparation — screening hits for validation docking")
    parser.add_argument("--config", "-c", type=str, required=True)
    parser.add_argument("--campaigns", type=str, required=True)
    parser.add_argument("--output", "-o", type=str, default=None)
    parser.add_argument("--log-level", type=str, default=None, choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    args = parser.parse_args()

    cc = load_yaml(args.campaigns)
    campaign_dir = Path(args.campaigns).parent
    campaign_id = cc.get("campaign_id", campaign_dir.name)
    mc = load_yaml(args.config)
    params = mc.get("parameters", {})

    # Resolve paths
    output_subdir = mc.get("outputs", {}).get("subdir", "00a_ligand_preparation")
    output_dir = args.output or str(Path("05_results") / campaign_id / output_subdir)

    # Ligand directory
    ligands_dir_key = cc.get("ligands_dir", "ligands/")
    ligand_dir = campaign_dir / ligands_dir_key
    if not ligand_dir.exists():
        logger.error(f"Ligands dir not found: {ligand_dir}")
        logger.error("Set ligands_dir in campaign_config.yaml")
        return 1

    docking_ph = cc.get("docking_ph", 6.3)
    atom_type = params.get("atom_type", "sybyl")
    charge_method = params.get("charge_method", "bcc")
    timeout = params.get("timeout_per_molecule", 300)
    protonation_tool = params.get("protonation_tool", "obabel")

    log_level = args.log_level or params.get("log_level", "INFO")
    logging.getLogger().setLevel(getattr(logging, log_level.upper()))

    logger.info("=" * 60)
    logger.info("  HIT_VALIDATION - Module 00a: Ligand Preparation")
    logger.info("=" * 60)
    logger.info(f"Campaign:    {campaign_id}")
    logger.info(f"Ligand dir:  {ligand_dir}")
    logger.info(f"Atom type:   {atom_type}")
    logger.info(f"Charges:     {charge_method}")
    logger.info(f"pH:          {docking_ph}")
    logger.info(f"Output:      {output_dir}")

    result = run_ligand_preparation(
        ligand_dir=str(ligand_dir),
        output_dir=output_dir,
        atom_type=atom_type,
        charge_method=charge_method,
        docking_ph=docking_ph,
        timeout_per_molecule=timeout,
        protonation_tool=protonation_tool,
    )

    if not result.get("success"):
        logger.error(f"Failed: {result.get('error')}")
        return 1

    logger.info("00a complete.")
    return 0


if __name__ == "__main__":
    sys.exit(main())