#!/usr/bin/env python3
"""
07b Residue Comparison Report - CLI (Hit Validation)
======================================================
Residue-by-residue comparison of each hit vs crystallographic reference.
07a decides purchase/template. 07b explains WHY.

Usage:
    python 02_scripts/07b_residue_comparison.py \
        --config 03_configs/07b_residue_comparison.yaml \
        --campaigns 04_data/campaigns/druglikeness_1/campaign_config.yaml

Project: hit_validation
Module: 07b
Version: 2.0 (2026-04-02)
"""

import argparse
import logging
import sys
import yaml
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-8s | %(message)s")
sys.path.insert(0, str(Path(__file__).parent.parent / "01_src"))

from hit_validation.m07_decision_report.residue_comparison import run_residue_comparison

logger = logging.getLogger(__name__)


def load_yaml(path):
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def _parse_zones(cc: dict) -> dict:
    """Parse zone definitions from campaign config into internal format."""
    raw_zones = cc.get("zones", {})
    if not raw_zones:
        return {}
    zones = {}
    for zone_id, zdef in raw_zones.items():
        zones[zone_id] = {
            "residues": set(zdef.get("residues", [])),
            "label": zdef.get("label", zone_id),
            "druglike": zdef.get("druglike", True),
            "description": zdef.get("description", ""),
            "color": zdef.get("color"),
        }
    return zones


def main():
    parser = argparse.ArgumentParser(
        description="07b Residue Comparison Report — hit vs reference per residue",
    )
    parser.add_argument("--config", "-c", type=str, required=True)
    parser.add_argument("--campaigns", type=str, required=True)
    parser.add_argument("--output", "-o", type=str, default=None)
    parser.add_argument("--log-level", type=str, default=None,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    args = parser.parse_args()

    cc = load_yaml(args.campaigns)
    campaign_dir = Path(args.campaigns).parent
    campaign_id = cc.get("campaign_id", campaign_dir.name)
    mc = load_yaml(args.config)
    params = mc.get("parameters", {})

    # --- Zones ---
    zones = _parse_zones(cc)

    # --- Output ---
    output_subdir = mc.get("outputs", {}).get("subdir", "07b_residue_comparison")
    output_dir = args.output or str(Path("05_results") / campaign_id / output_subdir)

    # --- Results base ---
    results_base = Path("05_results") / campaign_id

    # --- Footprint CSV (obligatory) ---
    footprint_csv = str(
        results_base / "04_dock6_analysis" / "04b_footprint_analysis"
        / "footprint_per_molecule.csv"
    )

    # --- Optional inputs ---
    # PLIP analysis dir
    plip_analysis_dir = str(results_base / "03a_plip_analysis")
    if not Path(plip_analysis_dir).exists():
        plip_analysis_dir = None

    # 07a decision summary
    decision_summary_csv = str(
        results_base / "07a_decision_report" / "decision_summary.csv"
    )
    if not Path(decision_summary_csv).exists():
        decision_summary_csv = None

    # --- Reference label ---
    ref_ctx = cc.get("reference_context", {})
    reference_label = ref_ctx.get("label", "Reference")

    # --- Module parameters ---
    energy_cutoff = params.get("energy_cutoff", 0.5)
    delta_better = params.get("delta_better_threshold", -1.0)
    delta_worse = params.get("delta_worse_threshold", 1.0)
    include_plip = params.get("include_plip", True)

    log_level = args.log_level or params.get("log_level", "INFO")
    logging.getLogger().setLevel(getattr(logging, log_level.upper()))

    logger.info("=" * 60)
    logger.info("  HIT_VALIDATION - Module 07b: Residue Comparison")
    logger.info("=" * 60)
    logger.info(f"Campaign:      {campaign_id}")
    logger.info(f"Footprint:     {footprint_csv}")
    logger.info(f"Reference:     {reference_label}")
    logger.info(f"PLIP dir:      {plip_analysis_dir or 'not found'}")
    logger.info(f"07a summary:   {decision_summary_csv or 'not found'}")
    logger.info(f"Output:        {output_dir}")

    result = run_residue_comparison(
        footprint_csv=footprint_csv,
        output_dir=output_dir,
        zones=zones,
        campaign_id=campaign_id,
        energy_cutoff=energy_cutoff,
        delta_better_threshold=delta_better,
        delta_worse_threshold=delta_worse,
        include_plip=include_plip,
        plip_analysis_dir=plip_analysis_dir,
        decision_summary_csv=decision_summary_csv,
        reference_label=reference_label,
    )

    if not result.get("success"):
        logger.error(f"Failed: {result.get('error')}")
        return 1

    logger.info("07b complete.")
    logger.info(f"Next: open {output_dir}/residue_comparison.html")
    return 0


if __name__ == "__main__":
    sys.exit(main())