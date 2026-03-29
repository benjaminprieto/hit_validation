#!/usr/bin/env python3
"""
07a Decision Report - CLI (Hit Validation)
=============================================
Generate consolidated validation evidence report with multi-criteria
pose selection and mol2 export for screening hits.

Usage:
    python 02_scripts/07a_decision_report.py \
        --config 03_configs/07a_decision_report.yaml \
        --campaigns 04_data/campaigns/NAMIKI_top20_pH63/campaign_config.yaml

Project: hit_validation
Module: 07a
Version: 2.0 (2026-03-29)
"""

import argparse
import logging
import sys
import yaml
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-8s | %(message)s")
sys.path.insert(0, str(Path(__file__).parent.parent / "01_src"))

from hit_validation.m07_decision_report.decision_report import run_decision_report

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
    parser = argparse.ArgumentParser(description="07a Decision Report — pose selection + validation evidence")
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

    results_base = Path("05_results") / campaign_id

    # Resolve paths
    output_subdir = mc.get("outputs", {}).get("subdir", "07a_decision_report")
    output_dir = args.output or str(results_base / output_subdir)

    # Input paths
    scores_csv = str(results_base / "01e_score_collection" / "dock6_scores.csv")
    all_poses_csv = str(results_base / "01e_score_collection" / "dock6_all_poses.csv")
    plip_dir = str(results_base / "03a_plip_analysis")
    footprint_dir = str(results_base / "04_dock6_analysis" / "04b_footprint_analysis")
    docking_dir = str(results_base / "01c_dock6_run")

    # Reference context — configurable reference data (informational only)
    ref_context = cc.get("reference_context", {})
    ref_scores_path = ref_context.get("scores_csv", params.get("reference_scores_path"))
    ref_plip_path = ref_context.get("plip_json", params.get("reference_plip_path"))
    ref_footprint_path = ref_context.get("footprint_csv", params.get("reference_footprint_path"))
    ref_label = ref_context.get("label", "Reference")

    # MMPBSA analysis (from 01h)
    mmpbsa_analysis_dir = str(results_base / "01h_mmpbsa_analysis")
    if not Path(mmpbsa_analysis_dir).exists():
        mmpbsa_analysis_dir = None

    # Zones from campaign config
    zones = _parse_zones(cc)

    # Pose selection weights
    pose_weights_cfg = params.get("pose_selection", {})
    pose_weights = pose_weights_cfg.get("weights")

    # Decision thresholds
    zone_cutoff = params.get("zone_energy_cutoff", -0.5)
    min_zones = params.get("min_zones_for_template", 2)
    grid_threshold = params.get("grid_score_purchase_threshold", -30.0)
    mmpbsa_dg_strong = params.get("mmpbsa_dg_strong_threshold", -20.0)
    mmpbsa_dg_moderate = params.get("mmpbsa_dg_moderate_threshold", -10.0)

    log_level = args.log_level or params.get("log_level", "INFO")
    logging.getLogger().setLevel(getattr(logging, log_level.upper()))

    logger.info("=" * 60)
    logger.info("  HIT_VALIDATION - Module 07a: Decision Report")
    logger.info("=" * 60)
    logger.info(f"Campaign:  {campaign_id}")
    logger.info(f"Scores:    {scores_csv}")
    logger.info(f"Zones:     {len(zones)} defined" if zones else "Zones:     none (zone analysis skipped)")
    logger.info(f"Output:    {output_dir}")

    result = run_decision_report(
        scores_csv=scores_csv,
        output_dir=output_dir,
        plip_dir=plip_dir,
        footprint_dir=footprint_dir,
        reference_scores_path=ref_scores_path,
        reference_plip_path=ref_plip_path,
        reference_footprint_path=ref_footprint_path,
        zone_energy_cutoff=zone_cutoff,
        min_zones_for_template=min_zones,
        grid_score_purchase_threshold=grid_threshold,
        campaign_id=campaign_id,
        reference_label=ref_label,
        mmpbsa_analysis_dir=mmpbsa_analysis_dir,
        mmpbsa_dg_strong_threshold=mmpbsa_dg_strong,
        mmpbsa_dg_moderate_threshold=mmpbsa_dg_moderate,
        zones=zones,
        docking_dir=docking_dir,
        all_poses_csv=all_poses_csv,
        pose_selection_weights=pose_weights,
    )

    if not result.get("success"):
        logger.error(f"Failed: {result.get('error')}")
        return 1

    logger.info("07a complete.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
