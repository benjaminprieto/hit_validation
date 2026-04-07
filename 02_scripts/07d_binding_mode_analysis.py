#!/usr/bin/env python3
"""
07d Binding Mode Analysis - CLI (Hit Validation)
==================================================
Compares molecules by binding mode. Clusters by IEV similarity,
detects emergent sub-pockets from co-contact network.

Input:  04b footprint_per_molecule.csv, optionally 01i ProLIF
Output: 05_results/{campaign}/07d_binding_mode_analysis/

Usage:
    python 02_scripts/07d_binding_mode_analysis.py \
        --config 03_configs/07d_binding_mode_analysis.yaml \
        --campaigns 04_data/campaigns/druglikeness_1/campaign_config.yaml

    # Skip ProLIF layer:
    python 02_scripts/07d_binding_mode_analysis.py \
        --config 03_configs/07d_binding_mode_analysis.yaml \
        --campaigns 04_data/campaigns/druglikeness_1/campaign_config.yaml \
        --no-prolif

    # Skip community detection:
    python 02_scripts/07d_binding_mode_analysis.py \
        --config 03_configs/07d_binding_mode_analysis.yaml \
        --campaigns 04_data/campaigns/druglikeness_1/campaign_config.yaml \
        --no-communities

Project: hit_validation
Module: 07d
Version: 1.0
"""

import argparse
import logging
import sys
import yaml
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-8s | %(message)s")
sys.path.insert(0, str(Path(__file__).parent.parent / "01_src"))

from hit_validation.m07_decision_report.binding_mode_analysis import run_binding_mode_analysis

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
        description="07d Binding Mode Analysis -- IEV clustering + emergent sub-pockets",
    )
    parser.add_argument("--config", "-c", type=str, required=True, help="Module YAML")
    parser.add_argument("--campaigns", type=str, required=True, help="Campaign YAML")
    parser.add_argument("--output", "-o", type=str, default=None)
    parser.add_argument("--no-prolif", action="store_true",
                        help="Skip ProLIF layer (Layer 2)")
    parser.add_argument("--no-communities", action="store_true",
                        help="Skip community detection (Layer 3)")
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
    output_subdir = mc.get("outputs", {}).get("subdir", "07d_binding_mode_analysis")
    output_dir = Path(args.output) if args.output else results_base / output_subdir

    footprint_csv = results_base / "04_dock6_analysis" / "04b_footprint_analysis" / "footprint_per_molecule.csv"
    if not footprint_csv.exists():
        logger.error(f"Required: {footprint_csv}")
        return 1

    # Optional: ProLIF from 01i
    md_dir = results_base / "01i_trajectory_analysis"
    md_dir = str(md_dir) if md_dir.exists() else None

    # =========================================================================
    # MERGE PARAMETERS
    # =========================================================================
    energy_cutoff = params.get("energy_cutoff", 0.5)
    similarity_method = params.get("similarity_method", "pearson")
    cocontact_min_weight = params.get("cocontact_min_weight", 2)
    community_resolution = params.get("community_resolution", 1.0)

    do_prolif = not args.no_prolif
    do_communities = not args.no_communities

    log_level = args.log_level or params.get("log_level", "INFO")
    logging.getLogger().setLevel(getattr(logging, log_level.upper()))

    # =========================================================================
    # SETUP LOGGING
    # =========================================================================
    setup_log_file(output_dir / "07d_binding_mode_analysis.log", log_level)

    logger.info("=" * 60)
    logger.info("  HIT_VALIDATION - Module 07d: Binding Mode Analysis")
    logger.info("=" * 60)
    logger.info(f"Campaign:       {campaign_id}")
    logger.info(f"Footprint:      {footprint_csv}")
    logger.info(f"ProLIF (01i):   {md_dir or 'not found'}")
    logger.info(f"Layer 2 (ProLIF): {'yes' if do_prolif else 'skip'}")
    logger.info(f"Layer 3 (communities): {'yes' if do_communities else 'skip'}")
    logger.info(f"Output:         {output_dir}")

    # =========================================================================
    # EXECUTE
    # =========================================================================
    result = run_binding_mode_analysis(
        campaign_id=campaign_id,
        footprint_csv=str(footprint_csv),
        output_dir=str(output_dir),
        campaign_config=cc,
        md_dir=md_dir,
        energy_cutoff=energy_cutoff,
        similarity_method=similarity_method,
        cocontact_min_weight=cocontact_min_weight,
        community_resolution=community_resolution,
        run_communities=do_communities,
        run_prolif=do_prolif,
    )

    if not result.get("success"):
        logger.error(f"Error: {result.get('error')}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())