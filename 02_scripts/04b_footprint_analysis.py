#!/usr/bin/env python3
"""
04b DOCK6 Footprint Analysis - CLI
======================================
Per-residue vdW + ES energy decomposition via DOCK6 footprint scoring.

Phase 2 only: Parse per-residue footprint data from 01d TXT files.
Builds mol2→PDB residue mapping so footprint and contacts share numbering.

Input:
    01d_footprint_rescore/{name}/{name}_fps_footprint_scored.txt
    00b_receptor_preparation/rec_charged.mol2  (for residue mapping)
    00b_receptor_preparation/rec_noH.pdb       (for residue mapping)

Output:
    05_results/{campaigns}/04b_footprint_analysis/

Project: hit_validation
Module: 04b (DOCK6 analysis)
Version: 1.2 (2026-03-23) — adds residue remapping mol2→PDB
"""

import argparse
import logging
import sys
import yaml
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-8s | %(message)s")
sys.path.insert(0, str(Path(__file__).parent.parent / "01_src"))

from hit_validation.m04_dock6_analysis.footprint_analysis import run_footprint_analysis

logger = logging.getLogger(__name__)


def load_yaml(path):
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def setup_log_file(log_path, log_level="INFO"):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(str(log_path), mode="w", encoding="utf-8")
    fh.setLevel(getattr(logging, log_level.upper(), logging.INFO))
    fh.setFormatter(logging.Formatter("%(asctime)s | %(levelname)-8s | %(message)s"))
    logging.getLogger().addHandler(fh)


def run_consolidation_mode(args) -> int:
    """Consolidate 04b footprint outputs across N replicas."""
    from hit_validation.utils.replica_consolidation_helpers import discover_replicas
    from hit_validation.m04_dock6_analysis.footprint_analysis import (
        consolidate_footprint_analysis_replicas,
    )

    cc = load_yaml(args.campaigns)
    campaign_id = cc.get("campaign_id", Path(args.campaigns).parent.name)
    shared_base = Path("05_results") / campaign_id

    replicas = discover_replicas(shared_base, args.n_replicas)
    if len(replicas) < 2:
        logger.error(f"Need at least 2 replicas, found {len(replicas)}")
        return 1

    output_dir = shared_base / "04b_footprint_analysis" / "consolidated"
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Consolidating 04b from {len(replicas)} replicas -> {output_dir}")

    result = consolidate_footprint_analysis_replicas(replicas, output_dir)
    if not result.get("success"):
        logger.error(f"Consolidation failed: {result.get('error')}")
        return 1
    logger.info(
        f"04b consolidation done: {result['n_molecules']} molecules x "
        f"{result['n_residues']} residues, n_replicas={result['n_replicas']}"
    )
    return 0


def main():
    parser = argparse.ArgumentParser(
        description="04b DOCK6 Footprint Analysis — per-residue energy with PDB numbering",
    )
    parser.add_argument("--config", "-c", type=str, help="Module config YAML")
    parser.add_argument("--campaigns", type=str, help="Campaign config YAML")
    parser.add_argument("--output", "-o", type=str, default=None)
    parser.add_argument("--log-level", type=str, default=None)
    parser.add_argument("--replica-id", type=int, default=None,
                        help="Replica ID (1-indexed). None=legacy non-replicated mode.")
    parser.add_argument("--consolidate-replicas", action="store_true",
                        help="Consolidate footprint outputs across N replicas into "
                             "04b_footprint_analysis/consolidated/ with mean +/- std.")
    parser.add_argument("--n-replicas", type=int, default=1,
                        help="Number of replicas to consolidate. Used only with --consolidate-replicas.")

    args = parser.parse_args()

    if args.consolidate_replicas:
        if args.n_replicas < 2:
            parser.error("--consolidate-replicas requires --n-replicas >= 2")
        if args.replica_id is not None:
            parser.error("--consolidate-replicas cannot be combined with --replica-id")
        return run_consolidation_mode(args)

    # Defaults
    rescore_dir = None
    analysis_dir = None
    receptor_mol2 = None
    receptor_pdb = None
    campaign_id = "direct"
    pharmacophore_threshold = 0.8
    energy_cutoff = -0.5
    log_level = "INFO"
    reference_context_csv = None
    reference_context_label = "Reference"
    zones = {}

    # --- Campaign config ---
    if args.campaigns:
        cc = load_yaml(args.campaigns)
        campaign_dir = Path(args.campaigns).parent
        campaign_id = cc.get("campaign_id", campaign_dir.name)

        # --- Path resolution ---
        shared_base = Path("05_results") / campaign_id
        results_base = shared_base
        if args.replica_id is not None:
            results_base = shared_base / f"replica_{args.replica_id}"

        # Input: 01d footprint rescore output (REPLICATED)
        rescore_dir = str(results_base / "01d_footprint_rescore")

        # Output: 04b_footprint_analysis directly under results_base (REPLICATED).
        # Legacy 04_dock6_analysis/04b_footprint_analysis/ subdir removed 2026-04-26.
        analysis_dir = str(results_base / "04b_footprint_analysis")

        # Receptor files for residue mapping (SHARED)
        rec_mol2 = shared_base / "00b_receptor_preparation" / "rec_charged.mol2"
        rec_pdb = shared_base / "00b_receptor_preparation" / "rec_noH.pdb"
        if rec_mol2.exists():
            receptor_mol2 = str(rec_mol2)
        if rec_pdb.exists():
            receptor_pdb = str(rec_pdb)

        # Reference context — imported pre-computed data (configurable)
        ref_context = cc.get("reference_context", {})
        reference_context_csv = ref_context.get("footprint_csv")
        reference_context_label = ref_context.get("label", "Reference")

        # Zones from campaign config
        raw_zones = cc.get("zones", {})
        for zone_id, zdef in raw_zones.items():
            zones[zone_id] = {
                "residues": set(zdef.get("residues", [])),
                "label": zdef.get("label", zone_id),
                "druglike": zdef.get("druglike", True),
                "description": zdef.get("description", ""),
                "color": zdef.get("color"),
            }

    # --- Module config ---
    if args.config:
        mc = load_yaml(args.config)
        params = mc.get("parameters", {})
        pharmacophore_threshold = params.get("pharmacophore_threshold", pharmacophore_threshold)
        energy_cutoff = params.get("energy_cutoff", energy_cutoff)
        log_level = params.get("log_level", log_level)

    # --- CLI overrides ---
    if args.output:
        analysis_dir = args.output
    if args.log_level:
        log_level = args.log_level

    if not rescore_dir or not Path(rescore_dir).exists():
        logger.error(f"Footprint rescore dir not found: {rescore_dir}")
        logger.error("Run module 01d (footprint re-scoring) first.")
        return 1

    if not analysis_dir:
        analysis_dir = "05_results/04b_footprint_analysis"

    logging.getLogger().setLevel(getattr(logging, log_level.upper(), logging.INFO))
    setup_log_file(Path(analysis_dir) / "04b_footprint_analysis.log", log_level)

    logger.info("=" * 60)
    logger.info("  HIT_VALIDATION - Module 04b: DOCK6 Footprint Analysis")
    logger.info("=" * 60)
    logger.info(f"  Campaign:      {campaign_id}")
    logger.info(f"  Footprint dir: {rescore_dir} (from 01d)")
    logger.info(f"  Receptor mol2: {receptor_mol2 or 'not found'}")
    logger.info(f"  Receptor PDB:  {receptor_pdb or 'not found'}")
    logger.info(f"  Output:        {analysis_dir}")

    result = run_footprint_analysis(
        footprint_dir=rescore_dir,
        output_dir=analysis_dir,
        receptor_mol2=receptor_mol2,
        receptor_pdb=receptor_pdb,
        pharmacophore_threshold=pharmacophore_threshold,
        energy_cutoff=energy_cutoff,
        campaign_id=campaign_id,
        reference_context_csv=reference_context_csv,
        reference_context_label=reference_context_label,
        zones=zones if zones else None,
    )

    if not result.get("success"):
        logger.error(f"Analysis failed: {result.get('error')}")
        return 1

    logger.info(f"\nNext: python 02_scripts/04c_binding_modes.py "
                f"--config 03_configs/04c_binding_modes.yaml "
                f"--campaigns {args.campaigns or '<campaign_config.yaml>'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())