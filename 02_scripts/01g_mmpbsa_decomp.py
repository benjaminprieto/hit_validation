#!/usr/bin/env python3
"""
01g MMPBSA Per-Residue Decomposition - CLI (Hit Validation)
==============================================================
Batch MMPBSA per-residue decomposition for all validated screening hits.
Runs sequentially (1 molecule at a time on GPU).

Two modes:
  single_point — 1-frame MMPBSA per molecule (fast, no error bars)
  md           — OpenMM MD → trajectory MMPBSA per molecule (mean ± std)

Input:  05_results/{campaign}/01c_dock6_run/{name}/{name}_scored.mol2
Output: 05_results/{campaign}/01g_mmpbsa_decomp/{name}/

Usage:
    # Full batch (all molecules, MD mode):
    python 02_scripts/01g_mmpbsa_decomp.py \\
        --config 03_configs/01g_mmpbsa_decomp.yaml \\
        --campaigns 04_data/campaigns/NAMIKI_top20_pH63/campaign_config.yaml

    # Single molecule only:
    python 02_scripts/01g_mmpbsa_decomp.py \\
        --config 03_configs/01g_mmpbsa_decomp.yaml \\
        --campaigns 04_data/campaigns/NAMIKI_top20_pH63/campaign_config.yaml \\
        --molecule MOL_001

    # Fast mode (single_point):
    python 02_scripts/01g_mmpbsa_decomp.py \\
        --config 03_configs/01g_mmpbsa_decomp.yaml \\
        --campaigns 04_data/campaigns/NAMIKI_top20_pH63/campaign_config.yaml \\
        --mode single_point

Project: hit_validation
Module: 01g
Version: 1.0 (2026-03-28)
"""

import argparse
import json
import logging
import secrets
import sys
import yaml
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-8s | %(message)s")
sys.path.insert(0, str(Path(__file__).parent.parent / "01_src"))

from hit_validation.m01_docking.mmpbsa_decomp import run_mmpbsa_batch, run_mmpbsa_decomp

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
        description="01g MMPBSA Per-Residue Decomposition — batch processing",
    )
    parser.add_argument("--config", "-c", type=str, required=True, help="Module YAML")
    parser.add_argument("--campaigns", type=str, required=True, help="Campaign YAML")
    parser.add_argument("--output", "-o", type=str, default=None)
    parser.add_argument("--molecule", type=str, default=None,
                        help="Process single molecule (default: all)")
    parser.add_argument("--mode", type=str, default=None,
                        choices=["single_point", "md"],
                        help="Override mode (default from YAML)")
    parser.add_argument("--production-ns", type=float, default=None,
                        help="MD production length in ns (overrides YAML)")
    parser.add_argument("--log-level", type=str, default=None,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    parser.add_argument("--replica-id", type=int, default=None,
                        help="Replica ID (1-indexed). None=legacy non-replicated mode.")
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
    # RESOLVE PATHS (shared_base vs results_base)
    # =========================================================================
    shared_base = Path("05_results") / campaign_id
    results_base = shared_base
    if args.replica_id is not None:
        results_base = shared_base / f"replica_{args.replica_id}"

    output_subdir = mc.get("outputs", {}).get("subdir", "01g_mmpbsa_decomp")
    output_dir = Path(args.output) if args.output else results_base / output_subdir

    # 01c is REPLICATED, receptor is SHARED
    dock6_dir = results_base / "01c_dock6_run"
    receptor_dir = shared_base / "00b_receptor_preparation"

    # =========================================================================
    # MERGE PARAMETERS (YAML + CLI overrides)
    # =========================================================================
    mode = args.mode or params.get("mode", "md")
    pose_selection = params.get("pose_selection", "best_score")
    pose_index = params.get("pose_index", 1)
    receptor_ff = params.get("receptor_ff", "ff14SB")
    ligand_ff = params.get("ligand_ff", "gaff2")
    charge_method = params.get("charge_method", "bcc")
    antechamber_timeout = params.get("antechamber_timeout", 600)

    mmpbsa_cfg = params.get("mmpbsa", {})

    # Campaign-level toggles (campaign_config.yaml `mmpbsa_decomp:` block).
    # Precedence: campaign_config > module YAML > code default.
    # Friendly key names: `mmgbsa`/`mmpbsa` (booleans) + `timeout` (seconds).
    cc_decomp = cc.get("mmpbsa_decomp", {}) or {}

    mmpbsa_idecomp = mmpbsa_cfg.get("idecomp", 1)
    mmpbsa_igb = mmpbsa_cfg.get("igb", 2)
    mmpbsa_saltcon = mmpbsa_cfg.get("saltcon", 0.15)
    mmpbsa_istrng = mmpbsa_cfg.get("istrng", 0.150)
    mmpbsa_fillratio = mmpbsa_cfg.get("fillratio", 4.0)
    mmpbsa_radiopt = mmpbsa_cfg.get("radiopt", 0)

    mmpbsa_run_gb = cc_decomp.get("mmgbsa", mmpbsa_cfg.get("mmgbsa", True))
    mmpbsa_run_pb = cc_decomp.get("mmpbsa", mmpbsa_cfg.get("mmpbsa", True))
    mmpbsa_timeout = cc_decomp.get("timeout", mmpbsa_cfg.get("timeout", 86400))

    if not (mmpbsa_run_gb or mmpbsa_run_pb):
        logger.error(
            "Both mmgbsa and mmpbsa are disabled in mmpbsa_decomp. "
            "At least one must be true. Aborting."
        )
        return 1

    # Audit which sources won (visible in run log)
    def _src(key):
        if key in cc_decomp:
            return "campaign_config"
        # The module YAML uses the same friendly keys, so we check there too.
        if key in mmpbsa_cfg:
            return "module_config"
        return "default"
    logger.info(
        f"MMPBSA config: mmgbsa={mmpbsa_run_gb} (source: {_src('mmgbsa')}), "
        f"mmpbsa={mmpbsa_run_pb} (source: {_src('mmpbsa')}), "
        f"timeout={mmpbsa_timeout}s ({mmpbsa_timeout / 3600:.1f}h, source: {_src('timeout')})"
    )

    md_params = params.get("md", {})
    if args.production_ns is not None:
        md_params["production_ns"] = args.production_ns

    # ---- Explicit OpenMM MD seed (Fase D / Sección 5) ----
    # Always generate an explicit, persisted seed when not pinned in the YAML.
    # Persisted to results_base/.seeds.json (shared with 01c entry).
    if mode == "md":
        if md_params.get("random_seed") in (None, 0, -1):
            md_params["random_seed"] = secrets.randbits(31)

        seeds_path = Path(results_base) / ".seeds.json"
        seeds_path.parent.mkdir(parents=True, exist_ok=True)
        existing_seeds: dict = {}
        if seeds_path.exists():
            try:
                existing_seeds = json.loads(seeds_path.read_text())
            except Exception:
                existing_seeds = {}
        existing_seeds["openmm_md_random_seed"] = int(md_params["random_seed"])
        seeds_path.write_text(json.dumps(existing_seeds, indent=2))
        logger.info(
            f"OpenMM MD random_seed = {md_params['random_seed']} "
            f"(persisted to {seeds_path})"
        )

    log_level = args.log_level or params.get("log_level", "INFO")
    logging.getLogger().setLevel(getattr(logging, log_level.upper()))

    # =========================================================================
    # SETUP LOGGING
    # =========================================================================
    setup_log_file(output_dir / "01g_mmpbsa_decomp.log", log_level)

    logger.info("=" * 60)
    logger.info("  HIT_VALIDATION - Module 01g: MMPBSA Decomposition")
    logger.info("=" * 60)
    logger.info(f"Campaign:   {campaign_id}")
    logger.info(f"Mode:       {mode}")
    logger.info(f"Molecule:   {args.molecule or 'ALL'}")
    logger.info(f"Output:     {output_dir}")

    # =========================================================================
    # EXECUTE
    # =========================================================================
    mmpbsa_params = dict(
        receptor_ff=receptor_ff,
        ligand_ff=ligand_ff,
        charge_method=charge_method,
        mmpbsa_idecomp=mmpbsa_idecomp,
        mmpbsa_igb=mmpbsa_igb,
        mmpbsa_saltcon=mmpbsa_saltcon,
        mmpbsa_istrng=mmpbsa_istrng,
        mmpbsa_fillratio=mmpbsa_fillratio,
        mmpbsa_radiopt=mmpbsa_radiopt,
        mmpbsa_timeout=mmpbsa_timeout,
        mmpbsa_run_gb=mmpbsa_run_gb,
        mmpbsa_run_pb=mmpbsa_run_pb,
        md_params=md_params,
        antechamber_timeout=antechamber_timeout,
    )

    if args.molecule:
        # Single molecule mode
        scored_mol2 = dock6_dir / args.molecule / f"{args.molecule}_scored.mol2"
        if not scored_mol2.exists():
            logger.error(f"Scored mol2 not found: {scored_mol2}")
            return 1

        rec_mol2 = receptor_dir / "rec_charged.mol2"
        rec_pdb = receptor_dir / "receptor_protonated.pdb"
        if not rec_pdb.exists():
            rec_pdb = receptor_dir / "rec_noH.pdb"

        mol_output = output_dir / args.molecule
        result = run_mmpbsa_decomp(
            scored_mol2=str(scored_mol2),
            receptor_mol2=str(rec_mol2),
            receptor_pdb=str(rec_pdb),
            output_dir=str(mol_output),
            mode=mode,
            pose_selection=pose_selection,
            pose_index=pose_index,
            **mmpbsa_params,
        )
    else:
        # Batch mode — pass shared/replicated paths explicitly
        molecules = None
        result = run_mmpbsa_batch(
            campaign_dir=str(campaign_dir),
            results_base=str(results_base),
            output_dir=str(output_dir),
            mode=mode,
            pose_selection=pose_selection,
            molecules=molecules,
            dock6_dir=str(dock6_dir),
            receptor_dir=str(receptor_dir),
            **mmpbsa_params,
        )

    if not result.get("success"):
        logger.error(f"Error: {result.get('error')}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
