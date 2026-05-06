#!/usr/bin/env python3
"""
01c DOCK6 Run - CLI (Rigid Re-Docking)
=========================================
Ejecuta DOCK6 rigid re-docking del ligando cristalográfico contra los grids.

Reads campaign_config.yaml:
    - grids.grid_dir             -> directorio con grids
    - grids.spheres_file         -> spheres_ligand.sph
    - grids.energy_grid          -> ligand.nrg
    - grids.bump_grid            -> ligand.bmp
    - ligand_mol2                -> path al mol2 del ligando cristalográfico

Reads module YAML (03_configs/01c_dock6_run.yaml):
    - search_method, max_orientations, num_scored_conformers, etc.

Grid routing:
    1. Busca grids en 05_results/{campaign_id}/01b_grid_generation/
    2. Si no existen, busca en campaign_config.grids.grid_dir

Ligand routing (reference docking):
    1. campaign_config.ligand_mol2 (crystallographic ligand mol2)
    2. Fallback: campaign_dir/ligands/*.mol2

Output: 05_results/{campaign_id}/01c_dock6_run/

Usage:
    python 02_scripts/01c_dock6_run.py --config 03_configs/01c_dock6_run.yaml \
        --campaigns 04_data/campaigns/NAMIKI_top20_pH63/campaign_config.yaml

Project: hit_validation
Module: 01c (DOCK6 rigid re-docking)
Version: 4.0 — adapted for reference docking (2026-03-25)
"""

import argparse
import json
import logging
import secrets
import sys
import yaml
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-8s | %(message)s",
)

sys.path.insert(0, str(Path(__file__).parent.parent / "01_src"))

from hit_validation.m01_docking.dock6_run import run_dock6_batch, resolve_grid_prefix

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


def run_consolidation_mode(args) -> int:
    """Consolidate DOCK6 per-mol outputs across N replicas via symlinks."""
    from hit_validation.utils.replica_consolidation_helpers import (
        discover_replicas, link_or_copy, load_representative_per_mol,
    )

    cc = load_yaml(args.campaigns)
    campaign_dir = Path(args.campaigns).parent
    campaign_id = cc.get("campaign_id", campaign_dir.name)
    shared_base = Path("05_results") / campaign_id

    replicas = discover_replicas(shared_base, args.n_replicas)
    if len(replicas) < 2:
        logger.error(f"Need at least 2 replicas, found {len(replicas)}")
        return 1

    rep_csv = shared_base / "_replica_metadata" / "representative_per_mol.csv"
    rep_per_mol = load_representative_per_mol(rep_csv)

    output_dir = shared_base / "01c_dock6_run" / "consolidated"
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Consolidating DOCK6 from {len(replicas)} replicas -> {output_dir}")

    n_linked = 0
    for mol, rep_id in rep_per_mol.items():
        src_dir = shared_base / f"replica_{rep_id}" / "01c_dock6_run" / mol
        if not src_dir.exists():
            logger.warning(f"Missing DOCK6 output for {mol} in replica_{rep_id}: {src_dir}")
            continue
        dst_dir = output_dir / mol
        link_or_copy(src_dir, dst_dir)
        n_linked += 1

    logger.info(f"DOCK6 consolidation: linked {n_linked}/{len(rep_per_mol)} molecules")
    return 0


def main():
    parser = argparse.ArgumentParser(
        description="01c DOCK6 Run — rigid re-docking of crystallographic ligand",
    )
    parser.add_argument("--config", "-c", type=str, required=True,
                        help="Module config YAML")
    parser.add_argument("--campaigns", type=str, required=True,
                        help="Campaign config YAML")
    parser.add_argument("--output", "-o", type=str, default=None)
    parser.add_argument("--name", type=str, default=None,
                        help="Dock only this molecule")
    parser.add_argument("--method", type=str, choices=["flex", "rigid"], default=None)
    parser.add_argument("--orientations", type=int, default=None)
    parser.add_argument("--timeout", type=int, default=None)
    parser.add_argument("--dry-run", action="store_true",
                        help="Generate input files only")
    parser.add_argument("--log-level", type=str, default=None,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    parser.add_argument("--replica-id", type=int, default=None,
                        help="Replica ID (1-indexed). None=legacy non-replicated mode.")
    parser.add_argument("--consolidate-replicas", action="store_true",
                        help="Consolidate per-mol DOCK6 outputs from N replicas into "
                             "01c_dock6_run/consolidated/ via symlinks to representative.")
    parser.add_argument("--n-replicas", type=int, default=1,
                        help="Number of replicas to consolidate. Used only with --consolidate-replicas.")

    args = parser.parse_args()

    if args.consolidate_replicas:
        if args.n_replicas < 2:
            parser.error("--consolidate-replicas requires --n-replicas >= 2")
        if args.replica_id is not None:
            parser.error("--consolidate-replicas cannot be combined with --replica-id")
        return run_consolidation_mode(args)

    # =========================================================================
    # RESOLVE PARAMETERS
    # =========================================================================

    cc = load_yaml(args.campaigns)
    campaign_dir = Path(args.campaigns).parent
    campaign_id = cc.get("campaign_id", campaign_dir.name)

    # --- Path resolution: shared_base vs results_base ---
    shared_base = Path("05_results") / campaign_id
    results_base = shared_base
    if args.replica_id is not None:
        results_base = shared_base / f"replica_{args.replica_id}"

    # --- Grid routing (SHARED input) ---
    gc = cc.get("grids", {})
    grid_dir_01b = shared_base / "01b_grid_generation"

    if (grid_dir_01b / gc.get("spheres_file", "spheres_ligand.sph")).exists():
        grid_path = grid_dir_01b
        logger.info(f"Using grids from 01b: {grid_path}")
    else:
        grid_dir = gc.get("grid_dir", "grids/")
        grid_path = Path(grid_dir) if Path(grid_dir).is_absolute() else campaign_dir / grid_dir
        if (grid_path / gc.get("spheres_file", "spheres_ligand.sph")).exists():
            logger.info(f"Using grids from campaigns: {grid_path}")
        else:
            logger.warning(f"Grids not found in 01b ({grid_dir_01b}) or campaigns ({grid_path})")

    spheres_file = str(grid_path / gc.get("spheres_file", "spheres_ligand.sph"))
    grid_prefix = resolve_grid_prefix(
        str(grid_path), gc.get("energy_grid", "ligand.nrg"),
    )

    # --- Ligand routing (SHARED input — 00a output) ---
    # Priority: 00a output > campaign_config.ligand_mol2 > campaign_dir/ligands/
    ligand_dir = None

    # 1. Output from 00a_ligand_preparation (primary source in hit_validation)
    ligand_dir_00a = shared_base / "00a_ligand_preparation"
    if ligand_dir_00a.exists() and list(ligand_dir_00a.glob("*/*.mol2")):
        ligand_dir = str(ligand_dir_00a)
        logger.info(f"Using ligands from 00a: {ligand_dir}")

    # 2. campaign_config.ligand_mol2
    if not ligand_dir:
        ligand_mol2_cfg = cc.get("ligand_mol2")
        if ligand_mol2_cfg:
            ligand_mol2_path = Path(ligand_mol2_cfg) if Path(ligand_mol2_cfg).is_absolute() else campaign_dir / ligand_mol2_cfg
            if ligand_mol2_path.is_dir():
                ligand_dir = str(ligand_mol2_path)
            elif ligand_mol2_path.is_file():
                ligand_dir = str(ligand_mol2_path.parent)

    # 3. campaign_dir/ligands/ or campaign_dir/reference/
    if not ligand_dir:
        for candidate in [campaign_dir / "ligands", campaign_dir / "reference"]:
            if candidate.exists() and list(candidate.glob("*.mol2")):
                ligand_dir = str(candidate)
                break

    if not ligand_dir:
        ligand_dir = str(ligand_dir_00a)  # Will fail with clear error below

    # --- Output (REPLICATED) ---
    output_dir = args.output or str(results_base / "01c_dock6_run")

    # --- Module config ---
    mc = load_yaml(args.config)
    params = mc.get("parameters", {})

    search_method = params.get("search_method", "rigid")
    max_orientations = params.get("max_orientations", 1000)
    num_scored_conformers = params.get("num_scored_conformers", 20)
    minimize = params.get("minimize", True)
    simplex_max_iterations = params.get("simplex_max_iterations", 500)
    timeout_per_molecule = params.get("timeout_per_molecule", 600)
    log_level = params.get("log_level", "INFO")

    # Extra DOCK6 params (grid scoring + minimization only, no secondary scores)
    extra_params = {}
    for key in ["min_anchor_size", "pruning_max_orients", "pruning_clustering_cutoff",
                "pruning_conformer_score_cutoff", "simplex_max_cycles",
                "simplex_score_converge", "simplex_cycle_converge",
                "simplex_trans_step", "simplex_rot_step", "simplex_tors_step",
                "simplex_random_seed",
                "num_final_scored_poses", "num_preclustered_conformers",
                "write_orientations"]:
        if key in params:
            extra_params[key] = params[key]

    # ---- Explicit DOCK6 seed (Fase D / Sección 5) ----
    # Always generate an explicit, persisted seed when not pinned in the YAML.
    # This makes the simplex non-deterministic across replicas (good) but
    # *reproducible* per replica (recoverable from .seeds.json).
    if "simplex_random_seed" not in extra_params or extra_params["simplex_random_seed"] in (0, -1):
        extra_params["simplex_random_seed"] = secrets.randbits(31)

    seeds_path = Path(results_base) / ".seeds.json"
    seeds_path.parent.mkdir(parents=True, exist_ok=True)
    existing_seeds: dict = {}
    if seeds_path.exists():
        try:
            existing_seeds = json.loads(seeds_path.read_text())
        except Exception:
            existing_seeds = {}
    existing_seeds["dock6_simplex_random_seed"] = int(extra_params["simplex_random_seed"])
    seeds_path.write_text(json.dumps(existing_seeds, indent=2))
    logger.info(
        f"DOCK6 simplex_random_seed = {extra_params['simplex_random_seed']} "
        f"(persisted to {seeds_path})"
    )

    # CLI overrides
    if args.method:
        search_method = args.method
    if args.orientations:
        max_orientations = args.orientations
    if args.timeout:
        timeout_per_molecule = args.timeout
    if args.log_level:
        log_level = args.log_level

    molecule_filter = [args.name] if args.name else None

    # =========================================================================
    # VALIDATE
    # =========================================================================

    if not Path(ligand_dir).exists():
        logger.error(f"Ligand directory not found: {ligand_dir}")
        logger.error("Run 00a_ligand_preparation first, or place mol2 files in campaign ligands/ dir.")
        return 1

    # =========================================================================
    # SETUP LOGGING
    # =========================================================================

    if log_level:
        logging.getLogger().setLevel(getattr(logging, log_level.upper(), logging.INFO))

    log_path = Path(output_dir) / "01c_dock6_run.log"
    setup_log_file(log_path, log_level)

    # =========================================================================
    # EXECUTE
    # =========================================================================

    logger.info("=" * 60)
    logger.info("  HIT_VALIDATION - Module 01c: DOCK6 Rigid Re-Docking")
    logger.info("=" * 60)
    logger.info(f"Campaign:      {campaign_id}")
    logger.info(f"Ligands:       {ligand_dir}")
    logger.info(f"Grid prefix:   {grid_prefix}")
    logger.info(f"Spheres:       {spheres_file}")
    logger.info(f"Method:        {search_method}")
    logger.info(f"Orientations:  {max_orientations}")
    logger.info(f"Minimize:      {minimize} (max iter: {simplex_max_iterations})")
    logger.info(f"Poses:         {num_scored_conformers}")
    logger.info(f"Timeout:       {timeout_per_molecule}s per molecule")
    logger.info("Scoring:       grid_score primary only (01d=footprint, 01f=GBSA rescore later)")
    if molecule_filter:
        logger.info(f"Filter:        {molecule_filter}")
    if args.dry_run:
        logger.info("*** DRY RUN ***")

    result = run_dock6_batch(
        ligand_mol2_dir=ligand_dir,
        spheres_file=spheres_file,
        grid_prefix=grid_prefix,
        output_dir=output_dir,
        search_method=search_method,
        max_orientations=max_orientations,
        num_scored_conformers=num_scored_conformers,
        minimize=minimize,
        simplex_max_iterations=simplex_max_iterations,
        timeout_per_molecule=timeout_per_molecule,
        molecule_filter=molecule_filter,
        dry_run=args.dry_run,
        **extra_params,
    )

    if result.get("error"):
        logger.error(f"Error: {result['error']}")
        return 1

    logger.info("")
    logger.info(f"{'=' * 60}")
    if args.dry_run:
        logger.info(f"  DRY RUN: {result['n_dry_run']} input files generated")
    else:
        logger.info(f"  {result['n_ok']}/{result['n_total']} dockings completed "
                     f"({result['total_runtime_sec']:.0f}s)")
    logger.info(f"{'=' * 60}")

    logger.info(f"Next: python 02_scripts/01d_footprint_rescore.py "
                f"--config 03_configs/01d_footprint_rescore.yaml "
                f"--campaigns {args.campaigns}")

    return 0 if result["n_failed"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
