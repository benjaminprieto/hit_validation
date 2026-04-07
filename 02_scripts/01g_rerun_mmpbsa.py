#!/usr/bin/env python3
"""
01g Rerun MMPBSA with keep_files=2 — Auxiliary Script
======================================================
Re-runs MMPBSA.py on existing topologies and trajectories to preserve
per-frame decomposition files (_MMPBSA_*). Does NOT re-run MD or
regenerate topologies.

Usage:
    # All molecules in campaign:
    python 02_scripts/01g_rerun_mmpbsa.py \
        --campaigns 04_data/campaigns/druglikeness_1/campaign_config.yaml

    # Specific molecules only:
    python 02_scripts/01g_rerun_mmpbsa.py \
        --campaigns 04_data/campaigns/druglikeness_1/campaign_config.yaml \
        --molecules PubChem-146072367 PubChem-59068465

Project: hit_validation
Module: 01g (auxiliary)
"""

import argparse
import logging
import subprocess
import sys
import time
import yaml
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-8s | %(message)s")
logger = logging.getLogger(__name__)

MMPBSA_INPUT = """\
&general
  startframe=1, endframe=99999, interval=1,
  verbose=2, keep_files=2,
/
&gb
  igb={igb}, saltcon={saltcon},
/
&decomp
  idecomp={idecomp},
  dec_verbose=3,
  print_res="all"
/
"""


def rerun_mmpbsa_molecule(mol_dir: Path, igb: int = 2, saltcon: float = 0.15,
                          idecomp: int = 1) -> dict:
    """Re-run MMPBSA.py for a single molecule using existing topologies/trajectory."""
    topo_dir = mol_dir / "topologies"
    traj_dir = mol_dir / "trajectory"
    mmpbsa_dir = mol_dir / "mmpbsa"

    complex_prmtop = topo_dir / "complex.prmtop"
    receptor_prmtop = topo_dir / "receptor.prmtop"
    ligand_prmtop = topo_dir / "ligand.prmtop"
    trajectory = traj_dir / "dry_trajectory.mdcrd"

    # Validate required files exist
    required = [complex_prmtop, receptor_prmtop, ligand_prmtop, trajectory]
    for f in required:
        if not f.exists():
            return {"success": False, "error": f"Missing: {f}"}

    mmpbsa_dir.mkdir(parents=True, exist_ok=True)

    # Write mmpbsa.in with keep_files=2
    mmpbsa_in = mmpbsa_dir / "mmpbsa_keepfiles.in"
    mmpbsa_in.write_text(MMPBSA_INPUT.format(igb=igb, saltcon=saltcon, idecomp=idecomp))

    # Check for solvated topology
    solvated_prmtop = topo_dir / "solvated_complex.prmtop"

    cmd = [
        "MMPBSA.py", "-O",
        "-i", str(mmpbsa_in),
        "-o", str(mmpbsa_dir / "FINAL_RESULTS_MMPBSA.dat"),
        "-do", str(mmpbsa_dir / "FINAL_DECOMP_MMPBSA.dat"),
        "-cp", str(complex_prmtop.resolve()),
        "-rp", str(receptor_prmtop.resolve()),
        "-lp", str(ligand_prmtop.resolve()),
        "-y", str(trajectory.resolve()),
    ]

    if solvated_prmtop.exists():
        cmd.extend(["-sp", str(solvated_prmtop.resolve())])

    logger.info(f"  MMPBSA.py -O: re-running with keep_files=2")
    logger.debug(f"    cmd: {' '.join(cmd)}")

    t0 = time.time()
    proc = subprocess.run(
        cmd, capture_output=True, text=True,
        timeout=7200,
        cwd=str(mmpbsa_dir),
    )
    runtime = time.time() - t0

    (mmpbsa_dir / "rerun_stdout.log").write_text(proc.stdout)
    (mmpbsa_dir / "rerun_stderr.log").write_text(proc.stderr)

    if proc.returncode != 0:
        error_lines = [l for l in proc.stderr.split("\n")
                       if "error" in l.lower() or "fatal" in l.lower()]
        error_msg = "; ".join(error_lines[:3]) if error_lines else proc.stderr[-300:]
        return {"success": False, "error": f"MMPBSA.py failed ({runtime:.1f}s): {error_msg}",
                "runtime_sec": round(runtime, 1)}

    decomp_file = mmpbsa_dir / "FINAL_DECOMP_MMPBSA.dat"
    if not decomp_file.exists():
        return {"success": False, "error": "FINAL_DECOMP_MMPBSA.dat not found",
                "runtime_sec": round(runtime, 1)}

    logger.info(f"  MMPBSA.py completed in {runtime:.1f}s")
    return {"success": True, "runtime_sec": round(runtime, 1)}


def main():
    parser = argparse.ArgumentParser(
        description="Re-run MMPBSA.py with keep_files=2 to preserve per-frame files",
    )
    parser.add_argument("--campaigns", type=str, required=True, help="Campaign YAML")
    parser.add_argument("--molecules", nargs="*", default=None,
                        help="Specific molecules (default: all)")
    parser.add_argument("--igb", type=int, default=2)
    parser.add_argument("--saltcon", type=float, default=0.15)
    parser.add_argument("--idecomp", type=int, default=1)
    args = parser.parse_args()

    with open(args.campaigns, "r", encoding="utf-8") as f:
        cc = yaml.safe_load(f)

    campaign_id = cc.get("campaign_id", Path(args.campaigns).parent.name)
    decomp_dir = Path("05_results") / campaign_id / "01g_mmpbsa_decomp"

    if not decomp_dir.exists():
        logger.error(f"MMPBSA output dir not found: {decomp_dir}")
        return 1

    # Discover molecules
    if args.molecules:
        mol_dirs = [decomp_dir / m for m in args.molecules]
    else:
        mol_dirs = sorted([d for d in decomp_dir.iterdir() if d.is_dir()])

    logger.info("=" * 60)
    logger.info("  01g RERUN MMPBSA — keep_files=2")
    logger.info("=" * 60)
    logger.info(f"Campaign:   {campaign_id}")
    logger.info(f"Molecules:  {len(mol_dirs)}")

    successes = 0
    failures = 0

    for mol_dir in mol_dirs:
        name = mol_dir.name
        if not mol_dir.exists():
            logger.warning(f"[SKIP] {name}: directory not found")
            failures += 1
            continue

        logger.info(f"\n{'─' * 40}")
        logger.info(f"Processing: {name}")

        result = rerun_mmpbsa_molecule(
            mol_dir, igb=args.igb, saltcon=args.saltcon, idecomp=args.idecomp
        )

        if result["success"]:
            logger.info(f"[OK] {name} ({result['runtime_sec']}s)")
            successes += 1
        else:
            logger.error(f"[FAIL] {name}: {result['error']}")
            failures += 1

    logger.info(f"\n{'=' * 60}")
    logger.info(f"Done: {successes} OK, {failures} failed")
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    sys.exit(main())