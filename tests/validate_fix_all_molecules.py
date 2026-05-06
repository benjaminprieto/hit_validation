#!/usr/bin/env python
"""
validate_fix_all_molecules.py

Valida el fix (clean-slate openbabel + -nc explícito) sobre TODAS las moléculas
de TODAS las campañas. Reporta:
  - Cuántas pasan Tier 1 con el fix nuevo (vs. fallaban antes).
  - Cuántas siguen fallando (y por qué).
  - Comparación con el preparation_status.csv original.

Output: validation_report.csv con una fila por molécula y un resumen agregado.

Uso:
    python validate_fix_all_molecules.py [--ph 6.3] [--max-workers 8] [--timeout 300]
"""
import argparse
import json
import multiprocessing as mp
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Dict, Any, Optional

import pandas as pd

CAMPAIGNS_DIR = Path("/home/bprieto/projects/hit_validation/04_data/campaigns")
RESULTS_DIR = Path("/home/bprieto/projects/hit_validation/05_results")


def run_step(cmd: list, timeout: int, cwd: Optional[str] = None) -> Dict[str, Any]:
    """Run a subprocess with timeout. Return dict with success, stdout, stderr."""
    try:
        r = subprocess.run(
            cmd, capture_output=True, text=True, timeout=timeout, cwd=cwd,
        )
        return {
            "success": r.returncode == 0,
            "returncode": r.returncode,
            "stdout": r.stdout,
            "stderr": r.stderr,
        }
    except subprocess.TimeoutExpired:
        return {"success": False, "returncode": -1, "stdout": "", "stderr": f"TIMEOUT after {timeout}s"}
    except FileNotFoundError as e:
        return {"success": False, "returncode": -1, "stdout": "", "stderr": str(e)}


def test_single_molecule(args_tuple) -> Dict[str, Any]:
    """Test the full fix pipeline on a single molecule."""
    sdf_path, campaign, ph, timeout = args_tuple
    sdf_path = Path(sdf_path)
    name = sdf_path.stem
    t0 = time.time()

    result = {
        "campaign": campaign,
        "ligand": name,
        "input_format": sdf_path.suffix.lower().lstrip("."),
        "tier1_with_fix": False,
        "formal_charge_pybel": None,
        "bcc_charge_sum": None,
        "n_atoms_input": None,
        "elapsed_s": None,
        "step": None,  # which step failed
        "error": None,
    }

    try:
        with tempfile.TemporaryDirectory() as td:
            wd = Path(td)

            # Step 1: Convert input to mol2 if needed
            input_mol2 = wd / "step1.mol2"
            if sdf_path.suffix.lower() == ".sdf":
                r = run_step(
                    ["obabel", f"-i{sdf_path.suffix.lstrip('.')}", str(sdf_path),
                     "-omol2", "-O", str(input_mol2)],
                    timeout=60,
                )
                if not r["success"] or not input_mol2.exists():
                    result["step"] = "obabel_convert"
                    result["error"] = (r["stderr"] or "no output")[:200]
                    return result
            else:
                shutil.copy2(sdf_path, input_mol2)

            # Step 2: Clean-slate protonation via pybel API
            try:
                from openbabel import pybel
                mol = next(pybel.readfile("mol2", str(input_mol2)))
                mol.OBMol.DeleteHydrogens()
                mol.OBMol.AddHydrogens(False, True, ph)
                charge = mol.OBMol.GetTotalCharge()
                result["formal_charge_pybel"] = charge
                result["n_atoms_input"] = mol.OBMol.NumAtoms()
                step2_mol2 = wd / "step2_clean_slate.mol2"
                mol.write("mol2", str(step2_mol2), overwrite=True)
            except Exception as e:
                result["step"] = "pybel_clean_slate"
                result["error"] = str(e)[:200]
                return result

            if not step2_mol2.exists() or step2_mol2.stat().st_size == 0:
                result["step"] = "pybel_write"
                result["error"] = "Empty mol2 after pybel"
                return result

            # Step 3: Antechamber Tier 1 with -nc explicit
            output_mol2 = wd / "output.mol2"
            r = run_step(
                ["antechamber",
                 "-fi", "mol2", "-fo", "mol2",
                 "-i", str(step2_mol2), "-o", str(output_mol2),
                 "-c", "bcc", "-at", "sybyl",
                 "-pf", "y", "-dr", "no",
                 "-nc", str(charge)],
                timeout=timeout,
                cwd=str(wd),
            )

            if not r["success"] or not output_mol2.exists() or output_mol2.stat().st_size == 0:
                result["step"] = "antechamber_tier1"
                result["error"] = (r["stderr"] or "")[:300]
                return result

            # Step 4: Verify BCC charge sum
            sum_q = 0.0
            atom_block = False
            for line in output_mol2.read_text().splitlines():
                if line.startswith("@<TRIPOS>ATOM"):
                    atom_block = True
                    continue
                if line.startswith("@<TRIPOS>") and atom_block:
                    break
                if atom_block:
                    parts = line.split()
                    if len(parts) >= 9:
                        try:
                            sum_q += float(parts[8])
                        except ValueError:
                            pass
            result["bcc_charge_sum"] = round(sum_q, 4)
            result["tier1_with_fix"] = True

    except Exception as e:
        result["step"] = "outer_exception"
        result["error"] = str(e)[:200]
    finally:
        result["elapsed_s"] = round(time.time() - t0, 2)

    return result


def get_original_status(campaign: str, ligand: str) -> Optional[str]:
    """Read tier from preparation_status.csv if available."""
    status_csv = RESULTS_DIR / campaign / "00a_ligand_preparation" / "preparation_status.csv"
    if not status_csv.exists():
        return None
    try:
        df = pd.read_csv(status_csv)
        # Match by Name column (with or without -01 suffix variations)
        for _, row in df.iterrows():
            if row.get("Name") == ligand or str(row.get("Name", "")).startswith(ligand):
                return f"tier{row.get('tier', '?')}_{row.get('status', '?')}"
        return "not_in_status"
    except Exception:
        return None


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ph", type=float, default=6.3)
    p.add_argument("--max-workers", type=int, default=8,
                   help="Parallel processes (default 8). Set to 1 for serial.")
    p.add_argument("--timeout", type=int, default=300,
                   help="Per-molecule antechamber timeout in seconds.")
    p.add_argument("--out-csv", default="validation_report.csv")
    p.add_argument("--campaigns", nargs="*", default=None,
                   help="Specific campaigns to test (default: all).")
    args = p.parse_args()

    # Build job list
    jobs = []
    if args.campaigns:
        camp_dirs = [CAMPAIGNS_DIR / c for c in args.campaigns]
    else:
        camp_dirs = sorted([d for d in CAMPAIGNS_DIR.iterdir() if d.is_dir()])

    for camp_dir in camp_dirs:
        ligands_dir = camp_dir / "ligands"
        if not ligands_dir.is_dir():
            continue
        for sdf in sorted(ligands_dir.glob("*.sdf")):
            jobs.append((str(sdf), camp_dir.name, args.ph, args.timeout))
        for mol2 in sorted(ligands_dir.glob("*.mol2")):
            jobs.append((str(mol2), camp_dir.name, args.ph, args.timeout))

    print(f"Total jobs: {len(jobs)}")
    print(f"Workers: {args.max_workers}, timeout per molecule: {args.timeout}s")
    print(f"Estimated max time: {len(jobs) * args.timeout / args.max_workers / 60:.1f} min (worst case)")
    print()

    t0 = time.time()
    if args.max_workers > 1:
        with mp.Pool(args.max_workers) as pool:
            results = []
            for i, r in enumerate(pool.imap_unordered(test_single_molecule, jobs), 1):
                results.append(r)
                if i % 10 == 0 or i == len(jobs):
                    elapsed = time.time() - t0
                    rate = i / elapsed
                    remaining = (len(jobs) - i) / rate if rate > 0 else 0
                    print(f"  [{i}/{len(jobs)}] elapsed {elapsed:.0f}s, ETA {remaining:.0f}s", flush=True)
    else:
        results = []
        for i, job in enumerate(jobs, 1):
            print(f"  [{i}/{len(jobs)}] {Path(job[0]).stem}", flush=True)
            results.append(test_single_molecule(job))

    # Cross-reference with original preparation_status.csv
    for r in results:
        r["original_status"] = get_original_status(r["campaign"], r["ligand"])

    df = pd.DataFrame(results)
    df.to_csv(args.out_csv, index=False)

    # SUMMARY
    print()
    print("=" * 70)
    print("RESUMEN GLOBAL")
    print("=" * 70)
    n_total = len(df)
    n_pass = df["tier1_with_fix"].sum()
    n_fail = n_total - n_pass
    print(f"Total moléculas testeadas: {n_total}")
    print(f"  Tier 1 BCC con fix: {n_pass}  ({100*n_pass/n_total:.1f}%)")
    print(f"  Falla aún con fix:  {n_fail}  ({100*n_fail/n_total:.1f}%)")

    print()
    print("--- Distribución de carga formal detectada por pybel ---")
    print(df["formal_charge_pybel"].value_counts(dropna=False).sort_index().to_string())

    print()
    print("--- Por campaña: éxito del fix ---")
    by_camp = df.groupby("campaign")["tier1_with_fix"].agg(["sum", "count"])
    by_camp["pct"] = (100 * by_camp["sum"] / by_camp["count"]).round(1)
    print(by_camp.to_string())

    print()
    print("--- Comparación con status original (si disponible) ---")
    if "original_status" in df.columns:
        ct = pd.crosstab(
            df["original_status"].fillna("unknown"),
            df["tier1_with_fix"].map({True: "PASS_with_fix", False: "FAIL_with_fix"}),
        )
        print(ct.to_string())

    if n_fail > 0:
        print()
        print("--- Moléculas que fallan AÚN con fix (top 20) ---")
        fail_df = df[~df["tier1_with_fix"]].head(20)
        cols = ["campaign", "ligand", "step", "error", "elapsed_s"]
        print(fail_df[cols].to_string(index=False))

    print()
    print(f"CSV completo: {args.out_csv}")
    print(f"Tiempo total: {(time.time() - t0)/60:.1f} min")


if __name__ == "__main__":
    main()
