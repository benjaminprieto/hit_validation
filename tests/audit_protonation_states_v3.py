#!/usr/bin/env python
"""
audit_protonation_states_v3.py

Versión corregida. Usa antechamber como ground truth para paridad
electrónica en vez de RDKit (que tiene bugs leyendo cargas formales
de mol2 producidos por obabel).

Pipeline replicado:
  SDF → obabel (mol2) → obabel -p <pH> (protonated mol2)
       → antechamber dry-run → captura "Total electrons" y "net charge"

Uso:
    python audit_protonation_states_v3.py <ligand_dir> [--ph 6.3] [--out audit.csv]
"""
import argparse
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

ELEC_RE  = re.compile(r"Total number of electrons:\s*(\d+);\s*net charge:\s*(-?\d+)")
ODDFLAG_RE = re.compile(r"number of electrons is odd", re.IGNORECASE)

def run_obabel_protonate(sdf_path: Path, ph: float, workdir: Path) -> Path | None:
    step1 = workdir / "step1_input.mol2"
    r1 = subprocess.run(
        ["obabel", "-isdf", str(sdf_path), "-omol2", "-O", str(step1)],
        capture_output=True, text=True
    )
    if r1.returncode != 0 or not step1.exists():
        return None
    step2 = workdir / "step2_protonated.mol2"
    r2 = subprocess.run(
        ["obabel", str(step1), "-O", str(step2), "-p", str(ph)],
        capture_output=True, text=True
    )
    if r2.returncode != 0 or not step2.exists():
        return None
    return step2

def antechamber_paridad(mol2_path: Path, workdir: Path) -> dict:
    """Corre antechamber Tier 1 y captura el reporte de paridad."""
    out = subprocess.run(
        ["antechamber", "-fi", "mol2", "-fo", "mol2",
         "-i", str(mol2_path), "-o", str(workdir / "out.mol2"),
         "-c", "bcc", "-at", "sybyl", "-dr", "no"],
        capture_output=True, text=True, cwd=str(workdir)
    )
    text = out.stdout + "\n" + out.stderr

    info = {"antechamber_succeeded": (out.returncode == 0)}
    m = ELEC_RE.search(text)
    if m:
        info["n_electrons"] = int(m.group(1))
        info["net_charge"] = int(m.group(2))
        info["parity_ok"] = (info["n_electrons"] % 2 == 0)
    else:
        info["n_electrons"] = None
        info["net_charge"] = None
        info["parity_ok"] = None
    info["odd_warning"] = bool(ODDFLAG_RE.search(text))

    # extraer carga neta real desde el mol2 protonado (más confiable que el regex)
    try:
        sum_q = 0.0
        atom_block = False
        for line in mol2_path.read_text().splitlines():
            if line.startswith("@<TRIPOS>ATOM"):
                atom_block = True
                continue
            if line.startswith("@<TRIPOS>") and atom_block:
                break
            if atom_block:
                parts = line.split()
                if len(parts) >= 9:
                    sum_q += float(parts[8])
        info["mol2_charge_sum"] = round(sum_q, 4)
    except Exception:
        info["mol2_charge_sum"] = None

    return info

def audit_molecule(sdf_path: Path, ph: float) -> dict:
    row = {"ligand": sdf_path.stem}
    with tempfile.TemporaryDirectory() as td:
        wd = Path(td)
        protonated = run_obabel_protonate(sdf_path, ph, wd)
        if protonated is None:
            row["error"] = "obabel failed"
            return row
        info = antechamber_paridad(protonated, wd)
        row.update(info)
    return row

def main():
    p = argparse.ArgumentParser()
    p.add_argument("ligand_dir", type=Path)
    p.add_argument("--ph", type=float, default=6.3)
    p.add_argument("--out", type=Path, default=Path("audit_v3.csv"))
    args = p.parse_args()

    if not args.ligand_dir.is_dir():
        sys.exit(f"ERROR: {args.ligand_dir} no es directorio")
    sdfs = sorted(args.ligand_dir.glob("*.sdf"))
    if not sdfs:
        sys.exit(f"ERROR: no .sdf en {args.ligand_dir}")
    if shutil.which("obabel") is None or shutil.which("antechamber") is None:
        sys.exit("ERROR: obabel y/o antechamber no en PATH. Activá hit_validation_env.")

    print(f"Auditando {len(sdfs)} ligandos a pH {args.ph}")
    print(f"(usando antechamber como ground truth para paridad)\n")

    results = []
    for i, sdf in enumerate(sdfs, 1):
        print(f"  [{i}/{len(sdfs)}] {sdf.stem}", flush=True)
        results.append(audit_molecule(sdf, args.ph))

    df = pd.DataFrame(results)
    df.to_csv(args.out, index=False)

    print("\n" + "=" * 70)
    print("RESUMEN (paridad reportada por antechamber)")
    print("=" * 70)

    if "error" in df.columns and df["error"].notna().any():
        n_err = df["error"].notna().sum()
        print(f"\n[!] {n_err} ligandos fallaron obabel:")
        print(df[df["error"].notna()][["ligand", "error"]].to_string(index=False))

    valid = df[df.get("error", pd.Series([None]*len(df))).isna()] if "error" in df.columns else df

    # Filas donde antechamber pudo reportar paridad
    has_parity = valid[valid["parity_ok"].notna()]
    no_parity = valid[valid["parity_ok"].isna()]

    if len(no_parity):
        print(f"\n[?] {len(no_parity)} ligandos: antechamber no reportó electrones (raro):")
        print(no_parity[["ligand", "antechamber_succeeded"]].to_string(index=False))

    print(f"\n--- Paridad post-obabel (según antechamber) ---")
    print(f"  Par (sqm OK):     {has_parity['parity_ok'].sum()} / {len(has_parity)}")
    print(f"  Impar (sqm fail): {(~has_parity['parity_ok']).sum()} / {len(has_parity)}")

    will_fail = has_parity[~has_parity["parity_ok"]]
    if len(will_fail):
        print(f"\n[!] Moléculas que VAN A FALLAR Tier 1:")
        cols = ["ligand", "n_electrons", "net_charge", "mol2_charge_sum"]
        print(will_fail[cols].to_string(index=False))

    # Distribución de carga neta de los que sí pasan
    if len(has_parity):
        print(f"\n--- Distribución de carga neta (post-obabel) ---")
        print(has_parity["net_charge"].value_counts().sort_index().to_string())

    print(f"\nCSV: {args.out}")

if __name__ == "__main__":
    main()