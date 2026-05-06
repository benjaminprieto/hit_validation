#!/usr/bin/env python
"""
audit_protonation_states_v2.py

Versión que reproduce el pipeline real: SDF → obabel -p <pH> → mol2 → audit.
Necesario porque el problema de paridad NO está en el SDF crudo (que viene
neutro con todos los OHs), sino en lo que obabel produce después de protonar
a pH operativo.

Uso:
    python audit_protonation_states_v2.py <ligand_dir> [--ph 6.3] [--out audit.csv]
"""
import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

IONIZABLE_GROUPS = {
    "phosphate_monoester_OH":  ("[OX2H1][PX4](=[OX1])([OX1,OX2])[OX2][#6]", 6.7),
    "phosphate_monoester_O-":  ("[OX1-][PX4](=[OX1])([OX1-,OX2])[OX2][#6]", 6.7),
    "phosphonate_OH":          ("[OX2H1][PX4](=[OX1])([OX1,OX2])[#6]",      7.0),
    "phosphonate_O-":          ("[OX1-][PX4](=[OX1])([OX1-,OX2])[#6]",       7.0),
    "carboxylic_acid":         ("[CX3](=O)[OX2H1]",                          4.5),
    "carboxylate":             ("[CX3](=O)[OX1-]",                           4.5),
    "imidazole_neutral":       ("[nX2H0]1[cX3H][nX3H1][cX3H][cX3H]1",        6.0),
    "imidazolium":             ("[nX3H1+]1[cX3H][nX3H1][cX3H][cX3H]1",       6.0),
    "primary_amine":           ("[NX3;H2;!$(NC=O);!$(N=*);!$([N+])]",        9.5),
    "secondary_amine":         ("[NX3;H1;!$(NC=O);!$(N=*);!$([N+])]",       10.5),
    "ammonium":                ("[NX4+;!$([N+]=*)]",                         9.5),
    "pyridine":                ("[nX2;R5,R6;!$(n[#1])]",                     5.2),
    "pyridinium":              ("[nX3H1+;R5,R6]",                            5.2),
    "tetrazole_neutral":       ("c1[nH]nnn1",                                4.9),
    "tetrazolate":             ("c1[n-]nnn1",                                4.9),
    "thiol":                   ("[SX2H1][#6]",                               9.0),
    "thiolate":                ("[SX1-][#6]",                                9.0),
}

def count_electrons(mol: Chem.Mol) -> int:
    z = sum(a.GetAtomicNum() for a in mol.GetAtoms())
    h_implicit = sum(a.GetTotalNumHs() - a.GetNumExplicitHs() for a in mol.GetAtoms())
    return z + h_implicit - Chem.GetFormalCharge(mol)

def run_obabel_protonate(sdf_path: Path, ph: float, workdir: Path) -> Path:
    """Reproduce el paso de obabel -p del pipeline real."""
    # Paso 1: SDF → mol2 (sin protonación todavía)
    step1 = workdir / "step1_input.mol2"
    r1 = subprocess.run(
        ["obabel", "-isdf", str(sdf_path), "-omol2", "-O", str(step1)],
        capture_output=True, text=True
    )
    if r1.returncode != 0 or not step1.exists():
        raise RuntimeError(f"obabel step 1 failed: {r1.stderr}")

    # Paso 2: protonar a pH operativo
    step2 = workdir / "step2_protonated.mol2"
    r2 = subprocess.run(
        ["obabel", str(step1), "-O", str(step2), "-p", str(ph)],
        capture_output=True, text=True
    )
    if r2.returncode != 0 or not step2.exists():
        raise RuntimeError(f"obabel protonation failed: {r2.stderr}")

    return step2

def audit_molecule(sdf_path: Path, ph: float, borderline_window: float = 1.5) -> dict:
    row = {"ligand": sdf_path.stem}

    with tempfile.TemporaryDirectory() as td:
        workdir = Path(td)
        try:
            mol2_path = run_obabel_protonate(sdf_path, ph, workdir)
        except Exception as e:
            row["error"] = f"obabel failed: {e}"
            return row

        # Audit del mol2 protonado
        mol = Chem.MolFromMol2File(str(mol2_path), removeHs=False, sanitize=False)
        if mol is None:
            row["error"] = "RDKit could not parse protonated mol2"
            return row
        try:
            Chem.SanitizeMol(mol)
        except Exception as e:
            row["error"] = f"Sanitize failed: {e}"
            return row

        row["smiles_protonated"] = Chem.MolToSmiles(mol)
        row["formal_charge"] = Chem.GetFormalCharge(mol)
        row["n_electrons"] = count_electrons(mol)
        row["n_heavy_atoms"] = mol.GetNumHeavyAtoms()
        row["parity_ok"] = (row["n_electrons"] % 2 == 0)

        borderline_count = 0
        borderline_groups = []
        for name, (smarts, pka) in IONIZABLE_GROUPS.items():
            patt = Chem.MolFromSmarts(smarts)
            if patt is None:
                continue
            n = len(mol.GetSubstructMatches(patt))
            row[name] = n
            if n > 0 and abs(ph - pka) < borderline_window:
                borderline_count += n
                borderline_groups.append(f"{name}({n}x,pKa={pka})")

        row["borderline_total"] = borderline_count
        row["borderline_groups"] = "; ".join(borderline_groups) if borderline_groups else ""
    return row

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("ligand_dir", type=Path)
    parser.add_argument("--ph", type=float, default=6.3)
    parser.add_argument("--window", type=float, default=1.5)
    parser.add_argument("--out", type=Path, default=Path("ionizable_audit_v2.csv"))
    args = parser.parse_args()

    if not args.ligand_dir.is_dir():
        sys.exit(f"ERROR: {args.ligand_dir} no es directorio")

    sdfs = sorted(args.ligand_dir.glob("*.sdf"))
    if not sdfs:
        sys.exit(f"ERROR: no .sdf en {args.ligand_dir}")

    if shutil.which("obabel") is None:
        sys.exit("ERROR: obabel no está en PATH. Activá hit_validation_env.")

    print(f"Auditando {len(sdfs)} ligandos POST-obabel a pH {args.ph}\n")
    results = [audit_molecule(s, args.ph, args.window) for s in sdfs]
    df = pd.DataFrame(results)
    df.to_csv(args.out, index=False)

    print("=" * 70)
    print("RESUMEN (estado POST-obabel, lo que va a antechamber)")
    print("=" * 70)

    if "error" in df.columns:
        n_err = df["error"].notna().sum()
        if n_err:
            print(f"\n[!] {n_err} ligandos fallaron el paso obabel:")
            print(df[df["error"].notna()][["ligand", "error"]].to_string(index=False))

    valid = df[df.get("error", pd.Series([None]*len(df))).isna()] if "error" in df.columns else df

    print(f"\n--- Paridad electrónica POST-obabel ---")
    print(f"  Par (sqm OK):     {valid['parity_ok'].sum()} / {len(valid)}")
    print(f"  Impar (sqm fail): {(~valid['parity_ok']).sum()} / {len(valid)}")

    if (~valid["parity_ok"]).any():
        print(f"\n[!] Moléculas que VAN A FALLAR Tier 1 (paridad impar post-obabel):")
        cols = ["ligand", "formal_charge", "n_electrons", "borderline_groups"]
        print(valid[~valid["parity_ok"]][cols].to_string(index=False))

    print(f"\n--- Distribución de grupos borderline ---")
    print(valid["borderline_total"].value_counts().sort_index().to_string())

    multi = valid[valid["borderline_total"] > 1]
    if len(multi):
        print(f"\n[!] >1 grupo borderline (caso combinatorio):")
        print(multi[["ligand", "borderline_total", "borderline_groups"]].to_string(index=False))
    else:
        print(f"\n  ✓ Ninguna molécula tiene >1 grupo borderline → fix simple alcanza")

    print(f"\nCSV: {args.out}")

if __name__ == "__main__":
    main()