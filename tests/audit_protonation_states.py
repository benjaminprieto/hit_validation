#!/usr/bin/env python
"""
audit_protonation_states.py

Audita los grupos ionizables presentes en un set de SDFs candidatos.
Identifica grupos cuyo pKa está cerca del pH operativo (6.3 = Golgi)
para predecir cuáles moléculas tendrán problemas de paridad electrónica
en sqm/antechamber Tier 1 (AM1-BCC).

Uso:
    python audit_protonation_states.py <ligand_dir> [--ph 6.3] [--out audit.csv]

Ejemplo:
    python audit_protonation_states.py \
        04_data/campaigns/SD1_druglikeness_2_PChem_136488320_2n/ligands \
        --ph 6.3 --out ionizable_audit.csv
"""
import argparse
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")  # callate las warnings de RDKit

# SMARTS para grupos ionizables comunes en small molecules drug-like
# pKa values: valores típicos tabulados; el pKa real depende del sustituyente
IONIZABLE_GROUPS = {
    # name: (SMARTS, pKa_aprox, ionization_at_low_pH, ionization_at_high_pH)
    "phosphate_monoester_OH":  ("[OX2H1][PX4](=[OX1])([OX1,OX2])[OX2][#6]", 6.7, "neutral_OH", "anion_O-"),
    "phosphate_monoester_O-":  ("[OX1-][PX4](=[OX1])([OX1-])[OX2][#6]",     6.7, "anion_-1",   "dianion_-2"),
    "phosphonate_OH":          ("[OX2H1][PX4](=[OX1])([OX1,OX2])[#6]",      7.0, "neutral",    "anion"),
    "carboxylic_acid":         ("[CX3](=O)[OX2H1]",                          4.5, "neutral",    "anion"),
    "imidazole_neutral":       ("[nX2H0]1[cX3H][nX3H1][cX3H][cX3H]1",        6.0, "cation",     "neutral"),
    "imidazolium":             ("[nX3H1+]1[cX3H][nX3H1][cX3H][cX3H]1",       6.0, "cation",     "neutral"),
    "primary_amine":           ("[NX3;H2;!$(NC=O);!$(N=*);!$([N+])]",        9.5, "cation",     "neutral"),
    "secondary_amine":         ("[NX3;H1;!$(NC=O);!$(N=*);!$([N+])]",       10.5, "cation",     "neutral"),
    "pyridine":                ("[nX2;R5,R6;!$(n[#1])]",                     5.2, "cation",     "neutral"),
    "tetrazole":               ("c1[nH]nnn1",                                4.9, "neutral",    "anion"),
    "thiol":                   ("[SX2H1][#6]",                               9.0, "neutral",    "anion"),
}

def count_electrons(mol: Chem.Mol) -> int:
    """Total de electrones = sum(Z) - carga formal neta (con signo)."""
    z_total = sum(a.GetAtomicNum() for a in mol.GetAtoms())
    # incluir Hs implícitos (cada H aporta 1 electrón)
    h_implicit = sum(a.GetTotalNumHs() - a.GetNumExplicitHs() for a in mol.GetAtoms())
    z_total += h_implicit
    charge = Chem.GetFormalCharge(mol)
    return z_total - charge

def audit_molecule(sdf_path: Path, ph: float, borderline_window: float = 1.5) -> dict:
    """Audita una molécula. Devuelve dict con conteos y flags."""
    mol = Chem.MolFromMolFile(str(sdf_path), removeHs=False)
    if mol is None:
        # intentar con sanitize=False
        mol = Chem.MolFromMolFile(str(sdf_path), removeHs=False, sanitize=False)
        if mol is None:
            return {"ligand": sdf_path.stem, "error": "RDKit parse failed"}
        try:
            Chem.SanitizeMol(mol)
        except Exception as e:
            return {"ligand": sdf_path.stem, "error": f"Sanitize failed: {e}"}

    row = {
        "ligand": sdf_path.stem,
        "smiles": Chem.MolToSmiles(mol),
        "formal_charge": Chem.GetFormalCharge(mol),
        "n_electrons": count_electrons(mol),
        "n_heavy_atoms": mol.GetNumHeavyAtoms(),
    }
    row["parity_ok"] = (row["n_electrons"] % 2 == 0)

    borderline_count = 0
    borderline_groups = []
    for name, (smarts, pka, low_state, high_state) in IONIZABLE_GROUPS.items():
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            print(f"  [WARN] SMARTS inválido para {name}: {smarts}", file=sys.stderr)
            continue
        n_matches = len(mol.GetSubstructMatches(patt))
        row[name] = n_matches
        if n_matches > 0 and abs(ph - pka) < borderline_window:
            borderline_count += n_matches
            borderline_groups.append(f"{name}({n_matches}x,pKa={pka})")

    row["borderline_total"] = borderline_count
    row["borderline_groups"] = "; ".join(borderline_groups) if borderline_groups else ""
    return row

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("ligand_dir", type=Path, help="Directorio con SDFs candidatos")
    parser.add_argument("--ph", type=float, default=6.3, help="pH operativo (default: 6.3, Golgi)")
    parser.add_argument("--window", type=float, default=1.5, help="Ventana borderline en unidades de pKa (default: 1.5)")
    parser.add_argument("--out", type=Path, default=Path("ionizable_audit.csv"), help="Output CSV")
    args = parser.parse_args()

    if not args.ligand_dir.is_dir():
        sys.exit(f"ERROR: {args.ligand_dir} no es un directorio válido")

    sdfs = sorted(args.ligand_dir.glob("*.sdf"))
    if not sdfs:
        sys.exit(f"ERROR: no se encontraron .sdf en {args.ligand_dir}")

    print(f"Auditando {len(sdfs)} ligandos en {args.ligand_dir}")
    print(f"pH = {args.ph}, ventana borderline = ±{args.window} unidades de pKa\n")

    results = [audit_molecule(s, args.ph, args.window) for s in sdfs]
    df = pd.DataFrame(results)
    df.to_csv(args.out, index=False)

    # Resúmenes
    print("=" * 70)
    print("RESUMEN")
    print("=" * 70)

    if "error" in df.columns:
        n_err = df["error"].notna().sum()
        if n_err:
            print(f"\n[!] {n_err} ligandos no pudieron parsearse:")
            print(df[df["error"].notna()][["ligand", "error"]].to_string(index=False))

    valid = df[df.get("error", pd.Series([None]*len(df))).isna()] if "error" in df.columns else df

    print(f"\n--- Paridad electrónica ---")
    print(f"  Par (sqm OK):    {valid['parity_ok'].sum()} / {len(valid)}")
    print(f"  Impar (sqm fail): {(~valid['parity_ok']).sum()} / {len(valid)}")

    if (~valid["parity_ok"]).any():
        print(f"\n[!] Moléculas con paridad impar (van a fallar Tier 1):")
        cols = ["ligand", "formal_charge", "n_electrons", "borderline_groups"]
        print(valid[~valid["parity_ok"]][cols].to_string(index=False))

    print(f"\n--- Distribución de grupos borderline por molécula ---")
    print(valid["borderline_total"].value_counts().sort_index().to_string())

    multi = valid[valid["borderline_total"] > 1]
    if len(multi):
        print(f"\n[!] Moléculas con >1 grupo borderline (caso combinatorio):")
        print(multi[["ligand", "borderline_total", "borderline_groups"]].to_string(index=False))
    else:
        print(f"\n  ✓ Ninguna molécula tiene >1 grupo borderline → fix simple alcanza")

    print(f"\n--- Tipos de grupos borderline encontrados ---")
    group_cols = [c for c in valid.columns if c in IONIZABLE_GROUPS]
    for col in group_cols:
        pka = IONIZABLE_GROUPS[col][1]
        if abs(args.ph - pka) < args.window:
            n_mols = (valid[col] > 0).sum()
            if n_mols:
                print(f"  {col} (pKa={pka}): presente en {n_mols} molécula(s)")

    print(f"\nCSV completo: {args.out}")

if __name__ == "__main__":
    main()