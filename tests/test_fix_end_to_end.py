#!/usr/bin/env python
"""
test_fix_end_to_end.py

Prueba real del fix aplicado al flujo completo de prepare_single_molecule().
Replica exactamente lo que hace ligand_preparation.py en producción, pero
con dos cambios:

  1. protonate_at_ph()  →  clean-slate via pybel API
  2. run_antechamber_tier1()  →  con -nc <charge> explícito

NO toca el código de producción. Solo lo replica con el fix aplicado.
Compara el output con el mol2 actualmente en producción.

Uso:
    python test_fix_end_to_end.py <campaign_name> [<molecule_name>]

Ejemplo:
    # Test sobre toda la campaña corriente (2 moléculas)
    python test_fix_end_to_end.py SD1_druglikeness_2_PChem_136488320_2n

    # O sobre una molécula específica
    python test_fix_end_to_end.py SD1_druglikeness_1 PubChem-136488320
"""
import argparse
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Dict, Any, Optional

CAMPAIGNS_DIR = Path("/home/bprieto/projects/hit_validation/04_data/campaigns")
RESULTS_DIR = Path("/home/bprieto/projects/hit_validation/05_results")


def convert_sdf_to_mol2(sdf_path: str, output_mol2: str) -> bool:
    """Replica exactamente convert_sdf_to_mol2() de ligand_preparation.py."""
    try:
        result = subprocess.run(
            ["obabel", "-isdf", sdf_path, "-omol2", "-O", output_mol2],
            capture_output=True, text=True, timeout=60,
        )
        return result.returncode == 0 and Path(output_mol2).exists()
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def protonate_at_ph_clean_slate(input_mol2: str, output_mol2: str, ph: float = 6.3) -> Optional[int]:
    """
    NUEVA versión: clean-slate via pybel API.
    Devuelve la carga formal total (para pasar a antechamber como -nc).
    Devuelve None si falla.
    """
    try:
        from openbabel import pybel
        mol = next(pybel.readfile("mol2", input_mol2))
        mol.OBMol.DeleteHydrogens()
        mol.OBMol.AddHydrogens(False, True, ph)
        charge = mol.OBMol.GetTotalCharge()
        mol.write("mol2", output_mol2, overwrite=True)
        if Path(output_mol2).exists() and Path(output_mol2).stat().st_size > 0:
            return charge
        return None
    except Exception as e:
        print(f"    [pybel error] {e}")
        return None


def run_antechamber_tier1_with_nc(input_mol2: str, output_mol2: str,
                                  formal_charge: int,
                                  timeout: int = 1200) -> Dict[str, Any]:
    """
    Replica run_antechamber_tier1() pero con -nc <formal_charge> explícito.
    Timeout default 1200s para coincidir con configuración de producción.
    """
    abs_input = str(Path(input_mol2).resolve())
    abs_output = str(Path(output_mol2).resolve())
    out_dir = str(Path(abs_output).parent)
    cmd = [
        "antechamber",
        "-fi", "mol2",
        "-fo", "mol2",
        "-i", abs_input,
        "-o", abs_output,
        "-c", "bcc",
        "-at", "sybyl",
        "-pf", "y",
        "-dr", "no",
        "-nc", str(formal_charge),  # ← CAMBIO CLAVE
    ]
    t0 = time.time()
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=timeout, cwd=out_dir,
        )
        elapsed = time.time() - t0
        success = (result.returncode == 0
                   and Path(output_mol2).exists()
                   and Path(output_mol2).stat().st_size > 0)
        return {
            "success": success,
            "tier": 1,
            "method": "AM1-BCC + Sybyl + nc explicit",
            "elapsed_s": round(elapsed, 1),
            "error": result.stderr[:300] if not success else None,
        }
    except subprocess.TimeoutExpired:
        return {"success": False, "tier": 1, "method": "AM1-BCC + Sybyl + nc explicit",
                "elapsed_s": timeout, "error": f"Timeout after {timeout}s"}


def prepare_single_molecule_with_fix(input_path: str, output_dir: Path,
                                     docking_ph: float = 6.3, timeout: int = 1200) -> Dict[str, Any]:
    """
    Replica prepare_single_molecule() de producción, pero:
      - usa protonate_at_ph_clean_slate (clean-slate)
      - usa run_antechamber_tier1_with_nc (-nc explícito)
      - NO tiene fallback a Tier 2 (queremos ver si Tier 1 funciona limpio)
    """
    inp = Path(input_path).resolve()
    out = output_dir / f"{inp.stem}.mol2"
    out.parent.mkdir(parents=True, exist_ok=True)

    result = {
        "ligand": inp.stem,
        "input_format": inp.suffix.lower().lstrip("."),
        "stages": {},
    }

    # Step 1: Convert SDF to mol2 if needed
    work_mol2 = output_dir / f"{inp.stem}_input.mol2"
    if inp.suffix.lower() in (".sdf", ".mol"):
        if not convert_sdf_to_mol2(str(inp), str(work_mol2)):
            result["success"] = False
            result["fail_at"] = "sdf_conversion"
            return result
        result["stages"]["sdf_to_mol2"] = "OK"
    else:
        shutil.copy2(str(inp), str(work_mol2))
        result["stages"]["sdf_to_mol2"] = "skipped (already mol2)"

    # Step 2: Clean-slate protonation (NEW)
    protonated_mol2 = output_dir / f"{inp.stem}_protonated.mol2"
    formal_charge = protonate_at_ph_clean_slate(str(work_mol2), str(protonated_mol2), docking_ph)
    if formal_charge is None:
        result["success"] = False
        result["fail_at"] = "protonation"
        return result
    result["stages"]["clean_slate_protonation"] = "OK"
    result["formal_charge"] = formal_charge

    # Step 3: Antechamber Tier 1 with -nc explicit (NEW)
    tier1_result = run_antechamber_tier1_with_nc(
        str(protonated_mol2), str(out), formal_charge, timeout
    )
    result["stages"]["tier1_with_nc"] = "OK" if tier1_result["success"] else "FAIL"
    result["tier1_elapsed_s"] = tier1_result["elapsed_s"]
    result["tier1_error"] = tier1_result.get("error")

    if not tier1_result["success"]:
        result["success"] = False
        result["fail_at"] = "tier1_with_nc"
        return result

    # Read final BCC charge sum
    sum_q = 0.0
    atom_block = False
    for line in out.read_text().splitlines():
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
    result["output_mol2"] = str(out)
    result["success"] = True
    return result


def find_existing_mol2(campaign: str, ligand: str) -> Optional[Path]:
    """Buscar el mol2 actual en producción (de la corrida original)."""
    candidate = RESULTS_DIR / campaign / "00a_ligand_preparation" / ligand / f"{ligand}.mol2"
    if candidate.exists():
        return candidate
    return None


def compare_mol2(fix_path: Path, original_path: Path) -> Dict[str, Any]:
    """Comparar el mol2 con fix vs el mol2 actual en producción."""
    cmp = {"original_path": str(original_path), "fix_path": str(fix_path)}

    def parse_mol2(p):
        info = {"header": "", "n_atoms": 0, "sum_charges": 0.0, "charges": []}
        atom_block = False
        for i, line in enumerate(p.read_text().splitlines()):
            if i == 5:  # línea de header con tipo de carga
                info["header"] = line.strip()
            if line.startswith("@<TRIPOS>ATOM"):
                atom_block = True
                continue
            if line.startswith("@<TRIPOS>") and atom_block:
                break
            if atom_block:
                parts = line.split()
                if len(parts) >= 9:
                    try:
                        q = float(parts[8])
                        info["charges"].append(q)
                        info["sum_charges"] += q
                        info["n_atoms"] += 1
                    except ValueError:
                        pass
        info["sum_charges"] = round(info["sum_charges"], 4)
        info["n_unique_charges"] = len(set(round(q, 3) for q in info["charges"]))
        return info

    fix_info = parse_mol2(fix_path)
    orig_info = parse_mol2(original_path)
    cmp["fix"] = fix_info
    cmp["original"] = orig_info

    # Diferencias clave
    cmp["diff_charge_method"] = (fix_info["header"], orig_info["header"])
    cmp["diff_sum_charges"] = (fix_info["sum_charges"], orig_info["sum_charges"])
    cmp["diff_n_unique"] = (fix_info["n_unique_charges"], orig_info["n_unique_charges"])
    cmp["same_n_atoms"] = (fix_info["n_atoms"] == orig_info["n_atoms"])

    return cmp


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("campaign", help="Nombre de la campaña (e.g., SD1_druglikeness_2_PChem_136488320_2n)")
    parser.add_argument("molecule", nargs="?", help="Nombre de molécula específica (opcional, default todas)")
    parser.add_argument("--ph", type=float, default=6.3)
    parser.add_argument("--timeout", type=int, default=1200,
                        help="Timeout antechamber Tier 1 (default 1200s, igual que producción)")
    args = parser.parse_args()

    campaign_dir = CAMPAIGNS_DIR / args.campaign / "ligands"
    if not campaign_dir.is_dir():
        sys.exit(f"ERROR: campaña no existe: {campaign_dir}")

    # Build list of molecules to test
    if args.molecule:
        # Buscar el archivo específico
        candidates = list(campaign_dir.glob(f"{args.molecule}.*"))
        candidates = [c for c in candidates if c.suffix.lower() in (".sdf", ".mol", ".mol2")]
        if not candidates:
            sys.exit(f"ERROR: molécula '{args.molecule}' no encontrada en {campaign_dir}")
        ligands = [candidates[0]]
    else:
        ligands = sorted(list(campaign_dir.glob("*.sdf")) + list(campaign_dir.glob("*.mol2")))

    print("=" * 78)
    print("PRUEBA REAL DEL FIX END-TO-END")
    print("=" * 78)
    print(f"Campaign:      {args.campaign}")
    print(f"Molecules:     {len(ligands)}")
    print(f"pH:            {args.ph}")
    print(f"Timeout T1:    {args.timeout}s")
    print()

    with tempfile.TemporaryDirectory(prefix="test_fix_") as td:
        td = Path(td)
        all_results = []

        for i, lig in enumerate(ligands, 1):
            print(f"[{i}/{len(ligands)}] {lig.stem}")
            print("-" * 78)

            mol_dir = td / lig.stem
            mol_dir.mkdir()

            t0 = time.time()
            res = prepare_single_molecule_with_fix(
                str(lig), mol_dir, docking_ph=args.ph, timeout=args.timeout
            )
            total = time.time() - t0
            res["total_elapsed_s"] = round(total, 1)

            print(f"  Carga formal pybel:    {res.get('formal_charge', '?')}")
            print(f"  Stages:                {res['stages']}")
            print(f"  Total time:            {total:.1f}s")
            if res.get("tier1_elapsed_s") is not None:
                print(f"  Tier 1 antechamber:    {res['tier1_elapsed_s']}s")

            if res.get("success"):
                print(f"  RESULT: ✓ PASS — BCC sum = {res['bcc_charge_sum']}")
                # Compare with original
                original = find_existing_mol2(args.campaign, lig.stem)
                if original:
                    cmp = compare_mol2(Path(res["output_mol2"]), original)
                    print()
                    print(f"  COMPARACIÓN vs producción actual:")
                    print(f"    header (fix vs orig):    {cmp['diff_charge_method']}")
                    print(f"    sum charges (fix/orig):  {cmp['diff_sum_charges']}")
                    print(f"    unique charges (f/o):    {cmp['diff_n_unique']}")
                    print(f"    same n_atoms:            {cmp['same_n_atoms']}")
                    res["comparison"] = cmp
                else:
                    print(f"  (Sin mol2 existente para comparar)")
            else:
                print(f"  RESULT: ✗ FAIL en stage '{res.get('fail_at')}'")
                if res.get("tier1_error"):
                    print(f"  Error: {res['tier1_error'][:200]}")

            all_results.append(res)
            print()

        # SUMMARY
        print("=" * 78)
        print("RESUMEN")
        print("=" * 78)
        n_pass = sum(1 for r in all_results if r.get("success"))
        print(f"PASS: {n_pass}/{len(all_results)}")
        if n_pass < len(all_results):
            print()
            print("Fallos:")
            for r in all_results:
                if not r.get("success"):
                    print(f"  {r['ligand']}: failed at {r.get('fail_at')}")

        print()
        print("Tabla compacta:")
        print(f"{'Ligand':<40s} {'Q':>4s} {'T1(s)':>8s} {'BCC sum':>10s} {'Status':<10s}")
        print("-" * 78)
        for r in all_results:
            q = r.get("formal_charge", "?")
            t1 = r.get("tier1_elapsed_s", "?")
            bcc = r.get("bcc_charge_sum", "?")
            status = "PASS" if r.get("success") else f"FAIL@{r.get('fail_at', '?')}"
            print(f"{r['ligand']:<40s} {str(q):>4s} {str(t1):>8s} {str(bcc):>10s} {status:<10s}")


if __name__ == "__main__":
    main()
