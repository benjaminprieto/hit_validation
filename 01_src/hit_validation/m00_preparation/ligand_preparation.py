"""
Ligand Preparation - Core Module (00a)
========================================
Prepares screening hit ligands for DOCK6 validation docking.

Unlike reference_docking (which preserves crystal coordinates),
hit_validation uses antechamber for AM1-BCC charges + Sybyl atom types.
Input molecules are screening hits from mol2/SDF files.

3-tier strategy with explicit warnings on degradation:
  Tier 1: clean-slate openbabel protonation (pybel API)
        + antechamber -c bcc -at sybyl -nc <formal_charge>
        Default timeout: 1800s (handles ~80-atom molecules).
  Tier 2: antechamber -c gas -at sybyl  (Gasteiger fallback)
        DEGRADED CHARGES -- emits warning, flags degraded_charges=True
  Tier 3: obabel conversion             (last-resort)
        DEGRADED CHARGES + atom typing -- emits warning

Pipeline per molecule:
    1. If SDF -> convert to mol2 with OpenBabel
    2. Clean-slate protonation (DeleteHydrogens + AddHydrogens(ph))
       Returns formal_charge from GetTotalCharge() -- REQUIRED for Tier 1.
    3. Run antechamber Tier 1 with -nc <formal_charge>
    4. If Tier 1 fails: warn, try Tier 2 (Gasteiger). degraded_charges=True.
    5. If Tier 2 fails: warn, try Tier 3 (OpenBabel). degraded_charges=True.
    6. Validate output mol2 (atom count, charges, Sybyl types)
    7. Log final tier + degraded_charges flag

Critical constraints:
  - atom_type MUST be sybyl. GAFF2 causes DOCK6 to silently fall back
    to rigid docking with zero molecular descriptors.
  - Molecules with degraded_charges=True propagate to docking but DOCK6
    ranking may be biased. MD/MMPBSA in 01g regenerates BCC from scratch,
    so final dG is unaffected.
  - Do NOT reuse reference_docking's 00a -- it preserves crystal coords.

Input:
  - campaign_dir/ligands/ (mol2 or SDF files)

Output:
  - 05_results/{campaign}/00a_ligand_preparation/
    - {name}/{name}.mol2          (DOCK6-ready, AM1-BCC charges, Sybyl types)
    - preparation_status.csv      (with degraded_charges column)
    - preparation_summary.txt     (with degraded molecules block if any)

Location: 01_src/hit_validation/m00_preparation/ligand_preparation.py
Project: hit_validation
Module: 00a (core)
Version: 1.1 (clean-slate + explicit -nc + degraded_charges flag)
"""

import logging
import shutil
import subprocess
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any, Union

import pandas as pd

logger = logging.getLogger(__name__)


# =============================================================================
# SDF -> MOL2 CONVERSION
# =============================================================================

def convert_sdf_to_mol2(sdf_path: str, output_mol2: str) -> bool:
    """Convert SDF to mol2 using OpenBabel."""
    try:
        result = subprocess.run(
            ["obabel", "-isdf", sdf_path, "-omol2", "-O", output_mol2],
            capture_output=True, text=True, timeout=60,
        )
        return result.returncode == 0 and Path(output_mol2).exists()
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


# =============================================================================
# PH-AWARE PROTONATION
# =============================================================================

def protonate_at_ph(input_mol2: str, output_mol2: str, ph: float = 6.3,
                    tool: str = "obabel") -> Optional[int]:
    """
    Protonate molecule at target pH using clean-slate openbabel via pybel API.

    Returns the total formal charge of the protonated molecule (to pass to
    antechamber as -nc). Returns None on failure.

    The clean-slate approach (DeleteHydrogens + AddHydrogens(ph)) gives correct
    protonation states for borderline pKa cases, where obabel CLI -p produces
    inconsistent monoanionic states for phosphate monoester at pH 6.3.

    Args:
        input_mol2: Input mol2 path
        output_mol2: Output mol2 path
        ph: Target pH (default 6.3 = Golgi pH for XT1)
        tool: Kept for backward compatibility but only "obabel" is supported.
              "ionization" tool is deprecated (use ionprofile package separately
              for auditing, not runtime).

    Returns:
        int: total formal charge (e.g., -2 for phosphate dianion at pH 6.3)
        None: if protonation failed
    """
    if tool not in ("obabel",):
        logger.warning(f"  protonation_tool='{tool}' not supported, "
                       f"using obabel clean-slate (only supported value)")

    try:
        from openbabel import pybel
        mol = next(pybel.readfile("mol2", input_mol2))
        # Clean-slate: remove existing H, then add fresh at target pH
        mol.OBMol.DeleteHydrogens()
        mol.OBMol.AddHydrogens(False, True, ph)
        formal_charge = mol.OBMol.GetTotalCharge()
        mol.write("mol2", output_mol2, overwrite=True)
        if Path(output_mol2).exists() and Path(output_mol2).stat().st_size > 0:
            return formal_charge
        return None
    except ImportError:
        logger.error("  openbabel/pybel not available -- cannot protonate")
        return None
    except Exception as e:
        logger.warning(f"  Protonation failed: {e}")
        return None


# =============================================================================
# ANTECHAMBER PREPARATION (3-TIER FALLBACK)
# =============================================================================

def run_antechamber_tier1(input_mol2: str, output_mol2: str,
                          formal_charge: int,
                          timeout: int = 1800) -> Dict[str, Any]:
    """
    Tier 1: Full AM1-BCC charges + Sybyl atom types with explicit net charge.

    The -nc flag is critical: without it, antechamber infers charge from the
    sum of partial charges in the input mol2, which is unreliable for mol2
    files produced by pybel (sums to non-integer values like -0.99 instead
    of -2). Passing -nc explicitly avoids odd-electron rejection by sqm.

    Args:
        input_mol2: Input mol2 (clean-slate protonated)
        output_mol2: Output mol2 path
        formal_charge: Net formal charge from pybel.GetTotalCharge() (REQUIRED)
        timeout: Max seconds for AM1 SCF convergence (default 1800s = 30min,
                 enough for molecules up to ~80 heavy atoms)
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
        "-c", "bcc",       # AM1-BCC charges
        "-at", "sybyl",    # Sybyl atom types (MUST for DOCK6)
        "-pf", "y",        # Remove intermediate files
        "-dr", "no",       # No default residue check
        "-nc", str(formal_charge),  # CRITICAL: explicit net charge
    ]
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=timeout,
            cwd=out_dir,
        )
        success = (result.returncode == 0
                    and Path(output_mol2).exists()
                    and Path(output_mol2).stat().st_size > 0)
        return {
            "success": success,
            "tier": 1,
            "method": "AM1-BCC + Sybyl",
            "formal_charge": formal_charge,
            "error": result.stderr[:300] if not success else None,
        }
    except subprocess.TimeoutExpired:
        return {"success": False, "tier": 1, "method": "AM1-BCC + Sybyl",
                "formal_charge": formal_charge,
                "error": f"Timeout after {timeout}s"}
    except FileNotFoundError:
        return {"success": False, "tier": 1, "method": "AM1-BCC + Sybyl",
                "formal_charge": formal_charge,
                "error": "antechamber not found in PATH"}


def run_antechamber_tier2(input_mol2: str, output_mol2: str,
                          timeout: int = 120) -> Dict[str, Any]:
    """
    Tier 2: Gasteiger charges + Sybyl atom types.

    Fast fallback -- no QM calculation. Solves timeouts on large molecules.
    Charge quality is lower than AM1-BCC but sufficient for docking.
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
        "-c", "gas",        # Gasteiger charges (fast, no QM)
        "-at", "sybyl",     # Sybyl atom types (MUST for DOCK6)
        "-pf", "y",
        "-dr", "no",
    ]
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=timeout,
            cwd=out_dir,
        )
        success = (result.returncode == 0
                    and Path(output_mol2).exists()
                    and Path(output_mol2).stat().st_size > 0)
        return {
            "success": success,
            "tier": 2,
            "method": "Gasteiger + Sybyl",
            "error": result.stderr[:300] if not success else None,
        }
    except subprocess.TimeoutExpired:
        return {"success": False, "tier": 2, "method": "Gasteiger + Sybyl",
                "error": f"Timeout after {timeout}s"}
    except FileNotFoundError:
        return {"success": False, "tier": 2, "method": "Gasteiger + Sybyl",
                "error": "antechamber not found in PATH"}


def run_obabel_tier3(input_mol2: str, output_mol2: str) -> Dict[str, Any]:
    """
    Tier 3: OpenBabel conversion (last resort).

    Preserves connectivity but atom typing may be degraded.
    Generates Gasteiger charges via obabel --partialcharge gasteiger.
    """
    try:
        result = subprocess.run(
            ["obabel", input_mol2, "-O", output_mol2,
             "--partialcharge", "gasteiger"],
            capture_output=True, text=True, timeout=60,
        )
        success = (result.returncode == 0
                    and Path(output_mol2).exists()
                    and Path(output_mol2).stat().st_size > 0)
        return {
            "success": success,
            "tier": 3,
            "method": "OpenBabel (degraded)",
            "error": result.stderr[:300] if not success else None,
        }
    except (FileNotFoundError, subprocess.TimeoutExpired) as e:
        return {"success": False, "tier": 3, "method": "OpenBabel (degraded)",
                "error": str(e)}


def prepare_single_molecule(
        input_path: str,
        output_mol2: str,
        atom_type: str = "sybyl",
        charge_method: str = "bcc",
        docking_ph: float = 6.3,
        timeout: int = 1800,
        protonation_tool: str = "obabel",
) -> Dict[str, Any]:
    """
    Prepare a single molecule with 3-tier antechamber fallback.

    Tier 1 uses clean-slate openbabel protonation (pybel API) and explicit
    -nc <charge> flag to antechamber. This fixes a silent bug where obabel
    CLI -p produces monoanionic states for borderline-pKa groups, causing
    sqm to reject due to odd electron parity.

    If Tier 1 fails, fallback proceeds through Tier 2 (Gasteiger) and Tier 3
    (obabel-only). These fallbacks are NOT silent: they emit prominent
    warnings and the result dict carries degraded_charges=True so downstream
    analysis can flag affected molecules.

    Note: MD/MMPBSA in 01g re-parametrizes with AM1-BCC from scratch, so
    final dG is unaffected by Tier 2/3 fallback. However, DOCK6 ranking
    (Grid_Score, footprint per-residue) IS affected because it uses the
    mol2 produced here directly. Molecules with degraded_charges=True
    should be reviewed before any purchase decision based on DOCK6 ranking.

    Args:
        input_path: Path to input mol2 or SDF
        output_mol2: Path for output DOCK6-ready mol2
        atom_type: Must be "sybyl" for DOCK6
        charge_method: "bcc" (AM1-BCC) -- Tier 1; fallback methods are automatic.
        docking_ph: pH for protonation (default 6.3 = Golgi pH)
        timeout: Per-molecule Tier 1 timeout (default 1800s = 30min)
        protonation_tool: Only "obabel" supported (clean-slate via pybel API).

    Returns:
        Dict with: success, tier (1/2/3 or 0 if all failed), method,
        formal_charge, degraded_charges (bool), warnings, error.
    """
    inp = Path(input_path).resolve()
    out = Path(output_mol2).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)

    warnings = []

    if atom_type != "sybyl":
        warnings.append(f"atom_type={atom_type} overridden to sybyl (required for DOCK6)")
        atom_type = "sybyl"

    # Step 1: Convert SDF to mol2 if needed
    work_mol2 = str(out.parent / f"{inp.stem}_input.mol2")
    if inp.suffix.lower() in (".sdf", ".mol"):
        if not convert_sdf_to_mol2(str(inp), work_mol2):
            return {"success": False, "tier": 0, "method": "sdf_conversion",
                    "error": f"Failed to convert SDF to mol2: {inp.name}",
                    "degraded_charges": False, "warnings": warnings}
    else:
        shutil.copy2(str(inp), work_mol2)

    # Step 2: Clean-slate protonation (returns formal_charge or None)
    protonated_mol2 = str(out.parent / f"{inp.stem}_protonated.mol2")
    formal_charge = protonate_at_ph(work_mol2, protonated_mol2, docking_ph,
                                     tool=protonation_tool)
    if formal_charge is None:
        return {"success": False, "tier": 0, "method": "protonation",
                "error": "Clean-slate protonation failed (pybel)",
                "degraded_charges": False, "warnings": warnings}

    warnings.append(f"Clean-slate protonation: charge={formal_charge:+d} at pH {docking_ph}")

    # Step 3a: Tier 1 (AM1-BCC + Sybyl) with -nc explicit
    result = run_antechamber_tier1(protonated_mol2, str(out), formal_charge, timeout)
    if result["success"]:
        result["degraded_charges"] = False
        _cleanup_intermediates(out.parent, inp.stem)
        result["warnings"] = warnings
        return result

    # Tier 1 failed -- emit prominent warning and try Tier 2
    logger.warning(
        f"    [WARN] Tier 1 (AM1-BCC) FAILED for {inp.stem} "
        f"(charge={formal_charge:+d}): {result.get('error', '')[:120]}"
    )
    logger.warning(
        f"    [WARN] Falling back to Tier 2 (Gasteiger). "
        f"DOCK6 ranking for this molecule will use degraded charges."
    )
    warnings.append(f"Tier 1 (AM1-BCC) failed: {result.get('error', '')[:150]}")

    # Step 3b: Tier 2 (Gasteiger + Sybyl)
    result = run_antechamber_tier2(protonated_mol2, str(out), min(timeout, 120))
    if result["success"]:
        warnings.append("DEGRADED CHARGES: Tier 2 (Gasteiger) used because Tier 1 failed")
        result["degraded_charges"] = True
        _cleanup_intermediates(out.parent, inp.stem)
        result["warnings"] = warnings
        return result

    logger.warning(
        f"    [WARN] Tier 2 (Gasteiger) also FAILED for {inp.stem}: "
        f"{result.get('error', '')[:120]}"
    )
    warnings.append(f"Tier 2 (Gasteiger) failed: {result.get('error', '')[:150]}")

    # Step 3c: Tier 3 (OpenBabel direct)
    logger.warning(
        f"    [WARN] Falling back to Tier 3 (OpenBabel) for {inp.stem}. "
        f"Atom typing will be degraded."
    )
    result = run_obabel_tier3(protonated_mol2, str(out))
    if result["success"]:
        warnings.append("DEGRADED CHARGES + TYPING: Tier 3 (OpenBabel) used")
        result["degraded_charges"] = True
    else:
        result["degraded_charges"] = False  # No mol2 produced, so flag is moot
    _cleanup_intermediates(out.parent, inp.stem)
    result["warnings"] = warnings
    return result


def _cleanup_intermediates(parent: Path, stem: str):
    """Remove intermediate files from antechamber."""
    for pattern in [f"{stem}_input.mol2", f"{stem}_protonated.mol2",
                    "ANTECHAMBER_*", "ATOMTYPE.INF", "sqm.*", "leap.log"]:
        for f in parent.glob(pattern):
            try:
                f.unlink()
            except OSError:
                pass


# =============================================================================
# VALIDATION
# =============================================================================

def validate_mol2(mol2_path: str) -> Dict[str, Any]:
    """
    Validate a prepared mol2 file for DOCK6 compatibility.

    Checks:
      - File exists and is non-empty
      - Has ATOM section with atoms
      - Atoms have charges (column 9)
      - Atom types look like Sybyl (contain '.')
    """
    path = Path(mol2_path)
    if not path.exists() or path.stat().st_size == 0:
        return {"valid": False, "error": "File empty or missing"}

    n_atoms = 0
    n_with_charge = 0
    n_sybyl_types = 0
    in_atom = False

    with open(mol2_path) as f:
        for line in f:
            if "@<TRIPOS>ATOM" in line:
                in_atom = True
                continue
            if line.startswith("@<TRIPOS>") and in_atom:
                break
            if in_atom and line.strip():
                parts = line.split()
                if len(parts) >= 6:
                    n_atoms += 1
                    atom_type = parts[5] if len(parts) > 5 else ""
                    if "." in atom_type:
                        n_sybyl_types += 1
                    if len(parts) >= 9:
                        try:
                            charge = float(parts[8])
                            if charge != 0.0:
                                n_with_charge += 1
                        except ValueError:
                            pass

    if n_atoms == 0:
        return {"valid": False, "error": "No atoms in ATOM section"}

    return {
        "valid": True,
        "n_atoms": n_atoms,
        "n_with_charge": n_with_charge,
        "n_sybyl_types": n_sybyl_types,
        "frac_charged": round(n_with_charge / n_atoms, 3),
        "frac_sybyl": round(n_sybyl_types / n_atoms, 3),
    }


# =============================================================================
# MAIN PIPELINE
# =============================================================================

def run_ligand_preparation(
        ligand_dir: Union[str, Path],
        output_dir: Union[str, Path],
        atom_type: str = "sybyl",
        charge_method: str = "bcc",
        docking_ph: float = 6.3,
        timeout_per_molecule: int = 1800,
        molecule_filter: Optional[List[str]] = None,
        protonation_tool: str = "obabel",
) -> Dict[str, Any]:
    """
    Prepare N screening hit ligands for DOCK6 validation docking.

    Scans ligand_dir for *.mol2 and *.sdf files. For each molecule,
    runs antechamber with AM1-BCC charges + Sybyl atom types using
    a 3-tier fallback strategy.

    Args:
        ligand_dir: Directory with input mol2/SDF files
        output_dir: Output directory for prepared mol2 files
        atom_type: Atom type system (must be "sybyl" for DOCK6)
        charge_method: Charge method ("bcc" or "gas")
        docking_ph: pH for protonation
        timeout_per_molecule: Timeout per molecule in seconds
        molecule_filter: Optional list of molecule names to process

    Returns:
        Dict with: success, n_total, n_ok, n_failed, results
    """
    ligand_dir = Path(ligand_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not ligand_dir.exists():
        return {"success": False, "error": f"Ligand dir not found: {ligand_dir}"}

    # Find input files (mol2 + SDF)
    input_files = sorted(
        list(ligand_dir.glob("*.mol2")) + list(ligand_dir.glob("*.sdf"))
    )

    if molecule_filter:
        filter_set = set(molecule_filter)
        input_files = [f for f in input_files if f.stem in filter_set]

    if not input_files:
        return {"success": False, "error": f"No mol2/SDF files in {ligand_dir}"}

    logger.info("=" * 60)
    logger.info("  Ligand Preparation (Hit Validation)")
    logger.info("=" * 60)
    logger.info(f"  Input dir:   {ligand_dir}")
    logger.info(f"  Molecules:   {len(input_files)}")
    logger.info(f"  Atom type:   {atom_type}")
    logger.info(f"  Charges:     {charge_method}")
    logger.info(f"  pH:          {docking_ph}")
    logger.info(f"  Protonation: {protonation_tool}")
    logger.info(f"  Timeout:     {timeout_per_molecule}s per molecule")

    results = []
    total_time = 0

    for idx, input_file in enumerate(input_files, 1):
        name = input_file.stem
        mol_out_dir = output_dir / name
        mol_out_dir.mkdir(parents=True, exist_ok=True)
        output_mol2 = str(mol_out_dir / f"{name}.mol2")

        logger.info(f"  [{idx}/{len(input_files)}] {name} ({input_file.suffix})")

        t0 = time.time()
        prep_result = prepare_single_molecule(
            input_path=str(input_file),
            output_mol2=output_mol2,
            atom_type=atom_type,
            charge_method=charge_method,
            docking_ph=docking_ph,
            timeout=timeout_per_molecule,
            protonation_tool=protonation_tool,
        )
        runtime = time.time() - t0
        total_time += runtime

        # Validate output
        validation = {}
        if prep_result["success"]:
            validation = validate_mol2(output_mol2)
            if not validation.get("valid"):
                prep_result["success"] = False
                prep_result["error"] = f"Validation failed: {validation.get('error')}"

        status = "OK" if prep_result["success"] else "FAILED"
        results.append({
            "Name": name,
            "status": status,
            "tier": prep_result.get("tier", 0),
            "method": prep_result.get("method", "unknown"),
            "output_mol2": output_mol2 if prep_result["success"] else None,
            "runtime_sec": round(runtime, 1),
            "n_atoms": validation.get("n_atoms", 0),
            "frac_charged": validation.get("frac_charged", 0),
            "frac_sybyl": validation.get("frac_sybyl", 0),
            "degraded_charges": prep_result.get("degraded_charges", False),
            "warnings": "; ".join(prep_result.get("warnings", [])),
            "error": prep_result.get("error"),
        })

        if status == "OK":
            tier_str = f"tier {prep_result['tier']}"
            degraded_str = " [DEGRADED CHARGES]" if prep_result.get("degraded_charges") else ""
            logger.info(f"    -> OK ({runtime:.1f}s, {tier_str}, "
                         f"{validation.get('n_atoms', 0)} atoms){degraded_str}")
        else:
            logger.warning(f"    -> FAILED: {prep_result.get('error', 'unknown')}")

    # --- Save status CSV ---
    df_status = pd.DataFrame(results)
    status_csv = output_dir / "preparation_status.csv"
    df_status.to_csv(status_csv, index=False, encoding="utf-8")

    # --- Save summary TXT ---
    n_ok = sum(1 for r in results if r["status"] == "OK")
    n_fail = len(results) - n_ok
    tier_counts = {}
    for r in results:
        if r["status"] == "OK":
            t = r["tier"]
            tier_counts[t] = tier_counts.get(t, 0) + 1

    degraded = [r for r in results if r["status"] == "OK" and r.get("degraded_charges")]
    n_degraded = len(degraded)

    summary_path = output_dir / "preparation_summary.txt"
    w = 70
    lines = [
        "=" * w,
        "00a LIGAND PREPARATION - SUMMARY (Hit Validation)",
        "=" * w,
        "",
        f"Date:              {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        f"Input dir:         {ligand_dir}",
        f"Molecules:         {len(results)}",
        f"Successful:        {n_ok}",
        f"Failed:            {n_fail}",
        f"Total time:        {total_time:.0f}s ({total_time / 60:.1f} min)",
        "",
        "Tier breakdown:",
    ]
    tier_labels = {
        1: "AM1-BCC + Sybyl",
        2: "Gasteiger + Sybyl   [WARN] degraded charges",
        3: "OpenBabel           [WARN][WARN] degraded charges + atom typing",
    }
    for t in sorted(tier_counts.keys()):
        lines.append(f"  Tier {t} ({tier_labels.get(t, 'unknown')}): {tier_counts[t]}")

    lines.extend([
        "",
        "-" * w,
        f"{'Name':<30} {'Status':>8} {'Tier':>5} {'Atoms':>6} {'Time(s)':>8} {'Deg':>4}",
        "-" * w,
    ])

    for r in results:
        deg = "Y" if r.get("degraded_charges") else "-"
        lines.append(
            f"{r['Name']:<30} {r['status']:>8} {r['tier']:>5} "
            f"{r['n_atoms']:>6} {r['runtime_sec']:>8.1f} {deg:>4}"
        )

    if n_fail > 0:
        lines.extend(["", "FAILURES:"])
        for r in results:
            if r["status"] == "FAILED":
                lines.append(f"  {r['Name']}: {r.get('error', 'unknown')}")

    if n_degraded > 0:
        lines.extend([
            "",
            "=" * w,
            f"[WARN] {n_degraded} molecule(s) use DEGRADED partial charges (Tier 2 or 3).",
            "=" * w,
            "These molecules entered the pipeline but their DOCK6 ranking",
            "(Grid_Score, footprint per-residue) may be biased toward charged species.",
            "MD/MMPBSA in 01g re-parametrizes with AM1-BCC from scratch, so final",
            "dG values are NOT contaminated. However, DOCK6-based decisions",
            "(top-N selection, purchase ranking) should be reviewed manually for",
            "these molecules.",
            "",
            "Affected molecules (see preparation_status.csv, column 'degraded_charges'):",
        ])
        for r in degraded:
            lines.append(f"  - {r['Name']} (Tier {r['tier']})")

    lines.extend(["", "=" * w])

    summary_path.write_text("\n".join(lines), encoding="utf-8")

    logger.info("")
    logger.info(f"{'=' * 60}")
    logger.info(f"  LIGAND PREP: {n_ok}/{len(results)} prepared "
                 f"({total_time:.0f}s total)")
    if n_fail > 0:
        logger.info(f"  FAILED: {n_fail}")
    if n_degraded > 0:
        logger.warning(
            f"  [WARN] {n_degraded} molecule(s) with DEGRADED charges (Tier 2/3). "
            f"Review DOCK6 ranking before purchase decisions. "
            f"See preparation_summary.txt for details."
        )
    logger.info(f"{'=' * 60}")

    return {
        "success": n_ok > 0,
        "n_total": len(results),
        "n_ok": n_ok,
        "n_failed": n_fail,
        "n_degraded": n_degraded,
        "total_runtime_sec": round(total_time, 1),
        "tier_counts": tier_counts,
        "results": results,
        "status_csv": str(status_csv),
        "summary_txt": str(summary_path),
        "output_dir": str(output_dir),
    }