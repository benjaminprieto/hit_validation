# INSTRUCTIVO: Fix Módulo 00a + Bugs Identificados en hit_validation

Este instructivo contiene tres cambios al pipeline de `hit_validation`, todos
validados con datos. Cada cambio está descrito con contexto, ubicación exacta,
código nuevo, y tests obligatorios antes de commit.

**NO mezclar con cambios a `molecular_docking`** — ese fix se aplicará en un
instructivo y PR separados, una vez validado este.

---

## CAMBIO 1: Fix del bug Tier 1 silencioso en `00a` (CRÍTICO)

### Contexto

El módulo 00a tiene un bug en la cadena protonación → antechamber Tier 1 que
causa que ~35% de las moléculas con grupos ionizables borderline (fosfato,
fosfonato, etc. cerca de pH 6.3) fallen Tier 1 y caigan silenciosamente a Tier 2
(Gasteiger). Esto contamina el ranking DOCK6/footprint con cargas degradadas en
la dimensión electrostática (~88% del Grid_Score para este pocket).

**Bug raíz validado en este orden:**

1. `obabel input.mol2 -O output.mol2 -p 6.3` (CLI) NO hace clean-slate. Para
   fosfato monoester a pH 6.3 deja la molécula como mono-anión (-1), generando
   paridad electrónica impar.
2. Pybel sí hace clean-slate (`DeleteHydrogens()` + `AddHydrogens(False, True, ph)`)
   y reporta carga total correcta (-2 para dianión).
3. Pero el writer de mol2 de pybel NO escribe la suma de cargas parciales con
   precisión integer; suma -0.99 cuando la carga real es -2.
4. Antechamber sin `-nc` infiere carga del mol2 redondeando la suma → calcula
   electrones impares → sqm rechaza pre-SCF en <1s.

**Validación con datos (no especulación):**

- Audit completo de 162 moléculas con timeout 1200s → 84% PASS con fix.
- Test end-to-end con timeout 1800s → 138/140 moléculas no-Namiki PASS (98.6%).
- Sin regresiones (las "regresiones" iniciales fueron falsos timeouts del test).
- 2/2 moléculas en SD1_druglikeness_2_PChem_136488320_2n (que cayeron a Tier 2
  silencioso) validadas individualmente: pasan Tier 1 BCC con suma -2.0000 exacta.

### Archivos a modificar

```
01_src/hit_validation/m00_preparation/ligand_preparation.py
03_configs/00a_ligand_preparation.yaml
02_scripts/00a_ligand_preparation.py    (revisar reporte/logging)
```

### 1.1. Modificar `protonate_at_ph()` para clean-slate via pybel API

**Ubicación:** `01_src/hit_validation/m00_preparation/ligand_preparation.py`,
función `protonate_at_ph()` (líneas ~74-103).

**Cambio de signatura:** la función ahora devuelve `Optional[int]` (carga formal
total del molécula protonada) en vez de `bool`. Devolver None indica fallo.

**NO usar el paquete `ionization` (`ionprofile`).** Es para auditoría batch
separada, no es dependencia runtime.

```python
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
        logger.error("  openbabel/pybel not available — cannot protonate")
        return None
    except Exception as e:
        logger.warning(f"  Protonation failed: {e}")
        return None
```

### 1.2. Modificar `run_antechamber_tier1()` para aceptar `formal_charge` y pasar `-nc`

**Ubicación:** misma función, líneas ~115-157.

**Cambio:** agregar parámetro `formal_charge: int` obligatorio, pasar `-nc <charge>`
explícitamente a antechamber. Subir el timeout default a 1800s (suficiente para
moléculas hasta ~80 átomos pesados; molécula PubChem-70087525 con 76 átomos y
+2 cargas requirió 45min de SCF con AM1).

```python
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
```

### 1.3. Mantener Tier 2 y Tier 3 con warnings explícitos + flag `degraded_charges`

**Decisión:** las funciones `run_antechamber_tier2()` y `run_obabel_tier3()`
SE MANTIENEN. La razón es que el módulo 01g re-parametriza con AM1-BCC desde
cero (validado: produce topología con suma de cargas -2.001 y 26 cargas únicas
distintas), por lo que la MD y MMPBSA finales NO se contaminan con Gasteiger.
Lo que sí se contamina es el ranking DOCK6 pre-MD (Grid_Score, footprint).

**Cambios necesarios:**

1. Agregar warning fuerte y visible cuando una molécula cae a Tier 2 o Tier 3.
2. Agregar columna `degraded_charges: bool` a `preparation_status.csv`.
3. Imprimir resumen al final del módulo indicando cuántas moléculas tienen
   cargas degradadas y advirtiendo del impacto en ranking DOCK6.

### 1.4. Modificar `prepare_single_molecule()` para integrar el fix + warnings

**Ubicación:** `prepare_single_molecule()`, líneas ~239-309.

```python
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
    final ΔG is unaffected by Tier 2/3 fallback. However, DOCK6 ranking
    (Grid_Score, footprint per-residue) IS affected because it uses the
    mol2 produced here directly. Molecules with degraded_charges=True
    should be reviewed before any purchase decision based on DOCK6 ranking.

    Args:
        input_path: Path to input mol2 or SDF
        output_mol2: Path for output DOCK6-ready mol2
        atom_type: Must be "sybyl" for DOCK6
        charge_method: "bcc" (AM1-BCC) — Tier 1; fallback methods are automatic.
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

    # Tier 1 failed — emit prominent warning and try Tier 2
    logger.warning(
        f"    ⚠️  Tier 1 (AM1-BCC) FAILED for {inp.stem} "
        f"(charge={formal_charge:+d}): {result.get('error', '')[:120]}"
    )
    logger.warning(
        f"    ⚠️  Falling back to Tier 2 (Gasteiger). "
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
        f"    ⚠️  Tier 2 (Gasteiger) also FAILED for {inp.stem}: "
        f"{result.get('error', '')[:120]}"
    )
    warnings.append(f"Tier 2 (Gasteiger) failed: {result.get('error', '')[:150]}")

    # Step 3c: Tier 3 (OpenBabel direct)
    logger.warning(
        f"    ⚠️  Falling back to Tier 3 (OpenBabel) for {inp.stem}. "
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
```

### 1.5. Actualizar el reporte (script orquestador)

**Ubicación:** `02_scripts/00a_ligand_preparation.py` (o wherever the orchestrator
lives that calls `prepare_single_molecule` and writes `preparation_status.csv` /
`preparation_summary.txt`).

**Cambios al CSV:** agregar columna `degraded_charges` (bool). El header ahora
debe ser:

```
Name,status,tier,method,output_mol2,runtime_sec,n_atoms,frac_charged,frac_sybyl,degraded_charges,warnings,error
```

**Cambios al summary.txt:** agregar sección clara al final indicando cuántas
moléculas tienen `degraded_charges=True` y el impacto:

```
======================================================================
00a LIGAND PREPARATION - SUMMARY
======================================================================

[... output existing ...]

Tier breakdown:
  Tier 1 (AM1-BCC + Sybyl):   X
  Tier 2 (Gasteiger + Sybyl): Y   ⚠ degraded charges
  Tier 3 (OpenBabel):         Z   ⚠⚠ degraded charges + atom typing
  FAILED:                     W

[... if Y+Z > 0, emit this prominent block: ...]

⚠️  WARNING: (Y+Z) molecules use degraded partial charges (Tier 2 or 3).
   These molecules entered the pipeline but their DOCK6 ranking
   (Grid_Score, footprint per-residue) may be biased toward charged species.
   MD/MMPBSA in 01g re-parametrizes with AM1-BCC from scratch, so final
   ΔG values are NOT contaminated. However, DOCK6-based decisions
   (top-N selection, purchase ranking) should be reviewed manually for
   these molecules.

   Affected molecules (see preparation_status.csv, column 'degraded_charges'):
   - mol_name_1 (Tier 2)
   - mol_name_2 (Tier 3)
   ...
```

### 1.6. Actualizar el config

**Ubicación:** `03_configs/00a_ligand_preparation.yaml`

```yaml
# Increased from 1200 to 1800: covers AM1-BCC SCF for molecules up to ~80 heavy
# atoms with multiple formal charges (e.g., dicationic 76-atom case took 45min).
timeout: 1800

# Only "obabel" is supported. The "ionization" option is deprecated and will
# emit a warning if set. The ionprofile package is for batch auditing, not runtime.
protonation_tool: "obabel"
```

### 1.7. Actualizar el header del módulo

Reemplazar el docstring de `ligand_preparation.py` (líneas 1-44) para reflejar
la nueva arquitectura:

```python
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
        ⚠ DEGRADED CHARGES — emits warning, flags degraded_charges=True
  Tier 3: obabel conversion             (last-resort)
        ⚠⚠ DEGRADED CHARGES + atom typing — emits warning

Pipeline per molecule:
    1. If SDF → convert to mol2 with OpenBabel
    2. Clean-slate protonation (DeleteHydrogens + AddHydrogens(ph))
       Returns formal_charge from GetTotalCharge() — REQUIRED for Tier 1.
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
    so final ΔG is unaffected.
  - Do NOT reuse reference_docking's 00a — it preserves crystal coords.

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
```

---

## CAMBIO 2: Aumentar frames MMPBSA en `01g` para análisis temporal robusto

### Contexto

Actualmente 01g produce 10 frames de MMPBSA sobre 50ns de MD (`stride 500`).
Esto es insuficiente para el análisis temporal del módulo 07e (Temporal MMPBSA
Analysis) que requiere:

- Detección de equilibración (pymbar/RCA): N≥30 puntos para confiabilidad.
- Block averaging para SEM: bloques de 5 con N=10 da ruido alto; con N=50
  los bloques de 10 dan SEM limpia.
- Detección de basin transitions: 25 frames es justo, 50 es el estándar.
- Cross-campaign KS tests: N≥30 para tests no-paramétricos confiables.

**Cambio:** aumentar a 50 frames por molécula. Coherente con `mmpbsa_interval_ps=1000`
ya configurado en el diseño de módulo 07e (memoria del proyecto: "50 frames/run").

### Archivos a modificar

```
03_configs/01g_mmpbsa_decomp.yaml
01_src/hit_validation/m01_validation/mmpbsa_decomp.py    (o equivalente)
```

### 2.1. Cambiar el stride de cpptraj

Buscar la línea donde cpptraj hace strip + stride. Actualmente probablemente:

```
strip :WAT,Cl-,Na+ outprefix dry
trajout dry_trajectory.mdcrd offset 500
```

Cambiar `offset 500` por `offset 100` (stride 100 sobre 5000 frames totales =
50 frames para MMPBSA).

### 2.2. Aumentar el timeout de MMPBSA.py

**Ubicación:** función que llama `subprocess.run(["MMPBSA.py", ...])` con
parámetro `timeout=`.

Cambiar `timeout=7200` por `timeout=28800` (8 horas). Justificación:

- 10 frames tomaron ~6400s en run anterior (Conformer3D_COMPOUND_CID_135966069).
- 50 frames extrapolando linealmente: ~32000s, pero con margen para receptores
  grandes y variabilidad: 28800s es conservador.
- Si una molécula no termina en 8h, hay otro problema (PB grid demasiado grande,
  receptor muy lento) que debe diagnosticarse, no acomodarse con más timeout.

### 2.3. Actualizar el config

```yaml
# 03_configs/01g_mmpbsa_decomp.yaml

# Frame stride for cpptraj strip — produces ~50 frames over 50ns MD
# (5000 total frames / stride 100 = 50 frames for MMPBSA).
# Required for temporal analysis in 07e: equilibration detection,
# block averaging, basin transitions, cross-campaign KS tests.
mmpbsa_stride: 100   # was 500

# MMPBSA.py subprocess timeout (seconds). Increased from 7200 to 28800
# (8 hours) to accommodate ~50 frames over the protein receptor (~11k atoms).
mmpbsa_timeout: 28800   # was 7200
```

### 2.4. Verificar coherencia con módulo 07e

El módulo 07e (`Temporal MMPBSA Analysis`, según memoria del proyecto) está
diseñado con `mmpbsa_interval_ps=1000`. Con `save_interval=10ps` (de la MD) y
stride 100 en cpptraj, esto da:

- Frame interval real: 10ps × 100 = 1000ps = 1ns entre frames de MMPBSA.
- 50ns / 1ns = 50 frames totales.

Ambos números coinciden. **NO modificar 07e** — solo verificar que su lectura
de los archivos `_MMPBSA_*.mdout.NN` siga funcionando con el N nuevo (50 en
vez de 10). Si 07e tiene un check explícito de "N == 10", removerlo.

---

## CAMBIO 3: Bug cosmético "0 atoms" en `01g`

### Contexto

El log de 01g muestra:

```
tleap: building gas-phase topologies
Gas-phase complex: 0 atoms
tleap: building solvated topology
Solvated complex: 0 atoms
```

Pero los `prmtop` reales tienen los átomos correctos (verificado con parmed:
complex.prmtop tiene 11345 atoms, ligando UNL con 34 atoms y suma de cargas
-2.001). El bug es del logging — probablemente una función que cuenta átomos
del prmtop antes de que tleap termine de escribir, o lee el archivo equivocado.

### Archivos a modificar

```
01_src/hit_validation/m01_validation/mmpbsa_decomp.py    (Step 3: Topology Construction)
```

### 3.1. Localizar y corregir el conteo

Buscar el bloque que loguea "Gas-phase complex" y "Solvated complex". Probablemente
hay un patrón tipo:

```python
n_atoms = some_function_that_returns_zero(prmtop_path)
logger.info(f"  Gas-phase complex: {n_atoms} atoms")
```

**Fix recomendado:** usar parmed (que ya está en deps) para contar átomos del
prmtop con confiabilidad:

```python
import parmed as pmd

def count_atoms_in_prmtop(prmtop_path: str) -> int:
    """Count atoms in a prmtop file using parmed."""
    try:
        p = pmd.load_file(str(prmtop_path))
        return len(p.atoms)
    except Exception as e:
        logger.warning(f"  Could not count atoms in {Path(prmtop_path).name}: {e}")
        return -1  # signal counting failure, NOT zero atoms
```

Y reemplazar las llamadas al conteo:

```python
# Antes:
logger.info(f"  Gas-phase complex: {n_atoms} atoms")

# Después:
n_atoms_gas = count_atoms_in_prmtop(complex_gas_prmtop)
if n_atoms_gas > 0:
    logger.info(f"  Gas-phase complex: {n_atoms_gas} atoms")
elif n_atoms_gas == 0:
    logger.error(f"  Gas-phase complex topology is EMPTY — tleap may have failed")
else:
    logger.warning(f"  Gas-phase complex: atom count unavailable")
```

Esto distingue tres casos: (1) éxito real, (2) tleap falló de verdad
(prmtop vacío, escenario crítico), (3) error de lectura (no crítico, log sigue).

### 3.2. Verificar que el problema no es timing

Si el bug viene de leer el prmtop antes de que tleap lo termine de escribir,
agregar una verificación explícita post-tleap:

```python
# Después de invocar tleap:
subprocess.run(["tleap", "-f", tleap_script], ...)

# Verificar que el archivo existe Y tiene tamaño no-trivial:
prmtop = Path(complex_prmtop_path)
if not prmtop.exists() or prmtop.stat().st_size < 1024:
    raise RuntimeError(f"tleap failed: {prmtop} empty or too small")
```

Esto es defensivo y captura el caso real (tleap falló) que el bug del logging
estaba enmascarando.

---

## TESTS DE VALIDACIÓN OBLIGATORIOS

Antes de hacer commit, correr estos tests sobre moléculas conocidas. NO hacer
commit si alguno falla.

### Test 1: No-regresión sobre moléculas que ya pasaban Tier 1

Ya validado en el audit completo, pero re-correr post-fix para confirmar:

```bash
# Esperado: tier=1, success=True, degraded_charges=False, suma BCC entera (±0.01)
python 02_scripts/00a_ligand_preparation.py \
    --config 03_configs/00a_ligand_preparation.yaml \
    --campaigns 04_data/campaigns/SD1_druglikeness_2_PChem_67908940/campaign_config.yaml
```

(Esa campaña tenía 7/7 PASS Tier 1 originalmente — todas deben seguir PASS.)

### Test 2: Recovery de moléculas que caían a Tier 2

```bash
python 02_scripts/00a_ligand_preparation.py \
    --config 03_configs/00a_ligand_preparation.yaml \
    --campaigns 04_data/campaigns/SD1_druglikeness_2_PChem_136488320_2n/campaign_config.yaml
```

Esperado: 2/2 moléculas con `tier=1`, `formal_charge=-2`, `degraded_charges=False`.

### Test 3: Failure explícito — Tier 2 fallback con warning

Tomar UNA molécula Namiki conocida (HTS1710-*.mol2 de campaña druglikeness_1)
y confirmar que:

- Tier 1 falla con warning visible.
- Cae a Tier 2 con warning visible adicional.
- Result final: `tier=2`, `degraded_charges=True`.
- `preparation_status.csv` muestra columna `degraded_charges=True`.

### Test 4: MMPBSA con 50 frames

Tomar la molécula que ya completó 01g exitosamente
(`SD1_druglikeness_2_PChem_136488320_2n/replica_1/`,
`Conformer3D_COMPOUND_CID_135966069`) y re-correr 01g para verificar:

- Cpptraj produce trayectoria con 50 frames (no 10).
- MMPBSA.py corre sin timeout (debería completar en ~5-7h).
- Output `_MMPBSA_*.mdout.50` (no `.10`) para indicar 50 frames procesados.

### Test 5: Logging "0 atoms" arreglado

Re-correr 01g y verificar que el log muestre:

```
Gas-phase complex: 11345 atoms       ← ya no 0
Solvated complex: 78752 atoms        ← ya no 0
```

---

## FOLLOW-UP ISSUES IDENTIFICADOS (NO incluir en este PR)

Estos bugs fueron identificados durante el debugging pero son separados del
fix actual. Crear issues separados, NO mezclarlos en este PR:

1. **`molecular_docking` tiene el mismo bug en su 00a.** Aplicar el mismo fix
   una vez validado este. PR separado en repo `molecular_docking`.

2. **cpptraj timeout en 01i (water bridges)**: 3600s no alcanza para procesar
   trayectoria solvatada con 50ns. Considerar reducir frecuencia de análisis
   o tamaño de trayectoria. Issue separado en módulo 01i.

3. **Validar PubChem-70087525**: SCF tarda 45min con AM1-BCC para esta
   molécula (76 átomos, +2). Considerar si tiene sentido incluirla en
   campañas o marcarla como caso especial. Decisión científica, no técnica.

4. **Re-evaluar campañas históricas con `degraded_charges=True`**: una vez
   implementado el flag, identificar moléculas pasadas que cayeron a Tier 2
   silenciosamente y considerar si vale re-correr 00a→01e con el fix
   nuevo para validar/corregir su ranking DOCK6 original. Decisión por
   campaña.

---

## RESUMEN DE ARCHIVOS A MODIFICAR

```
01_src/hit_validation/m00_preparation/ligand_preparation.py    (cambio 1)
03_configs/00a_ligand_preparation.yaml                         (cambio 1)
02_scripts/00a_ligand_preparation.py                           (cambio 1, reporte)
01_src/hit_validation/m01_validation/mmpbsa_decomp.py          (cambios 2 y 3)
03_configs/01g_mmpbsa_decomp.yaml                              (cambio 2)
```

## REGLAS NEGATIVAS

- NO usar el paquete `ionization` (`ionprofile`) como dependencia runtime.
- NO eliminar las funciones `run_antechamber_tier2` ni `run_obabel_tier3` —
  se mantienen como fallbacks pero con warnings y flag `degraded_charges`.
- NO incluir cambios a `molecular_docking` en este instructivo.
- NO usar paths absolutos hardcodeados a binarios (antechamber, obabel,
  cpptraj). Mantener uso del PATH como antes — es responsabilidad del entorno.
- NO modificar `convert_sdf_to_mol2()` ni `_cleanup_intermediates()` ni el
  flujo de `00b`/`00d`/`01b`/`01c`/`01d`/`01e`/`01f`/`01h`/`01i`/`03a`/`04b`.
- NO cambiar `mmpbsa_interval_ps` ni nada del módulo 07e — solo verificar
  que sigue siendo coherente con stride 100 / 50 frames.
