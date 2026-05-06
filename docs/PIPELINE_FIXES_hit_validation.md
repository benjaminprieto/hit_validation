# Pipeline Fixes — hit_validation

**Goal:** fix two real bugs in `run_pipeline.sh` and `01i_trajectory_analysis` that were exposed by the SD1 GB-only campaign log (2026-05-02).

**Scope:** This is NOT a refactor. Two surgical fixes only. Do not modify anything else.

---

## Topology

- Working directory: `/home/bprieto/projects/hit_validation`
- Conda env: `hit_validation_env` (Python 3.12.13). Do not recreate or downgrade.
- The `.git/` directory does NOT exist on the server — git lives only on Windows local. Do NOT run `git commit`, `git push`, `git tag`, or any GitHub command. The user commits from PyCharm Windows.
- DO NOT touch `~/projects/mmpbsa_analysis/` (separate project, sealed at v1.0.0).

---

## Context — what the SD1 log revealed

A campaign just finished `replica_1` (SD1_druglikeness_2_PChem_136488320_2n_GBonly). `replica_2` is still running on the server (do NOT kill it, do NOT modify anything that affects an in-progress run). The log revealed two bugs that degrade output quality but were not fatal.

### Bug 1 — Module execution order is wrong inside each replica

**Where:** `run_pipeline.sh`, line 112:

```bash
REPLICATED_MODULES=(01c 01d 01e 01f 01g 01h 01i 03a 04b)
```

The execution order is taken **from this array**. Comments on lines 230-231 explicitly confirm: *"The order of case clauses below is historical and does not affect execution."*

**Problem:**

- `04b_footprint_analysis` is the module that **generates** `residue_mapping.csv`.
- `01h_mmpbsa_analysis` and `01i_trajectory_analysis` **read** `residue_mapping.csv`.
- Currently 04b runs LAST → 01h and 01i can't find the mapping → degraded output.

**Evidence in the log:**

```
01h:  Footprint: not found
01i:  WARNING — No residue mapping found. Distance monitoring will be limited.
01i:  Monitored residues: 0
```

This means:
- 01h skips the MMPBSA-vs-DOCK6 footprint comparison.
- 01i does NOT compute per-frame distances to functional residues (E529, D494, R598, K599, W392) — only `d_min` is computed.

### Bug 2 — Water bridge analysis (cpptraj) times out at 3600 s

**Where:** `01_src/hit_validation/m01_docking/trajectory_analysis.py`, lines 144-156 (the `CPPTRAJ_SOLVATED_TEMPLATE` constant):

```
parm {solvated_prmtop}
trajin {production_dcd}        # ← reads ALL 5000 frames, no stride
strip :Na+,Cl-
hbond WB ...
```

**Problem:** for a 78,857-atom solvated system × 5000 frames, cpptraj exceeds the 3600 s timeout.

**Evidence in the log:**

```
WARNING — cpptraj timed out after 3600s
WARNING — Water bridge analysis failed: Timeout after 3600s
```

Both ligands of replica_1 lost their water bridge analysis.

---

## Fix 1 — Module order in `run_pipeline.sh`

### 1.1 Edit line 112

**Before:**
```bash
REPLICATED_MODULES=(01c 01d 01e 01f 01g 01h 01i 03a 04b)
```

**After:**
```bash
REPLICATED_MODULES=(01c 01d 01e 01f 01g 04b 01h 01i 03a)
```

**Rationale:**
- 04b moves to position 6 (right after 01g) → `residue_mapping.csv` exists when 01h and 01i run.
- 01h and 01i remain consecutive (positions 7-8), reading the mapping correctly.
- 03a (PLIP) moves to the end. PLIP is independent of `residue_mapping.csv`, so this is safe.

**Independence verified:** `grep -n "01h\|01i_trajectory\|mmpbsa" 01_src/hit_validation/m04_dock6_analysis/footprint_analysis.py` returned EMPTY. 04b reads only the receptor (mol2/PDB) and footprint files from 01d. It does NOT depend on 01h, 01i, or any MMPBSA output.

### 1.2 Update line 110 (`ALL_MODULES` array)

**Before:**
```bash
ALL_MODULES=(00a 00b 00d 01b 01c 01d 01e 01f 01g 01h 01i 03a 04b 06a 06b 07c 07a 07b 07d)
```

**After:**
```bash
ALL_MODULES=(00a 00b 00d 01b 01c 01d 01e 01f 01g 04b 01h 01i 03a 06a 06b 07c 07a 07b 07d)
```

### 1.3 Update line 30 (help message)

**Before:**
```bash
echo "Modules: 00a 00b 00d 01b 01c 01d 01f 01e 03a 04b 01g 01h 01i 06a 06b 07c 07a 07b 07d"
```

**After:**
```bash
echo "Modules: 00a 00b 00d 01b 01c 01d 01f 01e 01g 04b 01h 01i 03a 06a 06b 07c 07a 07b 07d"
```

### 1.4 Update lines 15-16 (header comment)

**Before:**
```
#   Phase 2 (REPLICATED, N×):    01c → 01d → 01e → 01f → 01g → 01h → 01i → 03a → 04b
```

**After:**
```
#   Phase 2 (REPLICATED, N×):    01c → 01d → 01e → 01f → 01g → 04b → 01h → 01i → 03a
```

### 1.5 Do NOT touch the consolidation phase

The "Phase 3 (CONSOLIDATION)" block (lines 270+) handles N>1 consolidation. Its order is correct as-is: `01h_consolidate → select_representative → 01c/01e/03a/04b/01i_consolidate`. **Do not modify Phase 3.**

---

## Fix 2 — Water bridge stride + configurable timeout

### 2.1 Add stride placeholder to `CPPTRAJ_SOLVATED_TEMPLATE`

**File:** `01_src/hit_validation/m01_docking/trajectory_analysis.py`, lines 144-156.

**Before:**
```python
CPPTRAJ_SOLVATED_TEMPLATE = """\
# 01i Water Bridge Analysis -- solvated trajectory
parm {solvated_prmtop}
trajin {production_dcd}

# Strip ions, keep water
strip :Na+,Cl-

# Water bridges: ...
hbond WB donormask ... \
  bridgeout {output_dir}/water_bridges.dat

run
quit
"""
```

**After:**
```python
CPPTRAJ_SOLVATED_TEMPLATE = """\
# 01i Water Bridge Analysis -- solvated trajectory
parm {solvated_prmtop}
trajin {production_dcd} 1 last {stride}

# Strip ions, keep water
strip :Na+,Cl-

# Water bridges: ...
hbond WB donormask ... \
  bridgeout {output_dir}/water_bridges.dat

run
quit
"""
```

**cpptraj syntax note:** `trajin file 1 last N` reads from frame 1 to the last frame, with stride N (every Nth frame). `N=1` is "every frame" (current behavior), `N=10` is "1 of every 10 frames".

### 2.2 Add `water_bridge_stride` parameter to the call chain

The function that builds the cpptraj solvated script and calls `.format()` on `CPPTRAJ_SOLVATED_TEMPLATE` (find it via `grep -n "CPPTRAJ_SOLVATED_TEMPLATE" 01_src/hit_validation/m01_docking/trajectory_analysis.py`) needs a new parameter.

The chain to update (search for `cpptraj_solvated_timeout` to find all callers):

1. The function that uses `CPPTRAJ_SOLVATED_TEMPLATE.format(...)` → add `water_bridge_stride: int = 10` to its signature, and include `stride=water_bridge_stride` in the `.format()` call.
2. The public function around line 686 that takes `cpptraj_solvated_timeout: int = 3600` → also add `water_bridge_stride: int = 10` and pass it down.
3. Any caller above that propagates the timeout kwarg.

**Default value: 10** (means 1 of every 10 frames, so 5000 → 500 frames typical).

### 2.3 Add YAML key

**File:** `03_configs/01i_trajectory_analysis.yaml`

After the existing `cpptraj_solvated_timeout: 3600` line (around line 28), add:

```yaml
  # Water bridge analysis stride. Reduces frames passed to cpptraj.
  # stride=10 means use 1 of every 10 frames (5000 → 500 frames typical).
  # stride=1 = every frame (slow on large solvated systems).
  water_bridge_stride: 10
```

### 2.4 Read YAML key in the CLI script

**File:** `02_scripts/01i_trajectory_analysis.py`

Wherever the script currently reads `cpptraj_solvated_timeout` from the module config, also read `water_bridge_stride` (with fallback default `10` if the key is missing) and pass it to the core trajectory_analysis function as `water_bridge_stride=...`.

---

## Verification gate

After both fixes, run all of these and confirm green:

```bash
# 1. Bash syntax check
bash -n run_pipeline.sh

# 2. Python imports (with PYTHONPATH set to 01_src/ as the project requires)
PYTHONPATH=01_src python -c "from hit_validation.m01_docking.trajectory_analysis import CPPTRAJ_SOLVATED_TEMPLATE; assert '{stride}' in CPPTRAJ_SOLVATED_TEMPLATE; print('import + stride placeholder OK')"

# 3. YAML parses + new key present
python -c "import yaml; cfg = yaml.safe_load(open('03_configs/01i_trajectory_analysis.yaml')); assert 'water_bridge_stride' in str(cfg), 'water_bridge_stride key missing'; print('yaml + new key OK')"

# 4. CLI --help still works
PYTHONPATH=01_src python 02_scripts/01i_trajectory_analysis.py --help

# 5. Confirm new module order
grep "^REPLICATED_MODULES" run_pipeline.sh
# Expected: REPLICATED_MODULES=(01c 01d 01e 01f 01g 04b 01h 01i 03a)

# 6. Confirm CPPTRAJ_SOLVATED_TEMPLATE has the stride placeholder
grep -A 1 "trajin {production_dcd}" 01_src/hit_validation/m01_docking/trajectory_analysis.py
# Expected output should include: trajin {production_dcd} 1 last {stride}
```

---

## Constraints (do not violate)

- DO NOT modify any module beyond what's listed above (`run_pipeline.sh`, `trajectory_analysis.py`, `01i_trajectory_analysis.py` script, `01i_trajectory_analysis.yaml`).
- DO NOT touch `04_data/campaigns/`.
- DO NOT delete anything from `05_results/`, `audit_outputs/`, `test_outputs/`, `logs/`, `.claude/`, etc.
- DO NOT kill the running pipeline processes (SD1 r2 and DIOS_pharmit_pH63_1).
- DO NOT run any `git` command.
- DO NOT touch `~/projects/mmpbsa_analysis/`.
- The fixes apply to FUTURE runs. The currently-running SD1 r2 will finish with the OLD order; selective re-run of `04b → 01h → 01i` for r1 and r2 is a SEPARATE task, not part of this one.

---

## Output expected

When verification gate passes, emit exactly this block:

```
=== PIPELINE FIXES COMPLETE ===

Files modified:
  - run_pipeline.sh
      • line 15-16: header comment updated
      • line 30: help message updated
      • line 110: ALL_MODULES order updated
      • line 112: REPLICATED_MODULES order updated
  - 01_src/hit_validation/m01_docking/trajectory_analysis.py
      • CPPTRAJ_SOLVATED_TEMPLATE: added {stride} placeholder
      • call chain: added water_bridge_stride parameter (default 10)
  - 02_scripts/01i_trajectory_analysis.py
      • reads water_bridge_stride from YAML, passes through to core
  - 03_configs/01i_trajectory_analysis.yaml
      • new key: water_bridge_stride (default 10, commented)

Verification gate: all green
  • bash -n run_pipeline.sh: clean
  • python imports: OK
  • yaml parses: OK
  • 01i --help: OK
  • REPLICATED_MODULES = (01c 01d 01e 01f 01g 04b 01h 01i 03a): verified
  • CPPTRAJ_SOLVATED_TEMPLATE has {stride}: verified

Suggested commit message (for user to use from Windows):
fix: pipeline module order (04b before 01h/01i) + water bridge cpptraj stride configurable

Awaiting user approval. No git operations performed.
```

---

## Stop conditions

If any of these happen, STOP and report — do not improvise:

- `bash -n run_pipeline.sh` returns syntax error.
- Python import fails.
- YAML doesn't parse.
- The call chain in `trajectory_analysis.py` is more complex than expected (more than 2 functions need to thread the new kwarg).

## Concurrency safety note

Two pipelines are currently running on the server (SD1 r2 and DIOS_pharmit_pH63_1). Editing files is generally safe because:

- `run_pipeline.sh`: bash loaded the script into memory at launch via `nohup bash run_pipeline.sh ...`. Editing the file on disk does NOT affect the running process. Future invocations get the new version.
- `trajectory_analysis.py`: SD1 r2 is currently in MD production (50ns, takes ~9-10h per molecule). It won't reach 01i until much later. By the time r2 invokes 01i, the file is already updated — so r2 will pick up the new version naturally.
- `01i_trajectory_analysis.py` script and YAML config: same as above. Only 01i invocation reads them.

If `lsof | grep trajectory_analysis.py` returns active reads from another bash process at the time of editing, that means the file is being read RIGHT NOW. In that case, wait a few minutes and retry. This is extremely unlikely given current MD progress.

Proceed.
