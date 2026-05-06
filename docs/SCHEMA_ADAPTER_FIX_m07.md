# Schema Adapter Fix — m07_decision_report

**Goal:** Fix the schema mismatch that breaks 07a/07b/07c/07d when `n_replicas > 1`.

**Scope:** Targeted fix only. Add a schema adapter helper that flattens consolidated columns (`*_mean`) to plain names (`total`, `vdw`, `es`, etc.) when present. Apply it consistently in the 4 modules. Re-run Phase 4 to validate.

**No rewrites. No new modules. No touching 07e.**

---

## Topology

- Working directory: `/home/bprieto/projects/hit_validation`
- Conda env: `hit_validation_env` (Python 3.12.13). Do not recreate.
- The `.git/` directory does NOT exist on the server. Do NOT run `git commit`, `git push`, `git tag`, or any GitHub command. The user commits from PyCharm Windows.
- DO NOT touch `~/projects/mmpbsa_analysis/`.
- DO NOT touch 07e (`temporal_analysis.py`, `07e_temporal_analysis.py`, etc.).
- The pipeline for SD1_druglikeness_2_PChem_136488320_2n is FINISHED. No live processes to worry about.

---

## Context

The pipeline finished with warnings:

```
PIPELINE COMPLETE — with warnings
Modules with errors: 06b 07c 07a 07b 07d
```

The root cause for the 07x failures (analyzed and confirmed) is a schema mismatch:

- **`replica_N/04b_footprint_analysis/footprint_per_molecule.csv`** has plain columns: `vdw, es, total, ref_vdw, ref_es, ref_total, delta_vdw, delta_es, delta_total, fps_score, hb_pose, hb_ref`.
- **`04b_footprint_analysis/consolidated/footprint_per_molecule.csv`** has only suffixed columns: `vdw_mean, vdw_std, vdw_sem, es_mean, ..., total_mean, ..., n_replicas_used`. The plain columns DO NOT exist.

Modules 07a, 07b, 07c, 07d read literal `df["total"]`, `df["vdw"]`, `df["es"]`, etc., and crash with `KeyError: 'total'`.

Notably:
- 07a (`decision_report.py`) ALREADY has the adapter pattern for `Grid_Score_mean` (lines 783-784) and `MMPBSA_dG_mean` (lines 650-654). The fix is to replicate the same pattern for footprint columns.
- `04b/consolidated/residue_consensus.csv` uses semantic names (`mean_total`, `mean_vdw`) instead of suffixes — that's why 06a (which reads it) works fine. We do NOT touch the 04b consolidator's output schema; we only adapt readers.

The fix is **targeted and additive**: a single shared helper function applied to 4 modules.

---

## What you must do

### Step 1 — Add a shared schema adapter helper

Create a new private utility (single function) used by all 4 modules. The cleanest place is to add it to `01_src/hit_validation/m07_decision_report/__init__.py` or a new file `01_src/hit_validation/m07_decision_report/_footprint_loader.py`.

**Recommended:** create `_footprint_loader.py` with a single function:

```python
"""
Footprint CSV loader with schema adapter.

Handles both single-replica (plain columns) and consolidated (*_mean suffix)
schemas, exposing a unified API to downstream m07 modules.
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterable

import pandas as pd

# Columns that the m07 consumers expect under their plain names.
# When the input CSV is from 04b/consolidated/, these come as <col>_mean.
_FOOTPRINT_NUMERIC_COLS: tuple[str, ...] = (
    "vdw", "es", "total",
    "ref_vdw", "ref_es", "ref_total",
    "delta_vdw", "delta_es", "delta_total",
    "fps_score", "hb_pose", "hb_ref",
)


def load_footprint_csv(
    path: Path | str,
    *,
    extra_required: Iterable[str] = (),
) -> pd.DataFrame:
    """
    Load a footprint CSV emitted by 04b. Transparently handles two schemas:

    1. Single-replica: columns are plain (`vdw`, `es`, `total`, ...).
    2. Consolidated (N>1): columns have `_mean` / `_std` / `_sem` suffix
       (e.g. `total_mean`). In that case we expose the plain names as
       aliases of the `_mean` columns. The `_std` and `_sem` columns
       are preserved as-is for downstream code that wants uncertainty.

    Args:
        path: Path to the CSV file.
        extra_required: Additional columns the caller needs. Raises if missing.

    Returns:
        A DataFrame with at least the columns in `_FOOTPRINT_NUMERIC_COLS`
        present (as plain names) plus whatever the original file already had.

    Raises:
        FileNotFoundError: if the path does not exist.
        ValueError: if neither plain nor `_mean`-suffixed columns are present.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Footprint CSV not found: {path}")

    df = pd.read_csv(path)

    # If the consolidated schema is present, alias <col>_mean → <col>.
    # Do NOT overwrite a column that already exists (single-replica wins).
    for col in _FOOTPRINT_NUMERIC_COLS:
        if col not in df.columns:
            mean_col = f"{col}_mean"
            if mean_col in df.columns:
                df[col] = df[mean_col]

    # Verify required columns are now present.
    missing = [c for c in _FOOTPRINT_NUMERIC_COLS if c not in df.columns]
    extra_missing = [c for c in extra_required if c not in df.columns]
    if missing or extra_missing:
        all_missing = sorted(set(missing) | set(extra_missing))
        raise ValueError(
            f"Missing columns in footprint CSV {path}: {all_missing}. "
            f"Available columns: {sorted(df.columns)}"
        )

    return df
```

Save it at: `01_src/hit_validation/m07_decision_report/_footprint_loader.py`

### Step 2 — Apply the adapter in the 4 modules

Find every place each of the 4 m07 modules calls `pd.read_csv(...)` on a footprint CSV and replace it with `load_footprint_csv(...)`.

Specifically (locations confirmed by previous diagnostic):

**`01_src/hit_validation/m07_decision_report/integrated_analysis.py`**
- Lines 261, 275-289: replace direct `df["total"]`, `df["vdw"]`, `df["es"]` access with the adapter pattern.
- Find the `pd.read_csv(footprint_csv)` call upstream and replace it with `load_footprint_csv(footprint_csv)`.

**`01_src/hit_validation/m07_decision_report/decision_report.py`**
- Line 750: `zone_rows["total"].sum()` is the trigger.
- Find the `pd.read_csv(footprint_csv)` upstream and replace it with `load_footprint_csv(footprint_csv)`.

**`01_src/hit_validation/m07_decision_report/residue_comparison.py`**
- Lines 47-49: the explicit `required = {Name, residue_id, vdw, es, total, ref_vdw, ...}` validation set.
- Replace the loader and remove the redundant `required` set check, since `load_footprint_csv` does it. Keep the `Name`, `residue_id` validation if present (those are not in the adapter's required list).

**`01_src/hit_validation/m07_decision_report/binding_mode_analysis.py`**
- Line 63: `df.groupby("residue_id")["total"]` triggers the crash.
- Lines 72, 100, 106 use `pivot(values="total" / "vdw" / "es")`.
- Lines 282, 385-393, 409: aggregations like `mol_df["total"]`.
- Replace the `pd.read_csv(footprint_csv)` upstream with `load_footprint_csv(footprint_csv)`. The downstream lines work as-is once the plain columns exist.

**Important conventions:**
- Use a top-of-file import: `from hit_validation.m07_decision_report._footprint_loader import load_footprint_csv`
- The function signature accepts `extra_required` for callers that need additional columns (e.g. `Name`, `residue_id`, `pose_id`). Use it where appropriate.
- Do NOT change the downstream code that uses `df["total"]`, `df["vdw"]`, etc. The adapter makes those columns exist.

### Step 3 — Verification gate

Run all of these and confirm green:

```bash
# 1. Bash syntax check (defensive, run_pipeline.sh shouldn't have changed)
bash -n run_pipeline.sh

# 2. New helper imports cleanly
PYTHONPATH=01_src /home/bprieto/mambaforge/envs/hit_validation_env/bin/python -c "
from hit_validation.m07_decision_report._footprint_loader import load_footprint_csv
import inspect
sig = inspect.signature(load_footprint_csv)
print('helper signature:', sig)
"

# 3. The 4 modules import cleanly
PYTHONPATH=01_src /home/bprieto/mambaforge/envs/hit_validation_env/bin/python -c "
from hit_validation.m07_decision_report import integrated_analysis, decision_report, residue_comparison, binding_mode_analysis
print('all 4 modules import OK')
"

# 4. CLI --help still works for all 4
for mod in 07a_decision_report 07b_residue_comparison 07c_integrated_analysis 07d_binding_mode_analysis; do
    echo "--- $mod ---"
    PYTHONPATH=01_src /home/bprieto/mambaforge/envs/hit_validation_env/bin/python 02_scripts/${mod}.py --help 2>&1 | tail -5 || echo "FAILED"
done

# 5. Real-world smoke test: load the consolidated CSV that crashed before
PYTHONPATH=01_src /home/bprieto/mambaforge/envs/hit_validation_env/bin/python -c "
from hit_validation.m07_decision_report._footprint_loader import load_footprint_csv
df = load_footprint_csv('05_results/SD1_druglikeness_2_PChem_136488320_2n/04b_footprint_analysis/consolidated/footprint_per_molecule.csv')
required = ['total', 'vdw', 'es', 'ref_total', 'delta_total']
for col in required:
    assert col in df.columns, f'MISSING: {col}'
    assert df[col].notna().any(), f'ALL NULL: {col}'
print(f'Consolidated CSV: {len(df)} rows, columns OK')
print(f'  Sample: total[0]={df[\"total\"].iloc[0]:.4f}, vdw[0]={df[\"vdw\"].iloc[0]:.4f}')
"

# 6. Real-world smoke test: load a single-replica CSV (must still work)
PYTHONPATH=01_src /home/bprieto/mambaforge/envs/hit_validation_env/bin/python -c "
from hit_validation.m07_decision_report._footprint_loader import load_footprint_csv
df = load_footprint_csv('05_results/SD1_druglikeness_2_PChem_136488320_2n/replica_1/04b_footprint_analysis/footprint_per_molecule.csv')
required = ['total', 'vdw', 'es', 'ref_total', 'delta_total']
for col in required:
    assert col in df.columns, f'MISSING: {col}'
    assert df[col].notna().any(), f'ALL NULL: {col}'
print(f'Single-replica CSV: {len(df)} rows, columns OK')
print(f'  Sample: total[0]={df[\"total\"].iloc[0]:.4f}, vdw[0]={df[\"vdw\"].iloc[0]:.4f}')
"
```

All 6 steps must report success. If any step fails, STOP and report.

### Step 4 — Re-run Phase 4 to validate end-to-end

After the verification gate passes, re-run Phase 4 (06a→07d) of the existing campaign to confirm the fix works in production:

```bash
cd ~/projects/hit_validation
nohup bash run_pipeline.sh \
  04_data/campaigns/SD1_druglikeness_2_PChem_136488320_2n/campaign_config.yaml \
  06a > /tmp/phase4_rerun.log 2>&1 &
echo "PID: $!"
disown

# Wait ~3 minutes, then check
sleep 180
tail -50 /tmp/phase4_rerun.log
```

Expected outcome:
- `06a` produces pharmacophore JSONs (already worked before — check it still does).
- `06b` may still fail (pharmit JSON not found is a separate bug, not in scope of this fix).
- `07c`, `07a`, `07b`, `07d` MUST complete without `KeyError: 'total'`.
- The pipeline emits `PIPELINE COMPLETE` (with or without warnings on 06b only).

If any of 07a/07b/07c/07d still fails with `KeyError: 'total'` or similar after the re-run → STOP, report which module and which line.

---

## Constraints (do not violate)

- DO NOT modify any module beyond:
  - new file `01_src/hit_validation/m07_decision_report/_footprint_loader.py`
  - existing files `integrated_analysis.py`, `decision_report.py`, `residue_comparison.py`, `binding_mode_analysis.py`
- DO NOT modify `04b_footprint_analysis.py` or its consolidator. The consolidated schema stays as-is.
- DO NOT modify any 07x CLI script (`02_scripts/07*_*.py`) unless absolutely necessary for an import path.
- DO NOT modify `run_pipeline.sh`.
- DO NOT touch `04_data/campaigns/`.
- DO NOT delete anything from `05_results/`, `audit_outputs/`, `test_outputs/`, `logs/`, `.claude/`, etc.
- DO NOT run any `git` command.
- DO NOT touch `~/projects/mmpbsa_analysis/`.
- DO NOT touch 07e (`02_scripts/07e_temporal_analysis.py`, `01_src/hit_validation/m07_decision_report/temporal_analysis.py`, `03_configs/07e_temporal_analysis.yaml`).

---

## Output expected

When verification gate passes AND Phase 4 re-runs successfully, emit exactly this block:

```
=== SCHEMA ADAPTER FIX COMPLETE ===

New file:
  - 01_src/hit_validation/m07_decision_report/_footprint_loader.py
      • load_footprint_csv(path, extra_required) helper

Files modified:
  - 01_src/hit_validation/m07_decision_report/integrated_analysis.py
      • replaced pd.read_csv with load_footprint_csv
  - 01_src/hit_validation/m07_decision_report/decision_report.py
      • replaced pd.read_csv with load_footprint_csv
  - 01_src/hit_validation/m07_decision_report/residue_comparison.py
      • replaced pd.read_csv with load_footprint_csv
      • removed redundant 'required' validation set
  - 01_src/hit_validation/m07_decision_report/binding_mode_analysis.py
      • replaced pd.read_csv with load_footprint_csv

Verification gate: all green
  • bash -n run_pipeline.sh: clean
  • _footprint_loader imports: OK
  • all 4 m07 modules import: OK
  • CLI --help (07a, 07b, 07c, 07d): OK
  • smoke test consolidated CSV: OK (rows + columns + non-null)
  • smoke test single-replica CSV: OK (rows + columns + non-null)

Phase 4 re-run: <STATUS>
  • 06a: <ok | fail>
  • 06b: <ok | fail (expected: pharmit JSON, separate bug)>
  • 07c: <ok | fail>
  • 07a: <ok | fail>
  • 07b: <ok | fail>
  • 07d: <ok | fail>
  • Final outputs available:
      - 07a_decision_report/{decision_report.html, decision_summary.csv, pose_selection_summary.csv, selected_poses/}
      - 07b_residue_comparison/{residue_comparison.html, residue_comparison_summary.csv, per_molecule/}
      - 07c_integrated_analysis/{integrated_report.html, integrated_summary.csv, unified_verdict.csv, consistency_matrix.csv, per_molecule/}
      - 07d_binding_mode_analysis/{binding_mode_report.html, iev_similarity_matrix.csv, molecule_clusters.csv, residue_cocontact_network.csv, emergent_subpockets.csv}

Suggested commit message (for user from Windows):
fix(m07): schema adapter for consolidated footprint CSV (unblocks 07a-d when n_replicas > 1)

Awaiting user approval. No git operations performed.
```

---

## Stop conditions

If any of these happen, STOP and report — do not improvise:

- `bash -n run_pipeline.sh` returns syntax error.
- New helper or modified module import fails.
- Smoke test shows null or missing columns after the adapter runs.
- Phase 4 re-run fails on 07a/07b/07c/07d with anything OTHER than `KeyError: 'total'` would indicate the fix didn't apply or there is a deeper issue.
- A module needs more than ~5 line changes to integrate the loader (suggests deeper coupling and warrants discussion before continuing).

Proceed.
