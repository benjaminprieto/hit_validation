#!/bin/bash
# =============================================================================
# run_pipeline.sh — Hit Validation Pipeline (con soporte de réplicas)
# =============================================================================
#
# Usage:
#   bash run_pipeline.sh <campaign_config.yaml> [start_from]
#
# Examples:
#   bash run_pipeline.sh 04_data/campaigns/hit_validation/campaign_config.yaml
#   bash run_pipeline.sh 04_data/campaigns/hit_validation/campaign_config.yaml 01c
#
# Pipeline phases:
#   Phase 1 (SHARED, 1×):        00a → 00b → 00d → 01b
#   Phase 2 (REPLICATED, N×):    01c → 01d → 01e → 01f → 01g → 04b → 01h → 01i → 03a
#   Phase 3 (CONSOLIDATION):     01h_consolidate → select_representative
#                                → 01c/01e/03a/04b/01i_consolidate    (only when N>1)
#   Phase 4 (POST-REPLICA):      06a → 06b → 07c → 07a → 07b → 07d
#       (07c emits unified_verdict.csv consumed by 07a, so 07c precedes 07a)
#
# n_replicas is read with precedence: campaign_config > pipeline_config > default 1
# =============================================================================
set -uo pipefail
# NOTE: no 'set -e' — we handle errors per-module to avoid aborting the
# entire pipeline when a batch module reports partial failures.

if [ $# -lt 1 ]; then
    echo "Usage: bash run_pipeline.sh <campaign_config.yaml> [start_from]"
    echo ""
    echo "Modules: 00a 00b 00d 01b 01c 01d 01f 01e 01g 04b 01h 01i 03a 06a 06b 07c 07a 07b 07d"
    echo ""
    echo "Examples:"
    echo "  bash run_pipeline.sh config.yaml        # full pipeline"
    echo "  bash run_pipeline.sh config.yaml 01c    # start from docking"
    echo "  bash run_pipeline.sh config.yaml 01g    # start from MMPBSA"
    echo "  bash run_pipeline.sh config.yaml 03a    # start from PLIP"
    exit 1
fi

CAMPAIGN="$1"
START="${2:-00a}"
CONFIGS="03_configs"
SCRIPTS="02_scripts"
PIPELINE_CONFIG="${CONFIGS}/pipeline_config.yaml"

if [ ! -f "$CAMPAIGN" ]; then
    echo "ERROR: Campaign config not found: $CAMPAIGN"
    exit 1
fi

if [ ! -f "$PIPELINE_CONFIG" ]; then
    echo "ERROR: Pipeline config not found: $PIPELINE_CONFIG"
    echo "Create it with: n_replicas: 1   (or N for replicated mode)"
    exit 1
fi

# Read n_replicas with precedence: campaign_config > pipeline_config > default 1
N_REPLICAS=$(python3 << PYEOF
import yaml
import sys

with open("${CAMPAIGN}") as f:
    campaign = yaml.safe_load(f) or {}
with open("${PIPELINE_CONFIG}") as f:
    pipeline = yaml.safe_load(f) or {}

n = campaign.get("n_replicas")
source = "campaign_config"
if n is None:
    n = pipeline.get("parameters", {}).get("n_replicas", 1)
    source = "pipeline_config (fallback)"

if not isinstance(n, int) or n < 1:
    print(f"ERROR: n_replicas must be int >= 1, got {n!r}", file=sys.stderr)
    sys.exit(1)

# Print as: <value>:<source>
print(f"{n}:{source}")
PYEOF
)

if [ -z "$N_REPLICAS" ]; then
    echo "ERROR: failed to determine n_replicas"
    exit 1
fi

# Parse "value:source" output
N_REPLICAS_SOURCE="${N_REPLICAS#*:}"
N_REPLICAS="${N_REPLICAS%%:*}"

if [ "$N_REPLICAS" -lt 1 ]; then
    echo "ERROR: n_replicas must be >= 1 (got $N_REPLICAS)"
    exit 1
fi

echo "[CONFIG] n_replicas = $N_REPLICAS (source: $N_REPLICAS_SOURCE)"

# Check if precomputed_grids_dir is set in campaign config (skips 00d/01b)
PRECOMPUTED_GRIDS=$(python3 -c "
import yaml
with open('$CAMPAIGN') as f:
    cc = yaml.safe_load(f)
grids = cc.get('grids', {})
d = grids.get('precomputed_grids_dir')
print(d if d else '')
" 2>/dev/null || echo "")

# All modules in canonical order — used for start-from logic.
# 07c precedes 07a because 07c emits unified_verdict.csv consumed by 07a.
ALL_MODULES=(00a 00b 00d 01b 01c 01d 01e 01f 01g 04b 01h 01i 03a 06a 06b 07c 07a 07b 07d)
SHARED_MODULES=(00a 00b 00d 01b)
REPLICATED_MODULES=(01c 01d 01e 01f 01g 04b 01h 01i 03a)

# Find start index
START_IDX=0
FOUND=0
for i in "${!ALL_MODULES[@]}"; do
    if [ "${ALL_MODULES[$i]}" = "$START" ]; then
        START_IDX=$i
        FOUND=1
        break
    fi
done
if [ "$FOUND" -eq 0 ]; then
    echo "ERROR: Unknown start module: $START"
    echo "Valid modules: ${ALL_MODULES[*]}"
    exit 1
fi

echo "============================================================"
echo "  HIT_VALIDATION — Pipeline (replicas-aware)"
echo "  Campaign:    $CAMPAIGN"
echo "  Start from:  $START"
echo "  Replicas:    $N_REPLICAS"
echo "  Started:     $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"
echo ""

FAILED_MODULES=()

# Helper: index of module in ALL_MODULES (for start-from filtering)
mod_index() {
    local target="$1"
    for i in "${!ALL_MODULES[@]}"; do
        if [ "${ALL_MODULES[$i]}" = "$target" ]; then
            echo "$i"
            return 0
        fi
    done
    echo "-1"
}

# Helper: should we run this module (based on START_IDX)?
should_run() {
    local mod="$1"
    local idx
    idx=$(mod_index "$mod")
    [ "$idx" -ge "$START_IDX" ]
}

# Helper: run a module with optional --replica-id
run_module() {
    local mod="$1"
    local script="$2"
    local config="$3"
    local label="$4"
    local replica_id="${5:-}"

    local extra_args=()
    if [ -n "$replica_id" ]; then
        extra_args=(--replica-id "$replica_id")
        echo "[$mod r$replica_id] $label — $(date '+%H:%M:%S')"
    else
        echo "[$mod] $label — $(date '+%H:%M:%S')"
    fi

    if python "$SCRIPTS/$script" \
        --config "$CONFIGS/$config" \
        --campaigns "$CAMPAIGN" \
        "${extra_args[@]}"; then
        echo ""
    else
        local rc=$?
        echo "[$mod] WARNING: $label exited with code $rc — continuing pipeline"
        echo ""
        FAILED_MODULES+=("${mod}${replica_id:+_r$replica_id}")
    fi
}

# ─────────────────────────────────────────────
# PHASE 1: SHARED (always 1×)
# ─────────────────────────────────────────────
for mod in "${SHARED_MODULES[@]}"; do
    should_run "$mod" || continue

    case "$mod" in
        00a) run_module 00a 00a_ligand_preparation.py 00a_ligand_preparation.yaml "Ligand Preparation" ;;
        00b) run_module 00b 00b_receptor_preparation.py 00b_receptor_preparation.yaml "Receptor Preparation" ;;
        00d)
            if [ -n "$PRECOMPUTED_GRIDS" ]; then
                echo "[00d] Binding Site Definition — SKIPPED (precomputed_grids_dir set)"
                echo ""
            else
                run_module 00d 00d_binding_site_definition.py 00d_binding_site_definition.yaml "Binding Site Definition"
            fi
            ;;
        01b)
            if [ -n "$PRECOMPUTED_GRIDS" ]; then
                echo "[01b] Grid Generation — SKIPPED (precomputed_grids_dir set)"
                echo ""
            else
                run_module 01b 01b_grid_generation.py 01b_grid_generation.yaml "Grid Generation"
            fi
            ;;
    esac
done

# ─────────────────────────────────────────────
# PHASE 2: REPLICATED
# ─────────────────────────────────────────────
run_replicated_modules() {
    local replica_id="$1"  # empty string = legacy, otherwise replica index

    for mod in "${REPLICATED_MODULES[@]}"; do
        should_run "$mod" || continue

        # Note: order of execution comes from REPLICATED_MODULES variable above.
        # The order of case clauses below is historical and does not affect execution.
        case "$mod" in
            01c) run_module 01c 01c_dock6_run.py 01c_dock6_run.yaml "DOCK6 Docking" "$replica_id" ;;
            01d) run_module 01d 01d_footprint_rescore.py 01d_footprint_rescore.yaml "Footprint Rescore" "$replica_id" ;;
            01f) run_module 01f 01f_gbsa_rescore.py 01f_gbsa_rescore.yaml "GB/SA Hawkins Rescore" "$replica_id" ;;
            01e) run_module 01e 01e_score_collection.py 01e_score_collection.yaml "Score Collection" "$replica_id" ;;
            03a) run_module 03a 03a_plip_interaction_analysis.py 03a_plip_interaction_analysis.yaml "PLIP Analysis" "$replica_id" ;;
            04b) run_module 04b 04b_footprint_analysis.py 04b_footprint_analysis.yaml "Footprint Analysis" "$replica_id" ;;
            01g) run_module 01g 01g_mmpbsa_decomp.py 01g_mmpbsa_decomp.yaml "MMPBSA Decomposition" "$replica_id" ;;
            01h) run_module 01h 01h_mmpbsa_analysis.py 01h_mmpbsa_analysis.yaml "MMPBSA Analysis" "$replica_id" ;;
            01i) run_module 01i 01i_trajectory_analysis.py 01i_trajectory_analysis.yaml "Trajectory Analysis" "$replica_id" ;;
        esac
    done
}

if [ "$N_REPLICAS" -gt 1 ]; then
    for r in $(seq 1 "$N_REPLICAS"); do
        echo "============================================================"
        echo "  REPLICA $r / $N_REPLICAS"
        echo "============================================================"
        run_replicated_modules "$r"

        # ---- capture per-replica metadata (Fase E / Sección 6) ----
        echo "[capture_replica_metadata r$r] — $(date '+%H:%M:%S')"
        if python "$SCRIPTS/utils/capture_replica_metadata.py" \
            --campaigns "$CAMPAIGN" \
            --replica-id "$r"; then
            echo ""
        else
            rc=$?
            echo "[capture_replica_metadata r$r] WARNING: exited with code $rc"
            echo ""
            FAILED_MODULES+=("capture_replica_metadata_r$r")
        fi
    done
else
    # Legacy single-run mode (no replica_*/ subdirs)
    run_replicated_modules ""
fi

# ─────────────────────────────────────────────
# PHASE 3: CONSOLIDATION — distributed across modules (only when N > 1)
# Order: 01h consolidates first (no representative needed) → select_representative
# → remaining modules consolidate using representative_per_mol.csv.
# ─────────────────────────────────────────────
run_consolidate_module() {
    local mod="$1"
    local script="$2"
    local config="$3"
    local label="$4"

    echo "[$mod consolidate] $label — $(date '+%H:%M:%S')"
    if python "$SCRIPTS/$script" \
        --config "$CONFIGS/$config" \
        --campaigns "$CAMPAIGN" \
        --consolidate-replicas \
        --n-replicas "$N_REPLICAS"; then
        echo ""
    else
        local rc=$?
        echo "[$mod consolidate] WARNING: exited with code $rc"
        echo ""
        FAILED_MODULES+=("${mod}_consolidate")
    fi
}

if [ "$N_REPLICAS" -gt 1 ]; then
    echo "============================================================"
    echo "  PHASE 3: distributed consolidation across $N_REPLICAS replicas"
    echo "============================================================"

    # Step 1: 01h consolidates first (its consolidated_mmpbsa.csv is the input
    # to select_representative_replicas).
    run_consolidate_module 01h 01h_mmpbsa_analysis.py 01h_mmpbsa_analysis.yaml "MMPBSA Analysis"

    # Step 2: select representative replica per molecule (writes
    # _replica_metadata/representative_per_mol.csv).
    echo "[select_representative_replicas] — $(date '+%H:%M:%S')"
    if python "$SCRIPTS/utils/select_representative_replicas.py" \
        --campaigns "$CAMPAIGN" \
        --n-replicas "$N_REPLICAS"; then
        echo ""
    else
        rc=$?
        echo "[select_representative_replicas] WARNING: exited with code $rc"
        echo ""
        FAILED_MODULES+=("select_representative")
    fi

    # Step 3: remaining modules consolidate (depend on representative_per_mol.csv).
    run_consolidate_module 01c 01c_dock6_run.py 01c_dock6_run.yaml "DOCK6"
    run_consolidate_module 01e 01e_score_collection.py 01e_score_collection.yaml "Score Collection"
    run_consolidate_module 03a 03a_plip_interaction_analysis.py 03a_plip_interaction_analysis.yaml "PLIP"
    run_consolidate_module 04b 04b_footprint_analysis.py 04b_footprint_analysis.yaml "Footprint Analysis"
    run_consolidate_module 01i 01i_trajectory_analysis.py 01i_trajectory_analysis.yaml "Trajectory Analysis"

    echo "[CONSOLIDATION] Phase 3 done"
    echo ""
fi

# ─────────────────────────────────────────────
# PHASE 4: POST-REPLICA modules (06a, 06b, 07a, 07b, 07c, 07d)
# 06a/06b/07a are repointed to consume <module>/consolidated/ when N>1.
# 07b/07c/07d still receive --n-replicas; full repointing for them is
# scope of Fase D (next instructivo).
# ─────────────────────────────────────────────
run_post_replica_module() {
    local mod="$1"
    local script="$2"
    local config="$3"
    local label="$4"

    if [ "$N_REPLICAS" -gt 1 ]; then
        echo "[$mod] $label (replica-aware, N=$N_REPLICAS) — $(date '+%H:%M:%S')"
        if python "$SCRIPTS/$script" \
            --config "$CONFIGS/$config" \
            --campaigns "$CAMPAIGN" \
            --n-replicas "$N_REPLICAS"; then
            echo ""
        else
            local rc=$?
            echo "[$mod] WARNING: $label exited with code $rc"
            echo ""
            FAILED_MODULES+=("$mod")
        fi
    else
        run_module "$mod" "$script" "$config" "$label"
    fi
}

if should_run 06a; then
    run_post_replica_module 06a 06a_pharmit_pharmacophore.py 06a_pharmit_pharmacophore.yaml "Pharmit Pharmacophore"
fi
if should_run 06b; then
    run_post_replica_module 06b 06b_pharmit_zone_selector.py 06b_pharmit_zone_selector.yaml "Pharmit Zone Selector"
fi
# 07c must run before 07a — it emits unified_verdict.csv consumed by 07a.
if should_run 07c; then
    run_post_replica_module 07c 07c_integrated_analysis.py 07c_integrated_analysis.yaml "Integrated Analysis"
fi
if should_run 07a; then
    run_post_replica_module 07a 07a_decision_report.py 07a_decision_report.yaml "Decision Report"
fi
if should_run 07b; then
    run_post_replica_module 07b 07b_residue_comparison.py 07b_residue_comparison.yaml "Residue Comparison"
fi
if should_run 07d; then
    run_post_replica_module 07d 07d_binding_mode_analysis.py 07d_binding_mode_analysis.yaml "Binding Mode Analysis"
fi

# ─────────────────────────────────────────────
# PHASE 5: campaign manifest (Fase E / Sección 7)
# Persists the audit trail (configs + per-replica metadata + verdicts).
# ─────────────────────────────────────────────
echo "[build_campaign_manifest] — $(date '+%H:%M:%S')"
if python "$SCRIPTS/utils/build_campaign_manifest.py" \
    --campaigns "$CAMPAIGN" \
    --n-replicas "$N_REPLICAS"; then
    echo ""
else
    rc=$?
    echo "[build_campaign_manifest] WARNING: exited with code $rc"
    echo ""
    FAILED_MODULES+=("build_campaign_manifest")
fi

echo "============================================================"
if [ ${#FAILED_MODULES[@]} -eq 0 ]; then
    echo "  PIPELINE COMPLETE — all modules OK"
else
    echo "  PIPELINE COMPLETE — with warnings"
    echo "  Modules with errors: ${FAILED_MODULES[*]}"
fi
echo "  Finished: $(date '+%Y-%m-%d %H:%M:%S')"
echo "  Results:  05_results/"
echo "============================================================"
