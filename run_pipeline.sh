#!/bin/bash
# =============================================================================
# run_pipeline.sh — Hit Validation Pipeline
# =============================================================================
#
# Usage:
#   bash run_pipeline.sh <campaign_config.yaml> [start_from]
#
# Examples:
#   bash run_pipeline.sh 04_data/campaigns/hit_validation/campaign_config.yaml
#   bash run_pipeline.sh 04_data/campaigns/hit_validation/campaign_config.yaml 01c
#
#   # Background (survives logout):
#   nohup bash run_pipeline.sh 04_data/campaigns/hit_validation/campaign_config.yaml > pipeline.log 2>&1 &
#
# Pipeline:
#   00a  Ligand preparation    (antechamber: AM1-BCC + Sybyl types, batch)
#   00b  Receptor preparation  (ChimeraX DockPrep → rec_charged.mol2)
#   00d  Binding site definition (by residues/coords — skipped if precomputed)
#   01b  Grid generation       (DMS → spheres → box 5.0Å → grids)
#   01c  DOCK6 docking         (grid_score primary, flex search)
#   01d  Footprint rescore     (fps_primary, per-residue vdW+ES)
#   01f  GB/SA rescore         (gbsa_hawkins_primary, with minimization)
#   01e  Score collection      (parse all scores → CSV/Excel)
#   03a  PLIP analysis         (docked pose + receptor → interactions)
#   04b  Footprint analysis    (per-residue consensus + zone coverage)
#   07a  Decision report       (pose selection + per-molecule validation evidence)
# =============================================================================
set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: bash run_pipeline.sh <campaign_config.yaml> [start_from]"
    echo ""
    echo "Modules: 00a 00b 00d 01b 01c 01d 01f 01e 03a 04b 07a"
    echo ""
    echo "Examples:"
    echo "  bash run_pipeline.sh config.yaml        # full pipeline"
    echo "  bash run_pipeline.sh config.yaml 01c    # start from docking"
    echo "  bash run_pipeline.sh config.yaml 03a    # start from PLIP"
    exit 1
fi

CAMPAIGN="$1"
START="${2:-00a}"
CONFIGS="03_configs"
SCRIPTS="02_scripts"

if [ ! -f "$CAMPAIGN" ]; then
    echo "ERROR: Campaign config not found: $CAMPAIGN"
    exit 1
fi

# Check if precomputed_grids_dir is set in campaign config
PRECOMPUTED_GRIDS=$(python3 -c "
import yaml, sys
with open('$CAMPAIGN') as f:
    cc = yaml.safe_load(f)
grids = cc.get('grids', {})
d = grids.get('precomputed_grids_dir')
print(d if d else '')
" 2>/dev/null || echo "")

# Module order
MODULES=(00a 00b 00d 01b 01c 01d 01f 01e 03a 04b 07a)

# Find start index
START_IDX=0
for i in "${!MODULES[@]}"; do
    if [ "${MODULES[$i]}" = "$START" ]; then
        START_IDX=$i
        break
    fi
done

echo "============================================================"
echo "  HIT_VALIDATION — Pipeline v1.0"
echo "  Campaign: $CAMPAIGN"
echo "  Start:    $START"
echo "  Started:  $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"
echo ""

run_module() {
    local mod="$1"
    local script="$2"
    local config="$3"
    local label="$4"

    echo "[$mod] $label — $(date '+%H:%M:%S')"
    python "$SCRIPTS/$script" \
        --config "$CONFIGS/$config" \
        --campaigns "$CAMPAIGN"
    echo ""
}

for i in "${!MODULES[@]}"; do
    [ "$i" -lt "$START_IDX" ] && continue

    case "${MODULES[$i]}" in
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
        01c) run_module 01c 01c_dock6_run.py 01c_dock6_run.yaml "DOCK6 Docking" ;;
        01d) run_module 01d 01d_footprint_rescore.py 01d_footprint_rescore.yaml "Footprint Rescore" ;;
        01f) run_module 01f 01f_gbsa_rescore.py 01f_gbsa_rescore.yaml "GB/SA Hawkins Rescore" ;;
        01e) run_module 01e 01e_score_collection.py 01e_score_collection.yaml "Score Collection" ;;
        03a) run_module 03a 03a_plip_interaction_analysis.py 03a_plip_interaction_analysis.yaml "PLIP Analysis" ;;
        04b) run_module 04b 04b_footprint_analysis.py 04b_footprint_analysis.yaml "Footprint Analysis" ;;
        07a) run_module 07a 07a_decision_report.py 07a_decision_report.yaml "Decision Report" ;;
    esac
done

echo "============================================================"
echo "  PIPELINE COMPLETE"
echo "  Finished: $(date '+%Y-%m-%d %H:%M:%S')"
echo "  Results:  05_results/"
echo "============================================================"