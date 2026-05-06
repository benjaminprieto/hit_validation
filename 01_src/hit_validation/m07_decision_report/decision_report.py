"""
Decision Report - Core Module (07a)
=====================================
Generates consolidated validation evidence reports for screening hits.

For each candidate molecule, produces:
  - Multi-criteria pose selection (composite scoring)
  - Scoring comparison (with optional reference context)
  - Sub-pocket zone coverage analysis (absolute energy thresholds)
  - PLIP interaction summary
  - Auto-generated recommendation (template/purchase/weak)
  - Selected pose mol2 export for downstream pipelines

Input:
  - 01e_score_collection/dock6_scores.csv (and dock6_all_poses.csv)
  - 01c_dock6_run/{name}/{name}_scored.mol2 (multi-pose mol2 files)
  - 03a_plip_analysis/{name}/interactions.json
  - 04b_footprint_analysis/subpocket_coverage.csv
  - Optional reference context data (configurable)

Output:
  - decision_report.html         (one-page-per-molecule report)
  - decision_summary.csv         (machine-readable ranking)
  - pose_selection_summary.csv   (per-molecule pose selection details)
  - selected_poses/{name}.mol2   (best pose per molecule)

Location: 01_src/hit_validation/m07_decision_report/decision_report.py
Project: hit_validation
Module: 07a (core)
Version: 2.0 (2026-03-29)
"""

import json
import logging
import re
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import numpy as np
import pandas as pd

from hit_validation.m07_decision_report._footprint_loader import load_footprint_csv

logger = logging.getLogger(__name__)


# =============================================================================
# ZONE DEFINITIONS — loaded from campaign_config at runtime
# =============================================================================
# No hardcoded zones. Zones are passed via the `zones` parameter.
# If no zones are provided, zone-based analysis is skipped.

ZONE_DEFINITIONS = {}  # empty default — populated from campaign_config


# =============================================================================
# REFERENCE DATA LOADING (generic — works with any reference ligand)
# =============================================================================

def load_reference_scores(ref_scores_path: Optional[str]) -> Dict[str, float]:
    """Load reference scores from a prior pipeline run (e.g., reference_docking)."""
    if not ref_scores_path or not Path(ref_scores_path).exists():
        return {}
    try:
        df = pd.read_csv(ref_scores_path)
        if len(df) == 0:
            return {}
        row = df.iloc[0]
        return {
            "Grid_Score": float(row.get("Grid_Score", 0)),
            "Grid_vdw_energy": float(row.get("Grid_vdw_energy", 0)),
            "Grid_es_energy": float(row.get("Grid_es_energy", 0)),
        }
    except Exception as e:
        logger.warning(f"  Could not load reference scores: {e}")
        return {}


def load_reference_plip(ref_plip_path: Optional[str]) -> Dict[str, Any]:
    """Load reference PLIP interactions from a prior pipeline run."""
    if not ref_plip_path or not Path(ref_plip_path).exists():
        return {}
    try:
        with open(ref_plip_path) as f:
            return json.load(f)
    except Exception as e:
        logger.warning(f"  Could not load reference PLIP: {e}")
        return {}


def load_reference_footprint(ref_footprint_path: Optional[str]) -> pd.DataFrame:
    """Load reference footprint consensus from a prior pipeline run."""
    if not ref_footprint_path or not Path(ref_footprint_path).exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(ref_footprint_path)
    except Exception as e:
        logger.warning(f"  Could not load reference footprint: {e}")
        return pd.DataFrame()


# =============================================================================
# POSE SELECTION — multi-criteria composite scoring
# =============================================================================

def _normalize_min_max(values: pd.Series) -> pd.Series:
    """Normalize to [0, 1] where more negative = better = 1.0."""
    vmin, vmax = values.min(), values.max()
    if vmin == vmax:
        return pd.Series(1.0, index=values.index)
    return (vmax - values) / (vmax - vmin)


def select_best_pose(
        mol_poses: pd.DataFrame,
        weights: Optional[Dict[str, float]] = None,
) -> Dict[str, Any]:
    """
    Select best pose from multi-pose data using composite scoring.

    Args:
        mol_poses: DataFrame with one row per pose, columns include score fields
        weights: Dict mapping score column → weight (default: equal weight grid+gbsa)

    Returns:
        Dict with best_pose_idx, composite_score, per-criterion scores, top5 poses
    """
    if weights is None:
        weights = {"Grid_Score": 0.35, "GBSA_Score": 0.35,
                   "FPS_es_energy": 0.15, "FPS_vdw_energy": 0.15}

    if len(mol_poses) == 0:
        return {"best_pose_idx": 0, "composite_score": 0, "top5": []}

    if len(mol_poses) == 1:
        row = mol_poses.iloc[0]
        return {
            "best_pose_idx": int(row.get("Pose_Index", 0)),
            "composite_score": 1.0,
            "top5": [{"pose_idx": int(row.get("Pose_Index", 0)),
                      "composite": 1.0}],
        }

    # Filter to available columns and re-normalize weights
    available = {k: v for k, v in weights.items() if k in mol_poses.columns}
    if not available:
        # Fallback: use Grid_Score if available
        if "Grid_Score" in mol_poses.columns:
            available = {"Grid_Score": 1.0}
        else:
            best_idx = 0
            return {"best_pose_idx": int(mol_poses.iloc[best_idx].get("Pose_Index", 0)),
                    "composite_score": 0, "top5": []}

    total_w = sum(available.values())
    norm_weights = {k: v / total_w for k, v in available.items()}

    # Normalize each criterion and compute composite
    composite = pd.Series(0.0, index=mol_poses.index)
    criterion_scores = {}
    for col, w in norm_weights.items():
        normed = _normalize_min_max(mol_poses[col].astype(float))
        criterion_scores[col] = normed
        composite += w * normed

    best_local_idx = composite.idxmax()
    best_row = mol_poses.loc[best_local_idx]

    # Top 5 poses
    top5_indices = composite.nlargest(min(5, len(composite))).index
    top5 = []
    for idx in top5_indices:
        pose_info = {
            "pose_idx": int(mol_poses.loc[idx].get("Pose_Index", idx)),
            "composite": round(float(composite.loc[idx]), 4),
        }
        for col in available:
            pose_info[col] = round(float(mol_poses.loc[idx][col]), 2)
        top5.append(pose_info)

    return {
        "best_pose_idx": int(best_row.get("Pose_Index", 0)),
        "composite_score": round(float(composite.loc[best_local_idx]), 4),
        "top5": top5,
    }


# =============================================================================
# MOL2 EXTRACTION — extract single pose from multi-pose mol2
# =============================================================================

def extract_pose_mol2(scored_mol2_path: str, pose_index: int) -> Optional[str]:
    """
    Extract a single pose from a DOCK6 multi-pose scored mol2 file.

    Returns the mol2 text for the specified pose (0-indexed), or None if not found.
    """
    path = Path(scored_mol2_path)
    if not path.exists():
        return None

    try:
        with open(path, "r") as f:
            content = f.read()
    except Exception:
        return None

    # Split by @<TRIPOS>MOLECULE markers
    molecules = re.split(r'(?=@<TRIPOS>MOLECULE)', content)
    molecules = [m for m in molecules if m.strip()]

    if pose_index < len(molecules):
        return molecules[pose_index]
    elif molecules:
        return molecules[0]  # fallback to first pose
    return None


# =============================================================================
# DECISION LOGIC
# =============================================================================

def classify_molecule(
        scores: Dict[str, float],
        zone_coverage: Dict[str, bool],
        zone_energies: Dict[str, float],
        zone_energy_cutoff: float = -0.5,
        min_zones_for_template: int = 2,
        grid_score_purchase_threshold: float = -30.0,
        zones: Optional[Dict[str, Dict]] = None,
) -> Dict[str, Any]:
    """
    Auto-classify a molecule based on scoring and zone coverage.

    Uses absolute energy thresholds (not relative to a reference).
    Zone definitions come from campaign_config (no hardcoded zones).

    Returns:
        Dict with recommendation, reasoning, covered_zones
    """
    grid_score = scores.get("Grid_Score", 0)

    # Count zones with meaningful energy contribution (absolute threshold)
    zones_covered = []
    for zone_id, covered in zone_coverage.items():
        if not covered:
            continue
        hit_energy = zone_energies.get(zone_id, 0)
        if hit_energy < zone_energy_cutoff:
            zones_covered.append(zone_id)

    n_zones = len(zones_covered)

    # Get zone labels for readable output
    active_zones = zones or ZONE_DEFINITIONS
    zone_labels = [active_zones.get(z, {}).get("label", z) for z in zones_covered]

    if n_zones >= min_zones_for_template:
        recommendation = "Strong Pharmit template candidate"
        reasoning = (f"Covers {n_zones} zones with meaningful energy: "
                     f"{', '.join(zone_labels)}")
    elif grid_score <= grid_score_purchase_threshold and n_zones >= 1:
        recommendation = "Purchase candidate"
        reasoning = (f"Grid Score {grid_score:.1f} <= {grid_score_purchase_threshold} "
                     f"and covers {n_zones} zone(s)")
    else:
        recommendation = "Weak candidate"
        reasoning = f"Grid Score {grid_score:.1f}, zones covered: {n_zones}"

    return {
        "recommendation": recommendation,
        "reasoning": reasoning,
        "zones_covered": zones_covered,
        "n_zones_covered": n_zones,
        "grid_score": grid_score,
    }


# =============================================================================
# HTML REPORT GENERATION
# =============================================================================

def generate_decision_html(
        molecules: List[Dict[str, Any]],
        campaign_id: str = "",
        ref_scores: Optional[Dict[str, float]] = None,
        reference_label: str = "Reference",
        zones: Optional[Dict[str, Dict]] = None,
        replica_aware: bool = False,
        n_replicas: int = 1,
) -> str:
    """Generate a one-page-per-molecule HTML decision report."""
    ref_scores = ref_scores or {}
    active_zones = zones or ZONE_DEFINITIONS

    html = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8">
<title>Decision Report: {campaign_id}</title>
<style>
body{{font-family:'Segoe UI',Arial,sans-serif;max-width:950px;margin:0 auto;padding:20px;background:#fafafa;color:#333}}
h1{{color:#1a5276;border-bottom:3px solid #1a5276;padding-bottom:10px}}
h2{{color:#2c3e50;margin-top:30px;border-bottom:1px solid #ddd;padding-bottom:5px}}
h3{{color:#34495e;margin-top:20px}}
.mol-card{{border:1px solid #ddd;border-radius:8px;padding:20px;margin:20px 0;background:white;page-break-inside:avoid}}
.strong{{border-left:5px solid #27ae60}}
.purchase{{border-left:5px solid #f39c12}}
.weak{{border-left:5px solid #95a5a6}}
.discordant{{border-left:5px solid #e74c3c}}
.rec-discordant{{background:#f8d7da;color:#721c24}}
.badge-replica{{display:inline-block;padding:3px 8px;border-radius:10px;font-size:11px;font-weight:bold;margin-left:6px}}
.badge-consistent{{background:#d4edda;color:#155724}}
.badge-discordant{{background:#f8d7da;color:#721c24}}
table{{width:100%;border-collapse:collapse;margin:10px 0;font-size:13px}}
th{{background:#f8f9fa;padding:8px;text-align:left;border-bottom:2px solid #dee2e6}}
td{{padding:6px 8px;border-bottom:1px solid #f0f0f0}}
.mono{{font-family:monospace}}
.rec{{font-size:14px;font-weight:600;padding:8px 12px;border-radius:4px;display:inline-block;margin:5px 0}}
.rec-strong{{background:#d4edda;color:#155724}}
.rec-purchase{{background:#fff3cd;color:#856404}}
.rec-weak{{background:#f8d7da;color:#721c24}}
.zone-yes{{color:#27ae60;font-weight:bold}}
.zone-no{{color:#e74c3c}}
.meta{{font-size:12px;color:#777}}
.highlight{{background:#e8f5e9}}
.footer{{margin-top:40px;padding-top:10px;border-top:1px solid #ddd;font-size:11px;color:#999}}
</style></head><body>

<h1>Hit Validation Decision Report</h1>
<p class="meta">Campaign: {campaign_id} | Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')} | Molecules: {len(molecules)}{f' | Replicas: {n_replicas}' if replica_aware else ''}</p>
"""
    if replica_aware:
        html += ('<p class="meta">Mode: <b>replica-aware</b>. Scores are reported '
                 'as mean &plusmn; std across replicas. Pose-consistency verdict '
                 'is produced by 07c (integrated analysis), not here.</p>')

    # Summary table
    html += """<h2>Summary Ranking</h2>
<table>
<tr><th>Rank</th><th>Name</th><th>Grid Score</th><th>Composite</th><th>Pose</th><th>Zones</th><th>Recommendation</th></tr>
"""
    if ref_scores:
        html += f"""<tr style="background:#f0f0f0;font-style:italic">
<td>&mdash;</td>
<td>{reference_label} (ref.)</td>
<td class="mono">{ref_scores.get('Grid_Score', 0):.1f}</td>
<td>&mdash;</td>
<td>&mdash;</td>
<td>&mdash;</td>
<td><span class="meta">reference</span></td>
</tr>
"""
    for i, mol in enumerate(molecules, 1):
        rec = mol.get("recommendation", "")
        rec_class = "rec-strong" if "Strong" in rec else ("rec-purchase" if "Purchase" in rec else "rec-weak")
        html += f"""<tr>
<td>{i}</td>
<td><b>{mol['name']}</b></td>
<td class="mono">{mol.get('grid_score', 0):.1f}</td>
<td class="mono">{mol.get('composite_score', 0):.3f}</td>
<td>{mol.get('selected_pose_idx', 0)}</td>
<td>{mol.get('n_zones_covered', 0)}</td>
<td><span class="rec {rec_class}">{rec}</span></td>
</tr>
"""
    html += "</table>"

    # Per-molecule cards
    for mol in molecules:
        name = mol["name"]
        rec = mol.get("recommendation", "Weak candidate")
        if "Discordant" in rec:
            card_class = "discordant"
            rec_class = "rec-discordant"
        elif "Strong" in rec:
            card_class = "strong"
            rec_class = "rec-strong"
        elif "Purchase" in rec:
            card_class = "purchase"
            rec_class = "rec-purchase"
        else:
            card_class = "weak"
            rec_class = "rec-weak"

        # Pose-consistency badge (consistent/discordant) is determined by 07c
        # in Fase D and is not rendered here.
        replica_badge = ""

        html += f"""
<div class="mol-card {card_class}">
<h3>{name}{replica_badge}</h3>
<div class="rec {rec_class}">{rec}</div>
<p class="meta">{mol.get('reasoning', '')}</p>
"""

        # Pose Selection Card
        top5 = mol.get("top5_poses", [])
        if top5:
            html += """<h4>Pose Selection (top 5)</h4>
<table><tr><th>Pose</th><th>Composite</th>"""
            # Add column headers for available scores
            score_cols = [k for k in top5[0].keys() if k not in ("pose_idx", "composite")]
            for col in score_cols:
                html += f"<th>{col}</th>"
            html += "</tr>"
            for pose in top5:
                is_selected = pose["pose_idx"] == mol.get("selected_pose_idx", -1)
                row_class = ' class="highlight"' if is_selected else ""
                html += f"<tr{row_class}><td>{'<b>' if is_selected else ''}{pose['pose_idx']}{'</b>' if is_selected else ''}</td>"
                html += f"<td class='mono'>{pose['composite']:.3f}</td>"
                for col in score_cols:
                    html += f"<td class='mono'>{pose.get(col, 0):.2f}</td>"
                html += "</tr>"
            html += "</table>"

        # Scoring Summary
        replica_stats = mol.get("replica_stats", {}) if replica_aware else {}
        html += """<h4>Scoring Summary</h4>
<table>
<tr><th>Score</th><th>Hit Value</th>"""
        if replica_aware:
            html += "<th>Std</th>"
        if ref_scores:
            html += f"<th>{reference_label}</th><th>Delta</th>"
        html += "</tr>"

        for score_key in ["Grid_Score", "Grid_vdw_energy", "Grid_es_energy"]:
            hit_val = mol.get("scores", {}).get(score_key, 0)
            html += f"""<tr>
<td>{score_key}</td>
<td class="mono">{hit_val:.2f}</td>"""
            if replica_aware:
                std_key = f"{score_key}_std"
                if std_key in replica_stats:
                    html += f'<td class="mono">&plusmn;{replica_stats[std_key]:.2f}</td>'
                else:
                    html += '<td class="mono">&mdash;</td>'
            if ref_scores:
                ref_val = ref_scores.get(score_key, 0)
                delta = hit_val - ref_val
                html += f"""<td class="mono">{ref_val:.2f}</td>
<td class="mono">{delta:+.2f}</td>"""
            html += "</tr>"

        # MMPBSA row
        mmpbsa_dg = mol.get("mmpbsa_dg")
        if mmpbsa_dg is not None:
            html += f"""<tr style="border-top:2px solid #dee2e6">
<td><b>MMPBSA dG</b></td>
<td class="mono"><b>{mmpbsa_dg:.2f}</b></td>"""
            if replica_aware:
                std = replica_stats.get("delta_total_std")
                if std is not None:
                    html += f'<td class="mono"><b>&plusmn;{std:.2f}</b></td>'
                else:
                    html += '<td class="mono">&mdash;</td>'
            if ref_scores:
                html += '<td class="mono">&mdash;</td><td class="mono">&mdash;</td>'
            html += "</tr>"
        html += "</table>"

        # Zone coverage
        if active_zones:
            html += """<h4>Sub-pocket Coverage</h4>
<table><tr><th>Zone</th><th>Covered</th><th>Hit Energy</th>"""
            if mol.get("ref_zone_energies"):
                html += f"<th>{reference_label}</th>"
            html += "</tr>"

            for zone_id, zdef in active_zones.items():
                covered = mol.get("zone_coverage", {}).get(zone_id, False)
                hit_e = mol.get("zone_energies", {}).get(zone_id, 0)
                cov_str = '<span class="zone-yes">YES</span>' if covered else '<span class="zone-no">no</span>'
                html += f"""<tr>
<td>{zdef.get('label', zone_id)}</td>
<td>{cov_str}</td>
<td class="mono">{hit_e:.2f}</td>"""
                if mol.get("ref_zone_energies"):
                    ref_e = mol["ref_zone_energies"].get(zone_id, 0)
                    html += f'<td class="mono">{ref_e:.2f}</td>'
                html += "</tr>"
            html += "</table>"

        # PLIP interactions
        plip_ints = mol.get("plip_interactions", [])
        if plip_ints:
            html += """<h4>PLIP Interactions</h4>
<table><tr><th>Type</th><th>Residue</th><th>Distance</th></tr>
"""
            for ix in plip_ints[:15]:
                html += f"""<tr>
<td>{ix.get('interaction_type', '')}</td>
<td>{ix.get('residue', '')}</td>
<td class="mono">{ix.get('distance', 0):.1f} A</td>
</tr>
"""
            html += "</table>"

        html += "</div>"

    html += f"""
<div class="footer">
hit_validation | Module 07a | Decision Report v2.0 | {datetime.now().strftime('%Y-%m-%d %H:%M')}
</div>
</body></html>"""

    return html


# =============================================================================
# MAIN PIPELINE
# =============================================================================

def run_decision_report(
        scores_csv: Union[str, Path],
        output_dir: Union[str, Path],
        plip_dir: Optional[str] = None,
        footprint_dir: Optional[str] = None,
        reference_scores_path: Optional[str] = None,
        reference_plip_path: Optional[str] = None,
        reference_footprint_path: Optional[str] = None,
        zone_energy_cutoff: float = -0.5,
        min_zones_for_template: int = 2,
        grid_score_purchase_threshold: float = -30.0,
        campaign_id: str = "",
        reference_label: str = "Reference",
        mmpbsa_analysis_dir: Optional[str] = None,
        mmpbsa_dg_strong_threshold: float = -20.0,
        mmpbsa_dg_moderate_threshold: float = -10.0,
        zones: Optional[Dict[str, Dict]] = None,
        docking_dir: Optional[str] = None,
        all_poses_csv: Optional[str] = None,
        pose_selection_weights: Optional[Dict[str, float]] = None,
        replica_metadata_dir: Optional[Union[str, Path]] = None,
        n_replicas: int = 1,
        score_std_threshold: float = 5.0,
        unified_verdict_csv: Optional[Union[str, Path]] = None,
) -> Dict[str, Any]:
    """
    Generate decision report with multi-criteria pose selection.

    Args:
        scores_csv: Path to 01e dock6_scores.csv. When n_replicas > 1, this
                    must point to the consolidated CSV (has _mean/_std/_sem cols).
        output_dir: Output directory
        plip_dir: Path to 03a_plip_analysis/ (or its consolidated/ when N>1).
        footprint_dir: Path to 04b_footprint_analysis/ (or consolidated/).
        reference_scores_path: Path to reference scores CSV (optional)
        reference_plip_path: Path to reference PLIP JSON (optional)
        reference_footprint_path: Path to reference footprint CSV (optional)
        zone_energy_cutoff: Absolute energy cutoff for zone coverage (kcal/mol)
        min_zones_for_template: Min zones for "strong template" recommendation
        grid_score_purchase_threshold: Grid Score threshold for purchase
        campaign_id: Campaign identifier
        reference_label: Display label for reference data
        zones: Binding site zone definitions from campaign_config
        docking_dir: Path to 01c_dock6_run/ (for mol2 extraction)
        all_poses_csv: Path to dock6_all_poses.csv (for multi-pose selection)
        pose_selection_weights: Weights for composite pose scoring
        replica_metadata_dir: Path to <campaign>/_replica_metadata/. When set,
                              representative_per_mol.csv is read for display.
                              The actual mean/std data already lives inside
                              the consolidated scores_csv and mmpbsa CSV.
        n_replicas: Number of replicas (default 1 = legacy mode).
        score_std_threshold: Grid_Score std threshold (kept for future Fase D
                              discordant-verdict logic; currently informational).

    Returns:
        Dict with success, output paths
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    active_zones = zones or ZONE_DEFINITIONS

    replica_aware = n_replicas > 1 and replica_metadata_dir is not None
    representative_df = None
    if replica_aware:
        rep_path = Path(replica_metadata_dir) / "representative_per_mol.csv"
        if rep_path.exists():
            representative_df = pd.read_csv(rep_path)
        else:
            logger.warning(f"  representative_per_mol.csv not found at {rep_path}")
            replica_aware = False

    logger.info("=" * 60)
    logger.info("  07a Decision Report v2.0"
                + (f" (replica-aware, N={n_replicas})" if replica_aware else ""))
    logger.info("=" * 60)

    # --- Load hit scores ---
    if not Path(scores_csv).exists():
        return {"success": False, "error": f"Scores CSV not found: {scores_csv}"}

    df_scores = pd.read_csv(scores_csv)
    logger.info(f"  Loaded {len(df_scores)} molecules from scores CSV")

    # --- Load unified verdict from 07c (single source of truth) ---
    unified_verdicts: Dict[str, Dict[str, Any]] = {}
    if unified_verdict_csv and Path(unified_verdict_csv).exists():
        try:
            df_uv = pd.read_csv(unified_verdict_csv)
            for _, row in df_uv.iterrows():
                name = row.get("Name")
                if pd.isna(name):
                    continue
                unified_verdicts[str(name)] = {
                    k: (None if pd.isna(v) else v) for k, v in row.to_dict().items()
                }
            logger.info(
                f"  Loaded unified verdicts for {len(unified_verdicts)} molecules "
                f"from {unified_verdict_csv}"
            )
        except Exception as e:
            logger.warning(f"  Could not load unified_verdict.csv: {e}")

    # --- Load all poses for multi-criteria selection ---
    df_all_poses = None
    if all_poses_csv and Path(all_poses_csv).exists():
        df_all_poses = pd.read_csv(all_poses_csv)
        logger.info(f"  Loaded {len(df_all_poses)} poses from all_poses CSV")

    # --- Load reference data ---
    ref_scores = load_reference_scores(reference_scores_path)
    ref_plip = load_reference_plip(reference_plip_path)
    ref_footprint = load_reference_footprint(reference_footprint_path)

    # Build reference zone energies from footprint
    ref_zone_energies = {}
    if not ref_footprint.empty and "residue_id" in ref_footprint.columns and active_zones:
        for zone_id, zdef in active_zones.items():
            zone_rows = ref_footprint[
                ref_footprint["residue_id"].apply(
                    lambda rid: str(rid).split(".")[0] in zdef["residues"]
                )
            ]
            if not zone_rows.empty and "mean_total" in zone_rows.columns:
                ref_zone_energies[zone_id] = zone_rows["mean_total"].sum()
            else:
                ref_zone_energies[zone_id] = 0.0

    # --- Load MMPBSA data (from 01h) ---
    mmpbsa_data = {}
    if mmpbsa_analysis_dir and Path(mmpbsa_analysis_dir).exists():
        consolidated_csv = Path(mmpbsa_analysis_dir) / "consolidated_mmpbsa.csv"
        if consolidated_csv.exists():
            try:
                df_mmpbsa = pd.read_csv(consolidated_csv)
                for _, row in df_mmpbsa.iterrows():
                    # In replica mode, _mean columns are present; in legacy
                    # mode, raw columns are present. Prefer _mean when available.
                    def _val(col):
                        mc = f"{col}_mean"
                        if mc in row and pd.notna(row[mc]):
                            return row[mc]
                        return row.get(col)
                    entry = {
                        "dG_total": _val("MMPBSA_dG_total"),
                        "vdW": _val("MMPBSA_vdW"),
                        "EEL": _val("MMPBSA_EEL"),
                        "EGB": _val("MMPBSA_EGB"),
                        "ESURF": _val("MMPBSA_ESURF"),
                    }
                    # Capture mean/std/sem for HTML/CSV in replica mode
                    for col in ("MMPBSA_dG_total", "delta_total"):
                        for suffix in ("_mean", "_std", "_sem"):
                            key = f"{col}{suffix}"
                            if key in row and pd.notna(row[key]):
                                entry[key] = float(row[key])
                    mmpbsa_data[row["Name"]] = entry
                logger.info(f"  Loaded MMPBSA data for {len(mmpbsa_data)} molecules")
            except Exception as e:
                logger.warning(f"  Could not load MMPBSA data: {e}")

        for mol_dir in Path(mmpbsa_analysis_dir).iterdir():
            if not mol_dir.is_dir():
                continue
            zone_json = mol_dir / "zone_summary.json"
            if zone_json.exists():
                try:
                    with open(zone_json) as f:
                        zone_s = json.load(f)
                    if mol_dir.name in mmpbsa_data:
                        mmpbsa_data[mol_dir.name]["zone_summary"] = zone_s
                    else:
                        mmpbsa_data[mol_dir.name] = {"zone_summary": zone_s}
                except Exception:
                    pass

    # --- Load hit footprint data ---
    hit_footprint = pd.DataFrame()
    if footprint_dir:
        fp_csv = Path(footprint_dir) / "footprint_per_molecule.csv"
        if fp_csv.exists():
            hit_footprint = load_footprint_csv(fp_csv)

    # --- Pose selection + mol2 export setup ---
    selected_poses_dir = output_dir / "selected_poses"
    selected_poses_dir.mkdir(parents=True, exist_ok=True)
    pose_selection_rows = []

    # --- Build per-molecule report data ---
    molecules = []

    for _, row in df_scores.iterrows():
        name = row["Name"]
        scores = {col: float(row[col]) for col in row.index
                  if isinstance(row[col], (int, float)) and col not in ("Rank", "n_poses")}

        # --- Multi-criteria pose selection ---
        pose_result = {"best_pose_idx": 0, "composite_score": 0, "top5": []}
        if df_all_poses is not None:
            mol_poses = df_all_poses[df_all_poses["Name"] == name].copy()
            if not mol_poses.empty:
                # Add Pose_Index if not present
                if "Pose_Index" not in mol_poses.columns:
                    mol_poses = mol_poses.reset_index(drop=True)
                    mol_poses["Pose_Index"] = range(len(mol_poses))
                pose_result = select_best_pose(mol_poses, pose_selection_weights)

                # Update scores to use selected pose's scores
                best_pose_rows = mol_poses[
                    mol_poses["Pose_Index"] == pose_result["best_pose_idx"]
                ]
                if not best_pose_rows.empty:
                    bp = best_pose_rows.iloc[0]
                    for col in bp.index:
                        if isinstance(bp[col], (int, float)) and col not in ("Rank", "n_poses", "Pose_Index"):
                            scores[col] = float(bp[col])

        selected_pose_idx = pose_result["best_pose_idx"]

        # --- Export selected pose mol2 ---
        if docking_dir:
            scored_mol2 = Path(docking_dir) / name / f"{name}_scored.mol2"
            pose_text = extract_pose_mol2(str(scored_mol2), selected_pose_idx)
            if pose_text:
                mol2_out = selected_poses_dir / f"{name}.mol2"
                mol2_out.write_text(pose_text, encoding="utf-8")

        # Zone coverage from footprint
        zone_coverage = {}
        zone_energies = {}
        if not hit_footprint.empty and active_zones:
            mol_fp = hit_footprint[hit_footprint["Name"] == name]
            for zone_id, zdef in active_zones.items():
                zone_rows = mol_fp[
                    mol_fp["residue_id"].apply(
                        lambda rid: str(rid).split(".")[0] in zdef["residues"]
                    )
                ]
                zone_total = zone_rows["total"].sum() if not zone_rows.empty else 0
                zone_energies[zone_id] = zone_total
                zone_coverage[zone_id] = zone_total < -0.5
        else:
            for zone_id in active_zones:
                zone_coverage[zone_id] = False
                zone_energies[zone_id] = 0.0

        # PLIP interactions
        plip_interactions = []
        if plip_dir:
            plip_json = Path(plip_dir) / name / "interactions.json"
            if plip_json.exists():
                try:
                    with open(plip_json) as f:
                        plip_data = json.load(f)
                    plip_interactions = plip_data.get("interactions", [])
                except Exception:
                    pass

        # MMPBSA data
        mol_mmpbsa = mmpbsa_data.get(name, {})
        mmpbsa_dg = mol_mmpbsa.get("dG_total")

        # --- Replica-aware overrides ---
        # In replica mode the consolidated dock6_scores.csv has Grid_Score_mean/std/sem
        # columns (no plain "Grid_Score"); copy mean into Grid_Score for downstream
        # logic and capture mean/std/sem stats for HTML/CSV display.
        replica_stats: Dict[str, Any] = {}
        if replica_aware:
            score_row = df_scores[df_scores["Name"] == name]
            if not score_row.empty:
                sr = score_row.iloc[0]
                if "Grid_Score_mean" in sr and pd.notna(sr["Grid_Score_mean"]):
                    scores["Grid_Score"] = float(sr["Grid_Score_mean"])
                for col in ("Grid_Score_mean", "Grid_Score_std", "Grid_Score_sem",
                            "Grid_vdw_energy_mean", "Grid_vdw_energy_std",
                            "Grid_es_energy_mean", "Grid_es_energy_std",
                            "GBSA_Score_mean", "GBSA_Score_std", "GBSA_Score_sem",
                            "Hawkins_GBSA_Score_mean", "Hawkins_GBSA_Score_std"):
                    if col in sr and pd.notna(sr[col]):
                        replica_stats[col] = float(sr[col])
            # MMPBSA mean/std came from consolidated_mmpbsa.csv (loaded above).
            for col in ("MMPBSA_dG_total_mean", "MMPBSA_dG_total_std",
                        "MMPBSA_dG_total_sem", "delta_total_mean",
                        "delta_total_std", "delta_total_sem"):
                if col in mol_mmpbsa:
                    replica_stats[col] = float(mol_mmpbsa[col])
            # Pose consistency verdict / inter-replica RMSD live in 07c (Fase D).
            # Not computed in 07a anymore.

        # Classify
        decision = classify_molecule(
            scores=scores,
            zone_coverage=zone_coverage,
            zone_energies=zone_energies,
            zone_energy_cutoff=zone_energy_cutoff,
            min_zones_for_template=min_zones_for_template,
            grid_score_purchase_threshold=grid_score_purchase_threshold,
            zones=active_zones,
        )

        # Enhance classification with MMPBSA if available (skipped if 07c verdict will override below)
        if mmpbsa_dg is not None and name not in unified_verdicts:
            n_zones = decision["n_zones_covered"]
            grid_score = scores.get("Grid_Score", 0)
            if mmpbsa_dg <= mmpbsa_dg_strong_threshold and n_zones >= min_zones_for_template:
                decision["recommendation"] = "Strong Pharmit template candidate"
                decision["reasoning"] = (
                    f"MMPBSA dG {mmpbsa_dg:.1f} <= {mmpbsa_dg_strong_threshold} "
                    f"AND covers {n_zones} zones: {', '.join(decision['zones_covered'])}"
                )
            elif mmpbsa_dg <= mmpbsa_dg_moderate_threshold and grid_score <= grid_score_purchase_threshold:
                if decision["recommendation"] == "Weak candidate":
                    decision["recommendation"] = "Purchase candidate"
                    decision["reasoning"] = (
                        f"MMPBSA dG {mmpbsa_dg:.1f} <= {mmpbsa_dg_moderate_threshold} "
                        f"AND Grid Score {grid_score:.1f} <= {grid_score_purchase_threshold}"
                    )

        # Discordant verdict logic moved to 07c (Fase D). 07a stays informational
        # in replica mode: mean +/- std are displayed but no purchase/discordant
        # verdict is forced based on inter-replica consistency. score_std_threshold
        # is still received as a parameter for forward-compat but unused here.

        # If 07c emitted a unified verdict, prefer that as the recommendation.
        uv = unified_verdicts.get(name)
        if uv:
            decision["recommendation"] = uv.get("verdict") or decision["recommendation"]
            decision["reasoning"] = uv.get("reasons") or decision["reasoning"]
            decision["unified_category"] = uv.get("category")
            decision["pose_max_rmsd_inter_replica"] = uv.get("pose_max_rmsd_inter_replica")
            decision["sign_change"] = uv.get("sign_change")

        molecules.append({
            "name": name,
            "scores": scores,
            "grid_score": scores.get("Grid_Score", 0),
            "composite_score": pose_result["composite_score"],
            "selected_pose_idx": selected_pose_idx,
            "top5_poses": pose_result["top5"],
            "zone_coverage": zone_coverage,
            "zone_energies": zone_energies,
            "ref_zone_energies": ref_zone_energies if ref_zone_energies else None,
            "plip_interactions": plip_interactions,
            "n_plip_interactions": len(plip_interactions),
            "mmpbsa_dg": mmpbsa_dg,
            "mmpbsa_data": mol_mmpbsa,
            "replica_stats": replica_stats,
            **decision,
        })

        # Pose selection summary row
        pose_selection_rows.append({
            "Name": name,
            "Selected_Pose_Index": selected_pose_idx,
            "Composite_Score": pose_result["composite_score"],
            "Grid_Score": scores.get("Grid_Score", 0),
            "Recommendation": decision["recommendation"],
            "N_Zones_Covered": decision["n_zones_covered"],
        })

    # Sort by recommendation strength then grid score
    def _rec_rank(rec: str) -> int:
        if "Discordant" in rec:
            return 4
        if "Strong" in rec:
            return 0
        if "Purchase" in rec:
            return 1
        if "Weak" in rec:
            return 2
        return 3
    molecules.sort(key=lambda m: (_rec_rank(m["recommendation"]), m["grid_score"]))

    # --- Generate HTML report ---
    html = generate_decision_html(
        molecules, campaign_id, ref_scores, reference_label, active_zones,
        replica_aware=replica_aware, n_replicas=n_replicas,
    )
    html_path = output_dir / "decision_report.html"
    html_path.write_text(html, encoding="utf-8")
    logger.info(f"  Saved: {html_path}")

    # --- Generate pose selection summary CSV ---
    if pose_selection_rows:
        df_pose = pd.DataFrame(pose_selection_rows)
        pose_csv = output_dir / "pose_selection_summary.csv"
        df_pose.to_csv(pose_csv, index=False, encoding="utf-8")
        logger.info(f"  Saved: {pose_csv}")

    # --- Generate summary CSV ---
    summary_rows = []
    for i, mol in enumerate(molecules, 1):
        row = {
            "Rank": i,
            "Name": mol["name"],
            "Grid_Score": mol["grid_score"],
            "Composite_Score": mol["composite_score"],
            "Selected_Pose_Index": mol["selected_pose_idx"],
            "MMPBSA_dG": mol.get("mmpbsa_dg"),
            "Recommendation": mol["recommendation"],
            "Zones_Covered": mol["n_zones_covered"],
            "Zones_List": ",".join(mol["zones_covered"]),
            "N_PLIP_Interactions": mol["n_plip_interactions"],
            "Reasoning": mol["reasoning"],
        }
        if replica_aware:
            stats = mol.get("replica_stats", {})
            row["Grid_Score_Std"] = stats.get("Grid_Score_std")
            # delta_total_std mirrors MMPBSA_dG_total_std (alias from 01h consolidator).
            row["MMPBSA_dG_Std"] = stats.get(
                "delta_total_std", stats.get("MMPBSA_dG_total_std")
            )
            # Pose RMSD inter-replica + Pose_Consistency from 07c unified verdict.
            pose_rmsd = mol.get("pose_max_rmsd_inter_replica")
            row["Pose_Max_RMSD"] = (
                float(pose_rmsd) if pose_rmsd is not None else None
            )
            cat = mol.get("unified_category")
            row["Pose_Consistency"] = (
                "Discordant" if cat == "discordant" else
                ("OK" if cat is not None else None)
            )
        summary_rows.append(row)

    df_summary = pd.DataFrame(summary_rows)
    summary_csv = output_dir / "decision_summary.csv"
    df_summary.to_csv(summary_csv, index=False, encoding="utf-8")
    logger.info(f"  Saved: {summary_csv}")

    # --- Log summary ---
    n_discordant = sum(1 for m in molecules if "Discordant" in m["recommendation"])
    n_strong = sum(1 for m in molecules
                   if "Strong" in m["recommendation"] and "Discordant" not in m["recommendation"])
    n_purchase = sum(1 for m in molecules
                     if "Purchase" in m["recommendation"] and "Discordant" not in m["recommendation"])
    n_weak = sum(1 for m in molecules
                 if "Weak" in m["recommendation"] and "Discordant" not in m["recommendation"])

    logger.info("")
    logger.info(f"{'=' * 60}")
    logger.info(f"  DECISION REPORT: {len(molecules)} molecules")
    logger.info(f"    Strong template: {n_strong}")
    logger.info(f"    Purchase:        {n_purchase}")
    logger.info(f"    Weak:            {n_weak}")
    if replica_aware:
        logger.info(f"    Discordant:      {n_discordant}")
    logger.info(f"  Selected poses exported to: {selected_poses_dir}")
    logger.info(f"{'=' * 60}")

    return {
        "success": True,
        "n_molecules": len(molecules),
        "n_strong": n_strong,
        "n_purchase": n_purchase,
        "n_weak": n_weak,
        "n_discordant": n_discordant,
        "html_path": str(html_path),
        "summary_csv": str(summary_csv),
        "selected_poses_dir": str(selected_poses_dir),
    }
