"""
Decision Report - Core Module (07a)
=====================================
Generates consolidated validation evidence reports for screening hits.

For each candidate molecule, produces:
  - Scoring comparison (with optional reference context)
  - Sub-pocket coverage analysis (absolute energy thresholds)
  - PLIP interaction summary
  - Auto-generated recommendation (template/purchase/weak)

Input:
  - 01e_score_collection/dock6_scores.csv
  - 03a_plip_analysis/{name}/interactions.json
  - 04b_footprint_analysis/hit_vs_udx_comparison.csv
  - 04b_footprint_analysis/subpocket_coverage.csv
  - Reference data from reference_docking

Output:
  - decision_report.html  (one-page-per-molecule report)
  - decision_summary.csv  (machine-readable ranking)

Location: 01_src/hit_validation/m07_decision_report/decision_report.py
Project: hit_validation
Module: 07a (core)
Version: 1.0 (2026-03-27)
"""

import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import pandas as pd

logger = logging.getLogger(__name__)


# =============================================================================
# ZONE DEFINITIONS (from 04b — same for XT1/UDX system)
# =============================================================================

ZONE_DEFINITIONS = {
    "phosphate": {
        "residues": {"ARG598", "LYS599"},
        "label": "Phosphate (salt bridges)",
    },
    "xylose": {
        "residues": {"TRP392", "TRP495", "TYR565", "SER575"},
        "label": "Xylose pocket",
    },
    "uracil": {
        "residues": {"ASP361", "ARG363"},
        "label": "Uracil pocket",
    },
    "ribose": {
        "residues": {"HIS335", "VAL333", "THR390"},
        "label": "Ribose / base recognition",
    },
    "catalytic": {
        "residues": {"GLU529"},
        "label": "Catalytic (GLU529)",
    },
}

KEY_RESIDUES = ["TRP392", "HIS335", "GLU529", "ASP494", "ARG598"]


# =============================================================================
# REFERENCE DATA LOADING
# =============================================================================

def load_reference_scores(ref_scores_path: Optional[str]) -> Dict[str, float]:
    """Load UDX reference scores from reference_docking output."""
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
    """Load UDX reference PLIP interactions."""
    if not ref_plip_path or not Path(ref_plip_path).exists():
        return {}
    try:
        with open(ref_plip_path) as f:
            return json.load(f)
    except Exception as e:
        logger.warning(f"  Could not load reference PLIP: {e}")
        return {}


def load_reference_footprint(ref_footprint_path: Optional[str]) -> pd.DataFrame:
    """Load UDX reference footprint consensus."""
    if not ref_footprint_path or not Path(ref_footprint_path).exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(ref_footprint_path)
    except Exception as e:
        logger.warning(f"  Could not load reference footprint: {e}")
        return pd.DataFrame()


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
) -> Dict[str, Any]:
    """
    Auto-classify a molecule based on scoring and zone coverage.

    Uses absolute energy thresholds (not relative to a reference).

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
        # Zone is meaningfully covered if energy < -0.5 kcal/mol (absolute)
        if hit_energy < zone_energy_cutoff:
            zones_covered.append(zone_id)

    n_zones = len(zones_covered)

    if n_zones >= min_zones_for_template:
        recommendation = "Strong Pharmit template candidate"
        reasoning = (f"Covers {n_zones} drug-like zones with meaningful energy: "
                     f"{', '.join(zones_covered)}")
    elif grid_score <= grid_score_purchase_threshold and (
            zone_coverage.get("xylose", False) or zone_coverage.get("uracil", False)):
        recommendation = "Purchase candidate"
        reasoning = (f"Grid Score {grid_score:.1f} <= {grid_score_purchase_threshold} "
                     f"and covers xylose or uracil zone")
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
) -> str:
    """Generate a one-page-per-molecule HTML decision report."""
    ref_scores = ref_scores or {}

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
.footer{{margin-top:40px;padding-top:10px;border-top:1px solid #ddd;font-size:11px;color:#999}}
</style></head><body>

<h1>Hit Validation Decision Report</h1>
<p class="meta">Campaign: {campaign_id} | Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')} | Molecules: {len(molecules)}</p>
"""

    # Summary table
    html += """<h2>Summary Ranking</h2>
<table>
<tr><th>Rank</th><th>Name</th><th>Grid Score</th><th>Zones</th><th>Recommendation</th></tr>
"""
    if ref_scores:
        html += f"""<tr style="background:#f0f0f0;font-style:italic">
<td>&mdash;</td>
<td>{reference_label} (ref.)</td>
<td class="mono">{ref_scores.get('Grid_Score', 0):.1f}</td>
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
<td>{mol.get('n_zones_covered', 0)}</td>
<td><span class="rec {rec_class}">{rec}</span></td>
</tr>
"""
    html += "</table>"

    # Per-molecule cards
    for mol in molecules:
        name = mol["name"]
        rec = mol.get("recommendation", "Weak candidate")
        card_class = "strong" if "Strong" in rec else ("purchase" if "Purchase" in rec else "weak")
        rec_class = "rec-strong" if "Strong" in rec else ("rec-purchase" if "Purchase" in rec else "rec-weak")

        html += f"""
<div class="mol-card {card_class}">
<h3>{name}</h3>
<div class="rec {rec_class}">{rec}</div>
<p class="meta">{mol.get('reasoning', '')}</p>

<h4>Scoring Comparison</h4>
<table>
<tr><th>Score</th><th>Hit Value</th><th>Reference</th><th>Delta</th></tr>
"""
        for score_key in ["Grid_Score", "Grid_vdw_energy", "Grid_es_energy"]:
            hit_val = mol.get("scores", {}).get(score_key, 0)
            if ref_scores:
                ref_val = ref_scores.get(score_key, 0)
                delta = hit_val - ref_val
                ref_str = f'{ref_val:.2f}'
                delta_str = f'{delta:+.2f}'
            else:
                ref_str = '&mdash;'
                delta_str = '&mdash;'
            html += f"""<tr>
<td>{score_key}</td>
<td class="mono">{hit_val:.2f}</td>
<td class="mono">{ref_str}</td>
<td class="mono">{delta_str}</td>
</tr>
"""
        # MMPBSA ΔG row (if available)
        mmpbsa_dg = mol.get("mmpbsa_dg")
        if mmpbsa_dg is not None:
            html += f"""<tr style="border-top:2px solid #dee2e6">
<td><b>MMPBSA ΔG</b></td>
<td class="mono"><b>{mmpbsa_dg:.2f}</b></td>
<td class="mono">&mdash;</td>
<td class="mono">&mdash;</td>
</tr>
"""
        html += "</table>"

        # Zone coverage
        html += """<h4>Sub-pocket Coverage</h4>
<table><tr><th>Zone</th><th>Covered</th><th>Hit Energy</th><th>Ref Energy</th></tr>
"""
        for zone_id, zdef in ZONE_DEFINITIONS.items():
            covered = mol.get("zone_coverage", {}).get(zone_id, False)
            hit_e = mol.get("zone_energies", {}).get(zone_id, 0)
            ref_zone_e = mol.get("ref_zone_energies", {})
            if ref_zone_e:
                ref_e_str = f'{ref_zone_e.get(zone_id, 0):.2f}'
            else:
                ref_e_str = '&mdash;'
            cov_str = '<span class="zone-yes">YES</span>' if covered else '<span class="zone-no">no</span>'
            html += f"""<tr>
<td>{zdef['label']}</td>
<td>{cov_str}</td>
<td class="mono">{hit_e:.2f}</td>
<td class="mono">{ref_e_str}</td>
</tr>
"""
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
hit_validation | Module 07a | Decision Report v1.0 | {datetime.now().strftime('%Y-%m-%d %H:%M')}
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
) -> Dict[str, Any]:
    """
    Generate decision report for all validated hits.

    Args:
        scores_csv: Path to 01e dock6_scores.csv
        output_dir: Output directory
        plip_dir: Path to 03a_plip_analysis/ (per-molecule PLIP results)
        footprint_dir: Path to 04b_footprint_analysis/
        reference_scores_path: Path to reference_docking scores CSV
        reference_plip_path: Path to reference_docking PLIP JSON
        reference_footprint_path: Path to reference_docking footprint CSV
        zone_energy_cutoff: Absolute energy cutoff for zone coverage (kcal/mol)
        min_zones_for_template: Min zones for "strong template" recommendation
        grid_score_purchase_threshold: Grid Score threshold for purchase
        campaign_id: Campaign identifier
        reference_label: Display label for reference data in report

    Returns:
        Dict with success, output paths
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("  07a Decision Report")
    logger.info("=" * 60)

    # --- Load hit scores ---
    if not Path(scores_csv).exists():
        return {"success": False, "error": f"Scores CSV not found: {scores_csv}"}

    df_scores = pd.read_csv(scores_csv)
    logger.info(f"  Loaded {len(df_scores)} molecules from scores CSV")

    # --- Load reference data ---
    ref_scores = load_reference_scores(reference_scores_path)
    ref_plip = load_reference_plip(reference_plip_path)
    ref_footprint = load_reference_footprint(reference_footprint_path)

    # Build reference zone energies from footprint
    ref_zone_energies = {}
    if not ref_footprint.empty and "residue_id" in ref_footprint.columns:
        for zone_id, zdef in ZONE_DEFINITIONS.items():
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
    mmpbsa_data = {}  # mol_name → {dG_total, zone_summary, ...}
    if mmpbsa_analysis_dir and Path(mmpbsa_analysis_dir).exists():
        # Load consolidated energies
        consolidated_csv = Path(mmpbsa_analysis_dir) / "consolidated_mmpbsa.csv"
        if consolidated_csv.exists():
            try:
                df_mmpbsa = pd.read_csv(consolidated_csv)
                for _, row in df_mmpbsa.iterrows():
                    mmpbsa_data[row["Name"]] = {
                        "dG_total": row.get("MMPBSA_dG_total"),
                        "vdW": row.get("MMPBSA_vdW"),
                        "EEL": row.get("MMPBSA_EEL"),
                        "EGB": row.get("MMPBSA_EGB"),
                        "ESURF": row.get("MMPBSA_ESURF"),
                    }
                logger.info(f"  Loaded MMPBSA data for {len(mmpbsa_data)} molecules")
            except Exception as e:
                logger.warning(f"  Could not load MMPBSA data: {e}")

        # Load per-molecule zone summaries
        for mol_dir in Path(mmpbsa_analysis_dir).iterdir():
            if not mol_dir.is_dir():
                continue
            zone_json = mol_dir / "zone_summary.json"
            if zone_json.exists():
                try:
                    with open(zone_json) as f:
                        zones = json.load(f)
                    if mol_dir.name in mmpbsa_data:
                        mmpbsa_data[mol_dir.name]["zone_summary"] = zones
                    else:
                        mmpbsa_data[mol_dir.name] = {"zone_summary": zones}
                except Exception:
                    pass

    # --- Load hit footprint data ---
    hit_footprint = pd.DataFrame()
    if footprint_dir:
        fp_csv = Path(footprint_dir) / "footprint_per_molecule.csv"
        if fp_csv.exists():
            hit_footprint = pd.read_csv(fp_csv)

    # --- Build per-molecule report data ---
    molecules = []

    for _, row in df_scores.iterrows():
        name = row["Name"]
        scores = {col: float(row[col]) for col in row.index
                  if isinstance(row[col], (int, float)) and col not in ("Rank", "n_poses")}

        # Zone coverage from footprint
        zone_coverage = {}
        zone_energies = {}
        if not hit_footprint.empty:
            mol_fp = hit_footprint[hit_footprint["Name"] == name]
            for zone_id, zdef in ZONE_DEFINITIONS.items():
                zone_rows = mol_fp[
                    mol_fp["residue_id"].apply(
                        lambda rid: str(rid).split(".")[0] in zdef["residues"]
                    )
                ]
                zone_total = zone_rows["total"].sum() if not zone_rows.empty else 0
                zone_energies[zone_id] = zone_total
                zone_coverage[zone_id] = zone_total < -0.5
        else:
            for zone_id in ZONE_DEFINITIONS:
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

        # MMPBSA data for this molecule
        mol_mmpbsa = mmpbsa_data.get(name, {})
        mmpbsa_dg = mol_mmpbsa.get("dG_total")

        # Classify (with MMPBSA enhancement when available)
        decision = classify_molecule(
            scores=scores,
            zone_coverage=zone_coverage,
            zone_energies=zone_energies,
            zone_energy_cutoff=zone_energy_cutoff,
            min_zones_for_template=min_zones_for_template,
            grid_score_purchase_threshold=grid_score_purchase_threshold,
        )

        # Enhance classification with MMPBSA if available
        if mmpbsa_dg is not None:
            n_zones = decision["n_zones_covered"]
            grid_score = scores.get("Grid_Score", 0)
            if mmpbsa_dg <= mmpbsa_dg_strong_threshold and n_zones >= min_zones_for_template:
                decision["recommendation"] = "Strong Pharmit template candidate"
                decision["reasoning"] = (
                    f"MMPBSA ΔG {mmpbsa_dg:.1f} <= {mmpbsa_dg_strong_threshold} "
                    f"AND covers {n_zones} zones: {', '.join(decision['zones_covered'])}"
                )
            elif mmpbsa_dg <= mmpbsa_dg_moderate_threshold and grid_score <= grid_score_purchase_threshold:
                if decision["recommendation"] == "Weak candidate":
                    decision["recommendation"] = "Purchase candidate"
                    decision["reasoning"] = (
                        f"MMPBSA ΔG {mmpbsa_dg:.1f} <= {mmpbsa_dg_moderate_threshold} "
                        f"AND Grid Score {grid_score:.1f} <= {grid_score_purchase_threshold}"
                    )

        molecules.append({
            "name": name,
            "scores": scores,
            "grid_score": scores.get("Grid_Score", 0),
            "zone_coverage": zone_coverage,
            "zone_energies": zone_energies,
            "ref_zone_energies": ref_zone_energies,
            "plip_interactions": plip_interactions,
            "n_plip_interactions": len(plip_interactions),
            "mmpbsa_dg": mmpbsa_dg,
            "mmpbsa_data": mol_mmpbsa,
            **decision,
        })

    # Sort by recommendation strength then grid score
    rec_order = {"Strong Pharmit template candidate": 0, "Purchase candidate": 1, "Weak candidate": 2}
    molecules.sort(key=lambda m: (rec_order.get(m["recommendation"], 3), m["grid_score"]))

    # --- Generate HTML report ---
    html = generate_decision_html(molecules, campaign_id, ref_scores, reference_label)
    html_path = output_dir / "decision_report.html"
    html_path.write_text(html, encoding="utf-8")
    logger.info(f"  Saved: {html_path}")

    # --- Generate summary CSV ---
    summary_rows = []
    for i, mol in enumerate(molecules, 1):
        summary_rows.append({
            "Rank": i,
            "Name": mol["name"],
            "Grid_Score": mol["grid_score"],
            "MMPBSA_dG": mol.get("mmpbsa_dg"),
            "Recommendation": mol["recommendation"],
            "Zones_Covered": mol["n_zones_covered"],
            "Zones_List": ",".join(mol["zones_covered"]),
            "N_PLIP_Interactions": mol["n_plip_interactions"],
            "Reasoning": mol["reasoning"],
        })

    df_summary = pd.DataFrame(summary_rows)
    summary_csv = output_dir / "decision_summary.csv"
    df_summary.to_csv(summary_csv, index=False, encoding="utf-8")
    logger.info(f"  Saved: {summary_csv}")

    # --- Log summary ---
    n_strong = sum(1 for m in molecules if "Strong" in m["recommendation"])
    n_purchase = sum(1 for m in molecules if "Purchase" in m["recommendation"])
    n_weak = sum(1 for m in molecules if "Weak" in m["recommendation"])

    logger.info("")
    logger.info(f"{'=' * 60}")
    logger.info(f"  DECISION REPORT: {len(molecules)} molecules")
    logger.info(f"    Strong template: {n_strong}")
    logger.info(f"    Purchase:        {n_purchase}")
    logger.info(f"    Weak:            {n_weak}")
    logger.info(f"{'=' * 60}")

    return {
        "success": True,
        "n_molecules": len(molecules),
        "n_strong": n_strong,
        "n_purchase": n_purchase,
        "n_weak": n_weak,
        "html_path": str(html_path),
        "summary_csv": str(summary_csv),
    }