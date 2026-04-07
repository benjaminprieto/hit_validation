"""
Residue Comparison Report - Core Module (07b)
================================================
Residue-by-residue comparison of each hit vs crystallographic reference.
07a decides purchase/template. 07b explains WHY — showing per-residue
energy gains and losses vs reference, organized by binding site zone.

HTML structure: one section per zone, each with a table of residues x molecules.
Plus a functional scores summary and an "Others" section for residues outside
defined zones.

Input:
  - footprint_per_molecule.csv (from 04b, obligatory)
  - campaign_config zones (obligatory)
  - PLIP interactions.json per molecule (from 03a, optional)
  - decision_summary.csv (from 07a, optional)

Output:
  - residue_comparison.html              (full HTML report)
  - residue_comparison_summary.csv       (per-zone energy totals per molecule)
  - per_molecule/{name}_residues.csv     (residue detail per molecule)

Location: 01_src/hit_validation/m07_decision_report/residue_comparison.py
Project: hit_validation
Module: 07b (core)
Version: 2.0 (2026-04-02)
"""

import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

logger = logging.getLogger(__name__)


# =============================================================================
# DATA LOADING
# =============================================================================

def load_footprint_data(footprint_csv: str) -> pd.DataFrame:
    """Load footprint_per_molecule.csv from 04b."""
    df = pd.read_csv(footprint_csv)
    required = {"Name", "residue_id", "vdw", "es", "total",
                "ref_vdw", "ref_es", "ref_total",
                "delta_vdw", "delta_es", "delta_total"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns in footprint CSV: {missing}")
    return df


def load_plip_interactions(plip_json: str) -> List[Dict]:
    """Load PLIP interactions from 03a for a single molecule."""
    if not Path(plip_json).exists():
        return []
    try:
        with open(plip_json, "r", encoding="utf-8") as f:
            data = json.load(f)
        if isinstance(data, dict):
            return data.get("interactions", [])
        return data if isinstance(data, list) else []
    except Exception as e:
        logger.debug(f"  Could not load PLIP: {e}")
        return []


def load_decision_summary(decision_csv: str) -> pd.DataFrame:
    """Load decision_summary.csv from 07a."""
    if not Path(decision_csv).exists():
        return pd.DataFrame()
    return pd.read_csv(decision_csv)


# =============================================================================
# ZONE ASSIGNMENT
# =============================================================================

def assign_zone(residue_id: str, zones: Dict) -> str:
    """Assign a residue to its zone label. Returns zone label or 'Others'."""
    res_prefix = residue_id.split(".")[0]
    for zone_id, zdef in zones.items():
        if res_prefix in zdef.get("residues", set()):
            return zdef.get("label", zone_id)
    return "Others"


def assign_zone_id(residue_id: str, zones: Dict) -> str:
    """Assign a residue to its zone_id. Returns zone_id or 'others'."""
    res_prefix = residue_id.split(".")[0]
    for zone_id, zdef in zones.items():
        if res_prefix in zdef.get("residues", set()):
            return zone_id
    return "others"


# =============================================================================
# MOLECULE NAME SHORTENING
# =============================================================================

def shorten_name(name: str) -> str:
    """
    Shorten molecule name for table headers.

    HTS1710-01239291-01 -> 01239291 (penultimate segment)
    PubChem-11165043    -> 11165043 (last segment, only 2 parts)
    Other               -> truncated to 10 chars
    """
    parts = name.split("-")
    if len(parts) >= 3:
        return parts[-2]
    elif len(parts) == 2:
        return parts[-1]
    return name[:10]


# =============================================================================
# PLIP SUMMARY PER RESIDUE
# =============================================================================

def plip_by_residue(interactions: List[Dict]) -> Dict[str, str]:
    """
    Summarize PLIP interactions by residue_id.

    Returns dict: residue_id -> "hbond 2.9A, pi_stack 4.4A"
    """
    res_map: Dict[str, str] = {}
    for ix in interactions:
        itype = ix.get("interaction_type", ix.get("type", "unknown"))
        residue = ix.get("residue", "")
        resnr = ix.get("residue_number", "")
        chain = ix.get("chain", "A")
        dist = ix.get("distance", ix.get("dist", 0))

        if residue and resnr:
            rid = f"{residue}{resnr}.{chain}"
        elif residue:
            rid = residue
        else:
            continue

        desc = f"{itype}"
        if dist:
            desc += f" {float(dist):.1f}A"

        if rid in res_map:
            res_map[rid] += f", {desc}"
        else:
            res_map[rid] = desc
    return res_map


# =============================================================================
# RESIDUE FILTERING
# =============================================================================

def filter_relevant_residues(df_mol: pd.DataFrame,
                             energy_cutoff: float = 0.5) -> pd.DataFrame:
    """
    Filter residues where hit or reference has significant energy.

    Keeps residues where |total| >= cutoff OR |ref_total| >= cutoff.
    """
    mask = (df_mol["total"].abs() >= energy_cutoff) | \
           (df_mol["ref_total"].abs() >= energy_cutoff)
    return df_mol[mask].copy()


# =============================================================================
# PER-MOLECULE RESIDUE TABLE
# =============================================================================

def build_molecule_table(
        df_mol: pd.DataFrame,
        zones: Dict,
        energy_cutoff: float = 0.5,
        plip_map: Optional[Dict[str, str]] = None,
) -> pd.DataFrame:
    """
    Build the per-molecule residue comparison table.

    Returns DataFrame sorted by delta_total ascending (hit more favorable first).
    """
    filtered = filter_relevant_residues(df_mol, energy_cutoff)
    if filtered.empty:
        return pd.DataFrame()

    rows = []
    for _, r in filtered.iterrows():
        row = {
            "residue_id": r["residue_id"],
            "zone": assign_zone(r["residue_id"], zones),
            "zone_id": assign_zone_id(r["residue_id"], zones),
            "ref_vdw": round(r["ref_vdw"], 2),
            "ref_es": round(r["ref_es"], 2),
            "ref_total": round(r["ref_total"], 2),
            "hit_vdw": round(r["vdw"], 2),
            "hit_es": round(r["es"], 2),
            "hit_total": round(r["total"], 2),
            "delta_vdw": round(r["delta_vdw"], 2),
            "delta_es": round(r["delta_es"], 2),
            "delta_total": round(r["delta_total"], 2),
            "plip": plip_map.get(r["residue_id"], "") if plip_map else "",
        }
        rows.append(row)

    result = pd.DataFrame(rows)
    if not result.empty:
        result = result.sort_values("delta_total", ascending=True)
    return result


# =============================================================================
# ZONE-LEVEL AGGREGATION
# =============================================================================

def build_zone_tables(
        df_fp: pd.DataFrame,
        zones: Dict,
        molecules: List[str],
        energy_cutoff: float = 0.5,
) -> Dict[str, pd.DataFrame]:
    """
    Build one table per zone: rows = residues, columns = ref + each molecule.

    Returns dict: zone_label -> DataFrame with columns:
        residue_id, ref_total, mol1_total, mol1_delta, mol2_total, mol2_delta, ...
    """
    # Assign zone to every row
    df = df_fp.copy()
    df["zone_label"] = df["residue_id"].apply(lambda r: assign_zone(r, zones))

    # Collect all relevant residues per zone
    zone_residues: Dict[str, set] = {}
    for _, r in df.iterrows():
        rid = r["residue_id"]
        zlabel = r["zone_label"]
        if abs(r["total"]) >= energy_cutoff or abs(r["ref_total"]) >= energy_cutoff:
            zone_residues.setdefault(zlabel, set()).add(rid)

    # Build ordered zone list: defined zones first (in config order), then Others
    zone_order = []
    for zone_id, zdef in zones.items():
        label = zdef.get("label", zone_id)
        if label in zone_residues:
            zone_order.append(label)
    if "Others" in zone_residues:
        zone_order.append("Others")

    result = {}
    for zlabel in zone_order:
        residues = sorted(zone_residues[zlabel])
        rows = []
        for rid in residues:
            row = {"residue_id": rid}
            # Reference (same across all molecules, take first)
            ref_vals = df[df["residue_id"] == rid]["ref_total"]
            row["ref_total"] = round(ref_vals.iloc[0], 2) if not ref_vals.empty else 0.0

            for mol in molecules:
                mol_data = df[(df["Name"] == mol) & (df["residue_id"] == rid)]
                if not mol_data.empty:
                    row[f"{mol}_total"] = round(mol_data.iloc[0]["total"], 2)
                    row[f"{mol}_delta"] = round(mol_data.iloc[0]["delta_total"], 2)
                else:
                    row[f"{mol}_total"] = 0.0
                    row[f"{mol}_delta"] = 0.0 - row["ref_total"]
            rows.append(row)
        result[zlabel] = pd.DataFrame(rows)

    return result


def compute_zone_totals(
        df_fp: pd.DataFrame,
        zones: Dict,
        molecules: List[str],
        energy_cutoff: float = 0.5,
) -> pd.DataFrame:
    """
    Compute total energy per zone per molecule.

    Returns DataFrame: one row per molecule, columns = zone labels + Others.
    """
    df = df_fp.copy()
    df["zone_label"] = df["residue_id"].apply(lambda r: assign_zone(r, zones))

    # Determine zone order
    zone_order = []
    for zone_id, zdef in zones.items():
        zone_order.append(zdef.get("label", zone_id))
    zone_order.append("Others")

    rows = []
    for mol in molecules:
        df_mol = df[df["Name"] == mol]
        row = {"Name": mol}
        for zlabel in zone_order:
            zone_data = df_mol[df_mol["zone_label"] == zlabel]
            relevant = zone_data[
                (zone_data["total"].abs() >= energy_cutoff) |
                (zone_data["ref_total"].abs() >= energy_cutoff)
            ]
            row[zlabel] = round(relevant["total"].sum(), 2) if not relevant.empty else 0.0
        rows.append(row)

    # Also compute reference totals
    ref_row = {"Name": "__reference__"}
    # Get reference values (same for all molecules, pick first molecule)
    if molecules:
        df_first = df[df["Name"] == molecules[0]]
        for zlabel in zone_order:
            zone_data = df_first[df_first["zone_label"] == zlabel]
            relevant = zone_data[
                (zone_data["total"].abs() >= energy_cutoff) |
                (zone_data["ref_total"].abs() >= energy_cutoff)
            ]
            ref_row[zlabel] = round(relevant["ref_total"].sum(), 2) if not relevant.empty else 0.0
        rows.insert(0, ref_row)

    return pd.DataFrame(rows)


# =============================================================================
# HTML REPORT GENERATION
# =============================================================================

_CSS = """\
body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
       max-width: 1400px; margin: 0 auto; padding: 20px; background: #fafafa;
       color: #333; line-height: 1.5; }
h1 { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }
h2 { color: #2c3e50; margin-top: 30px; }
h3 { color: #34495e; margin-top: 20px; }
.meta { color: #666; font-size: 0.9em; margin-bottom: 20px; }
.ref-note { background: #fff3cd; border-left: 4px solid #ffc107;
            padding: 10px 15px; margin: 15px 0; font-size: 0.9em; }
table { border-collapse: collapse; width: 100%; margin: 15px 0;
        font-size: 0.85em; background: white; }
th { background: #2c3e50; color: white; padding: 8px 10px; text-align: right;
     font-weight: 600; }
th:first-child { text-align: left; }
td { padding: 6px 10px; border-bottom: 1px solid #eee; text-align: right; }
td:first-child { text-align: left; }
tr:hover { background: #f5f5f5; }
tr.subtotal { font-weight: 700; border-top: 2px solid #2c3e50;
              background: #ecf0f1 !important; }
.better { color: #27ae60; font-weight: 600; }
.worse { color: #e74c3c; font-weight: 600; }
.neutral { color: #999; }
.catalytic { background: #fff3cd !important; }
.zone-section { background: white; border-radius: 8px; padding: 20px;
                margin: 20px 0; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }
.zone-section h3 { margin-top: 0; }
.badge { display: inline-block; padding: 3px 10px; border-radius: 12px;
         font-size: 0.8em; font-weight: 600; margin-left: 8px; }
.badge-strong { background: #d4edda; color: #155724; }
.badge-purchase { background: #cce5ff; color: #004085; }
.badge-weak { background: #f8d7da; color: #721c24; }
.legend { background: #f8f9fa; padding: 10px 15px; border-radius: 6px;
          margin: 10px 0; font-size: 0.85em; }
"""


def _delta_class(val: float, better_thresh: float, worse_thresh: float) -> str:
    """Return CSS class for a delta value."""
    if val < better_thresh:
        return "better"
    elif val > worse_thresh:
        return "worse"
    return "neutral"


def _delta_arrow(val: float, better_thresh: float, worse_thresh: float) -> str:
    """Return arrow indicator."""
    if val < better_thresh:
        return " \u25bc"
    elif val > worse_thresh:
        return " \u25b2"
    return ""


def _fmt_energy(val: float) -> str:
    """Format energy value."""
    if val == 0.0:
        return "-"
    return f"{val:.2f}"


def _is_catalytic(residue_id: str, zones: Dict) -> bool:
    """Check if residue is in the catalytic zone."""
    res_prefix = residue_id.split(".")[0]
    for zdef in zones.values():
        label_lower = zdef.get("label", "").lower()
        if "catalytic" in label_lower or "catalyt" in label_lower:
            if res_prefix in zdef.get("residues", set()):
                return True
    return False


def _recommendation_badge(rec: str) -> str:
    """Generate HTML badge for recommendation."""
    if not rec:
        return ""
    rec_lower = rec.lower()
    if "strong" in rec_lower:
        return f'<span class="badge badge-strong">{rec}</span>'
    elif "purchase" in rec_lower:
        return f'<span class="badge badge-purchase">{rec}</span>'
    return f'<span class="badge badge-weak">{rec}</span>'


def generate_html_report(
        zone_tables: Dict[str, pd.DataFrame],
        zone_totals: pd.DataFrame,
        molecules: List[str],
        zones: Dict,
        campaign_id: str,
        reference_label: str = "Reference",
        delta_better: float = -1.0,
        delta_worse: float = 1.0,
        decision_df: Optional[pd.DataFrame] = None,
        plip_maps: Optional[Dict[str, Dict[str, str]]] = None,
) -> str:
    """Generate the full HTML report organized by zone."""

    ref_label = reference_label
    short_names = {mol: shorten_name(mol) for mol in molecules}

    html = []
    html.append("<!DOCTYPE html>")
    html.append("<html><head><meta charset='utf-8'>")
    html.append(f"<title>07b Residue Comparison — {campaign_id}</title>")
    html.append(f"<style>{_CSS}</style>")
    html.append("</head><body>")

    # --- Header ---
    html.append("<h1>Residue Comparison Report</h1>")
    html.append('<div class="meta">')
    html.append(f"Campaign: <strong>{campaign_id}</strong> | "
                f"Molecules: <strong>{len(molecules)}</strong> | "
                f"Reference: <strong>{ref_label}</strong> | "
                f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    html.append("</div>")

    html.append('<div class="ref-note">')
    html.append(f"<strong>Reference source:</strong> campaign (flex docking, "
                f"box 5.0 A). {ref_label} was re-docked in the same context "
                f"as hits for fair comparison. Crystal-pose energies may differ.")
    html.append("</div>")

    html.append('<div class="legend">')
    html.append(f'<span class="better">\u25bc Better than {ref_label} '
                f'(delta &lt; {delta_better})</span> &nbsp;|&nbsp; '
                f'<span class="worse">\u25b2 Worse than {ref_label} '
                f'(delta &gt; {delta_worse})</span> &nbsp;|&nbsp; '
                f'<span class="neutral">Similar (|delta| &le; {abs(delta_better)})</span>')
    html.append("</div>")

    # --- Section 1: Functional Scores (per-zone energy totals) ---
    html.append("<h2>1. Functional Scores (Energy by Zone)</h2>")

    if not zone_totals.empty:
        zone_cols = [c for c in zone_totals.columns if c != "Name"]

        # Get reference row and molecule rows
        ref_row = zone_totals[zone_totals["Name"] == "__reference__"]
        mol_rows = zone_totals[zone_totals["Name"] != "__reference__"]

        html.append("<table>")
        header = "<tr><th>Molecule</th>"
        if decision_df is not None and not decision_df.empty:
            header += "<th>Rec. 07a</th>"
        for zc in zone_cols:
            header += f"<th>{zc}</th>"
        header += "</tr>"
        html.append(header)

        # Reference row
        if not ref_row.empty:
            r = ref_row.iloc[0]
            html.append("<tr class='subtotal'>")
            html.append(f"<td><strong>{ref_label}</strong></td>")
            if decision_df is not None and not decision_df.empty:
                html.append("<td>-</td>")
            for zc in zone_cols:
                html.append(f"<td><strong>{_fmt_energy(r[zc])}</strong></td>")
            html.append("</tr>")

        # Sort by catalytic zone energy (zone_E) ascending
        catalytic_label = None
        for zdef in zones.values():
            if "catalytic" in zdef.get("label", "").lower():
                catalytic_label = zdef["label"]
                break
        if catalytic_label and catalytic_label in mol_rows.columns:
            mol_rows = mol_rows.sort_values(catalytic_label, ascending=True)

        for _, row in mol_rows.iterrows():
            mol = row["Name"]
            html.append("<tr>")
            html.append(f"<td>{short_names.get(mol, mol)}</td>")

            if decision_df is not None and not decision_df.empty:
                match = decision_df[decision_df["Name"] == mol]
                rec = str(match.iloc[0].get("Recommendation", "")) if not match.empty else ""
                html.append(f"<td>{_recommendation_badge(rec)}</td>")

            ref_vals = ref_row.iloc[0] if not ref_row.empty else None
            for zc in zone_cols:
                val = row[zc]
                ref_val = ref_vals[zc] if ref_vals is not None else 0.0
                delta = val - ref_val
                cls = _delta_class(delta, delta_better, delta_worse)
                arrow = _delta_arrow(delta, delta_better, delta_worse)
                html.append(f'<td class="{cls}">{_fmt_energy(val)}{arrow}</td>')
            html.append("</tr>")

        html.append("</table>")

    # --- Section 2+: One section per zone ---
    section_num = 2
    for zlabel, ztable in zone_tables.items():
        if ztable.empty:
            continue

        html.append(f'<div class="zone-section">')
        html.append(f"<h3>{section_num}. {zlabel}</h3>")

        html.append("<table>")
        # Header: Residue | Reference | mol1 | mol2 | ...
        header = f"<tr><th>Residue</th><th>{ref_label}</th>"
        for mol in molecules:
            header += f"<th>{short_names[mol]}</th>"
        header += "</tr>"
        html.append(header)

        # Residue rows
        for _, row in ztable.iterrows():
            rid = row["residue_id"]
            cat_class = ' class="catalytic"' if _is_catalytic(rid, zones) else ""
            html.append(f"<tr{cat_class}>")
            html.append(f"<td>{rid}</td>")
            html.append(f"<td><strong>{_fmt_energy(row['ref_total'])}</strong></td>")

            for mol in molecules:
                val = row.get(f"{mol}_total", 0.0)
                delta = row.get(f"{mol}_delta", 0.0)
                cls = _delta_class(delta, delta_better, delta_worse)
                arrow = _delta_arrow(delta, delta_better, delta_worse)
                html.append(f'<td class="{cls}">{_fmt_energy(val)}{arrow}</td>')
            html.append("</tr>")

        # Subtotal row
        html.append('<tr class="subtotal">')
        html.append("<td><strong>Subtotal</strong></td>")
        ref_sum = ztable["ref_total"].sum()
        html.append(f"<td><strong>{ref_sum:.2f}</strong></td>")
        for mol in molecules:
            col = f"{mol}_total"
            mol_sum = ztable[col].sum() if col in ztable.columns else 0.0
            delta_sum = mol_sum - ref_sum
            cls = _delta_class(delta_sum, delta_better, delta_worse)
            arrow = _delta_arrow(delta_sum, delta_better, delta_worse)
            html.append(f'<td class="{cls}"><strong>{mol_sum:.2f}{arrow}</strong></td>')
        html.append("</tr>")

        html.append("</table>")

        # PLIP annotations under the zone table
        if plip_maps:
            zone_residues = set(ztable["residue_id"].tolist())
            plip_notes = []
            for mol in molecules:
                pmap = plip_maps.get(mol, {})
                for rid in sorted(zone_residues):
                    if rid in pmap:
                        plip_notes.append(
                            f"<strong>{short_names[mol]}</strong> {rid}: {pmap[rid]}")
            if plip_notes:
                html.append(f'<p style="font-size:0.82em;color:#555;">'
                            f'PLIP: {" | ".join(plip_notes)}</p>')

        html.append("</div>")
        section_num += 1

    # Footer
    html.append("<hr>")
    html.append('<p style="color: #999; font-size: 0.8em;">'
                "Generated by hit_validation 07b | "
                f"{datetime.now().strftime('%Y-%m-%d %H:%M')}</p>")
    html.append("</body></html>")

    return "\n".join(html)


# =============================================================================
# MAIN PIPELINE FUNCTION
# =============================================================================

def run_residue_comparison(
        footprint_csv: str,
        output_dir: str,
        zones: Dict,
        campaign_id: str = "",
        energy_cutoff: float = 0.5,
        delta_better_threshold: float = -1.0,
        delta_worse_threshold: float = 1.0,
        include_plip: bool = True,
        plip_analysis_dir: Optional[str] = None,
        decision_summary_csv: Optional[str] = None,
        reference_label: str = "Reference",
) -> Dict[str, Any]:
    """
    Run the residue comparison analysis.

    Returns:
        Dict with: success, n_molecules, output files
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    (output_path / "per_molecule").mkdir(parents=True, exist_ok=True)

    # --- Load obligatory data ---
    if not Path(footprint_csv).exists():
        return {"success": False, "error": f"Footprint CSV not found: {footprint_csv}"}

    logger.info("=" * 60)
    logger.info("  Residue Comparison Report (07b)")
    logger.info("=" * 60)
    logger.info(f"  Campaign:    {campaign_id}")
    logger.info(f"  Footprint:   {footprint_csv}")
    logger.info(f"  Reference:   {reference_label}")
    logger.info(f"  Cutoff:      {energy_cutoff} kcal/mol")

    try:
        df_fp = load_footprint_data(footprint_csv)
    except Exception as e:
        return {"success": False, "error": f"Failed to load footprint data: {e}"}

    molecules = sorted(df_fp["Name"].unique())
    logger.info(f"  Molecules:   {len(molecules)}")

    # --- Load optional data ---
    decision_df = None
    if decision_summary_csv:
        decision_df = load_decision_summary(decision_summary_csv)
        if not decision_df.empty:
            logger.info(f"  07a summary: {len(decision_df)} molecules")

    # --- Load PLIP per molecule ---
    plip_maps: Dict[str, Dict[str, str]] = {}
    if include_plip and plip_analysis_dir:
        for mol in molecules:
            plip_path = Path(plip_analysis_dir) / mol / "interactions.json"
            interactions = load_plip_interactions(str(plip_path))
            if interactions:
                plip_maps[mol] = plip_by_residue(interactions)
        if plip_maps:
            logger.info(f"  PLIP loaded: {len(plip_maps)} molecules")

    # --- Per-molecule tables + CSVs ---
    logger.info("  Processing molecules...")
    for idx, mol in enumerate(molecules, 1):
        logger.info(f"  [{idx}/{len(molecules)}] {mol}")
        df_mol = df_fp[df_fp["Name"] == mol].copy()
        mol_table = build_molecule_table(
            df_mol=df_mol,
            zones=zones,
            energy_cutoff=energy_cutoff,
            plip_map=plip_maps.get(mol),
        )
        if not mol_table.empty:
            mol_csv = output_path / "per_molecule" / f"{mol}_residues.csv"
            mol_table.to_csv(mol_csv, index=False, encoding="utf-8")

    # --- Zone tables ---
    logger.info("  Building zone tables...")
    zone_tables = build_zone_tables(
        df_fp=df_fp,
        zones=zones,
        molecules=molecules,
        energy_cutoff=energy_cutoff,
    )

    # --- Zone totals (functional scores) ---
    logger.info("  Computing zone totals...")
    zone_totals = compute_zone_totals(
        df_fp=df_fp,
        zones=zones,
        molecules=molecules,
        energy_cutoff=energy_cutoff,
    )
    summary_csv = output_path / "residue_comparison_summary.csv"
    if not zone_totals.empty:
        # Save without the __reference__ marker in Name
        save_totals = zone_totals.copy()
        save_totals.loc[save_totals["Name"] == "__reference__", "Name"] = reference_label
        save_totals.to_csv(summary_csv, index=False, encoding="utf-8")
        logger.info(f"  Saved: {summary_csv}")

    # --- HTML Report ---
    logger.info("  Generating HTML report...")
    html_content = generate_html_report(
        zone_tables=zone_tables,
        zone_totals=zone_totals,
        molecules=molecules,
        zones=zones,
        campaign_id=campaign_id,
        reference_label=reference_label,
        delta_better=delta_better_threshold,
        delta_worse=delta_worse_threshold,
        decision_df=decision_df,
        plip_maps=plip_maps if include_plip else None,
    )

    html_path = output_path / "residue_comparison.html"
    html_path.write_text(html_content, encoding="utf-8")
    logger.info(f"  Saved: {html_path}")

    # --- Summary ---
    n_zones = len(zone_tables)
    n_residues = sum(len(zt) for zt in zone_tables.values())

    logger.info("")
    logger.info(f"{'=' * 60}")
    logger.info(f"  07b COMPLETE: {len(molecules)} molecules, {n_zones} zones, {n_residues} residues")
    logger.info(f"{'=' * 60}")

    return {
        "success": True,
        "n_molecules": len(molecules),
        "n_zones": n_zones,
        "n_residues": n_residues,
        "html_report": str(html_path),
        "summary_csv": str(summary_csv),
        "output_dir": str(output_dir),
    }