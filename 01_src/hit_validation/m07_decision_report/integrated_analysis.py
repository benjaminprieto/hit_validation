"""
Integrated Residue Analysis — Core Module (07c)
==================================================
Per-residue, per-molecule report answering 7 questions using all pipeline data.
No predefined zones. No forced comparison against UDX. Each molecule on its own merit.

Questions:
  Q1. Does it bind?          -> MMPBSA dG (01h)
  Q2. Is the pose stable?    -> RMSD from MD (01i)
  Q3. What residues?         -> Footprint per-residue (04b)
  Q4. What interaction types? -> Footprint VDW/ES + MMPBSA decomp + PLIP
  Q5. Do interactions persist? -> MD distances + H-bond occupancy (01i)
  Q6. Water-mediated?        -> Water bridges (01i)
  Q7. Functional residues?   -> Annotation from campaign_config

Works incrementally: runs with whatever data exists. Missing sources marked "pending".

Input:
  Required: 04b footprint_per_molecule.csv, 01h consolidated_mmpbsa.csv
  Optional: 01h per_residue_decomp (per molecule), 03a PLIP, 01i MD data

Output:
  integrated_report.html
  integrated_summary.csv
  per_molecule/{mol}_residue_analysis.csv
  consistency_matrix.csv

Location: 01_src/hit_validation/m07_decision_report/integrated_analysis.py
Project: hit_validation
Module: 07c (core)
Version: 1.0
"""

import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# =============================================================================
# DATA LOADERS
# =============================================================================

def load_footprint(footprint_csv: Union[str, Path]) -> pd.DataFrame:
    """Load 04b footprint_per_molecule.csv."""
    return pd.read_csv(footprint_csv)


def load_mmpbsa_global(consolidated_csv: Union[str, Path]) -> pd.DataFrame:
    """Load 01h consolidated_mmpbsa.csv (one row per molecule)."""
    return pd.read_csv(consolidated_csv)


def load_mmpbsa_decomp(decomp_csv: Union[str, Path]) -> Optional[pd.DataFrame]:
    """Load 01h per_residue_decomp.csv for one molecule."""
    path = Path(decomp_csv)
    if not path.exists():
        return None
    return pd.read_csv(decomp_csv)


def load_plip(plip_json: Union[str, Path]) -> Optional[List[Dict]]:
    """Load 03a interactions.json for one molecule."""
    path = Path(plip_json)
    if not path.exists():
        return None
    with open(path) as f:
        data = json.load(f)
    return data.get("interactions", [])


def load_md_rmsd(rmsd_csv: Union[str, Path]) -> Optional[Dict]:
    """Load 01i rmsd_ligand.csv, return summary stats."""
    path = Path(rmsd_csv)
    if not path.exists():
        return None
    df = pd.read_csv(rmsd_csv)
    if "value" not in df.columns:
        return None
    return {
        "rmsd_avg": round(df["value"].mean(), 2),
        "rmsd_max": round(df["value"].max(), 2),
        "rmsd_final": round(df["value"].iloc[-1], 2),
        "rmsd_std": round(df["value"].std(), 2),
    }


def load_md_distances(distances_csv: Union[str, Path]) -> Optional[pd.DataFrame]:
    """Load 01i distances.csv (per-frame distances to monitored residues)."""
    path = Path(distances_csv)
    if not path.exists():
        return None
    return pd.read_csv(distances_csv)


def load_md_hbonds(hbond_csv: Union[str, Path]) -> Optional[pd.DataFrame]:
    """Load 01i hbond_occupancy.csv."""
    path = Path(hbond_csv)
    if not path.exists():
        return None
    return pd.read_csv(hbond_csv)


def load_water_bridges(wb_csv: Union[str, Path]) -> Optional[pd.DataFrame]:
    """Load 01i water_bridges.csv."""
    path = Path(wb_csv)
    if not path.exists():
        return None
    return pd.read_csv(wb_csv)


def load_prolif_occupancy(prolif_csv: Union[str, Path]) -> Optional[pd.DataFrame]:
    """Load 01i prolif_occupancy.csv."""
    path = Path(prolif_csv)
    if not path.exists():
        return None
    return pd.read_csv(prolif_csv)


# =============================================================================
# Q1: DOES IT BIND?
# =============================================================================

def assess_binding(
        mmpbsa_global: pd.DataFrame,
        mol_name: str,
        reject_threshold: float = 0.0,
        weak_threshold: float = -15.0,
) -> Dict:
    """
    Q1: Assess binding from MMPBSA global energy.
    Returns dict with dG, components, verdict.
    """
    row = mmpbsa_global[mmpbsa_global["Name"] == mol_name]
    if row.empty:
        return {"verdict": "NO_DATA", "dG": None}

    row = row.iloc[0]
    dG = row.get("MMPBSA_dG_total", row.get("MMPBSA_dG", None))
    vdw = row.get("MMPBSA_vdW", None)
    eel = row.get("MMPBSA_EEL", None)
    egb = row.get("MMPBSA_EGB", None)
    esurf = row.get("MMPBSA_ESURF", None)

    net_es = (eel + egb) if eel is not None and egb is not None else None

    if dG is None:
        verdict = "NO_DATA"
    elif dG > reject_threshold:
        verdict = "REJECT"
    elif dG > weak_threshold:
        verdict = "WEAK"
    else:
        verdict = "PASS"

    return {
        "verdict": verdict,
        "dG": round(dG, 2) if dG is not None else None,
        "vdw": round(vdw, 2) if vdw is not None else None,
        "eel": round(eel, 2) if eel is not None else None,
        "egb": round(egb, 2) if egb is not None else None,
        "esurf": round(esurf, 2) if esurf is not None else None,
        "net_es": round(net_es, 2) if net_es is not None else None,
    }


# =============================================================================
# Q2: IS THE POSE STABLE?
# =============================================================================

def assess_stability(
        rmsd_data: Optional[Dict],
        unstable_threshold: float = 4.0,
        mobile_threshold: float = 2.0,
) -> Dict:
    """Q2: Assess pose stability from RMSD."""
    if rmsd_data is None:
        return {"verdict": "PENDING", "rmsd_avg": None}

    rmsd_avg = rmsd_data["rmsd_avg"]

    if rmsd_avg > unstable_threshold:
        verdict = "UNSTABLE"
    elif rmsd_avg > mobile_threshold:
        verdict = "MOBILE"
    else:
        verdict = "STABLE"

    return {
        "verdict": verdict,
        **rmsd_data,
    }


# =============================================================================
# Q3-Q7: PER-RESIDUE ANALYSIS
# =============================================================================

def _build_pdb_to_seq_mapping(residue_mapping_csv: Optional[str]) -> Dict[str, str]:
    """
    Build PDB->sequential number mapping from 04b residue_mapping.csv.

    Mapping CSV format: mol2_sequential (ASP243), pdb_original (ASP494.A)
    Returns: {"ASP494.A": "243", "GLU529.A": "278", ...}
    """
    if not residue_mapping_csv or not Path(residue_mapping_csv).exists():
        return {}
    import re
    df = pd.read_csv(residue_mapping_csv)
    mapping = {}
    for _, row in df.iterrows():
        pdb_id = row["pdb_original"]  # e.g. "ASP494.A"
        seq_id = row["mol2_sequential"]  # e.g. "ASP243"
        seq_num = re.search(r"(\d+)", seq_id)
        if seq_num:
            mapping[pdb_id] = seq_num.group(1)
    return mapping


def build_residue_table(
        mol_name: str,
        footprint_df: pd.DataFrame,
        mmpbsa_decomp: Optional[pd.DataFrame],
        plip_interactions: Optional[List[Dict]],
        md_distances: Optional[pd.DataFrame],
        md_hbonds: Optional[pd.DataFrame],
        water_bridges: Optional[pd.DataFrame],
        prolif_occ: Optional[pd.DataFrame],
        functional_residues: Dict,
        energy_cutoff: float = 0.5,
        pdb_to_seq: Optional[Dict[str, str]] = None,
) -> pd.DataFrame:
    """
    Build the central per-residue table integrating all data sources.
    Each row = one residue with all evidence columns.
    """
    # --- Q3: Footprint residues ---
    mol_fp = footprint_df[footprint_df["Name"] == mol_name].copy()

    if mol_fp.empty:
        logger.warning(f"  No footprint data for {mol_name}")
        return pd.DataFrame()

    # Filter to significant residues
    mol_fp = mol_fp[
        (mol_fp["total"].abs() >= energy_cutoff) |
        (mol_fp["ref_total"].abs() >= energy_cutoff)
    ].copy()

    # Build base table
    rows = []
    for _, fp_row in mol_fp.iterrows():
        res_id = fp_row["residue_id"]
        res_name = res_id.split(".")[0] if "." in res_id else res_id

        row = {
            "residue_id": res_id,
            "residue_name": res_name,
            # Q3: What residues
            "fp_vdw": round(fp_row["vdw"], 2),
            "fp_es": round(fp_row["es"], 2),
            "fp_total": round(fp_row["total"], 2),
        }

        # --- Q4: MMPBSA decomposition ---
        if mmpbsa_decomp is not None:
            mmpbsa_row = mmpbsa_decomp[mmpbsa_decomp["residue_id"] == res_id]
            if not mmpbsa_row.empty:
                mr = mmpbsa_row.iloc[0]
                row["mmpbsa_vdw"] = round(mr.get("vdw", 0), 2)
                row["mmpbsa_es"] = round(mr.get("es", 0), 2)
                row["mmpbsa_gb"] = round(mr.get("gb", 0), 2)
                row["mmpbsa_sa"] = round(mr.get("sa", 0), 4)
                row["mmpbsa_total"] = round(mr.get("total", 0), 2)
            else:
                row.update({k: None for k in ["mmpbsa_vdw", "mmpbsa_es", "mmpbsa_gb", "mmpbsa_sa", "mmpbsa_total"]})
        else:
            row.update({k: "pending" for k in ["mmpbsa_vdw", "mmpbsa_es", "mmpbsa_gb", "mmpbsa_sa", "mmpbsa_total"]})

        # --- Q4: PLIP interaction type ---
        # PLIP JSON format: {"residue": "SER575", "residue_number": 575,
        #   "interaction_type": "hbond", "distance": 2.8, ...}
        # res_name is e.g. "SER575" (from footprint residue_id "SER575.A")
        row["plip_type"] = None
        row["plip_distance"] = None
        if plip_interactions:
            for ix in plip_interactions:
                plip_res = ix.get("residue", "")  # e.g. "SER575"
                if plip_res == res_name or res_name.startswith(plip_res):
                    row["plip_type"] = ix.get("interaction_type", "")
                    dist_val = ix.get("distance", 0)
                    row["plip_distance"] = round(dist_val, 1) if dist_val else None
                    break

        # --- Q5: MD persistence ---
        # distances.csv columns: frame, ARG598.A, LYS599.A, ... (PDB IDs with chain)
        if md_distances is not None and res_id in md_distances.columns:
            dist_series = md_distances[res_id]
            n_contact = (dist_series < 4.0).sum()
            n_total = len(dist_series)
            row["md_contact_pct"] = round(n_contact / n_total * 100, 1) if n_total > 0 else None
            row["md_dist_avg"] = round(dist_series.mean(), 2)
            row["md_dist_std"] = round(dist_series.std(), 2)
        else:
            row["md_contact_pct"] = "pending"
            row["md_dist_avg"] = "pending"
            row["md_dist_std"] = "pending"

        # H-bond occupancy from MD
        # hbond_occupancy.csv: acceptor=ASP_243@OD2, donor=UNL_707@O2, occupancy=0.86
        # cpptraj uses sequential numbering: ASP_243 = sequential 243
        # Use pdb_to_seq mapping to find the sequential ID for this PDB residue
        seq_num = pdb_to_seq.get(res_id, "") if pdb_to_seq else ""
        # Build cpptraj-style prefix for matching: e.g. "ASP_243"
        import re
        res_type = re.match(r"([A-Z]{3})", res_name)
        res_type_3 = res_type.group(1) if res_type else res_name[:3]
        cpptraj_prefix = f"{res_type_3}_{seq_num}" if seq_num else ""

        row["md_hbond_pct"] = "pending"
        if md_hbonds is not None and not md_hbonds.empty:
            if cpptraj_prefix:
                hb_match = md_hbonds[
                    md_hbonds["donor"].str.startswith(cpptraj_prefix, na=False) |
                    md_hbonds["acceptor"].str.startswith(cpptraj_prefix, na=False)
                ]
            else:
                hb_match = pd.DataFrame()
            if not hb_match.empty:
                row["md_hbond_pct"] = round(hb_match["occupancy"].max(), 2)
            else:
                row["md_hbond_pct"] = 0.0

        # ProLIF occupancy
        row["prolif_pct"] = "pending"
        row["prolif_types"] = "pending"
        if prolif_occ is not None and not prolif_occ.empty:
            prolif_match = prolif_occ[
                prolif_occ["protein_residue"].str.contains(res_type_3, na=False)
            ]
            if not prolif_match.empty:
                row["prolif_pct"] = round(prolif_match["occupancy"].max(), 1)
                row["prolif_types"] = ", ".join(prolif_match["interaction_type"].unique())
            else:
                row["prolif_pct"] = 0.0
                row["prolif_types"] = ""

        # --- Q6: Water bridges ---
        # water_bridges.csv: protein_residues="110:ASP 139:THR", n_frames=151
        # cpptraj sequential numbering: "110:ASP" = sequential 110
        # Match by "seq_num:RES_TYPE" pattern
        wb_pattern = f"{seq_num}:{res_type_3}" if seq_num else ""

        row["water_bridge_pct"] = "pending"
        if water_bridges is not None and not water_bridges.empty:
            if wb_pattern:
                wb_match = water_bridges[
                    water_bridges["protein_residues"].str.contains(wb_pattern, na=False)
                ]
            else:
                wb_match = pd.DataFrame()
            if not wb_match.empty:
                row["water_bridge_n_frames"] = int(wb_match["n_frames"].max())
                row["water_bridge_pct"] = row["water_bridge_n_frames"]
            else:
                row["water_bridge_pct"] = 0
                row["water_bridge_n_frames"] = 0

        # --- Q7: Functional annotation ---
        row["function"] = ""
        row["evidence"] = ""
        for func_key, func_info in functional_residues.items():
            if res_name.startswith(func_key) or func_key in res_name:
                row["function"] = func_info.get("function", "")
                row["evidence"] = func_info.get("evidence", "")
                break

        # --- Desolvation trap detection ---
        row["desolvation_trap"] = False
        if (isinstance(row.get("mmpbsa_total"), (int, float)) and
                isinstance(row.get("fp_total"), (int, float))):
            if row["fp_total"] < -2.0 and row["mmpbsa_total"] > -0.5:
                row["desolvation_trap"] = True

        # --- Consistency scoring ---
        row["consistency"] = compute_consistency(row)

        rows.append(row)

    df = pd.DataFrame(rows)
    df = df.sort_values("fp_total", ascending=True).reset_index(drop=True)

    return df


def compute_consistency(row: Dict) -> str:
    """
    Compute consistency across evidence sources for one residue.
    Returns: "high", "moderate", "low", "conflicting", "pending"
    """
    signals = []

    # Footprint signal
    fp_total = row.get("fp_total")
    if isinstance(fp_total, (int, float)):
        signals.append("strong" if fp_total < -1.0 else ("weak" if fp_total < -0.5 else "absent"))

    # MMPBSA signal
    mmpbsa = row.get("mmpbsa_total")
    if isinstance(mmpbsa, (int, float)):
        signals.append("strong" if mmpbsa < -1.0 else ("weak" if mmpbsa < -0.5 else "absent"))
    elif mmpbsa == "pending":
        signals.append("pending")

    # PLIP signal
    plip = row.get("plip_type")
    if plip and plip != "None":
        signals.append("strong")
    elif plip is None:
        signals.append("absent")

    # MD signal
    md_contact = row.get("md_contact_pct")
    if isinstance(md_contact, (int, float)):
        signals.append("strong" if md_contact > 70 else ("weak" if md_contact > 30 else "absent"))
    elif md_contact == "pending":
        signals.append("pending")

    n_strong = signals.count("strong")
    n_absent = signals.count("absent")
    n_pending = signals.count("pending")

    if n_pending == len(signals):
        return "pending"

    if n_strong >= 3:
        return "high"
    elif n_strong >= 2:
        return "moderate"
    elif n_strong >= 1 and n_absent >= 2:
        return "conflicting"
    elif n_strong >= 1:
        return "low"
    else:
        return "low"


# =============================================================================
# HTML REPORT GENERATION
# =============================================================================

def generate_html_report(
        molecules: List[Dict],
        campaign_id: str,
) -> str:
    """Generate the complete HTML report."""

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")

    css = """
body { font-family: -apple-system, 'Segoe UI', sans-serif; max-width: 1400px; margin: 0 auto;
       padding: 20px; background: #fafafa; color: #333; line-height: 1.5; }
h1 { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }
h2 { color: #2c3e50; margin-top: 30px; }
h3 { color: #34495e; margin-top: 20px; }
.meta { color: #666; font-size: 0.9em; }
table { border-collapse: collapse; width: 100%; margin: 10px 0; font-size: 12px; background: white; }
th { background: #2c3e50; color: white; padding: 6px 8px; text-align: right; font-weight: 500; font-size: 11px; }
th:first-child { text-align: left; }
td { padding: 5px 8px; border-bottom: 1px solid #eee; text-align: right; font-family: monospace; font-size: 12px; }
td:first-child { text-align: left; font-family: sans-serif; font-weight: 500; }
tr:hover { background: #f5f5f5; }
.card { background: white; border-radius: 8px; padding: 20px; margin: 20px 0;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1); page-break-inside: avoid; }
.pass { color: #27ae60; font-weight: 600; }
.reject { color: #e74c3c; font-weight: 600; }
.weak { color: #f39c12; font-weight: 600; }
.pending { color: #95a5a6; font-style: italic; }
.stable { color: #27ae60; }
.mobile { color: #f39c12; }
.unstable { color: #e74c3c; }
.func-row { background: #fffde7 !important; }
.trap-icon { color: #e74c3c; font-weight: bold; }
.consistency-high { color: #27ae60; }
.consistency-moderate { color: #f39c12; }
.consistency-low { color: #95a5a6; }
.consistency-conflicting { color: #e74c3c; }
.summary-row-reject td { background: #fde8e8; }
.mono { font-family: monospace; }
"""

    # --- Summary table ---
    summary_rows = []
    for mol in molecules:
        q1 = mol["q1"]
        q2 = mol["q2"]
        n_res = mol.get("n_residues", 0)
        n_func = mol.get("n_functional", 0)
        n_persist = mol.get("n_persistent", "-")
        n_traps = mol.get("n_desolvation_traps", 0)

        q1_class = q1["verdict"].lower()
        q2_class = q2["verdict"].lower()
        row_class = ' class="summary-row-reject"' if q1["verdict"] == "REJECT" else ""

        summary_rows.append(
            f'<tr{row_class}>'
            f'<td>{mol["name"]}</td>'
            f'<td class="mono">{q1["dG"] if q1["dG"] is not None else "-"}</td>'
            f'<td class="mono">{q1.get("net_es", "-")}</td>'
            f'<td class="{q1_class}">{q1["verdict"]}</td>'
            f'<td class="mono">{q2.get("rmsd_avg", "-")}</td>'
            f'<td class="{q2_class}">{q2["verdict"]}</td>'
            f'<td>{n_res}</td>'
            f'<td>{n_func}</td>'
            f'<td>{n_persist}</td>'
            f'<td>{n_traps}</td>'
            f'</tr>'
        )

    summary_html = f"""
<h2>1. Summary</h2>
<table>
<tr><th>Molecule</th><th>dG</th><th>netES</th><th>Q1 Bind</th>
<th>RMSD</th><th>Q2 Stable</th>
<th>Residues</th><th>Functional</th><th>Persistent</th><th>Traps</th></tr>
{''.join(summary_rows)}
</table>
"""

    # --- Per-molecule cards ---
    mol_cards = []
    for mol in molecules:
        if mol["q1"]["verdict"] == "REJECT":
            mol_cards.append(
                f'<div class="card" style="border-left: 4px solid #e74c3c;">'
                f'<h3>{mol["name"]} <span class="reject">-- REJECTED (dG = {mol["q1"]["dG"]})</span></h3>'
                f'<p>MMPBSA dG is positive -- molecule does not bind in solvent.</p>'
                f'</div>'
            )
            continue

        card = generate_molecule_card(mol)
        mol_cards.append(card)

    html = f"""<!DOCTYPE html>
<html><head><meta charset='utf-8'>
<title>07c Integrated Analysis -- {campaign_id}</title>
<style>{css}</style>
</head><body>
<h1>Integrated Residue Analysis</h1>
<div class="meta">
Campaign: <strong>{campaign_id}</strong> | Molecules: <strong>{len(molecules)}</strong> | Generated: {timestamp}
</div>
<p class="meta">Per-residue integration of footprint (04b), MMPBSA (01h), PLIP (03a), and MD (01i).
No predefined zones. Each molecule evaluated on its own merit.</p>

{summary_html}

<h2>2. Per-Molecule Analysis</h2>
{''.join(mol_cards)}

<div style="margin-top:40px; padding-top:10px; border-top:1px solid #ddd; font-size:11px; color:#999;">
hit_validation | Module 07c | Integrated Analysis v1.0 | {timestamp}
</div>
</body></html>"""

    return html


def generate_molecule_card(mol: Dict) -> str:
    """Generate HTML card for one molecule."""
    name = mol["name"]
    q1 = mol["q1"]
    q2 = mol["q2"]
    df = mol.get("residue_table")

    if df is None or df.empty:
        return f'<div class="card"><h3>{name}</h3><p>No residue data available.</p></div>'

    q1_class = q1["verdict"].lower()
    q2_class = q2["verdict"].lower()

    header = (
        f'<h3>{name}</h3>'
        f'<p class="mono">'
        f'dG: {q1["dG"]} | vdW: {q1.get("vdw", "-")} | netES: {q1.get("net_es", "-")} '
        f'-> <span class="{q1_class}">{q1["verdict"]}</span> | '
        f'RMSD: {q2.get("rmsd_avg", "-")} A '
        f'-> <span class="{q2_class}">{q2["verdict"]}</span></p>'
    )

    table_rows = []
    for _, row in df.iterrows():
        func_class = ' class="func-row"' if row.get("function") else ""
        trap_icon = ' <span class="trap-icon">!</span>' if row.get("desolvation_trap") else ""

        # Format MD values
        md_contact = row.get("md_contact_pct")
        if isinstance(md_contact, (int, float)):
            if md_contact > 70:
                md_str = f'<span class="pass">{md_contact}%</span>'
            elif md_contact > 30:
                md_str = f'<span class="weak">{md_contact}%</span>'
            else:
                md_str = f'<span style="color:#999">{md_contact}%</span>'
        else:
            md_str = '<span class="pending">pend</span>'

        md_hb = row.get("md_hbond_pct")
        md_hb_str = f'{md_hb}%' if isinstance(md_hb, (int, float)) else '<span class="pending">pend</span>'

        wb = row.get("water_bridge_pct")
        wb_str = f'{wb}%' if isinstance(wb, (int, float)) else '<span class="pending">pend</span>'

        cons = row.get("consistency", "")
        cons_class = f"consistency-{cons}" if cons else ""

        plip_str = ""
        if row.get("plip_type") and row["plip_type"] != "None":
            plip_str = f'{row["plip_type"]} {row.get("plip_distance", "")}A'

        mmpbsa_str = str(row.get("mmpbsa_total", "pend"))
        gb_str = str(row.get("mmpbsa_gb", "pend"))

        func_str = row.get("function", "")

        table_rows.append(
            f'<tr{func_class}>'
            f'<td>{row["residue_id"]}{trap_icon}</td>'
            f'<td>{row["fp_vdw"]}</td>'
            f'<td>{row["fp_es"]}</td>'
            f'<td><b>{row["fp_total"]}</b></td>'
            f'<td>{mmpbsa_str}</td>'
            f'<td>{gb_str}</td>'
            f'<td>{plip_str}</td>'
            f'<td>{md_str}</td>'
            f'<td>{md_hb_str}</td>'
            f'<td>{wb_str}</td>'
            f'<td>{func_str}</td>'
            f'<td class="{cons_class}">{cons}</td>'
            f'</tr>'
        )

    table = f"""
<table>
<tr>
<th>Residue</th>
<th>FP vdw</th><th>FP es</th><th>FP total</th>
<th>MMPBSA</th><th>GB desolv</th>
<th>PLIP</th>
<th>MD contact</th><th>MD Hbond</th><th>Water br.</th>
<th>Function</th><th>Consist.</th>
</tr>
{''.join(table_rows)}
</table>"""

    return f'<div class="card">{header}{table}</div>'


# =============================================================================
# MAIN ORCHESTRATOR
# =============================================================================

def run_integrated_analysis(
        campaign_id: str,
        footprint_csv: Union[str, Path],
        mmpbsa_global_csv: Union[str, Path],
        output_dir: Union[str, Path],
        campaign_config: Dict,
        mmpbsa_decomp_dir: Optional[str] = None,
        plip_dir: Optional[str] = None,
        md_dir: Optional[str] = None,
        residue_mapping_csv: Optional[str] = None,
        energy_cutoff: float = 0.5,
        reject_threshold: float = 0.0,
        weak_threshold: float = -15.0,
        unstable_threshold: float = 4.0,
        mobile_threshold: float = 2.0,
) -> Dict:
    """Run integrated analysis for all molecules in a campaign."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    per_mol_dir = output_dir / "per_molecule"
    per_mol_dir.mkdir(exist_ok=True)

    logger.info("=" * 60)
    logger.info(f"  07c Integrated Residue Analysis")
    logger.info(f"  Campaign: {campaign_id}")
    logger.info("=" * 60)

    df_footprint = load_footprint(footprint_csv)
    df_mmpbsa_global = load_mmpbsa_global(mmpbsa_global_csv)

    mol_names = sorted(df_footprint["Name"].unique())
    logger.info(f"  Molecules: {len(mol_names)}")

    functional_residues = campaign_config.get("functional_residues", {})
    logger.info(f"  Functional residues: {len(functional_residues)}")

    # Load PDB->sequential mapping for H-bond/water bridge matching
    pdb_to_seq = _build_pdb_to_seq_mapping(residue_mapping_csv)
    if pdb_to_seq:
        logger.info(f"  Residue mapping: {len(pdb_to_seq)} residues (PDB->sequential)")
    else:
        logger.warning("  No residue mapping — H-bond/water bridge matching will be limited")

    all_molecules = []
    summary_rows = []
    consistency_rows = []

    for mol_name in mol_names:
        logger.info(f"  Processing: {mol_name}")

        # Q1: Binding
        q1 = assess_binding(df_mmpbsa_global, mol_name, reject_threshold, weak_threshold)

        # Q2: Stability (from 01i)
        rmsd_data = None
        if md_dir:
            rmsd_data = load_md_rmsd(Path(md_dir) / mol_name / "rmsd_ligand.csv")
        q2 = assess_stability(rmsd_data, unstable_threshold, mobile_threshold)

        # Load optional per-molecule data
        mmpbsa_decomp = None
        if mmpbsa_decomp_dir:
            mmpbsa_decomp = load_mmpbsa_decomp(
                Path(mmpbsa_decomp_dir) / mol_name / "per_residue_decomp.csv"
            )

        plip_interactions = None
        if plip_dir:
            plip_interactions = load_plip(Path(plip_dir) / mol_name / "interactions.json")

        md_distances = None
        md_hbonds = None
        water_bridges = None
        prolif_occ = None
        if md_dir:
            md_mol_dir = Path(md_dir) / mol_name
            md_distances = load_md_distances(md_mol_dir / "distances.csv")
            md_hbonds = load_md_hbonds(md_mol_dir / "hbond_occupancy.csv")
            water_bridges = load_water_bridges(md_mol_dir / "water_bridges.csv")
            prolif_occ = load_prolif_occupancy(md_mol_dir / "prolif" / "prolif_occupancy.csv")

        # Build per-residue table (Q3-Q7)
        df_residues = build_residue_table(
            mol_name=mol_name,
            footprint_df=df_footprint,
            mmpbsa_decomp=mmpbsa_decomp,
            plip_interactions=plip_interactions,
            md_distances=md_distances,
            md_hbonds=md_hbonds,
            water_bridges=water_bridges,
            prolif_occ=prolif_occ,
            functional_residues=functional_residues,
            energy_cutoff=energy_cutoff,
            pdb_to_seq=pdb_to_seq,
        )

        # Save per-molecule CSV
        if not df_residues.empty:
            df_residues.to_csv(
                per_mol_dir / f"{mol_name}_residue_analysis.csv",
                index=False, encoding="utf-8",
            )

        # Compute summary stats
        n_func = len(df_residues[df_residues["function"] != ""]) if not df_residues.empty else 0
        n_traps = len(df_residues[df_residues["desolvation_trap"] == True]) if not df_residues.empty else 0

        n_persistent = "-"
        if not df_residues.empty and "md_hbond_pct" in df_residues.columns:
            numeric_hbonds = pd.to_numeric(df_residues["md_hbond_pct"], errors="coerce")
            n_persistent = int((numeric_hbonds > 0.50).sum())

        mol_data = {
            "name": mol_name,
            "q1": q1,
            "q2": q2,
            "residue_table": df_residues,
            "n_residues": len(df_residues),
            "n_functional": n_func,
            "n_persistent": n_persistent,
            "n_desolvation_traps": n_traps,
        }
        all_molecules.append(mol_data)

        summary_rows.append({
            "Name": mol_name,
            "MMPBSA_dG": q1["dG"],
            "MMPBSA_netES": q1.get("net_es"),
            "Q1_verdict": q1["verdict"],
            "RMSD_avg": q2.get("rmsd_avg"),
            "Q2_verdict": q2["verdict"],
            "N_residues": len(df_residues),
            "N_functional": n_func,
            "N_persistent": n_persistent,
            "N_desolvation_traps": n_traps,
        })

        # Consistency matrix rows (top 10 residues)
        if not df_residues.empty:
            for _, row in df_residues.head(10).iterrows():
                consistency_rows.append({
                    "Name": mol_name,
                    "residue_id": row["residue_id"],
                    "fp_signal": "strong" if row["fp_total"] < -1 else "weak",
                    "mmpbsa_signal": (
                        "strong" if isinstance(row.get("mmpbsa_total"), (int, float)) and row["mmpbsa_total"] < -1
                        else ("pending" if row.get("mmpbsa_total") == "pending" else "weak")
                    ),
                    "plip_signal": "strong" if row.get("plip_type") and row["plip_type"] != "None" else "absent",
                    "md_signal": (
                        "strong" if isinstance(row.get("md_contact_pct"), (int, float)) and row["md_contact_pct"] > 70
                        else ("pending" if row.get("md_contact_pct") == "pending" else "weak")
                    ),
                    "consistency": row.get("consistency", ""),
                })

    # Sort: PASS first, then by dG
    all_molecules.sort(key=lambda m: (
        {"PASS": 0, "WEAK": 1, "REJECT": 2, "NO_DATA": 3}.get(m["q1"]["verdict"], 3),
        m["q1"]["dG"] if m["q1"]["dG"] is not None else 999,
    ))

    # Generate HTML
    html = generate_html_report(all_molecules, campaign_id)
    html_path = output_dir / "integrated_report.html"
    html_path.write_text(html, encoding="utf-8")
    logger.info(f"  Saved: {html_path}")

    # Save summary CSV
    df_summary = pd.DataFrame(summary_rows)
    summary_path = output_dir / "integrated_summary.csv"
    df_summary.to_csv(summary_path, index=False, encoding="utf-8")
    logger.info(f"  Saved: {summary_path}")

    # Save consistency matrix
    if consistency_rows:
        df_consistency = pd.DataFrame(consistency_rows)
        cons_path = output_dir / "consistency_matrix.csv"
        df_consistency.to_csv(cons_path, index=False, encoding="utf-8")
        logger.info(f"  Saved: {cons_path}")

    n_pass = sum(1 for m in all_molecules if m["q1"]["verdict"] == "PASS")
    n_weak = sum(1 for m in all_molecules if m["q1"]["verdict"] == "WEAK")
    n_reject = sum(1 for m in all_molecules if m["q1"]["verdict"] == "REJECT")

    logger.info("")
    logger.info(f"  Results: {n_pass} PASS, {n_weak} WEAK, {n_reject} REJECT")
    logger.info("=" * 60)

    return {
        "success": True,
        "n_molecules": len(all_molecules),
        "n_pass": n_pass,
        "n_weak": n_weak,
        "n_reject": n_reject,
        "html_path": str(html_path),
        "summary_csv": str(summary_path),
    }