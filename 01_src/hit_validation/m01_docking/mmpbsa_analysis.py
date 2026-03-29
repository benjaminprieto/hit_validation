"""
MMPBSA Analysis — Core Module (01h)
======================================
Parses MMPBSA.py decomposition output and compares with DOCK6 footprint.

Takes the raw output from 01g (FINAL_DECOMP_MMPBSA.dat, FINAL_RESULTS_MMPBSA.dat)
and produces:
  - per_residue_decomp.csv:             vdW, ES, GB, SA, total per residue (PDB numbering)
  - comparison_invacuo_vs_solvated.csv: merge with 01d/04b footprint
  - zone_summary.json:                  energy by binding site zone
  - global_binding_energy.json:         total ΔG components

Pipeline:
    6. Parse FINAL_DECOMP_MMPBSA.dat → per-residue DataFrame
    7. Map sequential → PDB residue numbering
    8. Parse FINAL_RESULTS_MMPBSA.dat → global binding energy
    9. Compare with 01d footprint (ranking delta, zone summary)

Location: 01_src/hit_validation/m01_docking/mmpbsa_analysis.py
Project: hit_validation
Module: 01h (core)
Version: 1.0 (2026-03-28)
"""

import json
import logging
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# =============================================================================
# ZONE DEFINITIONS
# =============================================================================
# Sub-pocket classification from UDX PLIP analysis (06a).
# Shared with 04b footprint_analysis.py.

ZONE_DEFINITIONS = {
    "phosphate": {
        "residues": {"ARG598", "LYS599"},
        "label": "Phosphate (salt bridges)",
        "description": "Strong in-vacuo ES, expect large GB desolvation penalty.",
    },
    "xylose": {
        "residues": {"TRP392", "TRP495", "TYR565", "SER575"},
        "label": "Xylose pocket",
        "description": "Hydrophobic/HBA interactions, expect minimal solvation change.",
    },
    "uracil": {
        "residues": {"ASP361", "THR390", "ARG363"},
        "label": "Uracil pocket",
        "description": "Mixed polar/charged, moderate solvation screening.",
    },
    "ribose_bridge": {
        "residues": {"ASP494", "GLU529", "HIS335"},
        "label": "Ribose bridge",
        "description": "ASP494 shows +2.74 ES repulsion in-vacuo. Water bridge in PLIP.",
    },
}


# =============================================================================
# GLOBAL RESULTS PARSING
# =============================================================================

def parse_mmpbsa_global(
        results_file: Union[str, Path],
) -> Dict[str, float]:
    """
    Parse global MMPBSA results (total binding energy components).

    Looks for the "Differences (Complex - Receptor - Ligand)" section
    in FINAL_RESULTS_MMPBSA.dat.

    Returns:
        Dict with delta_total, vdw, eel, egb, esurf (kcal/mol)
    """
    text = Path(results_file).read_text()
    results = {}

    in_delta = False
    for line in text.split("\n"):
        if "Differences (Complex - Receptor - Ligand)" in line:
            in_delta = True
            continue
        if not in_delta:
            continue

        stripped = line.strip()
        if not stripped or stripped.startswith("---") or stripped.startswith("Energy"):
            continue

        parts = stripped.split()
        if len(parts) < 2:
            continue

        if parts[0] == "DELTA" and len(parts) >= 3:
            key = parts[1]
            try:
                val = float(parts[2])
            except ValueError:
                continue
            if key == "TOTAL":
                results["delta_total"] = val
            elif key == "G":
                if len(parts) >= 4:
                    subkey = parts[2]
                    try:
                        results[f"delta_g_{subkey}"] = float(parts[3])
                    except ValueError:
                        pass
            continue

        key = parts[0]
        try:
            val = float(parts[1])
        except ValueError:
            continue

        if key == "VDWAALS":
            results["vdw"] = val
        elif key == "EEL":
            results["eel"] = val
        elif key == "EGB":
            results["egb"] = val
        elif key == "ESURF":
            results["esurf"] = val

    return results


# =============================================================================
# DECOMPOSITION PARSING
# =============================================================================

def parse_decomp_output(
        decomp_file: Union[str, Path],
        receptor_pdb: Union[str, Path],
        output_dir: Union[str, Path],
        is_single_frame: bool = True,
) -> Dict[str, Any]:
    """
    Parse FINAL_DECOMP_MMPBSA.dat → per-residue CSV with PDB numbering.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    text = Path(decomp_file).read_text()

    residue_mapping = _build_residue_mapping(receptor_pdb)
    rows = _parse_decomp_sections(text, is_single_frame)

    if not rows:
        return {"success": False, "error": "No residue data found in decomp file"}

    df = pd.DataFrame(rows)

    if residue_mapping:
        df = _apply_pdb_numbering(df, residue_mapping)

    if "total" in df.columns:
        df = df.sort_values("total", ascending=True).reset_index(drop=True)

    csv_path = output_dir / "per_residue_decomp.csv"
    df.to_csv(csv_path, index=False, encoding="utf-8")
    logger.info(f"  Saved: {csv_path} ({len(df)} residues)")

    n_favorable = len(df[df["total"] < -0.5]) if "total" in df.columns else 0
    n_unfavorable = len(df[df["total"] > 0.5]) if "total" in df.columns else 0

    return {
        "success": True,
        "csv_path": str(csv_path),
        "df": df,
        "n_residues": len(df),
        "n_favorable": n_favorable,
        "n_unfavorable": n_unfavorable,
    }


def _parse_decomp_sections(text: str, is_single_frame: bool) -> List[Dict]:
    """
    Parse per-residue data from FINAL_DECOMP_MMPBSA.dat.

    We want the DELTAS > Total Energy Decomposition section.
    """
    rows = []
    lines = text.split("\n")

    in_deltas = False
    in_total_decomp = False

    for i, line in enumerate(lines):
        line_stripped = line.strip()

        if line_stripped == "DELTAS:":
            in_deltas = True
            continue

        if not in_deltas:
            continue

        if "Total Energy Decomposition" in line_stripped:
            in_total_decomp = True
            continue

        if not in_total_decomp:
            continue

        if line_stripped in ("Sidechain Energy Decomposition:",
                            "Backbone Energy Decomposition:"):
            break

        if line_stripped.startswith("Residue,") or line_stripped.startswith(",,"):
            continue
        if not line_stripped:
            continue

        parsed = _parse_decomp_csv_line(line_stripped, is_single_frame)
        if parsed:
            rows.append(parsed)

    return rows


def _parse_decomp_csv_line(line: str, is_single_frame: bool) -> Optional[Dict]:
    """Parse a single CSV residue line from MMPBSA.py decomp output."""
    parts = line.split(",")
    if len(parts) < 18:
        return None

    res_id = parts[0].strip()
    location = parts[1].strip()

    # Skip ligand residues
    if location.startswith("L "):
        return None

    res_tokens = res_id.split()
    if len(res_tokens) < 2:
        return None

    res_name = res_tokens[0]
    try:
        res_num_seq = int(res_tokens[1])
    except ValueError:
        return None

    try:
        internal = float(parts[2])
        vdw = float(parts[5])
        es = float(parts[8])
        gb = float(parts[11])
        sa = float(parts[14])
        total = float(parts[17])
    except (ValueError, IndexError):
        return None

    result = {
        "residue_name": res_name,
        "residue_number_seq": res_num_seq,
        "chain_type": "R",
        "internal": internal,
        "vdw": vdw,
        "es": es,
        "gb": gb,
        "sa": sa,
        "total": total,
    }

    if not is_single_frame:
        try:
            result["vdw_std"] = float(parts[6])
            result["es_std"] = float(parts[9])
            result["gb_std"] = float(parts[12])
            result["sa_std"] = float(parts[15])
            result["total_std"] = float(parts[18])
        except (ValueError, IndexError):
            pass

    return result


# =============================================================================
# RESIDUE NUMBERING
# =============================================================================

def _build_residue_mapping(
        receptor_pdb: Union[str, Path],
) -> Dict[int, Dict[str, Any]]:
    """Build sequential→PDB residue mapping from receptor PDB."""
    pdb_path = Path(receptor_pdb)
    if not pdb_path.exists():
        logger.warning(f"  Receptor PDB not found for mapping: {pdb_path}")
        return {}

    mapping = {}
    seen_residues = []
    prev_resid = None

    for line in pdb_path.read_text().split("\n"):
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue

        res_name = line[17:20].strip()
        chain = line[21:22].strip() or "A"
        try:
            res_num = int(line[22:26].strip())
        except ValueError:
            continue

        resid = (res_name, res_num, chain)
        if resid != prev_resid:
            seen_residues.append(resid)
            prev_resid = resid

    for i, (name, number, chain) in enumerate(seen_residues, 1):
        mapping[i] = {
            "name": name,
            "number": number,
            "chain": chain,
            "residue_id": f"{name}{number}.{chain}",
        }

    logger.info(f"  Residue mapping: {len(mapping)} residues (PDB numbering)")
    return mapping


def _apply_pdb_numbering(
        df: pd.DataFrame,
        mapping: Dict[int, Dict],
) -> pd.DataFrame:
    """Apply PDB residue numbering to decomp DataFrame."""
    if "residue_number_seq" not in df.columns:
        return df

    pdb_names = []
    pdb_numbers = []
    pdb_chains = []
    pdb_ids = []

    for _, row in df.iterrows():
        seq_num = row.get("residue_number_seq", 0)
        if seq_num in mapping:
            m = mapping[seq_num]
            pdb_names.append(m["name"])
            pdb_numbers.append(m["number"])
            pdb_chains.append(m["chain"])
            pdb_ids.append(m["residue_id"])
        else:
            pdb_names.append(row.get("residue_name", "UNK"))
            pdb_numbers.append(seq_num)
            pdb_chains.append("?")
            pdb_ids.append(f"{row.get('residue_name', 'UNK')}{seq_num}.?")

    df["residue_name_pdb"] = pdb_names
    df["residue_number_pdb"] = pdb_numbers
    df["chain"] = pdb_chains
    df["residue_id"] = pdb_ids

    return df


# =============================================================================
# FOOTPRINT COMPARISON
# =============================================================================

def compare_with_footprint(
        decomp_df: pd.DataFrame,
        footprint_csv: Union[str, Path],
        output_dir: Union[str, Path],
) -> Dict[str, Any]:
    """
    Compare MMPBSA per-residue decomp with 01d/04b footprint in-vacuo.
    """
    output_dir = Path(output_dir)
    fp_path = Path(footprint_csv)

    if not fp_path.exists():
        return {"success": False, "error": f"Footprint CSV not found: {fp_path}"}

    df_fp = pd.read_csv(fp_path)
    logger.info(f"  Footprint: {len(df_fp)} residues from 04b")

    if "residue_id" not in decomp_df.columns or "residue_id" not in df_fp.columns:
        return {"success": False, "error": "Missing residue_id column for merge"}

    df_merged = pd.merge(
        decomp_df,
        df_fp[["residue_id", "residue_name", "mean_vdw", "mean_es", "mean_total"]],
        on="residue_id",
        how="outer",
        suffixes=("_mmpbsa", "_fp"),
    )

    rename_map = {
        "mean_vdw": "fp_vdw",
        "mean_es": "fp_es",
        "mean_total": "fp_total",
    }
    df_merged = df_merged.rename(columns=rename_map)

    if "vdw" in df_merged.columns and "fp_vdw" in df_merged.columns:
        df_merged["delta_vdw"] = df_merged["vdw"] - df_merged["fp_vdw"]
    if "es" in df_merged.columns and "fp_es" in df_merged.columns:
        df_merged["delta_es"] = df_merged["es"] - df_merged["fp_es"]

    if "vdw" in df_merged.columns and "es" in df_merged.columns:
        df_merged["mmpbsa_invacuo"] = df_merged["vdw"] + df_merged["es"]

    if "gb" in df_merged.columns and "sa" in df_merged.columns:
        df_merged["solvation"] = df_merged["gb"] + df_merged["sa"]

    if "total" in df_merged.columns:
        df_merged["rank_mmpbsa"] = df_merged["total"].rank(ascending=True)
    if "fp_total" in df_merged.columns:
        df_merged["rank_fp"] = df_merged["fp_total"].rank(ascending=True)
    if "rank_mmpbsa" in df_merged.columns and "rank_fp" in df_merged.columns:
        df_merged["rank_delta"] = df_merged["rank_fp"] - df_merged["rank_mmpbsa"]

    df_merged = df_merged.sort_values("total", ascending=True, na_position="last")

    comp_csv = output_dir / "comparison_invacuo_vs_solvated.csv"
    df_merged.to_csv(comp_csv, index=False, encoding="utf-8")
    logger.info(f"  Saved: {comp_csv} ({len(df_merged)} residues)")

    zone_summary = _compute_zone_summary(df_merged)

    return {
        "success": True,
        "comparison_csv": str(comp_csv),
        "df": df_merged,
        "zone_summary": zone_summary,
    }


def _compute_zone_summary(df: pd.DataFrame) -> Dict[str, Dict]:
    """Compute energy summary by binding site zone."""
    summary = {}

    for zone_name, zdef in ZONE_DEFINITIONS.items():
        mask = df["residue_id"].apply(
            lambda rid: str(rid).split(".")[0] in zdef["residues"]
            if pd.notna(rid) else False
        )
        zone_df = df[mask]

        if zone_df.empty:
            continue

        zone_info = {
            "label": zdef["label"],
            "n_residues": len(zone_df),
            "residues": zone_df["residue_id"].tolist(),
        }

        for col in ["vdw", "es", "gb", "sa", "total", "fp_total", "solvation"]:
            if col in zone_df.columns:
                vals = zone_df[col].dropna()
                if len(vals) > 0:
                    zone_info[f"sum_{col}"] = round(vals.sum(), 2)
                    zone_info[f"mean_{col}"] = round(vals.mean(), 2)

        if "solvation" in zone_df.columns:
            zone_info["solvation_impact"] = round(
                zone_df["solvation"].dropna().sum(), 2
            )

        summary[zone_name] = zone_info

    return summary


# =============================================================================
# HTML REPORT
# =============================================================================

def generate_mmpbsa_html(
        decomp_df: pd.DataFrame,
        global_results: Dict[str, float],
        zone_summary: Optional[Dict] = None,
        comparison_df: Optional[pd.DataFrame] = None,
        is_single_frame: bool = True,
        campaign_id: str = "",
        molecule_name: str = "",
) -> str:
    """Generate MMPBSA per-residue decomposition HTML report."""
    from datetime import datetime

    mode_label = "Single-point (1 frame)" if is_single_frame else "MD trajectory (mean ± std)"
    dg = global_results.get("delta_total", None)
    dg_str = f"{dg:+.2f}" if dg is not None else "N/A"

    zone_colors = {
        "phosphate": "#1D9E75", "xylose": "#378ADD",
        "uracil": "#BA7517", "ribose_bridge": "#7F77DD",
    }

    df_sig = decomp_df[decomp_df["total"].abs() > 0.5].copy() if "total" in decomp_df.columns else decomp_df.head(0)
    df_fav = df_sig[df_sig["total"] < 0].sort_values("total")
    df_unfav = df_sig[df_sig["total"] > 0].sort_values("total", ascending=False)

    max_abs = max(df_sig["total"].abs().max(), 1) if len(df_sig) > 0 else 1

    title_extra = f" — {molecule_name}" if molecule_name else ""
    html = f"""\
<style>
.tbl{{width:100%;border-collapse:collapse;font-size:13px}}
.tbl th{{text-align:left;padding:6px 8px;border-bottom:0.5px solid #dee2e6;color:#666;font-weight:500;font-size:12px}}
.tbl td{{padding:5px 8px;border-bottom:0.5px solid #f0f0f0;vertical-align:top}}
.m{{font-family:monospace;font-size:12px}}
.bar{{height:8px;border-radius:4px;display:inline-block;vertical-align:middle}}
.sect{{font-size:11px;font-weight:500;color:#666;padding:10px 8px 4px;text-transform:uppercase;letter-spacing:0.5px}}
.card{{background:#f8f9fa;border-radius:8px;padding:14px;margin:8px 4px 8px 0;display:inline-block;text-align:center;min-width:110px;vertical-align:top}}
.card-val{{font-size:20px;font-weight:bold}}
.card-lbl{{font-size:10px;color:#666;margin-top:2px}}
.neg{{color:#27ae60}}
.pos{{color:#e74c3c}}
.neut{{color:#999}}
</style>

<div style="padding:0.5rem 0">

<div style="font-size:11px;color:#666;margin-bottom:4px">
MMPBSA Per-Residue Decomposition{title_extra} — {campaign_id} | {mode_label} | {datetime.now().strftime('%Y-%m-%d %H:%M')}
</div>
"""

    vdw = global_results.get("vdw", 0)
    eel = global_results.get("eel", 0)
    egb = global_results.get("egb", 0)
    esurf = global_results.get("esurf", 0)

    def _val_class(v):
        return "neg" if v < -1 else "pos" if v > 1 else "neut"

    html += f"""
<div style="margin:12px 0">
  <div class="card"><div class="card-val {_val_class(dg if dg else 0)}">{dg_str}</div><div class="card-lbl">ΔG bind (kcal/mol)</div></div>
  <div class="card"><div class="card-val {_val_class(vdw)}">{vdw:+.1f}</div><div class="card-lbl">vdW</div></div>
  <div class="card"><div class="card-val {_val_class(eel)}">{eel:+.1f}</div><div class="card-lbl">Electrostatic</div></div>
  <div class="card"><div class="card-val {_val_class(egb)}">{egb:+.1f}</div><div class="card-lbl">GB solvation</div></div>
  <div class="card"><div class="card-val {_val_class(esurf)}">{esurf:+.1f}</div><div class="card-lbl">SA (non-polar)</div></div>
</div>
"""

    # Zone comparison
    if zone_summary:
        html += """
<table class="tbl" style="margin-top:16px">
<tr><th>Zone</th><th>MMPBSA total</th><th>Solvation (GB+SA)</th><th></th></tr>
"""
        for zone_id, info in sorted(zone_summary.items(), key=lambda x: x[1].get("sum_total", 0)):
            total = info.get("sum_total", 0)
            solv = info.get("solvation_impact", 0)
            label = info.get("label", zone_id)
            color = zone_colors.get(zone_id, "#999")
            n_res = info.get("n_residues", 0)

            bar_w = min(80, int(abs(total) / max_abs * 80))
            bar_color = "#2ecc71" if total < 0 else "#e74c3c"

            html += f"""<tr>
<td><span style="display:inline-block;width:10px;height:10px;border-radius:50%;background:{color};margin-right:4px;vertical-align:middle"></span><b>{label}</b> <span style="font-size:11px;color:#999">({n_res} res)</span></td>
<td class="m" style="font-weight:500">{total:+.2f}</td>
<td class="m">{solv:+.2f}</td>
<td><span class="bar" style="width:{bar_w}px;background:{bar_color}"></span></td>
</tr>"""

        html += "</table>"

    # Top favorable residues
    html += """
<table class="tbl" style="margin-top:16px">
<tr><th colspan="7" class="sect">Top favorable residues (binding)</th></tr>
<tr><th>#</th><th>Residue</th><th>vdW</th><th>ES</th><th>GB</th><th>SA</th><th>TOTAL</th></tr>
"""
    for rank, (_, r) in enumerate(df_fav.head(15).iterrows(), 1):
        rid = r.get("residue_id", f"{r['residue_name']}{r['residue_number_seq']}")
        bar_w = min(80, int(abs(r["total"]) / max_abs * 80))
        html += f"""<tr>
<td class="m">{rank}</td>
<td><b>{rid}</b></td>
<td class="m">{r['vdw']:+.2f}</td>
<td class="m">{r['es']:+.2f}</td>
<td class="m">{r['gb']:+.2f}</td>
<td class="m">{r['sa']:+.4f}</td>
<td class="m" style="font-weight:500">{r['total']:+.2f} <span class="bar" style="width:{bar_w}px;background:#2ecc71;margin-left:4px"></span></td>
</tr>"""

    html += "</table>"

    # Unfavorable residues
    if len(df_unfav) > 0:
        html += """
<table class="tbl" style="margin-top:16px">
<tr><th colspan="7" class="sect">Unfavorable residues (repulsive)</th></tr>
<tr><th>#</th><th>Residue</th><th>vdW</th><th>ES</th><th>GB</th><th>SA</th><th>TOTAL</th></tr>
"""
        for rank, (_, r) in enumerate(df_unfav.head(10).iterrows(), 1):
            rid = r.get("residue_id", f"{r['residue_name']}{r['residue_number_seq']}")
            bar_w = min(80, int(abs(r["total"]) / max_abs * 80))
            html += f"""<tr>
<td class="m">{rank}</td>
<td><b>{rid}</b></td>
<td class="m">{r['vdw']:+.2f}</td>
<td class="m">{r['es']:+.2f}</td>
<td class="m">{r['gb']:+.2f}</td>
<td class="m">{r['sa']:+.4f}</td>
<td class="m" style="font-weight:500;color:#e74c3c">{r['total']:+.2f} <span class="bar" style="width:{bar_w}px;background:#e74c3c;margin-left:4px"></span></td>
</tr>"""

        html += "</table>"

    html += f"""
<div style="margin-top:16px;font-size:11px;color:#999">
MMPBSA.py 14.0 | idecomp=1 (per-residue) | igb=2 (OBC-II) | saltcon=0.15 M | AMBER ff14SB + GAFF2 | {mode_label}
</div>
</div>"""

    return html


# =============================================================================
# SINGLE-MOLECULE ANALYSIS PIPELINE
# =============================================================================

def run_mmpbsa_analysis(
        decomp_file: Union[str, Path],
        results_file: Union[str, Path],
        receptor_pdb: Union[str, Path],
        output_dir: Union[str, Path],
        is_single_frame: bool = True,
        compare_footprint: bool = True,
        footprint_csv: Optional[str] = None,
        campaign_id: str = "",
        molecule_name: str = "",
) -> Dict[str, Any]:
    """
    Run complete MMPBSA analysis pipeline (steps 6-9) for a single molecule.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info(f"  MMPBSA Analysis (01h){' — ' + molecule_name if molecule_name else ''}")
    logger.info("=" * 60)
    logger.info(f"  Decomp file: {Path(decomp_file).name}")
    logger.info(f"  Mode:        {'single_point' if is_single_frame else 'md (multi-frame)'}")

    # ─── Step 6: Parse global results ───
    logger.info("")
    logger.info("─── Step 6: Parse Global MMPBSA Results ───")

    global_results = {}
    if Path(results_file).exists():
        global_results = parse_mmpbsa_global(results_file)
        if global_results:
            logger.info(f"  ΔG_bind:  {global_results.get('delta_total', 'N/A')} kcal/mol")
            logger.info(f"    vdW:    {global_results.get('vdw', 'N/A')}")
            logger.info(f"    EEL:    {global_results.get('eel', 'N/A')}")
            logger.info(f"    EGB:    {global_results.get('egb', 'N/A')}")
            logger.info(f"    ESURF:  {global_results.get('esurf', 'N/A')}")

    global_json = output_dir / "global_binding_energy.json"
    with open(global_json, "w") as f:
        json.dump(global_results, f, indent=2)

    # ─── Step 7: Parse per-residue decomposition ───
    logger.info("")
    logger.info("─── Step 7: Parse Per-Residue Decomposition ───")

    result_parse = parse_decomp_output(
        decomp_file=decomp_file,
        receptor_pdb=receptor_pdb,
        output_dir=output_dir,
        is_single_frame=is_single_frame,
    )
    if not result_parse["success"]:
        return {"success": False, "error": f"Parse: {result_parse['error']}"}

    logger.info(f"  Residues: {result_parse['n_residues']} total, "
                f"{result_parse['n_favorable']} favorable (<-0.5), "
                f"{result_parse['n_unfavorable']} unfavorable (>+0.5)")

    # ─── Step 8: Compare with footprint ───
    comparison_result = None
    if compare_footprint and footprint_csv:
        logger.info("")
        logger.info("─── Step 8: Footprint Comparison ───")

        comparison_result = compare_with_footprint(
            decomp_df=result_parse["df"],
            footprint_csv=footprint_csv,
            output_dir=output_dir,
        )
        if comparison_result["success"]:
            zs = comparison_result.get("zone_summary", {})
            if zs:
                logger.info("  Zone Summary (MMPBSA total | solvation impact):")
                for zone, info in zs.items():
                    total = info.get("sum_total", 0)
                    solv = info.get("solvation_impact", 0)
                    logger.info(f"    {info.get('label', zone):30s}: "
                                f"{total:8.2f} | {solv:+8.2f} kcal/mol")

            zone_json = output_dir / "zone_summary.json"
            with open(zone_json, "w") as f:
                json.dump(zs, f, indent=2, default=str)

    elif compare_footprint and not footprint_csv:
        logger.info("  Footprint comparison: skipped (no CSV path provided)")

    # ─── Step 9: Generate HTML report ───
    html_path = None
    try:
        logger.info("")
        logger.info("─── Step 9: HTML Report ───")

        html = generate_mmpbsa_html(
            decomp_df=result_parse["df"],
            global_results=global_results,
            zone_summary=(comparison_result.get("zone_summary")
                          if comparison_result and comparison_result.get("success")
                          else None),
            comparison_df=(comparison_result.get("df")
                           if comparison_result and comparison_result.get("success")
                           else None),
            is_single_frame=is_single_frame,
            campaign_id=campaign_id,
            molecule_name=molecule_name,
        )

        html_path = output_dir / "mmpbsa_report.html"
        html_path.write_text(html, encoding="utf-8")
        logger.info(f"  Report: {html_path}")
    except Exception as e:
        logger.warning(f"  HTML generation failed: {e}")

    logger.info("")
    logger.info("=" * 60)
    logger.info(f"  MMPBSA Analysis Complete")
    logger.info("=" * 60)

    return {
        "success": True,
        "per_residue_csv": result_parse["csv_path"],
        "comparison_csv": (comparison_result["comparison_csv"]
                           if comparison_result and comparison_result.get("success")
                           else None),
        "html_report": str(html_path) if html_path else None,
        "global_results": global_results,
        "zone_summary": (comparison_result.get("zone_summary")
                         if comparison_result and comparison_result.get("success")
                         else None),
        "n_residues": result_parse["n_residues"],
        "n_favorable": result_parse["n_favorable"],
        "n_unfavorable": result_parse["n_unfavorable"],
    }


# =============================================================================
# BATCH ANALYSIS PIPELINE (hit_validation specific)
# =============================================================================

def run_mmpbsa_batch_analysis(
        mmpbsa_results_dir: Union[str, Path],
        receptor_pdb: Union[str, Path],
        output_dir: Union[str, Path],
        reference_context: Optional[Dict] = None,
        compare_footprint: bool = True,
        footprint_dir: Optional[str] = None,
        campaign_id: str = "",
) -> Dict[str, Any]:
    """
    Parse and analyze MMPBSA results for all molecules.

    Produces per-molecule CSVs + consolidated comparison tables.

    Args:
        mmpbsa_results_dir: Path to 05_results/{campaign}/01g_mmpbsa_decomp/
        receptor_pdb:       Receptor PDB for residue mapping
        output_dir:         Output directory
        reference_context:  Optional UDX reference data from campaign_config
        compare_footprint:  Compare with 01d footprint
        footprint_dir:      Path to 04b footprint analysis output
        campaign_id:        Campaign identifier

    Returns:
        Dict with success, per-molecule results, consolidated outputs
    """
    mmpbsa_results_dir = Path(mmpbsa_results_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Discover molecules
    molecules = sorted([
        d.name for d in mmpbsa_results_dir.iterdir()
        if d.is_dir() and (d / "mmpbsa" / "FINAL_DECOMP_MMPBSA.dat").exists()
    ])

    if not molecules:
        return {"success": False, "error": "No MMPBSA results found"}

    logger.info("=" * 60)
    logger.info("  MMPBSA Batch Analysis (01h)")
    logger.info("=" * 60)
    logger.info(f"  Molecules:  {len(molecules)}")
    logger.info(f"  Output:     {output_dir}")

    all_results = {}
    all_global_energies = []
    all_zone_summaries = []

    for mol_name in molecules:
        decomp_file = mmpbsa_results_dir / mol_name / "mmpbsa" / "FINAL_DECOMP_MMPBSA.dat"
        results_file = mmpbsa_results_dir / mol_name / "mmpbsa" / "FINAL_RESULTS_MMPBSA.dat"

        if not decomp_file.exists():
            logger.warning(f"  Skipping {mol_name}: no DECOMP file")
            continue

        # Detect single_point vs md from pipeline log or DECOMP stddev
        is_single_frame = True
        pipeline_log_path = mmpbsa_results_dir / mol_name / "01g_pipeline_log.json"
        if pipeline_log_path.exists():
            with open(pipeline_log_path) as f:
                plog = json.load(f)
            is_single_frame = (plog.get("mode") == "single_point")

        # Auto-detect from stddev values
        if is_single_frame:
            with open(decomp_file) as f:
                decomp_text = f.read()
            in_deltas = False
            for line in decomp_text.split("\n"):
                if "DELTAS:" in line:
                    in_deltas = True
                    continue
                if in_deltas and line.strip() and not line.startswith(",") and not line.startswith("Residue"):
                    parts = line.split(",")
                    if len(parts) >= 4:
                        try:
                            std_dev = float(parts[3])
                            if std_dev > 0.001:
                                is_single_frame = False
                                break
                        except (ValueError, IndexError):
                            continue

        # Resolve per-molecule footprint CSV
        fp_csv = None
        if compare_footprint and footprint_dir:
            # Try per-molecule footprint from 01d
            for candidate in [
                Path(footprint_dir) / "residue_consensus.csv",
                Path(footprint_dir) / mol_name / "residue_consensus.csv",
            ]:
                if candidate.exists():
                    fp_csv = str(candidate)
                    break

        mol_output = output_dir / mol_name
        result = run_mmpbsa_analysis(
            decomp_file=str(decomp_file),
            results_file=str(results_file),
            receptor_pdb=str(receptor_pdb),
            output_dir=str(mol_output),
            is_single_frame=is_single_frame,
            compare_footprint=compare_footprint,
            footprint_csv=fp_csv,
            campaign_id=campaign_id,
            molecule_name=mol_name,
        )
        all_results[mol_name] = result

        # Collect global energies for consolidated output
        if result.get("success") and result.get("global_results"):
            gr = result["global_results"]
            all_global_energies.append({
                "Name": mol_name,
                "MMPBSA_dG_total": gr.get("delta_total"),
                "MMPBSA_vdW": gr.get("vdw"),
                "MMPBSA_EEL": gr.get("eel"),
                "MMPBSA_EGB": gr.get("egb"),
                "MMPBSA_ESURF": gr.get("esurf"),
            })

        # Collect zone summaries
        if result.get("success") and result.get("zone_summary"):
            zs = result["zone_summary"]
            zone_row = {"Name": mol_name}
            for zone_id, info in zs.items():
                zone_row[f"{zone_id}_total"] = info.get("sum_total", 0)
                zone_row[f"{zone_id}_solvation"] = info.get("solvation_impact", 0)
            all_zone_summaries.append(zone_row)

    # ─── Consolidated outputs ───

    # 1. consolidated_mmpbsa.csv
    consolidated_csv = None
    if all_global_energies:
        df_consolidated = pd.DataFrame(all_global_energies)
        df_consolidated = df_consolidated.sort_values(
            "MMPBSA_dG_total", ascending=True, na_position="last"
        ).reset_index(drop=True)
        consolidated_csv = output_dir / "consolidated_mmpbsa.csv"
        df_consolidated.to_csv(consolidated_csv, index=False, encoding="utf-8")
        logger.info(f"  Saved: {consolidated_csv}")

    # 2. zone_comparison.csv
    zone_comparison_csv = None
    if all_zone_summaries:
        df_zones = pd.DataFrame(all_zone_summaries)
        zone_comparison_csv = output_dir / "zone_comparison.csv"
        df_zones.to_csv(zone_comparison_csv, index=False, encoding="utf-8")
        logger.info(f"  Saved: {zone_comparison_csv}")

    # 3. mmpbsa_batch_report.html
    html_path = None
    try:
        html_parts = [
            "<!DOCTYPE html><html><head><meta charset='utf-8'>",
            f"<title>MMPBSA Batch Report: {campaign_id}</title>",
            "<style>body{font-family:'Segoe UI',Arial,sans-serif;max-width:950px;margin:0 auto;padding:20px;background:#fafafa;color:#333}",
            "h1{color:#1a5276;border-bottom:3px solid #1a5276;padding-bottom:10px}",
            ".mol-section{border:1px solid #ddd;border-radius:8px;padding:15px;margin:15px 0;background:white}",
            "table{width:100%;border-collapse:collapse;font-size:13px}",
            "th{background:#f8f9fa;padding:8px;text-align:left;border-bottom:2px solid #dee2e6}",
            "td{padding:6px 8px;border-bottom:1px solid #f0f0f0}",
            ".mono{font-family:monospace}</style></head><body>",
            f"<h1>MMPBSA Batch Analysis — {campaign_id}</h1>",
            f"<p>{len(molecules)} molecules analyzed</p>",
        ]

        # Summary table
        if all_global_energies:
            html_parts.append("<h2>Global Binding Energies</h2><table>")
            html_parts.append("<tr><th>Rank</th><th>Name</th><th>ΔG total</th><th>vdW</th><th>EEL</th><th>EGB</th><th>ESURF</th></tr>")
            sorted_energies = sorted(all_global_energies,
                                     key=lambda x: x.get("MMPBSA_dG_total") or 0)
            for rank, e in enumerate(sorted_energies, 1):
                dg = e.get("MMPBSA_dG_total")
                dg_str = f"{dg:+.2f}" if dg is not None else "N/A"
                html_parts.append(f"<tr><td>{rank}</td><td><b>{e['Name']}</b></td>"
                                  f"<td class='mono'>{dg_str}</td>"
                                  f"<td class='mono'>{e.get('MMPBSA_vdW', 0):+.1f}</td>"
                                  f"<td class='mono'>{e.get('MMPBSA_EEL', 0):+.1f}</td>"
                                  f"<td class='mono'>{e.get('MMPBSA_EGB', 0):+.1f}</td>"
                                  f"<td class='mono'>{e.get('MMPBSA_ESURF', 0):+.4f}</td></tr>")
            html_parts.append("</table>")

        html_parts.append("</body></html>")
        html_content = "\n".join(html_parts)
        html_path = output_dir / "mmpbsa_batch_report.html"
        html_path.write_text(html_content, encoding="utf-8")
        logger.info(f"  Saved: {html_path}")
    except Exception as e:
        logger.warning(f"  Batch HTML generation failed: {e}")

    n_success = sum(1 for r in all_results.values() if r.get("success"))
    n_failed = sum(1 for r in all_results.values() if not r.get("success"))

    logger.info("")
    logger.info("=" * 60)
    logger.info(f"  MMPBSA Batch Analysis Complete: {n_success}/{len(molecules)}")
    logger.info("=" * 60)

    return {
        "success": True,
        "results": all_results,
        "n_total": len(molecules),
        "n_success": n_success,
        "n_failed": n_failed,
        "consolidated_csv": str(consolidated_csv) if consolidated_csv else None,
        "zone_comparison_csv": str(zone_comparison_csv) if zone_comparison_csv else None,
        "html_report": str(html_path) if html_path else None,
    }
