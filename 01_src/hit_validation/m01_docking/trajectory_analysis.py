"""
Trajectory Analysis — Core Module (01i)
==========================================
Extracts per-frame data from MD trajectories using cpptraj and ProLIF.
Produces raw per-frame data that 07c and 07d consume.

Two analysis engines:
  cpptraj: RMSD, distances, H-bonds, water bridges, RMSF
  ProLIF:  interaction fingerprints per-frame (14 types x residues)

Input:
  01g output: production.dcd, solvated.prmtop, complex.prmtop, dry_trajectory.mdcrd
  campaign_config.yaml: zones (for targeted distance monitoring)
  04b output: residue_mapping.csv (PDB->sequential numbering for cpptraj)

Output per molecule:
  rmsd_ligand.csv           -- RMSD of ligand vs frame 1
  rmsd_receptor.csv         -- RMSD of receptor backbone
  distances.csv             -- min distance to each monitored residue per frame
  hbond_occupancy.csv       -- H-bond donor/acceptor/occupancy
  water_bridges.csv         -- water-mediated contacts
  rmsf_per_residue.csv      -- RMSF per residue (Ca)
  prolif_fingerprints.pkl   -- ProLIF fingerprint object (for 07d)
  prolif_occupancy.csv      -- per-residue per-interaction-type occupancy

Location: 01_src/hit_validation/m01_docking/trajectory_analysis.py
Project: hit_validation
Module: 01i (core)
Version: 1.0
"""

import logging
import re
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# =============================================================================
# LIGAND AUTO-DETECTION
# =============================================================================

def detect_ligand_resname(complex_prmtop: Union[str, Path]) -> Optional[str]:
    """
    Auto-detect ligand residue name from AMBER prmtop.

    The ligand is the last residue in the complex topology (after all protein
    residues). Uses parmed to read the topology directly.

    Returns resname (e.g. "UNL", "MOL") or None if detection fails.
    """
    try:
        import parmed
        parm = parmed.load_file(str(complex_prmtop))
        return parm.residues[-1].name
    except Exception as e:
        logger.warning(f"  Ligand detection failed: {e}")
        return None


# =============================================================================
# CPPTRAJ SCRIPT TEMPLATES
# =============================================================================

CPPTRAJ_STRIPPED_TEMPLATE = """\
# 01i Trajectory Analysis -- stripped trajectory (no waters)
parm {complex_prmtop}
trajin {dry_trajectory}

# 1. RMSD ligand (align by protein backbone, measure ligand heavy atoms)
rms fit_backbone @CA
rms RMSD_lig :{lig_mask}&!@H= out {output_dir}/rmsd_ligand.dat mass

# 2. RMSD receptor backbone
rms RMSD_rec @CA out {output_dir}/rmsd_receptor.dat mass

# 3. Key distances -- minimum distance between ligand and each monitored residue
{distance_commands}

# 4. H-bond analysis -- ligand to protein
hbond LigHB donormask :{lig_mask} acceptormask :1-{n_protein_res} out {output_dir}/hbond_timeseries.dat avgout {output_dir}/hbond_avg.dat
hbond ProtHB donormask :1-{n_protein_res} acceptormask :{lig_mask} out {output_dir}/hbond_timeseries_rev.dat avgout {output_dir}/hbond_avg_rev.dat

# 5. RMSF per residue (Ca)
atomicfluct out {output_dir}/rmsf_per_residue.dat @CA byres

run
quit
"""

CPPTRAJ_SOLVATED_TEMPLATE = """\
# 01i Water Bridge Analysis -- solvated trajectory
parm {solvated_prmtop}
trajin {production_dcd}

# Strip ions, keep water
strip :Na+,Cl-

# Water bridges: waters that simultaneously H-bond ligand and protein
hbond WB donormask :{lig_mask} acceptormask :1-{n_protein_res} \
  solventdonor :WAT solventacceptor :WAT \
  bridgeout {output_dir}/water_bridges.dat

run
quit
"""

DISTANCE_CMD_TEMPLATE = "distance d_{res_id} :{lig_mask}&!@H= :{seq_resnum}&!@H= out {output_dir}/dist_{res_id}.dat"


# =============================================================================
# RESIDUE MAPPING
# =============================================================================

def load_residue_mapping(mapping_csv: Union[str, Path]) -> Dict[str, str]:
    """
    Load PDB->sequential residue mapping from 04b output.

    Returns:
        Dict mapping PDB id (e.g., "GLU529.A") to sequential number (e.g., "278")
    """
    df = pd.read_csv(mapping_csv)
    mapping = {}
    for _, row in df.iterrows():
        pdb_id = row["pdb_original"]  # e.g., "GLU529.A"
        seq_id = row["mol2_sequential"]  # e.g., "GLH278"
        seq_num = re.search(r"(\d+)", seq_id)
        if seq_num:
            mapping[pdb_id] = seq_num.group(1)
    return mapping


def get_monitored_residues(
        campaign_config: Dict,
        residue_mapping: Dict[str, str],
) -> List[Dict]:
    """
    Build list of residues to monitor distances, from:
    1. All residues in zones defined in campaign_config
    2. All residues in functional_residues

    Returns list of {pdb_id, seq_num, label}
    """
    monitored = {}

    # From zones
    zones = campaign_config.get("zones", {})
    for zone_id, zone_def in zones.items():
        for res_name in zone_def.get("residues", []):
            for pdb_id, seq_num in residue_mapping.items():
                if pdb_id.startswith(res_name):
                    monitored[pdb_id] = {
                        "pdb_id": pdb_id,
                        "seq_num": seq_num,
                        "label": zone_def.get("label", zone_id),
                    }

    # From functional_residues
    func_res = campaign_config.get("functional_residues", {})
    for res_name, info in func_res.items():
        for pdb_id, seq_num in residue_mapping.items():
            if pdb_id.startswith(res_name):
                if pdb_id not in monitored:
                    monitored[pdb_id] = {
                        "pdb_id": pdb_id,
                        "seq_num": seq_num,
                        "label": info.get("function", "functional"),
                    }

    return list(monitored.values())


# =============================================================================
# CPPTRAJ EXECUTION
# =============================================================================

def generate_cpptraj_stripped_script(
        complex_prmtop: str,
        dry_trajectory: str,
        output_dir: str,
        monitored_residues: List[Dict],
        lig_mask: str = "LIG",
        n_protein_res: int = 706,
) -> str:
    """Generate cpptraj input script for stripped trajectory analysis."""
    dist_commands = []
    for res in monitored_residues:
        cmd = DISTANCE_CMD_TEMPLATE.format(
            lig_mask=lig_mask,
            seq_resnum=res["seq_num"],
            res_id=res["pdb_id"].replace(".", "_"),
            output_dir=output_dir,
        )
        dist_commands.append(cmd)

    script = CPPTRAJ_STRIPPED_TEMPLATE.format(
        complex_prmtop=complex_prmtop,
        dry_trajectory=dry_trajectory,
        output_dir=output_dir,
        lig_mask=lig_mask,
        n_protein_res=n_protein_res,
        distance_commands="\n".join(dist_commands),
    )
    return script


def generate_cpptraj_solvated_script(
        solvated_prmtop: str,
        production_dcd: str,
        output_dir: str,
        lig_mask: str = "LIG",
        n_protein_res: int = 706,
) -> str:
    """Generate cpptraj input script for water bridge analysis."""
    script = CPPTRAJ_SOLVATED_TEMPLATE.format(
        solvated_prmtop=solvated_prmtop,
        production_dcd=production_dcd,
        output_dir=output_dir,
        lig_mask=lig_mask,
        n_protein_res=n_protein_res,
    )
    return script


def run_cpptraj(script_content: str, script_path: Path, timeout: int = 600) -> Dict:
    """Execute a cpptraj script."""
    script_path.write_text(script_content)

    logger.info(f"  Running cpptraj: {script_path.name}")
    try:
        proc = subprocess.run(
            ["cpptraj", "-i", script_path.name],
            capture_output=True, text=True, timeout=timeout,
            cwd=str(script_path.parent),
        )
    except subprocess.TimeoutExpired:
        logger.warning(f"  cpptraj timed out after {timeout}s")
        return {"success": False, "error": f"Timeout after {timeout}s"}

    if proc.returncode != 0:
        logger.error(f"  cpptraj failed: {proc.stderr[:500]}")
        return {"success": False, "error": proc.stderr[:500]}

    return {"success": True, "stdout": proc.stdout}


# =============================================================================
# CPPTRAJ OUTPUT PARSING
# =============================================================================

def parse_cpptraj_dat(dat_file: Union[str, Path], value_col: int = 1) -> pd.DataFrame:
    """
    Parse a cpptraj .dat output file (space-separated, comment lines start with #).
    Returns DataFrame with 'frame' and 'value' columns.
    """
    path = Path(dat_file)
    if not path.exists():
        return pd.DataFrame(columns=["frame", "value"])

    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) > value_col:
                try:
                    frame = int(float(parts[0]))
                    value = float(parts[value_col])
                    rows.append({"frame": frame, "value": value})
                except (ValueError, IndexError):
                    continue

    return pd.DataFrame(rows)


def parse_hbond_avg(avg_file: Union[str, Path]) -> pd.DataFrame:
    """
    Parse cpptraj hbond avgout file.

    Format (7 columns, space-separated):
      #Acceptor  DonorH  Donor  Frames  Frac  AvgDist  AvgAng
      ASP_243@OD2  UNL_707@H6  UNL_707@O2  43  0.8600  2.6672  166.4432

    Frac is 0.0-1.0 (fraction of frames, NOT percentage).
    """
    path = Path(avg_file)
    if not path.exists():
        return pd.DataFrame()

    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 6:
                try:
                    rows.append({
                        "acceptor": parts[0],
                        "donor_h": parts[1],
                        "donor": parts[2],
                        "frames": int(parts[3]),
                        "occupancy": float(parts[4]),  # 0.0-1.0 fraction
                        "avg_distance": float(parts[5]),
                        "avg_angle": float(parts[6]) if len(parts) > 6 else None,
                    })
                except (ValueError, IndexError):
                    continue

    return pd.DataFrame(rows)


def parse_water_bridges(
        bridge_file: Union[str, Path],
        lig_resname: str = "UNL",
) -> pd.DataFrame:
    """
    Parse cpptraj bridgeout file for water-mediated interactions.

    Format:
      #Bridging Solute Residues:
      Bridge Res 110:ASP  139:THR  707:UNL , 151 frames.
      Bridge Res 243:ASP  707:UNL , 25 frames.

    Only lines containing lig_resname are ligand-relevant.
    Extracts the protein residues bridging to the ligand and frame count.
    """
    path = Path(bridge_file)
    if not path.exists():
        return pd.DataFrame()

    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if lig_resname not in line:
                continue
            # Parse: "Bridge Res 110:ASP  139:THR  707:UNL , 151 frames."
            try:
                # Split on comma to separate residues from frame count
                parts = line.split(",")
                if len(parts) < 2:
                    continue
                # Extract frame count: " 151 frames."
                frame_str = parts[-1].strip().replace("frames.", "").replace("frame.", "").strip()
                n_frames = int(frame_str)
                # Extract residues: "Bridge Res 110:ASP  139:THR  707:UNL"
                res_part = parts[0].replace("Bridge Res", "").strip()
                residues = res_part.split()
                # Filter out the ligand itself, keep protein residues
                protein_residues = [r for r in residues if lig_resname not in r]
                rows.append({
                    "protein_residues": " ".join(protein_residues),
                    "all_residues": " ".join(residues),
                    "n_frames": n_frames,
                })
            except (ValueError, IndexError):
                continue

    return pd.DataFrame(rows)


def collect_distances(
        output_dir: Path,
        monitored_residues: List[Dict],
) -> pd.DataFrame:
    """
    Collect all per-frame distance files into a single DataFrame.
    Columns: frame, res1_dist, res2_dist, ...
    """
    all_dists = {}

    for res in monitored_residues:
        res_id = res["pdb_id"].replace(".", "_")
        dat_file = output_dir / f"dist_{res_id}.dat"
        df = parse_cpptraj_dat(dat_file)
        if not df.empty:
            all_dists[res["pdb_id"]] = df.set_index("frame")["value"]

    if not all_dists:
        return pd.DataFrame()

    result = pd.DataFrame(all_dists)
    result.index.name = "frame"
    return result.reset_index()


# =============================================================================
# PROLIF ANALYSIS
# =============================================================================

def run_prolif_analysis(
        complex_prmtop: str,
        dry_trajectory: str,
        output_dir: Path,
        lig_resname: str = "UNL",
) -> Dict:
    """
    Run ProLIF interaction fingerprint analysis on the stripped trajectory.

    Produces:
    - prolif_fingerprints.pkl: raw ProLIF fingerprint object
    - prolif_occupancy.csv: per-residue per-interaction-type occupancy (% of frames)
    - prolif_per_frame.csv: binary interaction matrix per frame
    """
    try:
        import prolif as plf
        import MDAnalysis as mda
    except ImportError:
        logger.warning("  ProLIF or MDAnalysis not installed. Skipping fingerprint analysis.")
        logger.warning("  Install with: pip install prolif mdanalysis")
        return {"success": False, "error": "ProLIF not installed"}

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("  ProLIF: loading trajectory...")

    u = mda.Universe(complex_prmtop, dry_trajectory, format="MDCRD")

    protein = u.select_atoms("protein")
    ligand = u.select_atoms("not protein")

    if len(ligand) == 0:
        logger.error(f"  No non-protein atoms found in topology")
        return {"success": False, "error": "No ligand atoms found (not protein)"}

    logger.info(f"  Protein: {len(protein)} atoms, Ligand: {len(ligand)} atoms")
    logger.info(f"  Frames: {len(u.trajectory)}")

    logger.info("  ProLIF: computing fingerprints...")
    fp = plf.Fingerprint()
    fp.run(u.trajectory, ligand, protein)

    df_fp = fp.to_dataframe()
    n_frames = len(df_fp)

    logger.info(f"  ProLIF: {n_frames} frames processed")

    # Compute occupancy (% of frames each interaction is present)
    occupancy = df_fp.mean()

    occ_rows = []
    for col_tuple, occ_value in occupancy.items():
        if isinstance(col_tuple, tuple) and len(col_tuple) >= 3:
            lig_res, prot_res, int_type = col_tuple[0], col_tuple[1], col_tuple[2]
        else:
            continue

        occ_rows.append({
            "protein_residue": str(prot_res),
            "interaction_type": str(int_type),
            "occupancy": round(float(occ_value) * 100, 1),
            "n_frames": int(occ_value * n_frames),
        })

    df_occ = pd.DataFrame(occ_rows)
    if not df_occ.empty:
        df_occ = df_occ.sort_values("occupancy", ascending=False)

    occ_path = output_dir / "prolif_occupancy.csv"
    df_occ.to_csv(occ_path, index=False)

    frame_path = output_dir / "prolif_per_frame.csv"
    df_fp.to_csv(frame_path)

    logger.info(f"  ProLIF: {len(df_occ)} interaction-residue pairs detected")
    logger.info(f"  Saved: {occ_path}")
    logger.info(f"  Saved: {frame_path}")

    return {
        "success": True,
        "n_frames": n_frames,
        "n_interactions": len(df_occ),
        "occupancy_csv": str(occ_path),
        "per_frame_csv": str(frame_path),
    }


# =============================================================================
# MAIN ORCHESTRATOR
# =============================================================================

def run_trajectory_analysis(
        mol_name: str,
        mmpbsa_dir: Union[str, Path],
        output_dir: Union[str, Path],
        campaign_config: Dict,
        residue_mapping_csv: Optional[str] = None,
        lig_resname: str = "UNL",
        n_protein_res: int = 706,
        cpptraj_timeout: int = 600,
        cpptraj_solvated_timeout: int = 3600,
        run_prolif: bool = True,
        reparse_only: bool = False,
) -> Dict:
    """
    Run full trajectory analysis for one molecule.

    If reparse_only=True, skip cpptraj execution and only re-parse
    existing .dat files with corrected parsers to regenerate .csv outputs.

    Args:
        mol_name: molecule identifier
        mmpbsa_dir: path to 01g_mmpbsa_decomp/{mol_name}/
        output_dir: where to write results
        campaign_config: parsed campaign_config.yaml
        residue_mapping_csv: path to 04b residue_mapping.csv
        lig_resname: ligand residue name in topology
        n_protein_res: number of protein residues
        cpptraj_timeout: timeout for cpptraj in seconds
        run_prolif: whether to run ProLIF analysis
    """
    mmpbsa_dir = Path(mmpbsa_dir).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    # Resolve file paths (absolute for cpptraj compatibility)
    complex_prmtop = (mmpbsa_dir / "topologies" / "complex.prmtop").resolve()
    dry_trajectory = (mmpbsa_dir / "trajectory" / "dry_trajectory.mdcrd").resolve()
    solvated_prmtop = (mmpbsa_dir / "topologies" / "solvated.prmtop").resolve()
    production_dcd = (mmpbsa_dir / "md" / "production.dcd").resolve()

    for f, label in [
        (complex_prmtop, "complex.prmtop"),
        (dry_trajectory, "dry_trajectory.mdcrd"),
    ]:
        if not f.exists():
            return {"success": False, "error": f"Missing {label}: {f}"}

    logger.info("=" * 60)
    logger.info(f"  Trajectory Analysis (01i) -- {mol_name}")
    logger.info("=" * 60)
    logger.info(f"  Complex prmtop: {complex_prmtop}")
    logger.info(f"  Dry trajectory: {dry_trajectory}")
    logger.info(f"  Solvated prmtop: {solvated_prmtop} ({'exists' if solvated_prmtop.exists() else 'MISSING'})")
    logger.info(f"  Production DCD: {production_dcd} ({'exists' if production_dcd.exists() else 'MISSING'})")

    results = {"mol_name": mol_name}

    # --- Auto-detect ligand resname from prmtop ---
    detected_resname = detect_ligand_resname(complex_prmtop)
    if detected_resname:
        logger.info(f"  Ligand auto-detected: {detected_resname}")
        lig_resname = detected_resname
    else:
        logger.warning(f"  Ligand auto-detection failed, using fallback: {lig_resname}")

    # --- Load residue mapping ---
    residue_mapping = {}
    if residue_mapping_csv and Path(residue_mapping_csv).exists():
        residue_mapping = load_residue_mapping(residue_mapping_csv)
        logger.info(f"  Residue mapping: {len(residue_mapping)} residues")
        # Auto-detect n_protein_res from mapping (more reliable than config)
        if len(residue_mapping) > 0:
            n_protein_res = len(residue_mapping)
            logger.info(f"  n_protein_res auto-detected from mapping: {n_protein_res}")
    else:
        logger.warning("  No residue mapping found. Distance monitoring will be limited.")

    # --- Determine monitored residues ---
    monitored = get_monitored_residues(campaign_config, residue_mapping)
    logger.info(f"  Monitored residues: {len(monitored)}")
    for res in monitored:
        logger.info(f"    {res['pdb_id']} (seq {res['seq_num']}) -- {res['label']}")

    # --- Part 1: Stripped trajectory analysis (cpptraj) ---
    # --- Part 1: Stripped trajectory analysis ---
    logger.info("")
    if reparse_only:
        logger.info("--- Part 1: Re-parsing existing .dat files (--reparse-only) ---")
    else:
        logger.info("--- Part 1: Stripped Trajectory (cpptraj) ---")

        script = generate_cpptraj_stripped_script(
            complex_prmtop=str(complex_prmtop),
            dry_trajectory=str(dry_trajectory),
            output_dir=str(output_dir),
            monitored_residues=monitored,
            lig_mask=lig_resname,
            n_protein_res=n_protein_res,
        )

        script_path = output_dir / "cpptraj_stripped.in"
        result_cpp = run_cpptraj(script, script_path, timeout=cpptraj_timeout)

        if not result_cpp["success"]:
            logger.error(f"  cpptraj stripped failed: {result_cpp.get('error', 'unknown')}")

    # Parse .dat files (runs in both normal and reparse_only mode)
    dat_exists = (output_dir / "rmsd_ligand.dat").exists()
    if dat_exists or reparse_only:
        df_rmsd_lig = parse_cpptraj_dat(output_dir / "rmsd_ligand.dat")
        df_rmsd_rec = parse_cpptraj_dat(output_dir / "rmsd_receptor.dat")

        if not df_rmsd_lig.empty:
            rmsd_avg = df_rmsd_lig["value"].mean()
            rmsd_max = df_rmsd_lig["value"].max()
            rmsd_final = df_rmsd_lig["value"].iloc[-1] if len(df_rmsd_lig) > 0 else None
            logger.info(f"  RMSD ligand: avg={rmsd_avg:.2f} A, max={rmsd_max:.2f} A")
            results["rmsd_avg"] = round(rmsd_avg, 2)
            results["rmsd_max"] = round(rmsd_max, 2)
            results["rmsd_final"] = round(rmsd_final, 2) if rmsd_final else None

            df_rmsd_lig.to_csv(output_dir / "rmsd_ligand.csv", index=False)
            df_rmsd_rec.to_csv(output_dir / "rmsd_receptor.csv", index=False)

        # Collect distances
        df_distances = collect_distances(output_dir, monitored)
        if not df_distances.empty:
            df_distances.to_csv(output_dir / "distances.csv", index=False)

            contact_occ = {}
            for col in df_distances.columns:
                if col == "frame":
                    continue
                n_contact = (df_distances[col] < 4.0).sum()
                n_total = len(df_distances)
                contact_occ[col] = round(n_contact / n_total * 100, 1) if n_total > 0 else 0

            results["contact_occupancy"] = contact_occ
            logger.info(f"  Distances: {len(df_distances)} frames x {len(df_distances.columns)-1} residues")

        # Parse H-bonds
        df_hb1 = parse_hbond_avg(output_dir / "hbond_avg.dat")
        df_hb2 = parse_hbond_avg(output_dir / "hbond_avg_rev.dat")
        df_hbonds = pd.concat([df_hb1, df_hb2], ignore_index=True)
        if not df_hbonds.empty:
            df_hbonds.to_csv(output_dir / "hbond_occupancy.csv", index=False)
            n_persistent = len(df_hbonds[df_hbonds["occupancy"] > 0.50])
            logger.info(f"  H-bonds: {len(df_hbonds)} total, {n_persistent} persistent (>50%)")
            results["n_hbonds"] = len(df_hbonds)
            results["n_persistent_hbonds"] = n_persistent

        # Parse RMSF
        df_rmsf = parse_cpptraj_dat(output_dir / "rmsf_per_residue.dat")
        if not df_rmsf.empty:
            df_rmsf.columns = ["residue_seq", "rmsf"]
            df_rmsf.to_csv(output_dir / "rmsf_per_residue.csv", index=False)

    # --- Part 2: Water bridges ---
    wb_dat = output_dir / "water_bridges.dat"
    if reparse_only and wb_dat.exists():
        logger.info("")
        logger.info("--- Part 2: Re-parsing water bridges (--reparse-only) ---")
        df_wb = parse_water_bridges(wb_dat, lig_resname=lig_resname)
        if not df_wb.empty:
            df_wb.to_csv(output_dir / "water_bridges.csv", index=False)
            logger.info(f"  Water bridges involving {lig_resname}: {len(df_wb)} detected")
            results["n_water_bridges"] = len(df_wb)
        else:
            logger.info("  Water bridges: none involving ligand")
    elif not reparse_only and solvated_prmtop.exists() and production_dcd.exists():
        logger.info("")
        logger.info("--- Part 2: Water Bridges (solvated trajectory) ---")

        script_solv = generate_cpptraj_solvated_script(
            solvated_prmtop=str(solvated_prmtop),
            production_dcd=str(production_dcd),
            output_dir=str(output_dir),
            lig_mask=lig_resname,
            n_protein_res=n_protein_res,
        )

        script_path_solv = output_dir / "cpptraj_solvated.in"
        result_solv = run_cpptraj(script_solv, script_path_solv, timeout=cpptraj_solvated_timeout)

        if result_solv["success"]:
            df_wb = parse_water_bridges(output_dir / "water_bridges.dat", lig_resname=lig_resname)
            if not df_wb.empty:
                df_wb.to_csv(output_dir / "water_bridges.csv", index=False)
                logger.info(f"  Water bridges involving {lig_resname}: {len(df_wb)} detected")
                results["n_water_bridges"] = len(df_wb)
            else:
                logger.info("  Water bridges: none involving ligand")
        else:
            logger.warning(f"  Water bridge analysis failed: {result_solv.get('error', '')}")
    elif not reparse_only:
        logger.warning("  Solvated topology/trajectory not found. Skipping water bridges.")

    # --- Part 3: ProLIF fingerprints ---
    if run_prolif and not reparse_only:
        logger.info("")
        logger.info("--- Part 3: ProLIF Interaction Fingerprints ---")

        prolif_dir = output_dir / "prolif"
        result_prolif = run_prolif_analysis(
            complex_prmtop=str(complex_prmtop),
            dry_trajectory=str(dry_trajectory),
            output_dir=prolif_dir,
            lig_resname=lig_resname,
        )

        if result_prolif["success"]:
            results["prolif_n_frames"] = result_prolif["n_frames"]
            results["prolif_n_interactions"] = result_prolif["n_interactions"]
            pass  # ProLIF results stored in CSVs, not in results dict
        else:
            logger.warning(f"  ProLIF failed: {result_prolif.get('error', '')}")

    # --- Summary ---
    logger.info("")
    logger.info("=" * 60)
    logger.info(f"  01i Complete: {mol_name}")
    if "rmsd_avg" in results:
        logger.info(f"    RMSD ligand: {results['rmsd_avg']} A (avg), {results['rmsd_max']} A (max)")
    if "n_persistent_hbonds" in results:
        logger.info(f"    H-bonds: {results.get('n_hbonds', 0)} total, {results['n_persistent_hbonds']} persistent")
    if "n_water_bridges" in results:
        logger.info(f"    Water bridges: {results['n_water_bridges']}")
    logger.info("=" * 60)

    results["success"] = True
    return results


# =============================================================================
# BATCH RUNNER
# =============================================================================

def run_trajectory_analysis_batch(
        campaign_id: str,
        mmpbsa_base_dir: Union[str, Path],
        output_base_dir: Union[str, Path],
        campaign_config: Dict,
        residue_mapping_csv: Optional[str] = None,
        molecules: Optional[List[str]] = None,
        lig_resname: str = "UNL",
        n_protein_res: int = 706,
        run_prolif: bool = True,
        cpptraj_timeout: int = 600,
        cpptraj_solvated_timeout: int = 3600,
        reparse_only: bool = False,
) -> Dict:
    """
    Run trajectory analysis for all (or selected) molecules in a campaign.
    """
    mmpbsa_base = Path(mmpbsa_base_dir)
    output_base = Path(output_base_dir)

    if molecules:
        mol_dirs = [mmpbsa_base / m for m in molecules if (mmpbsa_base / m).exists()]
    elif reparse_only:
        # In reparse mode, discover from existing output dirs
        mol_dirs = sorted([
            mmpbsa_base / d.name for d in output_base.iterdir()
            if d.is_dir() and (d / "rmsd_ligand.dat").exists()
        ])
    else:
        mol_dirs = sorted([
            d for d in mmpbsa_base.iterdir()
            if d.is_dir() and (d / "md" / "production.dcd").exists()
        ])

    logger.info("=" * 60)
    logger.info(f"  01i Trajectory Analysis Batch")
    logger.info(f"  Campaign: {campaign_id}")
    logger.info(f"  Molecules: {len(mol_dirs)}")
    logger.info("=" * 60)

    all_results = {}
    for mol_dir in mol_dirs:
        mol_name = mol_dir.name
        mol_output = output_base / mol_name

        try:
            result = run_trajectory_analysis(
                mol_name=mol_name,
                mmpbsa_dir=mol_dir,
                output_dir=mol_output,
                campaign_config=campaign_config,
                residue_mapping_csv=residue_mapping_csv,
                lig_resname=lig_resname,
                n_protein_res=n_protein_res,
                run_prolif=run_prolif,
                cpptraj_timeout=cpptraj_timeout,
                cpptraj_solvated_timeout=cpptraj_solvated_timeout,
                reparse_only=reparse_only,
            )
            all_results[mol_name] = result
        except Exception as e:
            logger.error(f"  FAILED {mol_name}: {e}")
            all_results[mol_name] = {"success": False, "error": str(e)}

    n_ok = sum(1 for r in all_results.values() if r.get("success"))
    n_fail = len(all_results) - n_ok

    logger.info("")
    logger.info(f"  Batch complete: {n_ok} OK, {n_fail} failed")

    return {
        "success": n_ok > 0,
        "n_molecules": len(all_results),
        "n_ok": n_ok,
        "n_failed": n_fail,
        "results": all_results,
    }