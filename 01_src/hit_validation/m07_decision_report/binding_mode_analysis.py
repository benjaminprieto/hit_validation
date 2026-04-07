"""
Binding Mode Analysis — Core Module (07d)
============================================
Compares molecules by how they interact with the receptor.
Groups molecules by binding mode similarity. Lets sub-pocket definitions
emerge from the data instead of predefined zones.

Three layers of analysis:
  Layer 1: Interaction Energy Vectors (IEV) from footprint -> clustering
  Layer 2: ProLIF fingerprints from MD -> clustering (if available)
  Layer 3: Residue co-contact network -> community detection -> emergent sub-pockets

Input:
  Required: 04b footprint_per_molecule.csv
  Optional: 01i ProLIF fingerprints (prolif_occupancy.csv per molecule)
  Optional: 01h per_residue_decomp.csv per molecule (for MMPBSA layer)

Output:
  binding_mode_report.html
  molecule_clusters.csv           -- cluster assignment per molecule
  iev_similarity_matrix.csv       -- pairwise similarity from footprint
  emergent_subpockets.csv         -- data-driven sub-pocket definitions
  residue_cocontact_network.csv   -- co-contact edge list
  dendrogram.png                  -- hierarchical clustering visualization

Location: 01_src/hit_validation/m07_decision_report/binding_mode_analysis.py
Project: hit_validation
Module: 07d (core)
Version: 1.0
"""

import logging
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# =============================================================================
# LAYER 1: INTERACTION ENERGY VECTORS (IEV) FROM FOOTPRINT
# =============================================================================

def build_iev_matrix(
        footprint_csv: Union[str, Path],
        energy_cutoff: float = 0.5,
) -> Tuple[pd.DataFrame, List[str], List[str]]:
    """
    Build Interaction Energy Vector matrix from 04b footprint data.

    Filters to residues with |energy| >= cutoff in at least one molecule.
    Returns:
        df_matrix: molecules (rows) x residues (columns), values = total energy
        mol_names: list of molecule names
        residue_ids: list of residue ids (columns)
    """
    df = pd.read_csv(footprint_csv)

    # Identify residues significant in any molecule
    residue_max = df.groupby("residue_id")["total"].apply(lambda x: x.abs().max())
    significant_residues = residue_max[residue_max >= energy_cutoff].index.tolist()

    logger.info(f"  IEV: {len(significant_residues)} significant residues "
                f"(|energy| >= {energy_cutoff} in at least one molecule)")

    # Pivot to matrix
    df_sig = df[df["residue_id"].isin(significant_residues)]
    df_matrix = df_sig.pivot_table(
        index="Name", columns="residue_id", values="total", aggfunc="first"
    ).fillna(0)

    mol_names = df_matrix.index.tolist()
    residue_ids = df_matrix.columns.tolist()

    logger.info(f"  IEV matrix: {len(mol_names)} molecules x {len(residue_ids)} residues")

    return df_matrix, mol_names, residue_ids


def build_iev_matrix_vdw_es(
        footprint_csv: Union[str, Path],
        energy_cutoff: float = 0.5,
) -> pd.DataFrame:
    """
    Build IEV matrix with separate VDW and ES columns per residue.
    This gives 2x features and captures the character of each interaction.
    """
    df = pd.read_csv(footprint_csv)

    residue_max = df.groupby("residue_id")["total"].apply(lambda x: x.abs().max())
    significant_residues = residue_max[residue_max >= energy_cutoff].index.tolist()

    df_sig = df[df["residue_id"].isin(significant_residues)]

    # Pivot VDW
    df_vdw = df_sig.pivot_table(
        index="Name", columns="residue_id", values="vdw", aggfunc="first"
    ).fillna(0)
    df_vdw.columns = [f"{c}_vdw" for c in df_vdw.columns]

    # Pivot ES
    df_es = df_sig.pivot_table(
        index="Name", columns="residue_id", values="es", aggfunc="first"
    ).fillna(0)
    df_es.columns = [f"{c}_es" for c in df_es.columns]

    # Combine
    df_combined = pd.concat([df_vdw, df_es], axis=1)

    logger.info(f"  IEV (vdw+es): {df_combined.shape[0]} molecules x {df_combined.shape[1]} features")

    return df_combined


def compute_iev_similarity(
        df_matrix: pd.DataFrame,
        method: str = "pearson",
) -> pd.DataFrame:
    """
    Compute pairwise similarity between molecules from IEV matrix.

    Methods:
      pearson: Pearson correlation (recommended for continuous energy vectors)
      cosine: cosine similarity
      tanimoto_binary: Tanimoto on binarized contacts (contact yes/no)
    """
    mol_names = df_matrix.index.tolist()
    n = len(mol_names)

    if method == "pearson":
        sim_matrix = df_matrix.T.corr(method="pearson")

    elif method == "cosine":
        from numpy.linalg import norm
        mat = df_matrix.values
        norms = norm(mat, axis=1, keepdims=True)
        norms[norms == 0] = 1
        normalized = mat / norms
        sim_matrix = pd.DataFrame(
            normalized @ normalized.T,
            index=mol_names, columns=mol_names,
        )

    elif method == "tanimoto_binary":
        binary = (df_matrix.abs() > 0.5).astype(int)
        mat = binary.values
        sim = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                a_and_b = np.sum(mat[i] & mat[j])
                a_or_b = np.sum(mat[i] | mat[j])
                sim[i, j] = a_and_b / a_or_b if a_or_b > 0 else 0
        sim_matrix = pd.DataFrame(sim, index=mol_names, columns=mol_names)

    else:
        raise ValueError(f"Unknown method: {method}")

    return sim_matrix


def cluster_molecules(
        sim_matrix: pd.DataFrame,
        method: str = "ward",
) -> Tuple[Any, pd.DataFrame]:
    """
    Hierarchical clustering on similarity matrix.

    Returns:
        linkage: scipy linkage matrix (for dendrogram)
        df_clusters: DataFrame with molecule, cluster_id
    """
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform

    # Convert similarity to distance
    dist_values = (1 - sim_matrix.clip(lower=-1, upper=1)).values.copy()
    np.fill_diagonal(dist_values, 0)

    # Ensure symmetry
    dist_values = (dist_values + dist_values.T) / 2

    # Condensed distance matrix for scipy
    condensed = squareform(dist_values, checks=False)

    Z = linkage(condensed, method=method)

    # Select optimal K by silhouette score
    best_k = 2
    best_score = -1

    mol_names = sim_matrix.index.tolist()

    if len(mol_names) > 3:
        try:
            from sklearn.metrics import silhouette_score
            for k in range(2, min(7, len(mol_names))):
                labels = fcluster(Z, k, criterion="maxclust")
                if len(set(labels)) < 2:
                    continue
                score = silhouette_score(dist_values, labels, metric="precomputed")
                if score > best_score:
                    best_score = score
                    best_k = k
        except ImportError:
            best_k = 3

    labels = fcluster(Z, best_k, criterion="maxclust")

    df_clusters = pd.DataFrame({
        "Name": mol_names,
        "cluster": labels,
    })

    logger.info(f"  Clustering: {best_k} clusters (silhouette={best_score:.3f})")
    for k in range(1, best_k + 1):
        members = df_clusters[df_clusters["cluster"] == k]["Name"].tolist()
        logger.info(f"    Cluster {k}: {len(members)} molecules -- {', '.join(members[:5])}")

    return Z, df_clusters


# =============================================================================
# LAYER 2: PROLIF FINGERPRINTS (IF AVAILABLE)
# =============================================================================

def load_prolif_fingerprints(
        md_dir: Union[str, Path],
        mol_names: List[str],
) -> Optional[pd.DataFrame]:
    """
    Load ProLIF occupancy data for all molecules, build a combined matrix.
    Columns = (residue, interaction_type), values = occupancy %.
    """
    md_dir = Path(md_dir)

    all_data = {}
    for mol_name in mol_names:
        occ_csv = md_dir / mol_name / "prolif" / "prolif_occupancy.csv"
        if not occ_csv.exists():
            continue

        df = pd.read_csv(occ_csv)
        for _, row in df.iterrows():
            key = f"{row['protein_residue']}_{row['interaction_type']}"
            if mol_name not in all_data:
                all_data[mol_name] = {}
            all_data[mol_name][key] = row["occupancy"]

    if not all_data:
        return None

    df_prolif = pd.DataFrame(all_data).T.fillna(0)
    logger.info(f"  ProLIF matrix: {df_prolif.shape[0]} molecules x {df_prolif.shape[1]} interaction features")

    return df_prolif


# =============================================================================
# LAYER 3: RESIDUE CO-CONTACT NETWORK -> EMERGENT SUB-POCKETS
# =============================================================================

def build_cocontact_network(
        footprint_csv: Union[str, Path],
        energy_cutoff: float = 0.5,
) -> pd.DataFrame:
    """
    Build residue co-contact network.

    For each pair of residues, count how many molecules contact both simultaneously
    (|energy| >= cutoff for both). Edge weight = number of molecules.

    Returns edge list DataFrame: residue_A, residue_B, weight, molecules
    """
    df = pd.read_csv(footprint_csv)

    # Per molecule, get set of contacted residues
    contact_sets = {}
    for mol_name, mol_df in df.groupby("Name"):
        contacted = set(mol_df[mol_df["total"].abs() >= energy_cutoff]["residue_id"])
        contact_sets[mol_name] = contacted

    all_residues = sorted(set().union(*contact_sets.values()))
    n_res = len(all_residues)

    logger.info(f"  Co-contact network: {n_res} residues across {len(contact_sets)} molecules")

    edges = []
    for i in range(n_res):
        for j in range(i + 1, n_res):
            res_a = all_residues[i]
            res_b = all_residues[j]

            n_both = 0
            mol_list = []
            for mol_name, contacts in contact_sets.items():
                if res_a in contacts and res_b in contacts:
                    n_both += 1
                    mol_list.append(mol_name)

            if n_both > 0:
                edges.append({
                    "residue_A": res_a,
                    "residue_B": res_b,
                    "weight": n_both,
                    "fraction": round(n_both / len(contact_sets), 2),
                    "molecules": ",".join(mol_list),
                })

    df_edges = pd.DataFrame(edges)
    logger.info(f"  Co-contact edges: {len(df_edges)}")

    return df_edges


def detect_communities(
        df_edges: pd.DataFrame,
        min_weight: int = 2,
        resolution: float = 1.0,
) -> Tuple[Dict[str, int], Dict[int, List[str]]]:
    """
    Detect communities in the co-contact network using Louvain algorithm.

    Returns:
        residue_to_community: {residue_id: community_id}
        communities: {community_id: [list of residue_ids]}
    """
    try:
        import networkx as nx
        from networkx.algorithms.community import louvain_communities
    except ImportError:
        logger.warning("  networkx not available. Skipping community detection.")
        return {}, {}

    df_filtered = df_edges[df_edges["weight"] >= min_weight]

    if df_filtered.empty:
        logger.warning(f"  No edges with weight >= {min_weight}. Try lower threshold.")
        return {}, {}

    G = nx.Graph()
    for _, row in df_filtered.iterrows():
        G.add_edge(row["residue_A"], row["residue_B"], weight=row["weight"])

    logger.info(f"  Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    communities_sets = louvain_communities(G, weight="weight", resolution=resolution, seed=42)

    residue_to_community = {}
    communities = {}

    for i, community in enumerate(communities_sets):
        sorted_members = sorted(community)
        communities[i] = sorted_members
        for residue in sorted_members:
            residue_to_community[residue] = i

    logger.info(f"  Detected {len(communities)} communities:")
    for cid, members in communities.items():
        logger.info(f"    Community {cid}: {len(members)} residues -- {', '.join(members[:6])}")

    return residue_to_community, communities


def characterize_communities(
        communities: Dict[int, List[str]],
        footprint_csv: Union[str, Path],
        functional_residues: Dict,
) -> pd.DataFrame:
    """
    Characterize each emergent sub-pocket (community).

    For each community: which molecules contact it, total energy,
    dominant interaction type (vdw vs es), functional annotations.
    """
    df = pd.read_csv(footprint_csv)

    rows = []
    for cid, residues in communities.items():
        community_data = df[df["residue_id"].isin(residues)]

        mol_energies = community_data.groupby("Name").agg({
            "vdw": "sum",
            "es": "sum",
            "total": "sum",
        }).reset_index()

        significant = mol_energies[mol_energies["total"].abs() > 1.0]

        total_vdw = community_data["vdw"].sum()
        total_es = community_data["es"].sum()
        character = "VDW-dominated" if abs(total_vdw) > abs(total_es) else "ES-dominated"

        func_annotations = []
        for res in residues:
            res_name = res.split(".")[0] if "." in res else res
            for func_key, func_info in functional_residues.items():
                if res_name.startswith(func_key) or func_key in res_name:
                    func_annotations.append(f"{func_key}: {func_info.get('function', '')}")

        rows.append({
            "community_id": cid,
            "n_residues": len(residues),
            "residues": ", ".join(residues),
            "n_molecules_contacting": len(significant),
            "character": character,
            "mean_energy": round(significant["total"].mean(), 2) if not significant.empty else 0,
            "functional_annotations": "; ".join(func_annotations) if func_annotations else "none",
        })

    return pd.DataFrame(rows)


# =============================================================================
# HTML REPORT
# =============================================================================

def generate_html_report(
        df_clusters: pd.DataFrame,
        sim_matrix: pd.DataFrame,
        df_communities: pd.DataFrame,
        df_edges: pd.DataFrame,
        campaign_id: str,
        dendrogram_b64: Optional[str] = None,
) -> str:
    """Generate the binding mode analysis HTML report."""

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")

    css = """
body { font-family: -apple-system, 'Segoe UI', sans-serif; max-width: 1200px; margin: 0 auto;
       padding: 20px; background: #fafafa; color: #333; }
h1 { color: #2c3e50; border-bottom: 3px solid #9b59b6; padding-bottom: 10px; }
h2 { color: #2c3e50; margin-top: 30px; }
.meta { color: #666; font-size: 0.9em; }
table { border-collapse: collapse; width: 100%; margin: 10px 0; font-size: 13px; background: white; }
th { background: #2c3e50; color: white; padding: 6px 8px; text-align: left; font-weight: 500; }
td { padding: 5px 8px; border-bottom: 1px solid #eee; }
.card { background: white; border-radius: 8px; padding: 20px; margin: 15px 0;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1); }
.cluster-tag { display: inline-block; padding: 2px 10px; border-radius: 10px;
               font-size: 12px; font-weight: 600; margin: 2px; }
.mono { font-family: monospace; }
img { max-width: 100%; margin: 15px 0; border-radius: 4px; }
"""

    # Cluster table
    cluster_rows = []
    for _, row in df_clusters.iterrows():
        cluster_rows.append(
            f'<tr><td>{row["Name"]}</td><td>{row["cluster"]}</td></tr>'
        )

    # Community table
    community_rows = []
    for _, row in df_communities.iterrows():
        community_rows.append(
            f'<tr><td>{row["community_id"]}</td>'
            f'<td>{row["n_residues"]}</td>'
            f'<td>{row["residues"]}</td>'
            f'<td>{row["n_molecules_contacting"]}</td>'
            f'<td>{row["character"]}</td>'
            f'<td class="mono">{row["mean_energy"]}</td>'
            f'<td>{row["functional_annotations"]}</td></tr>'
        )

    dendrogram_html = ""
    if dendrogram_b64:
        dendrogram_html = f'<img src="data:image/png;base64,{dendrogram_b64}" alt="Dendrogram">'

    html = f"""<!DOCTYPE html>
<html><head><meta charset='utf-8'>
<title>07d Binding Mode Analysis -- {campaign_id}</title>
<style>{css}</style>
</head><body>
<h1>Binding Mode Analysis</h1>
<div class="meta">
Campaign: <strong>{campaign_id}</strong> | Generated: {timestamp}
</div>
<p class="meta">Molecules clustered by interaction fingerprint similarity.
Sub-pockets emerge from residue co-contact patterns across all molecules.</p>

<h2>1. Molecule Clustering (IEV Similarity)</h2>
<div class="card">
{dendrogram_html}
<table>
<tr><th>Molecule</th><th>Cluster</th></tr>
{''.join(cluster_rows)}
</table>
</div>

<h2>2. Emergent Sub-pockets (Co-contact Communities)</h2>
<div class="card">
<p>Residues grouped by co-contact patterns: residues frequently contacted by the same molecules
form a community. These are data-driven sub-pockets, not predefined zones.</p>
<table>
<tr><th>ID</th><th>Residues</th><th>Members</th><th>Molecules</th>
<th>Character</th><th>Mean Energy</th><th>Functional</th></tr>
{''.join(community_rows)}
</table>
</div>

<div style="margin-top:40px; padding-top:10px; border-top:1px solid #ddd; font-size:11px; color:#999;">
hit_validation | Module 07d | Binding Mode Analysis v1.0 | {timestamp}
</div>
</body></html>"""

    return html


# =============================================================================
# DENDROGRAM GENERATION
# =============================================================================

def generate_dendrogram_b64(
        linkage_matrix,
        labels: List[str],
) -> Optional[str]:
    """Generate dendrogram as base64 PNG string."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from scipy.cluster.hierarchy import dendrogram
        import base64
        from io import BytesIO

        fig, ax = plt.subplots(figsize=(12, 6))

        short_labels = []
        for label in labels:
            parts = label.split("-")
            if len(parts) >= 3:
                short_labels.append(parts[-2])
            elif label.startswith("PubChem-"):
                short_labels.append(label.replace("PubChem-", ""))
            else:
                short_labels.append(label[:12])

        dendrogram(linkage_matrix, labels=short_labels, ax=ax, leaf_rotation=45, leaf_font_size=10)
        ax.set_ylabel("Distance (1 - Pearson correlation)")
        ax.set_title("Binding Mode Clustering")
        plt.tight_layout()

        buf = BytesIO()
        fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
        plt.close(fig)
        buf.seek(0)

        return base64.b64encode(buf.read()).decode("utf-8")

    except ImportError:
        logger.warning("  matplotlib not available. Skipping dendrogram.")
        return None


# =============================================================================
# MAIN ORCHESTRATOR
# =============================================================================

def run_binding_mode_analysis(
        campaign_id: str,
        footprint_csv: Union[str, Path],
        output_dir: Union[str, Path],
        campaign_config: Dict,
        md_dir: Optional[str] = None,
        energy_cutoff: float = 0.5,
        similarity_method: str = "pearson",
        cocontact_min_weight: int = 2,
        community_resolution: float = 1.0,
        run_communities: bool = True,
        run_prolif: bool = True,
) -> Dict:
    """
    Run binding mode analysis for a campaign.

    Layer 1: IEV clustering from footprint
    Layer 2: ProLIF clustering (if md_dir provided)
    Layer 3: Emergent sub-pockets from co-contact network
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info(f"  07d Binding Mode Analysis")
    logger.info(f"  Campaign: {campaign_id}")
    logger.info("=" * 60)

    functional_residues = campaign_config.get("functional_residues", {})

    # === LAYER 1: IEV from footprint ===
    logger.info("")
    logger.info("--- Layer 1: Interaction Energy Vectors ---")

    df_matrix, mol_names, residue_ids = build_iev_matrix(footprint_csv, energy_cutoff)

    sim_matrix = compute_iev_similarity(df_matrix, method=similarity_method)
    sim_path = output_dir / "iev_similarity_matrix.csv"
    sim_matrix.to_csv(sim_path)
    logger.info(f"  Saved: {sim_path}")

    Z, df_clusters = cluster_molecules(sim_matrix)
    clusters_path = output_dir / "molecule_clusters.csv"
    df_clusters.to_csv(clusters_path, index=False)
    logger.info(f"  Saved: {clusters_path}")

    dendrogram_b64 = generate_dendrogram_b64(Z, mol_names)

    # === LAYER 2: ProLIF (if available) ===
    if run_prolif and md_dir:
        logger.info("")
        logger.info("--- Layer 2: ProLIF Fingerprints ---")

        df_prolif = load_prolif_fingerprints(md_dir, mol_names)
        if df_prolif is not None:
            prolif_binary = (df_prolif > 30).astype(int)
            n = len(prolif_binary)
            mat = prolif_binary.values
            prolif_sim = np.zeros((n, n))
            for i in range(n):
                for j in range(n):
                    a_and_b = np.sum(mat[i] & mat[j])
                    a_or_b = np.sum(mat[i] | mat[j])
                    prolif_sim[i, j] = a_and_b / a_or_b if a_or_b > 0 else 0

            df_prolif_sim = pd.DataFrame(
                prolif_sim,
                index=prolif_binary.index,
                columns=prolif_binary.index,
            )
            prolif_sim_path = output_dir / "prolif_similarity_matrix.csv"
            df_prolif_sim.to_csv(prolif_sim_path)
            logger.info(f"  Saved: {prolif_sim_path}")
        else:
            logger.info("  No ProLIF data found. Skipping Layer 2.")
    else:
        logger.info("  Layer 2 skipped (no ProLIF data or --no-prolif).")

    # === LAYER 3: Emergent sub-pockets ===
    df_edges = pd.DataFrame()
    df_communities = pd.DataFrame()

    if run_communities:
        logger.info("")
        logger.info("--- Layer 3: Emergent Sub-pockets ---")

        df_edges = build_cocontact_network(footprint_csv, energy_cutoff)
        edges_path = output_dir / "residue_cocontact_network.csv"
        df_edges.to_csv(edges_path, index=False)
        logger.info(f"  Saved: {edges_path}")

        residue_to_community, communities = detect_communities(
            df_edges, min_weight=cocontact_min_weight, resolution=community_resolution,
        )

        df_communities = characterize_communities(communities, footprint_csv, functional_residues)
        communities_path = output_dir / "emergent_subpockets.csv"
        df_communities.to_csv(communities_path, index=False)
        logger.info(f"  Saved: {communities_path}")
    else:
        logger.info("  Layer 3 skipped (--no-communities).")

    # === HTML Report ===
    logger.info("")
    logger.info("--- Generating Report ---")

    html = generate_html_report(
        df_clusters=df_clusters,
        sim_matrix=sim_matrix,
        df_communities=df_communities,
        df_edges=df_edges,
        campaign_id=campaign_id,
        dendrogram_b64=dendrogram_b64,
    )

    html_path = output_dir / "binding_mode_report.html"
    html_path.write_text(html, encoding="utf-8")
    logger.info(f"  Saved: {html_path}")

    n_communities = len(df_communities) if not df_communities.empty else 0

    logger.info("")
    logger.info("=" * 60)
    logger.info(f"  07d Complete")
    logger.info(f"    Molecules: {len(mol_names)}")
    logger.info(f"    Clusters: {df_clusters['cluster'].nunique()}")
    logger.info(f"    Emergent sub-pockets: {n_communities}")
    logger.info("=" * 60)

    return {
        "success": True,
        "n_molecules": len(mol_names),
        "n_clusters": int(df_clusters["cluster"].nunique()),
        "n_communities": n_communities,
        "html_path": str(html_path),
    }