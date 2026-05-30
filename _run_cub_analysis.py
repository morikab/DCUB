"""
Execute the CUB clustering + local optimization runs + analysis plots.
Run with:  .venv/bin/python3 _run_cub_analysis.py
"""
import json, os, sys, time

# ── Non-interactive matplotlib backend (must be first) ──────────────────────
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Fix sys.path for local module imports ───────────────────────────────────
REPO_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(REPO_ROOT, "app"))
sys.path.insert(0, REPO_ROOT)
os.chdir(REPO_ROOT)

print("=" * 60)
print("Step 1: shared imports and setup")
print("=" * 60)

import random
import pathlib
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.stats as stats
import seaborn as sns
import matplotlib.patheffects as pe
from Bio import AlignIO
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import squareform
from sklearn.decomposition import PCA
from itertools import groupby
import matplotlib.patches as mpatches

from modules.main import run_modules
from analysis.input_testing_data.generate_input_testing_data_for_modules import generate_testing_data

# ── Paths (local equivalents of notebook cell 13) ───────────────────────────
data_directory_path = os.path.join(Path(os.getcwd()).resolve(), "analysis", "example_data")
arabidopsis_data_path = os.path.join(data_directory_path, "arabidopsis")
genomes_gb_path = os.path.join(arabidopsis_data_path, "arabidopsis_microbiome")
zorA_file_path = os.path.join(data_directory_path, "zorA_anti_phage_defense.fasta")
aligned_16s_path = os.path.join(arabidopsis_data_path, "arabidopsis_aligned_16s.fasta")
arabidopsis_microbiome_processed_path = os.path.join(arabidopsis_data_path, "arabidopsis_microbiome_processed")
arabidopsis_mapping_path = os.path.join(arabidopsis_data_path, "sample_to_organism.json")
PROCESSED_DIR = os.path.join(REPO_ROOT, "analysis", "results", "arabidopsis", "processed")
os.makedirs(PROCESSED_DIR, exist_ok=True)

# ── Organism list (full file paths, cell 14) ─────────────────────────────────
genomes_gb_directory = Path(genomes_gb_path)
organisms_list = [str(f) for f in genomes_gb_directory.iterdir() if f.is_file()]
print(f"  {len(organisms_list)} genome files found")

# ── Sample → organism name mapping ──────────────────────────────────────────
with open(arabidopsis_mapping_path) as f:
    samples_to_organism_name_mapping = json.load(f)

print("=" * 60)
print("Step 2: K2P distance matrix + hierarchical clustering")
print("=" * 60)

def k2p(seq1, seq2):
    transitions  = {("A","G"),("G","A"),("C","T"),("T","C")}
    transversions = {("A","C"),("C","A"),("A","T"),("T","A"),
                     ("G","C"),("C","G"),("G","T"),("T","G")}
    s1, s2 = str(seq1).upper(), str(seq2).upper()
    ts = tv = valid = 0
    for a, b in zip(s1, s2):
        if a not in "ACGT" or b not in "ACGT":
            continue
        valid += 1
        if a != b:
            if (a,b) in transitions:   ts += 1
            elif (a,b) in transversions: tv += 1
    if valid == 0: return np.nan
    P, Q = ts / valid, tv / valid
    try:
        return -0.5*np.log(1 - 2*P - Q) - 0.25*np.log(1 - 2*Q)
    except ValueError:
        return np.inf

alignment = AlignIO.read(aligned_16s_path, "fasta")
seq_ids   = [r.id for r in alignment]
n         = len(alignment)
matrix    = np.zeros((n, n))
for i in range(n):
    for j in range(i, n):
        d = k2p(str(alignment[i].seq).upper(), str(alignment[j].seq).upper())
        matrix[i,j] = matrix[j,i] = d

new_index = [samples_to_organism_name_mapping[sid] for sid in seq_ids]
dist_df   = pd.DataFrame(matrix, index=new_index, columns=new_index)
print(f"  K2P matrix: {dist_df.shape}")

# ── DBI function ─────────────────────────────────────────────────────────────
def dbi_from_distance_matrix(dist_matrix, labels):
    unique_labels = np.unique(labels)
    s_values = {}
    for lab in unique_labels:
        mask   = labels == lab
        within = dist_matrix[mask][:, mask]
        if within.shape[0] > 1:
            idx = np.triu_indices_from(within, k=1)
            s_values[lab] = np.mean(within[idx])
        else:
            s_values[lab] = 0
    db_ratios = []
    for i, li in enumerate(unique_labels):
        max_ratio = 0
        for j, lj in enumerate(unique_labels):
            if i == j: continue
            mi, mj = labels == li, labels == lj
            d_ij   = np.mean(dist_matrix[mi][:, mj])
            ratio  = (s_values[li] + s_values[lj]) / d_ij if d_ij > 0 else np.inf
            if ratio > max_ratio: max_ratio = ratio
        db_ratios.append(max_ratio)
    return np.mean(db_ratios)

# K2P optimal clustering
condensed_dist = squareform(dist_df.values)
Z_k2p = linkage(condensed_dist, method="average")

k2p_optimal_k, k2p_min_dbi, cluster_lables = 2, None, None
for k in range(2, 10):
    lbl = fcluster(Z_k2p, k, criterion="maxclust")
    dbi = dbi_from_distance_matrix(dist_df.values, lbl)
    if k2p_min_dbi is None or dbi < k2p_min_dbi:
        k2p_min_dbi, k2p_optimal_k, cluster_lables = dbi, k, lbl

cluster_assignment_dict = dict(zip(new_index, cluster_lables))
organism_clusters = [
    [k for k, v in cluster_assignment_dict.items() if v == n]
    for n in sorted(set(cluster_assignment_dict.values()))
]
print(f"  K2P optimal k={k2p_optimal_k} (DBI={k2p_min_dbi:.4f})")
for i, cl in enumerate(organism_clusters, 1):
    print(f"    K2P Cluster {i}: {len(cl)} orgs")

# K2P cluster distances (pkl)
cluster_distances = pd.DataFrame(
    index=new_index, columns=["cluster_id","intra_cluster_dist","inter_cluster_dist"]
)
for org in new_index:
    cid = cluster_assignment_dict[org]
    cluster_distances.loc[org, "cluster_id"] = cid
    same  = [o for o in organism_clusters[cid-1] if o != org]
    other = [o for cl in organism_clusters for o in cl
             if organism_clusters.index(cl)+1 != cid]
    cluster_distances.loc[org, "intra_cluster_dist"] = (
        np.nanmean([dist_df.loc[org, o] for o in same]) if same else np.nan
    )
    cluster_distances.loc[org, "inter_cluster_dist"] = (
        np.nanmean([dist_df.loc[org, o] for o in other]) if other else np.nan
    )
cluster_distances.to_pickle(os.path.join(PROCESSED_DIR, "cluster_distances.pkl"))
print("  K2P cluster distances saved.")

print("=" * 60)
print("Step 3: CUB distance matrix + hierarchical clustering")
print("=" * 60)

# Load CAI profiles
cai_data, cai_organism_names = [], []
for fname in sorted(os.listdir(arabidopsis_microbiome_processed_path)):
    if fname.endswith(".json"):
        with open(os.path.join(arabidopsis_microbiome_processed_path, fname)) as f:
            content = json.load(f)
        cai_data.append(content.get("cai_weights", {}))
        cai_organism_names.append(os.path.splitext(fname)[0])

arabidopsis_cai_weights = pd.DataFrame(cai_data, index=cai_organism_names)
df = arabidopsis_cai_weights
print(f"  CAI profiles: {arabidopsis_cai_weights.shape}")

# CUB distance matrix (1 - Spearman ρ)
from scipy.stats import spearmanr as _spearmanr
_mat = arabidopsis_cai_weights.values
_n   = _mat.shape[0]
_cub = np.zeros((_n, _n))
for i in range(_n):
    for k in range(i+1, _n):
        r, _ = _spearmanr(_mat[i], _mat[k])
        _cub[i,k] = _cub[k,i] = 1.0 - r
cai_distance_df = pd.DataFrame(_cub, index=cai_organism_names, columns=cai_organism_names)
print(f"  CUB distance matrix: {cai_distance_df.shape}")

# CUB DBI clustering
condensed_cub = squareform(cai_distance_df.values, checks=False)
Z_cub = linkage(condensed_cub, method="average")

cub_optimal_k, cub_min_dbi, cub_optimal_labels = 2, None, None
print("  CUB DBI scan:")
for k in range(2, 10):
    lbl = fcluster(Z_cub, k, criterion="maxclust")
    dbi = dbi_from_distance_matrix(cai_distance_df.values, lbl)
    print(f"    k={k}: DBI={dbi:.4f}")
    if cub_min_dbi is None or dbi < cub_min_dbi:
        cub_min_dbi, cub_optimal_k, cub_optimal_labels = dbi, k, lbl

cub_labels_k4 = fcluster(Z_cub, 4, criterion="maxclust")

def _make_clusters(labels, names):
    d = dict(zip(names, labels))
    return d, [[k for k,v in d.items() if v==n] for n in sorted(set(d.values()))]

cub_assignment_optimal, cub_org_clusters_optimal = _make_clusters(cub_optimal_labels, cai_organism_names)
cub_assignment_k4,      cub_org_clusters_k4      = _make_clusters(cub_labels_k4,      cai_organism_names)

print(f"\n  Optimal CUB k={cub_optimal_k} (DBI={cub_min_dbi:.4f})")
for i, cl in enumerate(cub_org_clusters_optimal, 1):
    print(f"    CUB-opt Cluster {i}: {len(cl)} orgs  — {', '.join(sorted(cl))}")
print("  CUB k=4:")
for i, cl in enumerate(cub_org_clusters_k4, 1):
    print(f"    Cluster {i}: {len(cl)} orgs  — {', '.join(sorted(cl))}")

print("=" * 60)
print("Step 4: Generate CUB clustering figures")
print("=" * 60)

bg = "#f7f8f5"

def _cub_combined_figure(cai_dist_df, cub_lbl_arr, org_names, k, title_suffix, save_path):
    sns.set(style="whitegrid", context="paper")
    cluster_ids = sorted(set(cub_lbl_arr))
    base_pal    = sns.color_palette("tab10", len(cluster_ids))
    cl_colors   = {cid: base_pal[i] for i, cid in enumerate(cluster_ids)}
    cl_series   = pd.Series(dict(zip(org_names, cub_lbl_arr)), name="cub_cluster")

    pca_2d = PCA(n_components=2).fit_transform(arabidopsis_cai_weights.loc[org_names].values)
    pca_df = pd.DataFrame({"PC1": pca_2d[:,0], "PC2": pca_2d[:,1], "organism": org_names})
    pca_df = pca_df.join(cl_series, on="organism")

    _cond   = squareform(cai_dist_df.loc[org_names, org_names].values, checks=False)
    _Z      = linkage(_cond, method="average")
    _ddata  = dendrogram(_Z, no_plot=True)
    leaf_ord = _ddata["leaves"]
    ord_names = [org_names[i] for i in leaf_ord]
    ord_dist  = cai_dist_df.loc[ord_names, ord_names]
    ord_cl    = pd.Series(dict(zip(org_names, cub_lbl_arr))).loc[ord_names]
    row_colors = [cl_colors[c] for c in ord_cl]

    strip_blocks = []
    pos = 0
    for cid, grp in groupby(ord_cl.values):
        cnt = sum(1 for _ in grp)
        strip_blocks.append((cid, pos, pos+cnt-1))
        pos += cnt

    def _sx(x): return (x - 5.0) / 10.0 + 0.5

    from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
    n_org = len(org_names)
    fig   = plt.figure(figsize=(12, 15), facecolor=bg)
    gs    = GridSpec(3, 1, figure=fig, height_ratios=[1.0, 0.09, 1.6], hspace=0.22)

    ax_pca = fig.add_subplot(gs[0,0])
    ax_pca.set_facecolor(bg)
    ax_pca.grid(True, color="#d4dadb", linestyle="--", linewidth=1)
    ax_pca.axhline(0, color="grey", linewidth=2.5)
    ax_pca.axvline(0, color="grey", linewidth=2.5)
    ax_pca.scatter(pca_df["PC1"], pca_df["PC2"],
                   c=[cl_colors[r] for r in pca_df["cub_cluster"]],
                   s=70, edgecolor="black", alpha=0.9)
    ax_pca.set_title(f"PCA of CAI Profiles by CUB Cluster (k={k}, {title_suffix})",
                     fontweight="bold", pad=15)
    ax_pca.set_xlabel("Principal Component 1"); ax_pca.set_ylabel("Principal Component 2")
    ax_pca.text(-0.05, 1.02, "A.", transform=ax_pca.transAxes,
                fontsize=14, fontweight="bold", va="bottom", ha="left", clip_on=False)

    ax_leg = fig.add_subplot(gs[1,0])
    ax_leg.axis("off"); ax_leg.set_facecolor(bg)
    leg_handles = [mpatches.Patch(color=cl_colors[c], label=f"Cluster {c}") for c in cluster_ids]
    leg = ax_leg.legend(handles=leg_handles, loc="center", ncol=len(cluster_ids),
                        fontsize=9, frameon=True, framealpha=0.95,
                        title="CUB Clusters", title_fontsize=10)
    leg.get_frame().set_edgecolor("#cccccc")

    ax_b = fig.add_subplot(gs[2,0])
    ax_b.axis("off"); ax_b.set_facecolor("none")
    ax_b.text(-0.05, 1.01, "B.", transform=ax_b.transAxes,
              fontsize=14, fontweight="bold", va="bottom", ha="left", clip_on=False)
    bot = GridSpecFromSubplotSpec(2, 3, subplot_spec=gs[2,0],
                                  height_ratios=[0.32, 1.0], width_ratios=[0.10, 1.0, 0.04],
                                  hspace=0.04, wspace=0.03)
    ax_den = fig.add_subplot(bot[0,1]); ax_strip = fig.add_subplot(bot[1,0])
    ax_hm  = fig.add_subplot(bot[1,1]); ax_cb    = fig.add_subplot(bot[1,2])

    sns.heatmap(ord_dist, cmap="viridis", ax=ax_hm,
                xticklabels=ord_names, yticklabels=False,
                cbar=True, cbar_ax=ax_cb,
                cbar_kws={"label": "CUB distance (1 − Spearman ρ)"})
    ax_hm.set_facecolor(bg); ax_hm.set_xlabel("Organisms")
    ax_hm.tick_params(axis="x", labelrotation=90, labelsize=7)
    for _, start, _ in strip_blocks[1:]:
        ax_hm.axhline(start, color="white", linewidth=2.0)
        ax_hm.axvline(start, color="white", linewidth=2.0)

    strip_arr = np.array(row_colors).reshape(n_org, 1, 3)
    ax_strip.imshow(strip_arr, aspect="auto", origin="upper")
    ax_strip.set_xticks([])
    ax_strip.set_yticks(np.arange(n_org) + 0.5)
    ax_strip.set_yticklabels(ord_names, fontsize=7)
    ax_strip.tick_params(axis="y", length=0); ax_strip.yaxis.tick_left()
    ax_strip.set_facecolor(bg)
    for sp in ax_strip.spines.values(): sp.set_visible(False)
    for cid, start, end in strip_blocks:
        ax_strip.text(0, (start+end)/2, str(cid), ha="center", va="center",
                      fontsize=8, fontweight="bold", color="white",
                      path_effects=[pe.withStroke(linewidth=2.0, foreground="black")], clip_on=True)

    for xs, ys in zip(_ddata["icoord"], _ddata["dcoord"]):
        ax_den.plot([_sx(x) for x in xs], ys, color="grey", linewidth=1.6)
    ax_den.set_xlim(0, n_org)
    ax_den.set_ylim(0, max(max(y) for y in _ddata["dcoord"]) * 1.08)
    ax_den.set_xticks([]); ax_den.set_yticks([])
    ax_den.set_facecolor(bg)
    ax_den.set_title(f"CUB Distance Heatmap — Hierarchical Clustering (k={k}, {title_suffix})",
                     fontweight="bold", pad=10)
    for sp in ax_den.spines.values(): sp.set_visible(False)

    fig.subplots_adjust(left=0.08, right=0.96, top=0.96, bottom=0.05)
    fig.savefig(save_path, bbox_inches="tight", facecolor=bg)
    plt.close(fig)
    print(f"  Saved: {save_path}")

_cub_combined_figure(cai_distance_df, cub_optimal_labels, cai_organism_names,
                     k=cub_optimal_k, title_suffix="DBI-optimal",
                     save_path=os.path.join(REPO_ROOT, "combined_pca_cub_clustering_optimal.svg"))
_cub_combined_figure(cai_distance_df, cub_labels_k4, cai_organism_names,
                     k=4, title_suffix="k=4",
                     save_path=os.path.join(REPO_ROOT, "combined_pca_cub_clustering_k4.svg"))

# Comparison figure (K2P vs CUB-opt vs CUB-k4)
_common = sorted(set(new_index) & set(cai_organism_names))
_k2p_a  = dict(zip(new_index, cluster_lables))
comp_df = pd.DataFrame({
    "organism":   _common,
    "k2p_cluster":  [_k2p_a[o]               for o in _common],
    "cub_optimal":  [cub_assignment_optimal[o] for o in _common],
    "cub_k4":       [cub_assignment_k4[o]       for o in _common],
})
X_c   = arabidopsis_cai_weights.loc[_common].values
pca_c = PCA(n_components=2).fit_transform(X_c)
comp_df["PC1"], comp_df["PC2"] = pca_c[:,0], pca_c[:,1]

sns.set(style="whitegrid", context="paper")
fig_cmp, axes_cmp = plt.subplots(1, 3, figsize=(17, 5.5), facecolor=bg)
fig_cmp.suptitle("Cluster Compositions: K2P vs CUB", fontsize=13, fontweight="bold", y=1.01)
for ax, col, ttl in zip(axes_cmp,
    ["k2p_cluster", "cub_optimal", "cub_k4"],
    [f"K2P (k={k2p_optimal_k})", f"CUB DBI-optimal (k={cub_optimal_k})", "CUB k=4"]):
    cvals = sorted(comp_df[col].unique())
    pal   = sns.color_palette("tab10", len(cvals))
    cmap  = {c: pal[i] for i,c in enumerate(cvals)}
    ax.set_facecolor(bg)
    ax.grid(True, color="#d4dadb", linestyle="--", linewidth=0.8)
    ax.axhline(0, color="grey", linewidth=1.5); ax.axvline(0, color="grey", linewidth=1.5)
    ax.scatter(comp_df["PC1"], comp_df["PC2"],
               c=[cmap[v] for v in comp_df[col]], s=65, edgecolor="black", alpha=0.9)
    ax.set_title(ttl, fontweight="bold", pad=10)
    ax.set_xlabel("PC 1"); ax.set_ylabel("PC 2")
    handles = [mpatches.Patch(color=cmap[c], label=f"C{c}") for c in cvals]
    ax.legend(handles=handles, fontsize=8)
plt.tight_layout()
fig_cmp.savefig(os.path.join(REPO_ROOT, "cub_vs_k2p_crosstab.svg"), bbox_inches="tight", facecolor=bg)
plt.close(fig_cmp)
print("  Saved: cub_vs_k2p_crosstab.svg")

# Cross-tabs
print("\n  Cross-tab K2P vs CUB-optimal:")
print(pd.crosstab(comp_df["k2p_cluster"], comp_df["cub_optimal"], rownames=["K2P"], colnames=["CUB-opt"]).to_string())
print("\n  Cross-tab K2P vs CUB-k4:")
print(pd.crosstab(comp_df["k2p_cluster"], comp_df["cub_k4"], rownames=["K2P"], colnames=["CUB-k4"]).to_string())

# Save CUB cluster distances
def _cub_cluster_distances(assignment_dict, dist_df):
    rows = []
    for org, cid in assignment_dict.items():
        same  = [o for o,c in assignment_dict.items() if c == cid and o != org]
        diff  = [o for o,c in assignment_dict.items() if c != cid]
        intra = np.nanmean([dist_df.loc[org,o] for o in same]) if same else np.nan
        inter = np.nanmean([dist_df.loc[org,o] for o in diff]) if diff else np.nan
        rows.append({"organism": org, "cluster_id": cid,
                     "intra_cluster_dist": intra, "inter_cluster_dist": inter})
    return pd.DataFrame(rows).set_index("organism")

cub_distances_optimal = _cub_cluster_distances(cub_assignment_optimal, cai_distance_df)
cub_distances_k4      = _cub_cluster_distances(cub_assignment_k4,      cai_distance_df)
cub_distances_optimal.to_pickle(os.path.join(PROCESSED_DIR, "cub_cluster_distances_optimal.pkl"))
cub_distances_k4.to_pickle(     os.path.join(PROCESSED_DIR, "cub_cluster_distances_k4.pkl"))
print("  CUB cluster distance pkls saved.")

print("=" * 60)
print("Step 5: CUB cluster-vs-all optimization runs")
print("=" * 60)

_LOCAL_METHODS = [
    ("zscore_bulk_aa_ratio", "_zscore_bulk_ratio"),
    ("zscore_bulk_aa_diff",  "_zscore_bulk_diff"),
    ("single_codon_ratio",   "_single_ratio"),
    ("single_codon_diff",    "_single_diff"),
]

def _extract_result_record(wanted_paths, unwanted_paths, result):
    def _name(p): return os.path.splitext(os.path.basename(p))[0]
    evals = result.get("evaluation", [])
    if len(evals) == 1:
        orf = result.get("orf")
    else:
        final_ev = result.get("final_evaluation", {})
        orf = next(
            (result["orf"][i] for i,ev in enumerate(evals)
             if ev.get("average_distance_score") == final_ev.get("average_distance_score")),
            None
        )
    if orf is None:
        return {"error_message": "orf summary not found"}
    final_ev = result.get("final_evaluation", {})
    return {
        "wanted_organisms":  "|".join(sorted(_name(p) for p in wanted_paths)),
        "unwanted_organisms": "|".join(sorted(_name(p) for p in unwanted_paths)),
        "average_distance_score": final_ev.get("average_distance_score"),
        "average_distance_non_normalized_score": final_ev.get("average_distance_non_normalized_score"),
        "weakest_link_score": final_ev.get("weakest_link_score"),
        "ratio_score": final_ev.get("ratio_score"),
        "orf_optimization_cub_index": result.get("module_input", {}).get("orf_optimization_cub_index"),
        "run_time": orf.get("run_time"),
        "final_sequence": final_ev.get("final_sequence"),
    }


def run_cub_cluster_vs_all_local(org_clusters, label):
    path_cache = {}
    method_dfs = {}
    for method, suffix in _LOCAL_METHODS:
        print(f"  [{label}] {method} ...", flush=True)
        t0 = time.time()
        records = []
        for cluster_orgs in org_clusters:
            key = tuple(sorted(cluster_orgs))
            if key not in path_cache:
                _wanted_names = set(cluster_orgs)
                _wanteds  = [os.path.join(genomes_gb_path, o)+".gbff" for o in cluster_orgs]
                _unwanted = [f for f in organisms_list
                             if os.path.splitext(os.path.basename(f))[0] not in _wanted_names]
                path_cache[key] = (_wanteds, _unwanted)
            wanted_paths, unwanted_paths = path_cache[key]

            try:
                user_input = generate_testing_data(
                    orf_optimization_method=method,
                    orf_optimization_cub_index="CAI",
                    wanted_hosts=wanted_paths,
                    unwanted_hosts=unwanted_paths,
                    genome_path=genomes_gb_path,
                    sequence_file_path=zorA_file_path,
                    output_path=f"arabidopsis/cub_{label}_cluster_vs_all",
                )
                result = run_modules(user_input, should_run_output_module=False)
                record = _extract_result_record(wanted_paths, unwanted_paths, result)
            except Exception as e:
                print(f"    ERROR: {e}")
                record = {"error_message": str(e)}

            if "error_message" not in record:
                records.append(record)
            else:
                print(f"    Skipped cluster {cluster_orgs}: {record['error_message'][:80]}")

        df = pd.DataFrame(records)
        key_cols = {"wanted_organisms", "unwanted_organisms"}
        df = df.rename(columns={c: f"{c}{suffix}" for c in df.columns if c not in key_cols})
        method_dfs[method] = df
        print(f"    {len(records)} records, {time.time()-t0:.1f}s")

    keys = ["wanted_organisms", "unwanted_organisms"]
    bulk   = method_dfs["zscore_bulk_aa_ratio"].merge(method_dfs["zscore_bulk_aa_diff"], on=keys, how="inner")
    single = method_dfs["single_codon_ratio"].merge(method_dfs["single_codon_diff"], on=keys, how="inner")
    merged = bulk.merge(single, on=keys, how="inner")
    print(f"  [{label}] Merged: {len(merged)} rows")
    return merged


t_run = time.time()
print(f"Running CUB DBI-optimal k={cub_optimal_k} ...")
cub_cva_optimal = run_cub_cluster_vs_all_local(cub_org_clusters_optimal, f"optimal_k{cub_optimal_k}")
cub_cva_optimal.to_pickle(os.path.join(PROCESSED_DIR, "cub_cluster_vs_all_optimal.pkl"))
print(f"  Saved cub_cluster_vs_all_optimal.pkl  ({time.time()-t_run:.1f}s total)")

t_run = time.time()
print(f"Running CUB k=4 ...")
cub_cva_k4 = run_cub_cluster_vs_all_local(cub_org_clusters_k4, "k4")
cub_cva_k4.to_pickle(os.path.join(PROCESSED_DIR, "cub_cluster_vs_all_k4.pkl"))
print(f"  Saved cub_cluster_vs_all_k4.pkl  ({time.time()-t_run:.1f}s total)")

print("=" * 60)
print("Step 6: Analysis plots")
print("=" * 60)

# Load the one_vs_all baseline (pre-existing from K2P runs)
one_vs_all_path = os.path.join(PROCESSED_DIR, "one_vs_all.pkl")
one_vs_all = pd.read_pickle(one_vs_all_path)
print(f"  one_vs_all loaded: {len(one_vs_all)} rows")
print(f"  columns: {list(one_vs_all.columns)}")

# ── Helpers ──────────────────────────────────────────────────────────────────
_SCORE = "average_distance_score_zscore_bulk_ratio"

def _enrich_cva(cva_df, cub_dist_df):
    cva = cva_df.copy()
    cva["first_org"]  = cva["wanted_organisms"].str.split("|").str[0]
    org_to_cid = dict(zip(cub_dist_df.index, cub_dist_df["cluster_id"]))
    cva["cluster_id"] = cva["first_org"].map(org_to_cid)
    dist_agg = (cub_dist_df.groupby("cluster_id")[["intra_cluster_dist","inter_cluster_dist"]]
                .mean().reset_index())
    return cva.merge(dist_agg, on="cluster_id", how="left")

cub_cva_opt_enhanced = _enrich_cva(cub_cva_optimal, cub_distances_optimal)
cub_cva_k4_enhanced  = _enrich_cva(cub_cva_k4,      cub_distances_k4)

def _cub_strip_plot(cva_enhanced, org_clusters, one_vs_all_df, cluster_label, save_path):
    _org_to_cid = {org: cid for cid, cl in enumerate(org_clusters, 1) for org in cl}
    _clust_exp  = (cva_enhanced
                   .assign(_wl=lambda x: x["wanted_organisms"].str.split("|"))
                   .explode("_wl").rename(columns={"_wl": "organism"})
                   [["organism", _SCORE]].assign(source="Cluster"))
    _single_exp = (one_vs_all_df[["wanted_organisms", _SCORE]]
                   .rename(columns={"wanted_organisms": "organism"})
                   .assign(source="Single"))
    for _d in [_clust_exp, _single_exp]:
        _d["cluster_id"] = _d["organism"].map(_org_to_cid)

    _df_plot = pd.concat([_clust_exp, _single_exp], ignore_index=True).dropna(subset=["cluster_id", _SCORE])
    _cids    = sorted(_df_plot["cluster_id"].dropna().unique())
    _org_ord = [o for cid in _cids
                for o in sorted(_df_plot.loc[_df_plot["cluster_id"]==cid,"organism"].unique())]

    sns.set(style="whitegrid", context="paper", font_scale=1.1)
    _fig, _ax = plt.subplots(figsize=(14, 6))
    _fig.suptitle(f"Cluster vs Single Organism: {cluster_label}", fontweight="bold", y=1.01)
    sns.stripplot(data=_df_plot, x="organism", y=_SCORE, hue="source",
                  dodge=True, size=8, alpha=0.75, order=_org_ord,
                  palette={"Cluster":"#2E86AB","Single":"#A23B72"},
                  edgecolor="white", linewidth=0.8, ax=_ax)
    for _i, _cid in enumerate(_cids):
        _members = [o for o in _org_ord if _org_to_cid.get(o)==_cid]
        if _members:
            _s = _org_ord.index(_members[0]) - 0.5
            _e = _org_ord.index(_members[-1]) + 0.5
            _ax.axvspan(_s, _e, color="lightgrey", alpha=0.15 if _i%2==0 else 0.08, zorder=0)
            _ax.text((_s+_e)/2, _ax.get_ylim()[1]*0.98, f"C{int(_cid)}",
                     ha="center", va="top", fontsize=9, fontweight="bold", color="grey")
    _mc = _df_plot[_df_plot["source"]=="Cluster"][_SCORE].mean()
    _ms = _df_plot[_df_plot["source"]=="Single"][_SCORE].mean()
    _pct = 100*(_mc-_ms)/abs(_ms) if _ms else 0
    _ax.text(0.02, 0.98, f"Δ = {_mc:.3f} vs {_ms:.3f} ({_pct:+.1f}%)",
             transform=_ax.transAxes, va="top", ha="left", fontsize=11,
             bbox=dict(boxstyle="round,pad=0.4", facecolor="white", alpha=0.9))
    _ax.set_xlabel("Host organism", fontsize=12)
    _ax.set_ylabel("Optimization score", fontsize=12)
    _ax.tick_params(axis="x", rotation=45, labelsize=9)
    _ax.legend(title="Run type", frameon=True)
    sns.despine()
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close(_fig)
    print(f"  Saved: {save_path}")


def _cub_summary_plot(cva_enhanced, org_clusters, one_vs_all_df,
                      cluster_label, ax=None, save_path=None):
    _org_to_cid = {org: cid for cid, cl in enumerate(org_clusters, 1) for org in cl}
    _cluster_sc = (cva_enhanced[["wanted_organisms", _SCORE]]
                   .assign(_fst=lambda x: x["wanted_organisms"].str.split("|").str[0])
                   .assign(cluster_id=lambda x: x["_fst"].map(_org_to_cid))
                   .dropna(subset=["cluster_id"])
                   .groupby("cluster_id")[_SCORE].mean().reset_index()
                   .rename(columns={_SCORE: "cluster_score"}))
    _single_sc  = (one_vs_all_df[["wanted_organisms", _SCORE]]
                   .rename(columns={"wanted_organisms": "organism"})
                   .assign(cluster_id=lambda x: x["organism"].map(_org_to_cid))
                   .dropna(subset=["cluster_id", _SCORE])
                   .groupby("cluster_id")[_SCORE]
                   .agg(single_mean="mean", single_median="median").reset_index())
    _sum = _cluster_sc.merge(_single_sc, on="cluster_id", how="inner")
    _cids = sorted(_sum["cluster_id"].unique())
    _xs   = np.arange(len(_cids))
    _si   = _sum.set_index("cluster_id")

    standalone = ax is None
    if standalone:
        sns.set(style="whitegrid", context="paper", font_scale=1.1)
        _fig_s, ax = plt.subplots(figsize=(max(7, len(_cids)*2), 5))

    ax.plot(_xs, _si.loc[_cids,"cluster_score"], "s-", color="#1B5E20",
            lw=2.5, ms=8, label="Cluster")
    ax.plot(_xs, _si.loc[_cids,"single_median"], "o--", color="#FFB84D",
            lw=2.0, ms=6, label="Individual org (median)")
    for _i in range(len(_cids)):
        ax.axvspan(_i-0.4, _i+0.4, color="lightgrey",
                   alpha=0.15 if _i%2==0 else 0.08)
    _n = {cid: len(org_clusters[int(cid)-1]) for cid in _cids}
    ax.set_xticks(_xs)
    ax.set_xticklabels([f"C{int(c)}\n({_n[c]} orgs)" for c in _cids], fontsize=10, fontweight="bold")
    ax.set_xlabel("CUB cluster", fontsize=11); ax.set_ylabel("Optimization score", fontsize=11)
    ax.set_title(cluster_label, fontweight="bold")
    ax.legend(frameon=True, fontsize=10)
    sns.despine(ax=ax)

    if standalone:
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches="tight")
            plt.close(_fig_s)
            print(f"  Saved: {save_path}")
    print(f"\n  {cluster_label} summary:")
    print("  " + _sum[["cluster_id","cluster_score","single_median","single_mean"]].to_string(index=False))
    return _sum


_cub_strip_plot(cub_cva_opt_enhanced, cub_org_clusters_optimal, one_vs_all,
                f"CUB DBI-optimal (k={cub_optimal_k})",
                save_path=os.path.join(REPO_ROOT, f"graph_cub_k{cub_optimal_k}_cluster_vs_single.pdf"))

_cub_strip_plot(cub_cva_k4_enhanced, cub_org_clusters_k4, one_vs_all,
                "CUB k=4",
                save_path=os.path.join(REPO_ROOT, "graph_cub_k4_cluster_vs_single.pdf"))

sns.set(style="whitegrid", context="paper", font_scale=1.1)
fig_panel, (ax_opt, ax_k4) = plt.subplots(1, 2, figsize=(14, 5))
fig_panel.suptitle("CUB Cluster vs Single Organism — Summary", fontweight="bold")
_sum_opt = _cub_summary_plot(cub_cva_opt_enhanced, cub_org_clusters_optimal, one_vs_all,
                             f"DBI-optimal k={cub_optimal_k}", ax=ax_opt)
_sum_k4  = _cub_summary_plot(cub_cva_k4_enhanced, cub_org_clusters_k4, one_vs_all,
                             "k=4 (forced)", ax=ax_k4)
plt.tight_layout()
fig_panel.savefig(os.path.join(REPO_ROOT, "graph_cub_cluster_summary_panel.pdf"),
                  dpi=300, bbox_inches="tight")
plt.close(fig_panel)
print("  Saved: graph_cub_cluster_summary_panel.pdf")

# K2P cluster summary (for comparison panel — load existing K2P data)
cluster_vs_all_path = os.path.join(PROCESSED_DIR, "cluster_vs_all.pkl")
one_vs_all_enhanced_path = os.path.join(PROCESSED_DIR, "one_vs_all_enhanced.pkl")
cluster_vs_all = pd.read_pickle(cluster_vs_all_path)
one_vs_all_enhanced = pd.read_pickle(one_vs_all_enhanced_path)

def _k2p_summary(cva_enhanced, org_clst, one_vs_all_df):
    _org_to_cid = {o: cid for cid, cl in enumerate(org_clst, 1) for o in cl}
    _cluster_sc = (cva_enhanced[["wanted_organisms", _SCORE]]
                   .assign(_fst=lambda x: x["wanted_organisms"].str.split("|").str[0])
                   .assign(cluster_id=lambda x: x["_fst"].map(_org_to_cid))
                   .dropna(subset=["cluster_id"])
                   .groupby("cluster_id")[_SCORE].mean().reset_index()
                   .rename(columns={_SCORE:"cluster_score"}))
    _single_sc  = (one_vs_all_df[["wanted_organisms", _SCORE]]
                   .rename(columns={"wanted_organisms":"organism"})
                   .assign(cluster_id=lambda x: x["organism"].map(_org_to_cid))
                   .dropna(subset=["cluster_id", _SCORE])
                   .groupby("cluster_id")[_SCORE]
                   .agg(single_mean="mean", single_median="median").reset_index())
    return _cluster_sc.merge(_single_sc, on="cluster_id", how="inner")

# Build K2P cluster_vs_all_enhanced
cluster_vs_all_enhanced = cluster_vs_all.copy()
cluster_vs_all_enhanced["first_org"]  = cluster_vs_all_enhanced["wanted_organisms"].str.split("|").str[0]
cluster_vs_all_enhanced["cluster_id"] = cluster_vs_all_enhanced["first_org"].map(
    dict(zip(cluster_distances.index, cluster_distances["cluster_id"]))
)
dist_agg = (cluster_distances.groupby("cluster_id")[["intra_cluster_dist","inter_cluster_dist"]]
            .mean().reset_index())
cluster_vs_all_enhanced = cluster_vs_all_enhanced.merge(dist_agg, on="cluster_id", how="left")

_k2p_sum = _k2p_summary(cluster_vs_all_enhanced, organism_clusters, one_vs_all_enhanced)
print("\n  K2P summary:")
print("  " + _k2p_sum[["cluster_id","cluster_score","single_median"]].to_string(index=False))

sns.set(style="whitegrid", context="paper", font_scale=1.05)
fig3, axes3 = plt.subplots(1, 3, figsize=(18, 5))
fig3.suptitle("Cluster vs Individual Organism: K2P vs CUB Clustering", fontweight="bold")

_panels = [
    (_k2p_sum, organism_clusters,        axes3[0], f"K2P (k={k2p_optimal_k})"),
    (_sum_opt, cub_org_clusters_optimal, axes3[1], f"CUB DBI-optimal (k={cub_optimal_k})"),
    (_sum_k4,  cub_org_clusters_k4,      axes3[2], "CUB k=4"),
]
for _sum_p, _oc_p, _ax_p, _ttl_p in _panels:
    _cids_p = sorted(_sum_p["cluster_id"].unique())
    _xs_p   = np.arange(len(_cids_p))
    _si_p   = _sum_p.set_index("cluster_id")
    _n_p    = {c: len(_oc_p[int(c)-1]) for c in _cids_p}
    _ax_p.plot(_xs_p, _si_p.loc[_cids_p,"cluster_score"],
               "s-", color="#1B5E20", lw=2.5, ms=8, label="Cluster")
    _ax_p.plot(_xs_p, _si_p.loc[_cids_p,"single_median"],
               "o--", color="#FFB84D", lw=2.0, ms=6, label="Indiv. org (median)")
    for _i_p in range(len(_cids_p)):
        _ax_p.axvspan(_i_p-0.4, _i_p+0.4, color="lightgrey",
                      alpha=0.15 if _i_p%2==0 else 0.08)
    _ax_p.set_xticks(_xs_p)
    _ax_p.set_xticklabels([f"C{int(c)}\n({_n_p[c]})" for c in _cids_p], fontsize=10, fontweight="bold")
    _ax_p.set_xlabel("Cluster"); _ax_p.set_ylabel("Score")
    _ax_p.set_title(_ttl_p, fontweight="bold")
    _ax_p.legend(frameon=True, fontsize=8)
    sns.despine(ax=_ax_p)
plt.tight_layout()
fig3.savefig(os.path.join(REPO_ROOT, "graph_k2p_vs_cub_cluster_comparison.pdf"),
             dpi=300, bbox_inches="tight")
plt.close(fig3)
print("  Saved: graph_k2p_vs_cub_cluster_comparison.pdf")

print("\n" + "=" * 60)
print("ALL DONE")
print("=" * 60)
print("\nGenerated files:")
for f in [
    "combined_pca_cub_clustering_optimal.svg",
    "combined_pca_cub_clustering_k4.svg",
    "cub_vs_k2p_crosstab.svg",
    f"graph_cub_k{cub_optimal_k}_cluster_vs_single.pdf",
    "graph_cub_k4_cluster_vs_single.pdf",
    "graph_cub_cluster_summary_panel.pdf",
    "graph_k2p_vs_cub_cluster_comparison.pdf",
]:
    path = os.path.join(REPO_ROOT, f)
    exists = os.path.exists(path)
    size = os.path.getsize(path) if exists else 0
    print(f"  {'OK' if exists else 'MISSING':6s}  {f}  ({size:,} bytes)")
print()
for pkl in [
    "cluster_distances.pkl",
    "cub_cluster_distances_optimal.pkl",
    "cub_cluster_distances_k4.pkl",
    "cub_cluster_vs_all_optimal.pkl",
    "cub_cluster_vs_all_k4.pkl",
]:
    path = os.path.join(PROCESSED_DIR, pkl)
    exists = os.path.exists(path)
    size = os.path.getsize(path) if exists else 0
    print(f"  {'OK' if exists else 'MISSING':6s}  processed/{pkl}  ({size:,} bytes)")
