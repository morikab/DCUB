"""
Generate the "Clustered Optimization Performance" comparison figure.
  Panel A: K2P k=4 vs CUB k=4 — per-cluster aggregate score vs individual single-org scores
  Panel B: Mean silhouette score per cluster (K2P vs CUB, aggregated with ± std bars)

Run with:  .venv/bin/python3 _plot_comparison_figure.py
"""
import json, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd
from Bio import AlignIO
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr
from sklearn.metrics import silhouette_samples
import seaborn as sns

REPO_ROOT = os.path.dirname(os.path.abspath(__file__))
arab_dir  = os.path.join(REPO_ROOT, "analysis", "example_data", "arabidopsis")
proc_dir  = os.path.join(REPO_ROOT, "analysis", "results", "arabidopsis", "processed")

SCORE     = "average_distance_score_zscore_bulk_ratio"
BG        = "#f7f8f5"

# ── 1. Rebuild K2P distance matrix + k=4 clustering ──────────────────────────
print("Building K2P distance matrix …")
with open(os.path.join(arab_dir, "sample_to_organism.json")) as f:
    samples_to_org = json.load(f)

def _k2p_dist(s1, s2):
    ts_set = {("A","G"),("G","A"),("C","T"),("T","C")}
    tv_set = {("A","C"),("C","A"),("A","T"),("T","A"),
              ("G","C"),("C","G"),("G","T"),("T","G")}
    a, b = str(s1).upper(), str(s2).upper()
    ts = tv = valid = 0
    for x, y in zip(a, b):
        if x not in "ACGT" or y not in "ACGT": continue
        valid += 1
        if x != y:
            if (x, y) in ts_set: ts += 1
            elif (x, y) in tv_set: tv += 1
    if valid == 0: return np.nan
    P, Q = ts/valid, tv/valid
    try:    return -0.5*np.log(1 - 2*P - Q) - 0.25*np.log(1 - 2*Q)
    except: return np.inf

alignment = AlignIO.read(os.path.join(arab_dir, "arabidopsis_aligned_16s.fasta"), "fasta")
seq_ids   = [r.id for r in alignment]
k2p_names = [samples_to_org[sid] for sid in seq_ids]
n = len(alignment)
k2p_mat = np.zeros((n, n))
for i in range(n):
    for j in range(i, n):
        d = _k2p_dist(alignment[i].seq, alignment[j].seq)
        k2p_mat[i, j] = k2p_mat[j, i] = d
k2p_dist_df = pd.DataFrame(k2p_mat, index=k2p_names, columns=k2p_names)

Z_k2p      = linkage(squareform(k2p_mat), method="average")
k2p_lbl_k4 = fcluster(Z_k2p, 4, criterion="maxclust")
k2p_assign  = dict(zip(k2p_names, k2p_lbl_k4.astype(int)))
k2p_sizes   = pd.Series(k2p_lbl_k4).value_counts().sort_index().to_dict()
print(f"  K2P k=4 cluster sizes: {k2p_sizes}")

# ── 2. Rebuild CUB distance matrix + k=4 clustering ──────────────────────────
print("Building CUB distance matrix …")
cai_proc = os.path.join(arab_dir, "arabidopsis_microbiome_processed")
cai_data, cai_names = [], []
for fname in sorted(os.listdir(cai_proc)):
    if fname.endswith(".json"):
        with open(os.path.join(cai_proc, fname)) as f:
            c = json.load(f)
        cai_data.append(c.get("cai_weights", {}))
        cai_names.append(os.path.splitext(fname)[0])

cai_mat_df = pd.DataFrame(cai_data, index=cai_names)
_m = cai_mat_df.values
_nn = _m.shape[0]
cub_mat = np.zeros((_nn, _nn))
for i in range(_nn):
    for j in range(i+1, _nn):
        r, _ = spearmanr(_m[i], _m[j])
        cub_mat[i, j] = cub_mat[j, i] = 1.0 - r
cub_dist_df = pd.DataFrame(cub_mat, index=cai_names, columns=cai_names)

Z_cub      = linkage(squareform(cub_mat, checks=False), method="average")
cub_lbl_k4 = fcluster(Z_cub, 4, criterion="maxclust")
cub_assign  = dict(zip(cai_names, cub_lbl_k4.astype(int)))
cub_sizes   = pd.Series(cub_lbl_k4).value_counts().sort_index().to_dict()
print(f"  CUB k=4 cluster sizes: {cub_sizes}")

# ── 3. Load PKLs and compute per-cluster aggregate scores ────────────────────
ova     = pd.read_pickle(os.path.join(proc_dir, "one_vs_all_enhanced.pkl"))
cva_k2p = pd.read_pickle(os.path.join(proc_dir, "cluster_vs_all.pkl"))      # raw HPC 4-cluster runs
cva_cub = pd.read_pickle(os.path.join(proc_dir, "cub_cluster_vs_all_k4.pkl"))

# Map cluster IDs via first organism in the wanted list
def _add_cluster_id(df, assignment):
    df = df.copy()
    df["_first_org"] = df["wanted_organisms"].str.split("|").str[0]
    df["cluster_id"] = df["_first_org"].map(assignment)
    return df

cva_k2p = _add_cluster_id(cva_k2p, k2p_assign)
cva_cub = _add_cluster_id(cva_cub, cub_assign)

# One aggregate score per cluster
k2p_cl_scores = (cva_k2p.dropna(subset=["cluster_id"])
                 .groupby("cluster_id")[SCORE].mean().to_dict())
cub_cl_scores = (cva_cub.dropna(subset=["cluster_id"])
                 .groupby("cluster_id")[SCORE].mean().to_dict())

print(f"\nK2P cluster scores: { {k: round(v,3) for k,v in sorted(k2p_cl_scores.items())} }")
print(f"CUB cluster scores: { {k: round(v,3) for k,v in sorted(cub_cl_scores.items())} }")

# Add per-organism cluster assignments for OVA individual scores
ova["k2p_cluster"] = ova["wanted_organisms"].map(k2p_assign)
ova["cub_cluster"] = ova["wanted_organisms"].map(cub_assign)

# ── 4. Silhouette scores (aggregated per cluster, ± std) ─────────────────────
k2p_sil = silhouette_samples(k2p_dist_df.loc[k2p_names, k2p_names].values,
                              k2p_lbl_k4, metric="precomputed")
cub_sil = silhouette_samples(cub_dist_df.loc[cai_names, cai_names].values,
                              cub_lbl_k4, metric="precomputed")

k2p_sil_agg = (pd.DataFrame({"sil": k2p_sil, "cluster": k2p_lbl_k4})
               .groupby("cluster")["sil"].agg(mean="mean", std="std")
               .reindex(range(1, 5)))
cub_sil_agg = (pd.DataFrame({"sil": cub_sil, "cluster": cub_lbl_k4})
               .groupby("cluster")["sil"].agg(mean="mean", std="std")
               .reindex(range(1, 5)))

print(f"\nK2P silhouette:\n{k2p_sil_agg.round(3)}")
print(f"\nCUB silhouette:\n{cub_sil_agg.round(3)}")

# ── 5. Build figure ───────────────────────────────────────────────────────────
sns.set(style="whitegrid", context="paper", font_scale=1.1)
np.random.seed(42)

fig = plt.figure(figsize=(14, 11), facecolor=BG)
fig.suptitle("Clustered Optimization Performance", fontsize=16, fontweight="bold", y=0.997)

gs   = GridSpec(2, 2, figure=fig, height_ratios=[1.15, 1.0], hspace=0.48, wspace=0.36)
ax_a1 = fig.add_subplot(gs[0, 0])
ax_a2 = fig.add_subplot(gs[0, 1])
ax_b  = fig.add_subplot(gs[1, :])

# Cluster size lookup for axis labels
k2p_cl_sizes = {int(c): k2p_sizes.get(int(c), "?") for c in range(1, 5)}
cub_cl_sizes = {int(c): cub_sizes.get(int(c), "?") for c in range(1, 5)}

CL_PAL = ["#4E79A7", "#F28E2B", "#E15759", "#76B7B2"]

# ── Panel A helper ────────────────────────────────────────────────────────────
def _panel_a(ax, cl_scores, cl_sizes, ova_df, cluster_col, score_col, title, panel_label):
    cids      = sorted(int(c) for c in cl_scores)
    cl_colors = {cid: CL_PAL[i] for i, cid in enumerate(cids)}
    x_pos     = {cid: i for i, cid in enumerate(cids)}

    ax.set_facecolor(BG)
    ax.grid(True, axis="y", color="#d4dadb", linestyle="--", linewidth=0.8, zorder=0)

    for cid in cids:
        x   = x_pos[cid]
        col = cl_colors[cid]

        # Individual organism scores — jittered scatter
        org_sc = ova_df.loc[ova_df[cluster_col] == cid, score_col].dropna().values
        jitter = np.random.uniform(-0.22, 0.22, size=len(org_sc))
        ax.scatter(x + jitter, org_sc, color=col, s=45, alpha=0.55,
                   edgecolor="white", linewidth=0.4, zorder=3)

        # Cluster aggregate score — horizontal line + diamond
        cs = cl_scores[cid]
        ax.plot([x - 0.38, x + 0.38], [cs, cs], color=col,
                lw=3.0, zorder=4, alpha=0.92, solid_capstyle="round")
        ax.scatter(x, cs, color=col, s=180, marker="D",
                   edgecolor="black", linewidth=1.2, zorder=5)

    ax.set_xticks(list(x_pos.values()))
    ax.set_xticklabels(
        [f"Cluster {c}\n(n={cl_sizes.get(c,'?')})" for c in cids],
        fontsize=9.5
    )
    ax.set_xlabel("Cluster", fontsize=10)
    ax.set_ylabel("Avg. distance score (z-score)", fontsize=9)
    ax.set_title(title, fontweight="bold", pad=10, fontsize=11)
    ax.text(-0.10, 1.07, panel_label, transform=ax.transAxes,
            fontsize=14, fontweight="bold", va="bottom", ha="left", clip_on=False)

    legend_elems = [
        Line2D([0], [0], marker="D", color="w", markerfacecolor="#555",
               markeredgecolor="black", markersize=10, label="Cluster score"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#888",
               markersize=8, alpha=0.6, label="Single-org score"),
    ]
    ax.legend(handles=legend_elems, fontsize=8, loc="upper right", framealpha=0.9)
    sns.despine(ax=ax)


_panel_a(ax_a1, k2p_cl_scores, k2p_cl_sizes, ova, "k2p_cluster", SCORE,
         "K2P Clustering (k = 4)", "A.")
_panel_a(ax_a2, cub_cl_scores, cub_cl_sizes, ova, "cub_cluster", SCORE,
         "CUB Clustering (k = 4)", "")

# ── Panel B: grouped silhouette bars ─────────────────────────────────────────
ax_b.set_facecolor(BG)
ax_b.grid(True, axis="y", color="#d4dadb", linestyle="--", linewidth=0.8, zorder=0)

xs    = np.arange(4)
BAR_W = 0.34
K2P_C = "#4E79A7"
CUB_C = "#E15759"

k2p_means = k2p_sil_agg["mean"].values.astype(float)
k2p_stds  = k2p_sil_agg["std"].fillna(0).values.astype(float)
cub_means = cub_sil_agg["mean"].values.astype(float)
cub_stds  = cub_sil_agg["std"].fillna(0).values.astype(float)

bars_k2p = ax_b.bar(xs - BAR_W/2, k2p_means, BAR_W, yerr=k2p_stds, capsize=6,
                    color=K2P_C, alpha=0.82, edgecolor="black", linewidth=0.8,
                    label="K2P clusters", error_kw={"elinewidth": 1.6, "ecolor": "#333"},
                    zorder=3)
bars_cub = ax_b.bar(xs + BAR_W/2, cub_means, BAR_W, yerr=cub_stds, capsize=6,
                    color=CUB_C, alpha=0.82, edgecolor="black", linewidth=0.8,
                    label="CUB clusters", error_kw={"elinewidth": 1.6, "ecolor": "#333"},
                    zorder=3)

# Value labels above each bar
for bars, means in [(bars_k2p, k2p_means), (bars_cub, cub_means)]:
    for bar, val in zip(bars, means):
        if not np.isnan(val):
            ax_b.text(bar.get_x() + bar.get_width()/2,
                      bar.get_height() + max(k2p_stds.max(), cub_stds.max()) + 0.008,
                      f"{val:.2f}", ha="center", va="bottom",
                      fontsize=8.5, fontweight="bold", color="#222")

ax_b.axhline(0, color="#555", linewidth=1.0, linestyle="--", zorder=2)

# x-axis: show cluster sizes for both methods
xtick_labels = [
    f"Cluster {i+1}\nK2P n={k2p_cl_sizes.get(i+1,'?')} / CUB n={cub_cl_sizes.get(i+1,'?')}"
    for i in range(4)
]
ax_b.set_xticks(xs)
ax_b.set_xticklabels(xtick_labels, fontsize=9)
ax_b.set_xlabel("Cluster", fontsize=11)
ax_b.set_ylabel("Mean silhouette score  (± std)", fontsize=10.5)
ax_b.set_title("B.  Within-Cluster Cohesion: Mean Silhouette Score per Cluster",
               fontweight="bold", pad=10, fontsize=11, loc="left")
ax_b.legend(fontsize=10, framealpha=0.9)
sns.despine(ax=ax_b)

out = os.path.join(REPO_ROOT, "graph_comparison_k2p_vs_cub_performance.pdf")
fig.savefig(out, dpi=300, bbox_inches="tight", facecolor=BG)
plt.close(fig)
print(f"\nSaved: {out} ({os.path.getsize(out):,} bytes)")
