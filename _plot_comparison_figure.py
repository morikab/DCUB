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
import matplotlib.cm as cm
import matplotlib.colors as mcolors
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

# ── 4. Per-organism silhouette scores (for scatter in Panel B) ───────────────
k2p_sil = silhouette_samples(k2p_dist_df.loc[k2p_names, k2p_names].values,
                              k2p_lbl_k4, metric="precomputed")
cub_sil = silhouette_samples(cub_dist_df.loc[cai_names, cai_names].values,
                              cub_lbl_k4, metric="precomputed")

# Per-organism table: name | k2p_sil | k2p_cluster | cub_sil | cub_cluster
k2p_sil_df = pd.DataFrame({"org": k2p_names, "k2p_sil": k2p_sil, "k2p_cluster": k2p_lbl_k4})
cub_sil_df = pd.DataFrame({"org": cai_names,  "cub_sil": cub_sil, "cub_cluster": cub_lbl_k4})
sil_df = k2p_sil_df.merge(cub_sil_df, on="org")
sil_df["same_cluster"] = sil_df["k2p_cluster"] == sil_df["cub_cluster"]

print(f"\nSilhouette scatter data (n={len(sil_df)}):")
print(sil_df[["org","k2p_cluster","k2p_sil","cub_cluster","cub_sil"]].round(3).to_string())

# ── 5. Within-cluster dispersions for K2P and CUB clusters ───────────────────
def _mean_tri(dist_df, members):
    if len(members) < 2:
        return 0.0
    sub = dist_df.loc[members, members].values
    return sub[np.triu_indices(len(members), 1)].mean()

k2p_disp_rows = []
for cid in range(1, 5):
    members = [nm for nm, cl in k2p_assign.items() if cl == cid]
    k2p_disp_rows.append(dict(cluster=cid, n=len(members),
        k2p_disp=_mean_tri(k2p_dist_df, members),
        cub_disp=_mean_tri(cub_dist_df, members),
        score=k2p_cl_scores.get(cid, np.nan)))
k2p_disp_df = pd.DataFrame(k2p_disp_rows)

cub_disp_rows = []
for cid in range(1, 5):
    members = [nm for nm, cl in cub_assign.items() if cl == cid]
    cub_disp_rows.append(dict(cluster=cid, n=len(members),
        k2p_disp=_mean_tri(k2p_dist_df, members),
        cub_disp=_mean_tri(cub_dist_df, members),
        score=cub_cl_scores.get(cid, np.nan)))
cub_disp_df = pd.DataFrame(cub_disp_rows)

print(f"\nK2P cluster dispersions:\n{k2p_disp_df.round(4).to_string()}")
print(f"\nCUB cluster dispersions:\n{cub_disp_df.round(4).to_string()}")

# ── 6. Build figure ───────────────────────────────────────────────────────────
sns.set(style="whitegrid", context="paper", font_scale=1.1)
np.random.seed(42)

fig = plt.figure(figsize=(14, 11), facecolor=BG)
fig.suptitle("Clustered Optimization Performance", fontsize=16, fontweight="bold", y=0.97)
fig.subplots_adjust(top=0.90)

gs    = GridSpec(2, 2, figure=fig, height_ratios=[1.15, 1.0], hspace=0.52, wspace=0.42)
ax_a1 = fig.add_subplot(gs[0, 0])
ax_a2 = fig.add_subplot(gs[0, 1])
ax_b  = fig.add_subplot(gs[1, 0])   # Panel B: CUB compactness vs opt score
ax_c  = fig.add_subplot(gs[1, 1])   # Panel C: K2P vs CUB dispersion

k2p_cl_sizes = {int(c): k2p_sizes.get(int(c), "?") for c in range(1, 5)}
cub_cl_sizes = {int(c): cub_sizes.get(int(c), "?") for c in range(1, 5)}

CL_PAL      = ["#4E79A7", "#F28E2B", "#E15759", "#76B7B2"]
K2P_COLOR   = "#4E79A7"
CUB_COLOR   = "#E15759"
DARK_GREEN  = "#2d6a27"
ORANGE      = "#F28E2B"
cmap_score  = cm.RdYlGn

all_scores = pd.concat([k2p_disp_df["score"], cub_disp_df["score"]]).dropna()
score_norm = mcolors.Normalize(vmin=all_scores.min(), vmax=all_scores.max())

# Shared Y range for Panel A — cluster scores + per-cluster medians
k2p_meds = [ova.loc[ova["k2p_cluster"] == c, SCORE].dropna().median()
            for c in range(1, 5)]
cub_meds = [ova.loc[ova["cub_cluster"] == c, SCORE].dropna().median()
            for c in range(1, 5)]
all_a_vals = (list(k2p_cl_scores.values()) + list(cub_cl_scores.values()) +
              [v for v in k2p_meds + cub_meds if not np.isnan(v)])
a_ymin = min(all_a_vals) - 0.25
a_ymax = max(all_a_vals) + 0.25

# Bubble size for Panels B/C (proportional to n, no size legend)
all_n = list(k2p_disp_df["n"]) + list(cub_disp_df["n"])
def _bsz(n_val):
    return 80 + 400 * (n_val - min(all_n)) / max(max(all_n) - min(all_n), 1)

# ── Panel A helper: line chart — cluster score vs per-cluster median org score ──
def _panel_a(ax, cl_scores, cl_sizes, cluster_col, prefix, title, panel_label):
    cids    = sorted(int(c) for c in cl_scores)
    x_pos   = list(range(len(cids)))

    ax.set_facecolor(BG)
    ax.grid(True, color="#d4dadb", linestyle="-", linewidth=0.8, zorder=0)

    cl_vals  = [cl_scores[cid] for cid in cids]
    med_vals = [ova.loc[ova[cluster_col] == cid, SCORE].dropna().median()
                for cid in cids]

    ax.plot(x_pos, cl_vals, color=DARK_GREEN, linewidth=2.0,
            marker="s", markersize=8, zorder=4, label="Cluster")
    ax.plot(x_pos, med_vals, color=ORANGE, linewidth=1.8,
            marker="o", markersize=7, linestyle="--", zorder=3,
            label="Individual wanted organisms (median)")

    ax.set_xticks(x_pos)
    ax.set_xticklabels(
        [f"{prefix}-C{c}\n(n={cl_sizes.get(c,'?')})" for c in cids], fontsize=9.5)
    ax.set_ylim(a_ymin, a_ymax)
    ax.set_xlabel("Cluster", fontsize=10, fontweight="bold")
    ax.set_ylabel("Optimization Score", fontsize=10, fontweight="bold")
    ax.set_title(title, fontweight="bold", pad=10, fontsize=11)
    ax.text(-0.10, 1.07, panel_label, transform=ax.transAxes,
            fontsize=14, fontweight="bold", va="bottom", ha="left", clip_on=False)
    ax.legend(fontsize=8, loc="upper right", framealpha=0.9)
    sns.despine(ax=ax)

_panel_a(ax_a1, k2p_cl_scores, k2p_cl_sizes, "k2p_cluster", "K2P",
         "K2P Clustering (k = 4)", "A.")
_panel_a(ax_a2, cub_cl_scores, cub_cl_sizes, "cub_cluster", "CUB",
         "CUB Clustering (k = 4)", "")

# ── shared annotation helper ──────────────────────────────────────────────────
def _ann(ax, x, y, txt, ox, oy):
    ax.annotate(txt, xy=(x, y), xytext=(ox, oy), textcoords="offset points",
                fontsize=8.2, color="#222",
                bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                          edgecolor="none", alpha=0.88),
                zorder=6)

# ── Panel B: CUB compactness vs opt score — K2P (circles) + CUB (triangles) ──
ax_b.set_facecolor(BG)
ax_b.grid(True, color="#d4dadb", linestyle="--", linewidth=0.8, zorder=0)

# Offsets chosen from actual coordinates (see docstring above)
# K2P: C1(0.072,5.80) C2(0.011,1.68) C3(0.091,3.17) C4(0.182,2.64)
# CUB: C1(0.080,6.67) C2(0.064,3.03) C3(0.091,3.17)=same C4(0.000,3.82)
_lab_b_k2p = {1: (9, 6), 2: (9, 6), 3: (9, 14), 4: (9, -18)}
_lab_b_cub = {1: (9, 8), 2: (9, -18), 3: (9, -20), 4: (9, 6)}

for _, row in k2p_disp_df.iterrows():
    cid = int(row.cluster)
    ax_b.scatter(row.cub_disp, row.score, s=_bsz(row.n),
                 color=K2P_COLOR, marker="o", edgecolor="black", linewidth=1.1,
                 zorder=4, alpha=0.85)
    ox, oy = _lab_b_k2p[cid]
    _ann(ax_b, row.cub_disp, row.score, f"K2P-C{cid} (n={int(row.n)})", ox, oy)

for _, row in cub_disp_df.iterrows():
    cid = int(row.cluster)
    if pd.isna(row.score): continue
    ax_b.scatter(row.cub_disp, row.score, s=_bsz(row.n),
                 color=CUB_COLOR, marker="^", edgecolor="black", linewidth=1.1,
                 zorder=4, alpha=0.85)
    ox, oy = _lab_b_cub[cid]
    _ann(ax_b, row.cub_disp, row.score, f"CUB-C{cid} (n={int(row.n)})", ox, oy)

x_vals_b = list(k2p_disp_df["cub_disp"]) + list(cub_disp_df["cub_disp"])
ax_b.set_xlim(-0.005, max(x_vals_b) * 1.45)
ax_b.set_xlabel("Within-cluster CUB dispersion", fontsize=10, fontweight="bold")
ax_b.set_ylabel("Cluster optimization score", fontsize=10, fontweight="bold")
ax_b.set_title("CUB Compactness vs. Optimization Score",
               fontweight="bold", pad=8, fontsize=11)
ax_b.text(-0.14, 1.07, "B.", transform=ax_b.transAxes,
          fontsize=14, fontweight="bold", va="bottom", ha="left", clip_on=False)
ax_b.legend(handles=[
    Line2D([0],[0], marker="o", color="w", markerfacecolor=K2P_COLOR,
           markeredgecolor="black", markersize=9, label="K2P cluster"),
    Line2D([0],[0], marker="^", color="w", markerfacecolor=CUB_COLOR,
           markeredgecolor="black", markersize=9, label="CUB cluster"),
], fontsize=8.5, loc="upper right", framealpha=0.9)
sns.despine(ax=ax_b)

# ── Panel C: K2P vs CUB dispersion + identity line (K2P clusters only) ────────
ax_c.set_facecolor(BG)
ax_c.grid(True, color="#d4dadb", linestyle="--", linewidth=0.8, zorder=0)

# Equal-scale axes to make identity line meaningful
c_lim = max(k2p_disp_df["k2p_disp"].max(), k2p_disp_df["cub_disp"].max()) * 1.18
ax_c.plot([0, c_lim], [0, c_lim], color="#888", linewidth=1.3,
          linestyle="--", zorder=1)
ax_c.set_xlim(0, c_lim)
ax_c.set_ylim(0, c_lim)

# Offsets: C1(0.072,0.072) C2(0.011,0.011) C3(0.015,0.091) C4(0.026,0.182)
_lab_c = {1: (9, 6), 2: (9, -18), 3: (9, 6), 4: (9, 6)}

for _, row in k2p_disp_df.iterrows():
    cid = int(row.cluster)
    col = cmap_score(score_norm(row.score))
    ax_c.scatter(row.k2p_disp, row.cub_disp, s=_bsz(row.n),
                 color=col, marker="o", edgecolor="black", linewidth=1.1,
                 zorder=3, alpha=0.90)
    ox, oy = _lab_c[cid]
    _ann(ax_c, row.k2p_disp, row.cub_disp, f"K2P-C{cid} (n={int(row.n)})", ox, oy)

# Label the identity line
ax_c.text(c_lim * 0.62, c_lim * 0.55, "K2P = CUB",
          fontsize=8, color="#888", rotation=45, va="center", ha="center",
          rotation_mode="anchor")

ax_c.set_xlabel("Within-cluster K2P dispersion", fontsize=10, fontweight="bold")
ax_c.set_ylabel("Within-cluster CUB dispersion", fontsize=10, fontweight="bold")
ax_c.set_title("Cluster Compactness: K2P vs. CUB Space",
               fontweight="bold", pad=8, fontsize=11)
ax_c.text(-0.14, 1.07, "C.", transform=ax_c.transAxes,
          fontsize=14, fontweight="bold", va="bottom", ha="left", clip_on=False)

sm = cm.ScalarMappable(cmap=cmap_score, norm=score_norm)
sm.set_array([])
cb = fig.colorbar(sm, ax=ax_c, shrink=0.72, pad=0.03, aspect=18)
cb.set_label("Optimization score", fontsize=9, fontweight="bold")
sns.despine(ax=ax_c)

out = os.path.join(REPO_ROOT, "graph_comparison_k2p_vs_cub_performance.svg")
fig.savefig(out, bbox_inches="tight", facecolor=BG)
plt.close(fig)
print(f"\nSaved: {out} ({os.path.getsize(out):,} bytes)")
