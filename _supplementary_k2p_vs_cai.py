"""
Supplementary Figure: K2P clustering vs CAI-profile clustering
  Panel A — PCA of CAI profiles, coloured by K2P cluster
  Panel B — K2P distance heatmap with hierarchical-clustering dendrogram;
            left sidebars show K2P and CAI cluster assignments side-by-side;
            organism names as outer-left y-axis labels on the K2P strip;
            cluster scheme label in bold below each colour-coded column.

  CAI clustering = Spearman correlation of cai_weights profiles (same source
  as the DCUB SequenceFamilyModule).
"""
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as mgs
import matplotlib.colors as mcolors
import numpy as np, pandas as pd
from Bio import AlignIO
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.cluster.hierarchy import dendrogram as sp_dend
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
import seaborn as sns

REPO     = "/Users/shimka/PycharmProjects/Igem_TAU_2021"
arab_dir = os.path.join(REPO, "analysis", "example_data", "arabidopsis")
BG       = "#f7f8f5"
CL_PAL   = {1: "#4E79A7", 2: "#F28E2B", 3: "#E15759", 4: "#76B7B2"}

# ── K2P distance matrix + clustering ──────────────────────────────────────────
with open(os.path.join(arab_dir, "sample_to_organism.json")) as f:
    s2o = json.load(f)

def _k2p(s1, s2):
    ts = {("A","G"),("G","A"),("C","T"),("T","C")}
    tv = {("A","C"),("C","A"),("A","T"),("T","A"),
          ("G","C"),("C","G"),("G","T"),("T","G")}
    a, b = str(s1).upper(), str(s2).upper()
    P = Q = v = 0
    for x, y in zip(a, b):
        if x not in "ACGT" or y not in "ACGT": continue
        v += 1
        if x != y:
            if (x, y) in ts: P += 1
            elif (x, y) in tv: Q += 1
    if v == 0: return np.nan
    P, Q = P/v, Q/v
    try:    return -0.5*np.log(1-2*P-Q) - 0.25*np.log(1-2*Q)
    except: return np.inf

al        = AlignIO.read(os.path.join(arab_dir, "arabidopsis_aligned_16s.fasta"), "fasta")
k2p_names = [s2o[r.id] for r in al]
n         = len(al)
K         = np.zeros((n, n))
for i in range(n):
    for j in range(i, n): K[i,j] = K[j,i] = _k2p(al[i].seq, al[j].seq)

Z_k2p      = linkage(squareform(K), method="average")
k2p_lbl    = fcluster(Z_k2p, 4, "maxclust")
k2p_assign = dict(zip(k2p_names, k2p_lbl.astype(int)))
print(f"K2P cluster sizes: {pd.Series(k2p_lbl).value_counts().sort_index().to_dict()}")

# ── CAI clustering (Spearman on cai_weights — same source as DCUB pipeline) ───
cai_proc = os.path.join(arab_dir, "arabidopsis_microbiome_processed")
cai_data, cai_names = [], []
for fn in sorted(os.listdir(cai_proc)):
    if fn.endswith(".json"):
        with open(os.path.join(cai_proc, fn)) as fh: c = json.load(fh)
        cai_data.append(c.get("cai_weights", {}))
        cai_names.append(os.path.splitext(fn)[0])
M_df = pd.DataFrame(cai_data, index=cai_names)
nc = len(cai_names)
Cm = np.zeros((nc, nc))
for i in range(nc):
    for j in range(i+1, nc):
        r, _ = spearmanr(M_df.iloc[i].fillna(0), M_df.iloc[j].fillna(0))
        Cm[i,j] = Cm[j,i] = 1 - r
Z_cai      = linkage(squareform(Cm, checks=False), method="average")
cai_lbl    = fcluster(Z_cai, 4, "maxclust")
cai_assign = dict(zip(cai_names, cai_lbl.astype(int)))
print(f"CAI cluster sizes: {pd.Series(cai_lbl).value_counts().sort_index().to_dict()}")

# ── dendrogram leaf ordering (K2P) ────────────────────────────────────────────
dend_info  = sp_dend(Z_k2p, no_plot=True)
leaf_ord   = dend_info["leaves"]
ord_names  = [k2p_names[i] for i in leaf_ord]
K_ord      = K[np.ix_(leaf_ord, leaf_ord)]
max_dend_h = max(max(y) for y in dend_info["dcoord"])

# ── organisms whose K2P and CAI cluster-neighbour sets differ meaningfully ────
def _jaccard(a, b):
    u = len(a | b)
    return len(a & b) / u if u else 1.0

diff_orgs = set()
for nm in k2p_names:
    k2p_co = frozenset(x for x, cl in k2p_assign.items() if cl == k2p_assign[nm])
    cai_co = frozenset(x for x, cl in cai_assign.items() if cl == cai_assign.get(nm, -1))
    if _jaccard(k2p_co, cai_co) < 0.4:
        diff_orgs.add(nm)
print(f"K2P vs CAI meaningfully differing organisms ({len(diff_orgs)}): {diff_orgs}")

# ── figure ─────────────────────────────────────────────────────────────────────
sns.set(style="whitegrid", context="paper", font_scale=1.0)
fig = plt.figure(figsize=(13, 22), facecolor=BG)

# Three outer rows: PCA | legend strip | heatmap panel
outer = mgs.GridSpec(3, 1, figure=fig, height_ratios=[0.85, 0.07, 2.9],
                     hspace=0.08, left=0.13, right=0.96, top=0.97, bottom=0.02)

# ══ Panel A: PCA of CAI profiles ══════════════════════════════════════════════
ax_pca = fig.add_subplot(outer[0])
ax_pca.set_facecolor(BG)

M_pca  = M_df.reindex(k2p_names).fillna(0)
coords = PCA(n_components=2, random_state=42).fit_transform(M_pca.values)

for cid in sorted(CL_PAL):
    mask = np.array([k2p_assign[nm] == cid for nm in k2p_names])
    ax_pca.scatter(coords[mask, 0], coords[mask, 1],
                   color=CL_PAL[cid], s=65, alpha=0.88, edgecolor="white",
                   linewidth=0.5, zorder=3)

ax_pca.axhline(0, color="#ccc", lw=0.8, zorder=1)
ax_pca.axvline(0, color="#ccc", lw=0.8, zorder=1)
ax_pca.grid(True, color="#d4dadb", linestyle="--", linewidth=0.7, zorder=0)
ax_pca.set_xlabel("Principal Component 1", fontsize=11)
ax_pca.set_ylabel("Principal Component 2", fontsize=11)
ax_pca.set_title("PCA of CAI Profiles by K2P Phylogenetic Cluster",
                 fontsize=12, fontweight="bold", pad=8)
ax_pca.text(-0.07, 1.06, "A.", transform=ax_pca.transAxes,
            fontsize=14, fontweight="bold", va="bottom")
sns.despine(ax=ax_pca)

# ── Legend row — between A and B ───────────────────────────────────────────────
ax_leg = fig.add_subplot(outer[1])
ax_leg.axis("off")
ax_leg.set_facecolor(BG)
leg_handles = [mpatches.Patch(color=CL_PAL[cid], label=f"Cluster {cid}")
               for cid in sorted(CL_PAL)]
ax_leg.legend(handles=leg_handles, title="Phylogenetic Clusters",
              fontsize=9, title_fontsize=9,
              loc="center", framealpha=0.9, ncol=4)

# ══ Panel B: clustermap ═══════════════════════════════════════════════════════
# Width ratios: [k2p_sidebar(+names) | cai_sidebar | heatmap | gap | cbar]
# Dendrogram spans heatmap+gap+cbar (cols 2:5) so it fills the full row width;
# xlim is scaled to keep leaf positions aligned with the heatmap columns.
wr_k2p, wr_cai, wr_heat, wr_gap, wr_cbar = 0.030, 0.030, 1.0, 0.022, 0.052

inner = mgs.GridSpecFromSubplotSpec(
    2, 5, subplot_spec=outer[2],
    height_ratios=[0.42, 1],
    width_ratios=[wr_k2p, wr_cai, wr_heat, wr_gap, wr_cbar],
    hspace=0.0, wspace=0.01,
)
ax_dend    = fig.add_subplot(inner[0, 2])
ax_sid_k2p = fig.add_subplot(inner[1, 0])
ax_sid_cai = fig.add_subplot(inner[1, 1])
ax_heat    = fig.add_subplot(inner[1, 2])
ax_cbar    = fig.add_subplot(inner[1, 4])

# ── Dendrogram — manual plot from pre-computed icoord/dcoord (no clipping) ───
def _sx(x): return (x - 5.0) / 10.0

for xs, ys in zip(dend_info["icoord"], dend_info["dcoord"]):
    ax_dend.plot([_sx(x) for x in xs], ys, color="#333", linewidth=0.8)
ax_dend.set_xlim(-0.5, n - 0.5)
ax_dend.set_ylim(0, max_dend_h * 1.35)
ax_dend.set_facecolor(BG)
for sp in ax_dend.spines.values(): sp.set_visible(False)
ax_dend.set_xticks([])
ax_dend.set_yticks([])
ax_dend.text(-0.048, 1.10, "B.", transform=ax_dend.transAxes,
             fontsize=14, fontweight="bold", va="bottom")
ax_dend.set_title("K2P Distance Heatmap with Hierarchical Clustering",
                  fontsize=12, fontweight="bold", pad=6)

# ── Heatmap ───────────────────────────────────────────────────────────────────
vmax = np.percentile(K_ord[K_ord > 0], 99)
im = ax_heat.imshow(K_ord, cmap="viridis", aspect="auto",
                    vmin=0, vmax=vmax, interpolation="nearest")
ax_heat.set_xticks(range(n))
ax_heat.set_xticklabels(ord_names, rotation=90, fontsize=5, ha="right")
ax_heat.xaxis.set_tick_params(pad=1, length=1.5)
ax_heat.set_yticks([])
ax_heat.set_xlabel("Organisms", fontsize=9, labelpad=4)

prev_cl = None
for i, nm in enumerate(ord_names):
    cl = k2p_assign[nm]
    if cl != prev_cl and prev_cl is not None:
        ax_heat.axhline(i - 0.5, color="white", lw=1.5, zorder=5)
        ax_heat.axvline(i - 0.5, color="white", lw=1.5, zorder=5)
    prev_cl = cl

# ── Colorbar ──────────────────────────────────────────────────────────────────
cb = plt.colorbar(im, cax=ax_cbar)
cb.set_label("K2P distance", fontsize=9)
cb.ax.tick_params(labelsize=7)

# ── sidebar helper ────────────────────────────────────────────────────────────
def _draw_sidebar(ax, assign_dict, label, highlighted_orgs, show_org_names=False):
    rgb = np.array([[mcolors.to_rgb(CL_PAL[assign_dict.get(nm, 1)])] for nm in ord_names])
    ax.imshow(rgb, aspect="auto", interpolation="nearest")
    ax.set_xticks([])
    ax.set_xlabel(label, fontsize=6, labelpad=3, fontweight="bold")

    if show_org_names:
        ax.set_yticks(range(n))
        ax.set_yticklabels(ord_names, fontsize=5)
        ax.yaxis.tick_left()
        ax.yaxis.set_tick_params(length=0, pad=2)
    else:
        ax.set_yticks([])

    # cluster number labels inside the strip
    prev_cl = None; cl_start = 0
    for i, nm in enumerate(ord_names):
        cl = assign_dict.get(nm, 1)
        if cl != prev_cl:
            if prev_cl is not None:
                ax.text(0, (cl_start + i - 1) / 2, str(prev_cl),
                        ha="center", va="center", fontsize=6,
                        color="white", fontweight="bold")
            prev_cl = cl; cl_start = i
    ax.text(0, (cl_start + n - 1) / 2, str(prev_cl),
            ha="center", va="center", fontsize=6,
            color="white", fontweight="bold")

    for i, nm in enumerate(ord_names):
        if nm in highlighted_orgs:
            ax.add_patch(mpatches.Rectangle(
                (-0.5, i - 0.5), 1, 1,
                fill=False, edgecolor="white", linewidth=2.0, zorder=6))

# K2P sidebar carries the organism name labels on its left y-axis
_draw_sidebar(ax_sid_k2p, k2p_assign, "K2P", set(), show_org_names=True)
_draw_sidebar(ax_sid_cai, cai_assign,  "CAI", diff_orgs)

# ── save ──────────────────────────────────────────────────────────────────────
out = os.path.join(REPO, "_supplementary_k2p_vs_cai.svg")
fig.savefig(out, bbox_inches="tight", facecolor=BG)
plt.close(fig)
print(f"\nSaved: {out} ({os.path.getsize(out):,} bytes)")
