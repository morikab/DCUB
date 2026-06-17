"""
Regenerate the "Model Performance Across Phylogenetic Divergence and Host Group Sizes"
unified figure (v3) with:
  - Title embedded inside the figure (not clipped by viewport)
  - Font / style matching the homo_sapiens comparative-analysis figures
  - Panel A: density-scatter + KDE contour overlay (cleaner than plain scatter)
  - Panel B: scatter + LOWESS trend (unchanged)
  - Panel C: strip + mean-marker (unchanged)

Run with:  .venv/bin/python3 _plot_unified_figure.py
"""
import json, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns
from scipy.stats import spearmanr, pearsonr, gaussian_kde as _gkde
from scipy.stats import stats as scipy_stats  # noqa – not used directly
from scipy import stats
from Bio import AlignIO
import statsmodels.api as sm

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
arab_dir    = os.path.join(REPO_ROOT, "analysis", "example_data", "arabidopsis")
proc_dir    = os.path.join(REPO_ROOT, "analysis", "results", "arabidopsis", "processed")
figures_dir = os.path.join(REPO_ROOT, "analysis", "results", "arabidopsis", "figures")
os.makedirs(figures_dir, exist_ok=True)

# ── Font settings (matches homo_sapiens comparative-analysis figures) ─────────
BASE_FONT   = 18
TITLE_FONT  = 22   # sub-panel titles
SUPTITLE_FONT = 28 # main figure title  (≈ TITLE_FONT+6 from the homo-sapiens notebook)
LABEL_FONT  = 18
TICK_FONT   = 16
LEGEND_FONT = 15
ANNOT_FONT  = 15
BOLD        = "bold"

plt.rcParams.update({
    "font.size":             BASE_FONT,
    "axes.titlesize":        TITLE_FONT,
    "axes.labelsize":        LABEL_FONT,
    "xtick.labelsize":       TICK_FONT,
    "ytick.labelsize":       TICK_FONT,
    "legend.fontsize":       LEGEND_FONT,
    "font.weight":           "bold",
})
sns.set(style="whitegrid", context="talk",
        rc={"axes.facecolor": "#FAFAFA", "grid.color": "#E0E0E0"})

# ── Palette (Set2 accents, matching homo-sapiens figures) ─────────────────────
_SCATTER_CMAP  = "Blues"           # density-coloured scatter cmap
_SCATTER_ALPHA = 0.65
_TREND_COLOR   = "#1A6B5A"        # dark teal  (LOWESS line)
_STRIP_COLOR   = "#4393C3"        # mid-blue   (strip dots)
_MEAN_COLOR    = "#D95F02"        # Set2 orange (mean marker) – pops against blue strip
_ERR_COLOR     = "#7570B3"        # Set2 purple (error bars)

# ── 1. Build K2P pairwise distance matrix ─────────────────────────────────────
print("Building K2P distance matrix …")
with open(os.path.join(arab_dir, "sample_to_organism.json")) as f:
    sample_to_org = json.load(f)

def _k2p(seq1, seq2):
    ts_set = {("A","G"),("G","A"),("C","T"),("T","C")}
    tv_set = {("A","C"),("C","A"),("A","T"),("T","A"),
              ("G","C"),("C","G"),("G","T"),("T","G")}
    P = Q = valid = 0
    for a, b in zip(str(seq1).upper(), str(seq2).upper()):
        if a not in "ACGT" or b not in "ACGT":
            continue
        valid += 1
        if a != b:
            if (a, b) in ts_set: P += 1
            elif (a, b) in tv_set: Q += 1
    if valid == 0:
        return np.nan
    P /= valid; Q /= valid
    try:    return -0.5*np.log(1 - 2*P - Q) - 0.25*np.log(1 - 2*Q)
    except: return np.inf

alignment = AlignIO.read(os.path.join(arab_dir, "arabidopsis_aligned_16s.fasta"), "fasta")
seq_ids   = [r.id for r in alignment]
org_names = [sample_to_org[sid] for sid in seq_ids]
n = len(alignment)
k2p_mat = np.zeros((n, n))
for i in range(n):
    for j in range(i, n):
        d = _k2p(alignment[i].seq, alignment[j].seq)
        k2p_mat[i, j] = k2p_mat[j, i] = d
pairwise_distances_df = pd.DataFrame(k2p_mat, index=org_names, columns=org_names)
print(f"  K2P matrix: {n}×{n}")

# ── 2. Build CAI distance matrix ──────────────────────────────────────────────
print("Building CAI distance matrix …")
cai_proc = os.path.join(arab_dir, "arabidopsis_microbiome_processed")
cai_data, cai_names = [], []
for fname in sorted(os.listdir(cai_proc)):
    if fname.endswith(".json"):
        with open(os.path.join(cai_proc, fname)) as f:
            c = json.load(f)
        cai_data.append(c.get("cai_weights", {}))
        cai_names.append(os.path.splitext(fname)[0])

cai_weights_df = pd.DataFrame(cai_data, index=cai_names)
m = cai_weights_df.values
nn = m.shape[0]
cai_mat = np.zeros((nn, nn))
for i in range(nn):
    for j in range(i+1, nn):
        r, _ = spearmanr(m[i], m[j])
        cai_mat[i, j] = cai_mat[j, i] = 1.0 - r
cai_distance_df = pd.DataFrame(cai_mat, index=cai_names, columns=cai_names)
print(f"  CAI matrix: {nn}×{nn}")

# ── 3. d1_flat / d2_flat for Panel A ─────────────────────────────────────────
common_orgs  = cai_distance_df.index.intersection(pairwise_distances_df.index)
d1 = cai_distance_df.loc[common_orgs, common_orgs]
d2 = pairwise_distances_df.loc[common_orgs, common_orgs]
triu = np.triu_indices_from(d1, k=1)
d1_flat = d1.values[triu]
d2_flat = d2.values[triu]
print(f"  Panel-A pairs: {len(d1_flat)}")

# ── 4. Load pairwise results → Panel B ────────────────────────────────────────
res = pd.read_pickle(os.path.join(proc_dir, "pairwise.pkl"))
res["organism_distance"] = res.apply(
    lambda row: pairwise_distances_df[row.wanted_organisms][row.unwanted_organisms],
    axis=1
)
print(f"  Panel-B rows: {len(res)}")

# ── 5. Load size-splits results → Panel C ────────────────────────────────────
new_res = pd.read_pickle(os.path.join(proc_dir, "size_splits.pkl"))
print(f"  Panel-C rows: {len(new_res)}")

# ── 6. Build figure ────────────────────────────────────────────────────────────
print("Rendering figure …")
np.random.seed(42)

fig = plt.figure(figsize=(16, 15))
# Leave 10 % of height at the top so the suptitle sits inside the figure
fig.subplots_adjust(top=0.90, bottom=0.08, left=0.08, right=0.97,
                    hspace=0.58, wspace=0.42)

fig.suptitle(
    "Model Performance Across Phylogenetic Divergence and Host Group Sizes",
    fontsize=SUPTITLE_FONT, fontweight=BOLD,
    y=0.97,          # inside the figure (< 1.0) – always captured by the SVG viewport
)

gs   = GridSpec(2, 4, figure=fig, wspace=0.42, hspace=0.58)
ax_A = fig.add_subplot(gs[0, 0:2])
ax_B = fig.add_subplot(gs[0, 2:4])
ax_C = fig.add_subplot(gs[1, 1:3])

# ── Panel A: K2P vs CAI – density-coloured scatter ───────────────────────────
ax = ax_A
_xy = np.vstack([d1_flat, d2_flat])
_z  = _gkde(_xy)(_xy)
_ord = _z.argsort()
sc_A = ax.scatter(
    d1_flat[_ord], d2_flat[_ord],
    c=_z[_ord], cmap=_SCATTER_CMAP,
    s=28, alpha=_SCATTER_ALPHA, edgecolors="none",
)
cb = plt.colorbar(sc_A, ax=ax, label="Density", shrink=0.78, pad=0.02)
cb.ax.tick_params(labelsize=TICK_FONT - 2)

pear_A, _  = pearsonr(d1_flat, d2_flat)
spear_A, _ = spearmanr(d1_flat, d2_flat)
ax.annotate(
    f"Pearson  r = {pear_A:.2f}\nSpearman ρ = {spear_A:.2f}",
    xy=(0.05, 0.95), xycoords="axes fraction",
    ha="left", va="top", fontsize=ANNOT_FONT, fontweight=BOLD,
    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.88),
)
ax.set_xlabel("CAI distance",      fontsize=LABEL_FONT, fontweight=BOLD)
ax.set_ylabel("Kimura-2 distance", fontsize=LABEL_FONT, fontweight=BOLD)
ax.set_title("A.  K2P vs CAI pairwise distances",
             fontsize=TITLE_FONT, fontweight=BOLD, pad=12, loc="left")
ax.tick_params(axis="both", labelsize=TICK_FONT)

# ── Panel B: Evaluation Score vs K2P – scatter + LOWESS ──────────────────────
ax = ax_B
score_B  = "average_distance_score_zscore_bulk_ratio"
res_copy = res.copy()
r_B, _   = spearmanr(res_copy["organism_distance"], res_copy[score_B])

ax.scatter(res_copy["organism_distance"], res_copy[score_B],
           alpha=0.30, color=_STRIP_COLOR, s=25, edgecolors="none")

_lowess = sm.nonparametric.lowess(
    res_copy[score_B].values,
    res_copy["organism_distance"].values,
    frac=0.25, return_sorted=True,
)
ax.plot(_lowess[:, 0], _lowess[:, 1],
        color=_TREND_COLOR, linewidth=2.8, label="Trend")

ax.annotate(
    f"Spearman ρ = {r_B:.2f}",
    xy=(0.97, 0.05), xycoords="axes fraction",
    ha="right", va="bottom", fontsize=ANNOT_FONT, fontweight=BOLD,
    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.88),
)
ax.set_xlabel("Kimura-2 distance",    fontsize=LABEL_FONT, fontweight=BOLD)
ax.set_ylabel("Evaluation Score",     fontsize=LABEL_FONT, fontweight=BOLD)
ax.set_title("B.  Evaluation Score vs K2P distance",
             fontsize=TITLE_FONT, fontweight=BOLD, pad=12, loc="left")
ax.tick_params(axis="both", labelsize=TICK_FONT)
ax.legend(fontsize=LEGEND_FONT, framealpha=0.88, loc="upper left")

# ── Panel C: Evaluation Score vs group size – strip + mean marker ─────────────
ax = ax_C
score_col_C = "average_distance_score_zscore_bulk_diff"
df_plot_C   = new_res.copy()
df_plot_C["wanted_count"] = df_plot_C["wanted_organisms"].str.split("|").apply(len)
df_plot_C["total_count"]  = df_plot_C["wanted_count"] * 2

_rng = np.random.default_rng(42)
for _gc, _sub in df_plot_C.groupby("total_count"):
    _jitter = _rng.uniform(-0.18, 0.18, size=len(_sub))
    ax.scatter(
        np.full(len(_sub), _gc) + _jitter,
        _sub[score_col_C],
        s=22, alpha=0.35, color=_STRIP_COLOR, edgecolors="none", zorder=2,
    )

df_summary_C = df_plot_C.groupby("total_count")[score_col_C].agg(["mean", "std", "count"])
df_summary_C["se"] = df_summary_C["std"] / df_summary_C["count"] ** 0.5
ax.scatter(df_summary_C.index, df_summary_C["mean"],
           s=110, color=_MEAN_COLOR, zorder=4, label="Mean", marker="D")
ax.errorbar(df_summary_C.index, df_summary_C["mean"], yerr=df_summary_C["se"],
            fmt="none", color=_ERR_COLOR, capsize=6, linewidth=2, zorder=5)

ax.set_xticks(sorted(df_plot_C["total_count"].unique()))
ax.tick_params(axis="both", labelsize=TICK_FONT)
ax.set_xlabel("Number of organisms (wanted + unwanted)", fontsize=LABEL_FONT, fontweight=BOLD)
ax.set_ylabel("Evaluation Score",                        fontsize=LABEL_FONT, fontweight=BOLD)
ax.set_title("C.  Evaluation Score vs group size",
             fontsize=TITLE_FONT, fontweight=BOLD, pad=12, loc="left")
ax.legend(fontsize=LEGEND_FONT, framealpha=0.88, loc="upper right")
ax.yaxis.grid(True, alpha=0.35, zorder=0)
ax.set_axisbelow(True)

# ── Save ──────────────────────────────────────────────────────────────────────
out_svg = os.path.join(figures_dir, "unified_figure_ABC_v2.svg")
fig.savefig(out_svg, format="svg", bbox_inches="tight")
plt.close(fig)
print(f"\nSaved → {out_svg}  ({os.path.getsize(out_svg):,} bytes)")
