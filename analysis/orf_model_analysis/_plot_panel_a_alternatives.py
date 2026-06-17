"""
Generate three standalone Panel-A alternatives for evaluation.
Run with:  .venv/bin/python3 _plot_panel_a_alternatives.py
Outputs (PNG, 150 DPI) are written next to this script.
"""
import json, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr, pearsonr, gaussian_kde as _gkde
from Bio import AlignIO

REPO_ROOT = os.path.dirname(os.path.abspath(__file__))
arab_dir  = os.path.join(REPO_ROOT, "analysis", "example_data", "arabidopsis")

LABEL_FONT = 18
TITLE_FONT = 20
ANNOT_FONT = 15
TICK_FONT  = 15
BOLD       = "bold"
BG         = "#FAFAFA"

sns.set(style="whitegrid", context="talk",
        rc={"axes.facecolor": BG, "grid.color": "#E0E0E0"})

# ── Build K2P matrix ──────────────────────────────────────────────────────────
print("Building K2P matrix …")
with open(os.path.join(arab_dir, "sample_to_organism.json")) as f:
    sample_to_org = json.load(f)

def _k2p(seq1, seq2):
    ts = {("A","G"),("G","A"),("C","T"),("T","C")}
    tv = {("A","C"),("C","A"),("A","T"),("T","A"),
          ("G","C"),("C","G"),("G","T"),("T","G")}
    P = Q = valid = 0
    for a, b in zip(str(seq1).upper(), str(seq2).upper()):
        if a not in "ACGT" or b not in "ACGT": continue
        valid += 1
        if a != b:
            if (a,b) in ts: P += 1
            elif (a,b) in tv: Q += 1
    if valid == 0: return np.nan
    P /= valid; Q /= valid
    try:    return -0.5*np.log(1-2*P-Q) - 0.25*np.log(1-2*Q)
    except: return np.inf

alignment = AlignIO.read(os.path.join(arab_dir, "arabidopsis_aligned_16s.fasta"), "fasta")
seq_ids   = [r.id for r in alignment]
org_names = [sample_to_org[s] for s in seq_ids]
n = len(alignment)
k2p_mat = np.zeros((n, n))
for i in range(n):
    for j in range(i, n):
        d = _k2p(alignment[i].seq, alignment[j].seq)
        k2p_mat[i,j] = k2p_mat[j,i] = d
pairwise_df = pd.DataFrame(k2p_mat, index=org_names, columns=org_names)

# ── Build CAI matrix ──────────────────────────────────────────────────────────
print("Building CAI matrix …")
cai_proc = os.path.join(arab_dir, "arabidopsis_microbiome_processed")
cai_data, cai_names = [], []
for fname in sorted(os.listdir(cai_proc)):
    if fname.endswith(".json"):
        with open(os.path.join(cai_proc, fname)) as f:
            c = json.load(f)
        cai_data.append(c.get("cai_weights", {}))
        cai_names.append(os.path.splitext(fname)[0])
cai_weights = pd.DataFrame(cai_data, index=cai_names)
m  = cai_weights.values
nn = m.shape[0]
cai_mat = np.zeros((nn, nn))
for i in range(nn):
    for j in range(i+1, nn):
        r, _ = spearmanr(m[i], m[j])
        cai_mat[i,j] = cai_mat[j,i] = 1.0 - r
cai_df = pd.DataFrame(cai_mat, index=cai_names, columns=cai_names)

# ── Flatten upper triangle ────────────────────────────────────────────────────
common = cai_df.index.intersection(pairwise_df.index)
d1 = cai_df.loc[common, common]
d2 = pairwise_df.loc[common, common]
triu    = np.triu_indices_from(d1, k=1)
d1_flat = d1.values[triu]   # CAI distance
d2_flat = d2.values[triu]   # K2P distance
pear,  _ = pearsonr(d1_flat, d2_flat)
spear, _ = spearmanr(d1_flat, d2_flat)
stats_txt = f"Pearson r = {pear:.2f}\nSpearman ρ = {spear:.2f}"
print(f"  {len(d1_flat)} pairs ready")

ANNOT_KW = dict(
    xy=(0.05, 0.95), xycoords="axes fraction",
    ha="left", va="top", fontsize=ANNOT_FONT, fontweight=BOLD,
    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.88),
)

# ════════════════════════════════════════════════════════════════════════════════
# Option 1 – Hexbin
# ════════════════════════════════════════════════════════════════════════════════
print("Rendering Option 1: Hexbin …")
fig1, ax1 = plt.subplots(figsize=(7, 6), facecolor=BG)
ax1.set_facecolor(BG)
hb = ax1.hexbin(d1_flat, d2_flat, gridsize=28, cmap="Blues",
                mincnt=1, linewidths=0.3)
plt.colorbar(hb, ax=ax1, label="Count").ax.tick_params(labelsize=TICK_FONT - 2)
ax1.annotate(stats_txt, **ANNOT_KW)
ax1.set_xlabel("CAI distance",      fontsize=LABEL_FONT, fontweight=BOLD)
ax1.set_ylabel("Kimura-2 distance", fontsize=LABEL_FONT, fontweight=BOLD)
ax1.set_title("A.  K2P vs CAI — Hexbin density",
              fontsize=TITLE_FONT, fontweight=BOLD, pad=10, loc="left")
ax1.tick_params(labelsize=TICK_FONT)
sns.despine(ax=ax1)
fig1.tight_layout()
figures_dir = os.path.join(REPO_ROOT, "analysis", "results", "arabidopsis", "figures")
os.makedirs(figures_dir, exist_ok=True)
out1 = os.path.join(figures_dir, "panel_a_option1_hexbin.png")
fig1.savefig(out1, dpi=150, bbox_inches="tight")
plt.close(fig1)
print(f"  → {out1}")

# ════════════════════════════════════════════════════════════════════════════════
# Option 2 – KDE contour-only (no individual dots)
# ════════════════════════════════════════════════════════════════════════════════
print("Rendering Option 2: KDE contour-fill …")
fig2, ax2 = plt.subplots(figsize=(7, 6), facecolor=BG)
ax2.set_facecolor(BG)

xg = np.linspace(d1_flat.min(), d1_flat.max(), 120)
yg = np.linspace(d2_flat.min(), d2_flat.max(), 120)
XX, YY = np.meshgrid(xg, yg)
kde    = _gkde(np.vstack([d1_flat, d2_flat]))
ZZ     = kde(np.vstack([XX.ravel(), YY.ravel()])).reshape(XX.shape)

cf = ax2.contourf(XX, YY, ZZ, levels=12, cmap="Blues", alpha=0.85)
ax2.contour( XX, YY, ZZ, levels=12, colors="#1A3A5A", linewidths=0.7, alpha=0.55)
plt.colorbar(cf, ax=ax2, label="Density").ax.tick_params(labelsize=TICK_FONT - 2)
ax2.scatter(d1_flat, d2_flat, s=12, color="#333", alpha=0.18, edgecolors="none", zorder=3)
ax2.annotate(stats_txt, **ANNOT_KW)
ax2.set_xlabel("CAI distance",      fontsize=LABEL_FONT, fontweight=BOLD)
ax2.set_ylabel("Kimura-2 distance", fontsize=LABEL_FONT, fontweight=BOLD)
ax2.set_title("A.  K2P vs CAI — KDE contour fill",
              fontsize=TITLE_FONT, fontweight=BOLD, pad=10, loc="left")
ax2.tick_params(labelsize=TICK_FONT)
sns.despine(ax=ax2)
fig2.tight_layout()
out2 = os.path.join(figures_dir, "panel_a_option2_kde_contour.png")
fig2.savefig(out2, dpi=150, bbox_inches="tight")
plt.close(fig2)
print(f"  → {out2}")

# ════════════════════════════════════════════════════════════════════════════════
# Option 3 – Joint distribution (scatter + marginal histograms)
# ════════════════════════════════════════════════════════════════════════════════
print("Rendering Option 3: Joint distribution …")
_xy   = np.vstack([d1_flat, d2_flat])
_z    = _gkde(_xy)(_xy)
_ord  = _z.argsort()

g = sns.JointGrid(x=d1_flat, y=d2_flat, height=7, ratio=5)
g.figure.set_facecolor(BG)
g.ax_joint.set_facecolor(BG)

sc = g.ax_joint.scatter(
    d1_flat[_ord], d2_flat[_ord],
    c=_z[_ord], cmap="Blues", s=22, alpha=0.65, edgecolors="none",
)
plt.colorbar(sc, ax=g.ax_joint, label="Density",
             shrink=0.82, pad=0.02).ax.tick_params(labelsize=TICK_FONT - 2)

xg2 = np.linspace(d1_flat.min(), d1_flat.max(), 120)
yg2 = np.linspace(d2_flat.min(), d2_flat.max(), 120)
XX2, YY2 = np.meshgrid(xg2, yg2)
ZZ2 = _gkde(_xy)(np.vstack([XX2.ravel(), YY2.ravel()])).reshape(XX2.shape)
g.ax_joint.contour(XX2, YY2, ZZ2, levels=5, colors="#E6550D",
                   linewidths=0.9, alpha=0.6)

g.ax_marg_x.hist(d1_flat, bins=28, color="#4393C3", alpha=0.75, edgecolor="white", lw=0.4)
g.ax_marg_y.hist(d2_flat, bins=28, color="#4393C3", alpha=0.75, edgecolor="white", lw=0.4,
                 orientation="horizontal")

g.ax_joint.annotate(stats_txt, **ANNOT_KW)
g.ax_joint.set_xlabel("CAI distance",      fontsize=LABEL_FONT, fontweight=BOLD)
g.ax_joint.set_ylabel("Kimura-2 distance", fontsize=LABEL_FONT, fontweight=BOLD)
g.ax_joint.tick_params(labelsize=TICK_FONT)
g.ax_marg_x.set_ylabel("Count", fontsize=TICK_FONT - 2)
g.ax_marg_y.set_xlabel("Count", fontsize=TICK_FONT - 2)
g.figure.suptitle("A.  K2P vs CAI — Joint distribution",
                   fontsize=TITLE_FONT, fontweight=BOLD, y=1.01)
g.figure.tight_layout()
out3 = os.path.join(figures_dir, "panel_a_option3_joint.png")
g.figure.savefig(out3, dpi=150, bbox_inches="tight")
plt.close(g.figure)
print(f"  → {out3}")

print("\nDone. Three PNGs written to the repo root.")
