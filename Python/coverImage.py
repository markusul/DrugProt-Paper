# Cover-style striking image for the Drug-Prot paper.
# A cell after drug treatment, with protein expression at 6, 24 and 48 hours as three rings
# in blue, orange and green. The lines are the significant dependencies among the 300 proteins
# with the most dependencies. One protein is selected as in the application. Its dependencies
# are highlighted, and the panel shows which drugs change its expression and which drug pairs interact.
# The cell, the drugs and the order of the proteins around the rings are an illustration.
#
# The script reads the evidence files of the application directly (the folder with the Parquet files).
# Usage: python3 coverImage.py data/parquet figures/CoverImage --protein SPRR3
# The submitted image uses the protein SPRR3 and the default 600 dpi (--dpi).
# Without --protein, the script selects the protein with the most significant drug effects and interactions
# in the panel among the proteins with at least 15 dependencies in the image.
# It prints the name of the selected protein and other candidates for --protein.
# It writes <name>.png and <name>.tif, 7.5 by 5 in, and needs numpy, matplotlib, opencv-python-headless, pillow and duckdb.
#
# Caption for the submission, one sentence per line, <protein> is the name that the script prints
#   Drug-Prot follows the proteome of breast cancer cells after drug treatment and makes the results available in an interactive application.
#   Two drugs act on a cell, and the three rings show protein expression at 6, 24 and 48 hours in blue, orange and green.
#   The lines are the significant dependencies between consecutive time points among the 300 proteins with the most dependencies.
#   The protein <protein> is selected as in the application, and its dependencies are highlighted.
#   The panel shows the adjusted p-values of the drug effects on this protein for the twelve drugs that were also given in pairs, with single drugs on the diagonal, drug pairs off the diagonal, and small p-values in white.
#   The image is drawn from the estimated dependencies and drug effects by a Python script written with Claude (Anthropic).
import argparse
import os
import re
import numpy as np
import duckdb
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon
import cv2
from PIL import Image

parser = argparse.ArgumentParser()
parser.add_argument("evidence", help="folder with the Parquet files of the application")
parser.add_argument("out", help="name of the output files without extension")
parser.add_argument("--protein", default=None, help="name of the selected protein")
parser.add_argument("--dpi", type=int, default=600)
args = parser.parse_args()
out, dpi = args.out, args.dpi
W, H = 7.5, 5.0
N = 300
ALPHA = 0.05
rng = np.random.default_rng(7)

BLUE, ORANGE, GREEN = np.array([0.31, 0.62, 1.0]), np.array([1.0, 0.54, 0.30]), np.array([0.25, 0.88, 0.63])
TIME = [BLUE, ORANGE, GREEN]
PALE = np.array([0.80, 0.89, 1.0])


# ---------- the evidence ----------
con = duckdb.connect()
def table(name):
    return "read_parquet('%s')" % os.path.join(args.evidence, name + ".parquet").replace("'", "''")

protein_name = con.execute(f"SELECT name FROM {table('proteins')} ORDER BY protein").fetchnumpy()["name"]
n_prot = len(protein_name) + 1

# significant dependencies of the whole network, Benjamini-Hochberg over all tests as in R/pvalues/pValVis.R
# only p-values below ALPHA can be significant, so only those are ranked
n_tests = con.execute(f"SELECT count(*) FROM {table('protein_pvalue')}").fetchone()[0]
threshold = con.execute(f"""
    SELECT max(pvalue) FROM (SELECT pvalue, row_number() OVER (ORDER BY pvalue) AS k
                             FROM {table('protein_pvalue')} WHERE pvalue <= {ALPHA})
    WHERE pvalue <= k * {ALPHA} / {n_tests}""").fetchone()[0]
links = con.execute(f"SELECT source, target, time_idx FROM {table('protein_pvalue')} WHERE pvalue <= {threshold!r}").fetchnumpy()
src, tgt, tt = links["source"].astype(int), links["target"].astype(int), links["time_idx"].astype(int) - 1
print(f"{len(src)} significant dependencies among {n_tests} tests")

# ---------- the network ----------
degree = np.bincount(src, minlength=n_prot) + np.bincount(tgt, minlength=n_prot)
# the N proteins with the most dependencies, a protein that is selected by name is always among them
top = np.argsort(-degree, kind="stable")[:N]
if args.protein is not None:
    wanted = np.where(protein_name == args.protein)[0]
    if len(wanted) == 0:
        raise SystemExit(f"protein {args.protein} is not in the evidence files")
    if wanted[0] + 1 not in top:
        top[-1] = wanted[0] + 1
index = -np.ones(n_prot, int)
index[top] = np.arange(N)
keep = (index[src] >= 0) & (index[tgt] >= 0)
ei, ej, et = index[src[keep]], index[tgt[keep]], tt[keep]

# order of the proteins around the rings, so that dependent proteins are close to each other
# start from a spectral embedding and move each protein towards its neighbours a number of times
A = np.zeros((N, N))
cross = ei != ej
np.add.at(A, (ei[cross], ej[cross]), 1)
A = A + A.T
d = np.maximum(A.sum(1), 1)
vals, vecs = np.linalg.eigh(np.eye(N) - A / np.sqrt(d[:, None] * d[None, :]))
first = int(np.sum(vals < 1e-9))
f, g = vecs[:, first] / np.sqrt(d), vecs[:, first + 1] / np.sqrt(d)
pos = np.empty(N)
pos[np.argsort(np.arctan2(g - np.median(g), f - np.median(f)))] = np.arange(N)
for _ in range(200):
    a = 2 * np.pi * pos / N
    cx, cy = 0.8 * np.cos(a), 0.8 * np.sin(a)
    np.add.at(cx, ei, np.cos(a[ej])); np.add.at(cy, ei, np.sin(a[ej]))
    np.add.at(cx, ej, np.cos(a[ei])); np.add.at(cy, ej, np.sin(a[ei]))
    pos[np.argsort(np.arctan2(cy, cx))] = np.arange(N)

in_drawing = np.bincount(ei[cross], minlength=N) + np.bincount(ej[cross], minlength=N)


# ---------- drug effects on the proteins in the image ----------
# names of the drugs of a treatment, as in replace_drug_ids() of the application
lookup = dict(con.execute(f"SELECT id, name FROM {table('drug_lookup')}").fetchall())
def drugs_of(label):
    ids = [i for i in re.split("_|`", label) if i not in ("", "drug", ":")]
    return tuple(lookup.get(i, i) if "#" in i or i.replace("#", "").replace(".", "", 1).isdigit() else i for i in ids)

treatments = [drugs_of(r[0]) for r in con.execute(f"SELECT raw_label FROM {table('treatments')} ORDER BY treatment").fetchall()]
drug_order = [drugs_of(r[0])[0] for r in con.execute(f"SELECT raw_label FROM {table('drug_order')} ORDER BY ord").fetchall()]
# the panel shows the drugs that were also given in pairs, in the order of the application
in_pairs = {d_ for tr in treatments if len(tr) == 2 for d_ in tr}
panel_drugs = [d_ for d_ in drug_order if d_ in in_pairs]
K = len(panel_drugs)
# cell of the panel for each treatment, single drugs are on the diagonal
cells = [(k, panel_drugs.index(tr[0]), panel_drugs.index(tr[-1]))
         for k, tr in enumerate(treatments) if all(d_ in panel_drugs for d_ in tr)]

# p-values of the direct effects on each protein in the image
# Holm over all treatments and time points of the protein, then the minimum over time, as the application does for one protein
raw = con.execute(f"""SELECT protein, treatment, pvalue FROM {table('drug_pvalue')}
                      WHERE effect_type = 'direct' AND protein IN ({', '.join(str(i) for i in top)})""").fetchnumpy()
p_treat = np.full((N, len(treatments)), np.nan)
for k in range(N):
    rows = raw["protein"] == top[k]
    p_, tr_ = raw["pvalue"][rows], raw["treatment"][rows].astype(int) - 1
    o = np.argsort(p_, kind="stable")
    holm = np.empty(len(o))
    holm[o] = np.maximum.accumulate(np.minimum(1, (len(o) - np.arange(len(o))) * p_[o]))
    np.fmin.at(p_treat[k], tr_, holm)
n_single = np.array([sum(p_treat[k, c] < ALPHA for c, a_, b_ in cells if a_ == b_) for k in range(N)])
n_pair = np.array([sum(p_treat[k, c] < ALPHA for c, a_, b_ in cells if a_ != b_) for k in range(N)])


# ---------- the selected protein ----------
# without --protein, it is the protein with the most significant drug effects and interactions in the panel
# among the proteins with at least MIN_DEPENDENCIES dependencies in the image
# ties go to the protein with more dependencies
MIN_DEPENDENCIES = 15
ranking = np.lexsort((-in_drawing, -(n_single + n_pair) * (in_drawing >= MIN_DEPENDENCIES)))
selected = int(ranking[0]) if args.protein is None else int(index[wanted[0] + 1])
selected_id = int(top[selected])
describe = lambda k: (f"{protein_name[top[k] - 1]}: {in_drawing[k]} dependencies in the image, "
                      f"significant in the panel are {n_single[k]} single drugs and {n_pair[k]} pairs")
print("selected protein " + describe(selected))
print("other candidates for --protein")
for k in [k for k in ranking[:9] if k != selected][:8]:
    print("  " + describe(k))
p_mat = np.full((K, K), np.nan)
for c, a_, b_ in cells:
    p_mat[a_, b_] = p_mat[b_, a_] = p_treat[selected, c]
# the selected protein sits on the spoke
SPOKE = np.deg2rad(-4)
ang = SPOKE + 2 * np.pi * (pos - pos[selected]) / N


# ---------- geometry ----------
CX, CY, R, PH = 3.64, 2.5, 1.95, 0.9
RING = [0.36, 0.60, 0.84]

def membrane(a, scale=1.0, cx=None, r=R, ph=PH):
    return r * scale * (1 + 0.13 * np.cos(2 * a + ph) + 0.09 * np.cos(3 * a + 2 * ph) + 0.05 * np.cos(5 * a + 3 * ph))

def polar(a, f):
    r = membrane(a) * f
    return CX + r * np.cos(a), CY + r * np.sin(a)

def edge_path(t, i, j, n=48):
    a1, a2 = ang[i], ang[j]
    da = (a2 - a1 + np.pi) % (2 * np.pi) - np.pi
    u = np.linspace(0, 1, n)
    a = a1 + da * u * u * (3 - 2 * u)
    x, y = polar(a, RING[t] + (RING[t + 1] - RING[t]) * u)
    return np.c_[x, y], u, abs(da)

def bezier(p0, p1, p2, p3, n=40):
    u = np.linspace(0, 1, n)[:, None]
    return (1 - u) ** 3 * np.array(p0) + 3 * (1 - u) ** 2 * u * np.array(p1) + 3 * (1 - u) * u ** 2 * np.array(p2) + u ** 3 * np.array(p3)

def round_rect(x0, x1, y0, y1, r, n=12):
    a = np.linspace(0, np.pi / 2, n)
    return np.c_[np.r_[x1 - r + r * np.cos(a), x0 + r - r * np.sin(a), x0 + r - r * np.cos(a), x1 - r + r * np.sin(a)],
                 np.r_[y1 - r + r * np.sin(a), y1 - r + r * np.cos(a), y0 + r - r * np.sin(a), y0 + r - r * np.cos(a)]]


# ---------- drawing helpers ----------
def new_layer():
    fig = plt.figure(figsize=(W, H), dpi=dpi)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W); ax.set_ylim(0, H); ax.set_aspect("equal"); ax.axis("off")
    fig.patch.set_facecolor("black"); ax.set_facecolor("black")
    return fig, ax

def grab(fig):
    fig.canvas.draw()
    a = np.asarray(fig.canvas.buffer_rgba())[:, :, :3].astype(np.float32) / 255.0
    plt.close(fig)
    return a

def blur(a, sigma_in):
    return cv2.GaussianBlur(a, (0, 0), sigma_in * dpi)

def glow(layer, sharp, *halos):
    out_ = sharp * layer
    for weight, sigma in halos:
        out_ += weight * blur(layer, sigma)
    return out_

# a protein as a ribbon diagram, a sheet of three strands between two helices
def ribbon(ax, x, y, s, color, angle=-30, lw=1.1):
    t = np.linspace(np.pi, 9 * np.pi, 220)
    def half(end):
        return np.r_[bezier((-0.80, 0.95), (-0.94, 0.90), (-0.99, 0.72), (-0.99, 0.5)),
                     np.c_[-0.82 + 0.17 * np.cos(t), 0.5 - (t - np.pi) / (8 * np.pi) + 0.11 * np.sin(t)],
                     bezier((-0.99, -0.5), (-0.99, -0.84), (-0.36, end - 0.34), (-0.36, end))]
    chains = [half(-0.55), -half(-0.62),
              bezier((-0.36, 0.62), (-0.36, 0.88), (0, 0.88), (0, 0.55)),
              bezier((0, -0.62), (0, -0.88), (0.36, -0.88), (0.36, -0.55))]
    arrow = np.array([[-0.08, -0.55], [-0.08, 0.34], [-0.16, 0.34], [0, 0.62], [0.16, 0.34], [0.08, 0.34], [0.08, -0.55]])
    a_ = np.deg2rad(angle)
    rot = np.array([[np.cos(a_), -np.sin(a_)], [np.sin(a_), np.cos(a_)]])
    place = lambda q: (q * s) @ rot.T + [x, y]
    for c in chains:
        q = place(c)
        ax.plot(q[:, 0], q[:, 1], color=color, lw=lw, solid_capstyle="round", solid_joinstyle="round")
    for k in range(3):
        q = place(arrow * [1, (-1) ** k] + [0.36 * (k - 1), 0])
        ax.add_patch(Polygon(q, closed=True, facecolor=(0.07, 0.16, 0.32), edgecolor=color, lw=lw, joinstyle="miter"))

def capsule(ax, cx, cy, angle, length=0.62, r=0.26, fill=(0.95, 0.97, 1.0), lw=1.8):
    t = np.linspace(-np.pi / 2, np.pi / 2, 40)
    right = np.c_[length / 2 + r * np.cos(t), r * np.sin(t)]
    body = np.r_[right, -right]
    half = np.r_[[[0, -r]], right, [[0, r]]]
    rot = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    place = lambda q: q @ rot.T + [cx, cy]
    ax.add_patch(Polygon(place(body), closed=True, facecolor=(0.10, 0.19, 0.36), edgecolor="none"))
    ax.add_patch(Polygon(place(half), closed=True, facecolor=fill, edgecolor="none"))
    ax.add_patch(Polygon(place(body), closed=True, facecolor="none", edgecolor=PALE, lw=lw, joinstyle="round"))
    m = place(np.array([[0, -r], [0, r]]))
    ax.plot(m[:, 0], m[:, 1], color=PALE, lw=lw * 0.8)
    hl = place(np.c_[np.linspace(-length / 2, length / 2, 20), np.full(20, 0.55 * r)])
    ax.plot(hl[:, 0], hl[:, 1], color=(1, 1, 1, 0.55), lw=lw * 1.2, solid_capstyle="round")


# ---------- layers ----------
# dependencies, the colour runs from the earlier to the later time point
# long ones are fainter, the ones of the selected protein go to their own layer
fig, ax = new_layer()
fig_s, ax_s = new_layer()
segs, cols, segs_s, cols_s = [], [], [], []
for i, j, t in zip(ei, ej, et):
    path, u, span = edge_path(t, i, j)
    c = TIME[t] + (TIME[t + 1] - TIME[t]) * u[:-1, None]
    seg = np.stack([path[:-1], path[1:]], axis=1)
    if (i == selected or j == selected) and i != j:
        segs_s.append(seg); cols_s.append(np.c_[c, np.full(len(c), 0.95)])
    else:
        segs.append(seg); cols.append(np.c_[c, np.full(len(c), 0.36 * (0.3 + 0.7 * np.exp(-span / 0.5)))])
ax.add_collection(LineCollection(np.concatenate(segs), colors=np.concatenate(cols), linewidths=0.5, capstyle="butt"))
if segs_s:
    ax_s.add_collection(LineCollection(np.concatenate(segs_s), colors=np.concatenate(cols_s), linewidths=0.9, capstyle="butt"))
L_edges, L_sel = grab(fig), grab(fig_s)

# proteins, the partners of the selected protein are larger
fig, ax = new_layer()
deg_ring = np.zeros((3, N))
partner = np.zeros((3, N), bool)
for i, j, t in zip(ei, ej, et):
    deg_ring[t, i] += 1; deg_ring[t + 1, j] += 1
    if i != j and i == selected: partner[t + 1, j] = True
    if i != j and j == selected: partner[t, i] = True
for k in range(3):
    x, y = polar(ang, RING[k])
    ax.scatter(x, y, s=1.0 + 0.55 * np.sqrt(deg_ring[k]) + 5.0 * partner[k], color=TIME[k], linewidths=0, zorder=3)
L_nodes = grab(fig)

# the selected protein at the three time points, joined by a thread
fig, ax = new_layer()
fig_m, ax_m = new_layer()
for k in range(3):
    x_, y_ = polar(SPOKE, RING[k])
    ax_m.add_patch(plt.Circle((x_, y_), 0.235, color="white"))
    ribbon(ax, x_, y_, 0.17, TIME[k])
L_prot = grab(fig)
M_prot = blur(grab(fig_m), 0.02)[:, :, :1]
fig, ax = new_layer()
fr = np.linspace(RING[0], RING[2], 80)
pts = np.c_[polar(np.full(80, SPOKE), fr)]
u = (fr - RING[0]) / (RING[2] - RING[0])
col = [TIME[0] + (TIME[1] - TIME[0]) * min(2 * v, 1) + (TIME[2] - TIME[1]) * max(2 * v - 1, 0) for v in u[:-1]]
ax.add_collection(LineCollection([pts[k:k + 2] for k in range(79)], colors=col, linewidths=1.6, capstyle="butt"))
L_thread = grab(fig) * (1 - M_prot)

# the cell
a = np.linspace(0, 2 * np.pi, 800)
fig, ax = new_layer()
x, y = CX + membrane(a) * np.cos(a), CY + membrane(a) * np.sin(a)
ax.plot(x, y, color=PALE, lw=2.0, solid_joinstyle="round")
ax.plot(CX + membrane(a, 0.955) * np.cos(a), CY + membrane(a, 0.955) * np.sin(a), color=(*PALE, 0.35), lw=0.7)
L_cell = grab(fig)
fig, ax = new_layer()
ax.add_patch(Polygon(np.c_[x, y], closed=True, facecolor=(0.10, 0.22, 0.42), edgecolor="none"))
L_fill = grab(fig)
fig, ax = new_layer()
na = np.linspace(0, 2 * np.pi, 300)
nr = 0.23 * R * (1 + 0.06 * np.cos(2 * na + 1.0) + 0.04 * np.cos(3 * na))
for q in np.linspace(1, 0.05, 24):
    c = np.array([0.30, 0.38, 0.72]) + (1 - q) * np.array([0.34, 0.34, 0.22])
    ax.add_patch(Polygon(np.c_[CX + 0.05 * (1 - q) + q * nr * np.cos(na), CY + 0.04 * (1 - q) + q * nr * np.sin(na)],
                         closed=True, facecolor=c, edgecolor="none"))
ax.plot(CX + nr * np.cos(na), CY + nr * np.sin(na), color=(0.80, 0.86, 1.0), lw=0.9)
ax.add_patch(plt.Circle((CX + 0.10, CY + 0.07), 0.085, color=(0.90, 0.93, 1.0)))
L_nuc = grab(fig)

# two drugs and a trail of particles towards the cell
fig, ax = new_layer()
capsule(ax, 0.72, 4.0, np.deg2rad(-38))
capsule(ax, 1.3, 4.58, np.deg2rad(-20), length=0.5, r=0.21, fill=(1.0, 0.72, 0.50))
n = 46
u = rng.random(n) ** 0.8
ax.scatter(1.0 + u * 1.35 + rng.normal(0, 0.10, n), 3.8 - u * 0.7 + rng.normal(0, 0.15, n),
           s=1.5 + 5 * (1 - u) + rng.random(n) * 2, color=(0.96, 0.92, 0.85), linewidths=0, alpha=0.9)
L_drug = grab(fig)

# neighbouring cells, cut off by the frame
fig, ax = new_layer()
for cx, cy, rr, ph in ((0.15, -0.25, 1.6, 4.1), (7.75, 5.35, 1.25, 0.2), (7.9, -0.1, 1.3, 2.0)):
    r_ = membrane(a, r=rr, ph=ph)
    ax.add_patch(Polygon(np.c_[cx + r_ * np.cos(a), cy + r_ * np.sin(a)], closed=True,
                         facecolor=(0.035, 0.075, 0.155), edgecolor=(*PALE, 0.40), lw=1.2))
    ax.add_patch(plt.Circle((cx + 0.1 * rr * np.cos(ph + 2.5), cy + 0.1 * rr * np.sin(ph + 2.5)), 0.24 * rr, color=(0.07, 0.12, 0.24)))
L_far = grab(fig)

# the application: a panel with the drug effects on the selected protein
# single drugs are on the diagonal and drug pairs off the diagonal, as in the heatmap of the application
# the grey scale is the one of the application, small p-values are white and pairs without data are left out
def app_grey(p):
    return np.interp(p, [0, 0.01, 0.05, 0.10, 0.25, 0.50, 1], np.array([255, 217, 150, 115, 64, 26, 0]) / 255)

pitch = min(0.108, 1.3 / K)
PX1 = 7.38
PX0 = PX1 - K * pitch - 0.24
PY0 = 2.42 - (K * pitch + 0.48) / 2
PY1 = PY0 + K * pitch + 0.48
panel = round_rect(PX0, PX1, PY0, PY1, 0.09)
fig, ax = new_layer()
ax.add_patch(Polygon(panel, closed=True, facecolor="white", edgecolor="none"))
M_panel = grab(fig)[:, :, :1]
fig, ax = new_layer()
ax.add_patch(Polygon(panel, closed=True, facecolor="none", edgecolor=(*PALE, 0.9), lw=1.0, joinstyle="round"))
# a line from the selected protein to the panel
x_, y_ = polar(SPOKE, RING[2])
ax.plot([x_ + 0.25, PX0], [y_, y_], color=(*PALE, 0.8), lw=0.8)
ax.add_patch(plt.Circle((PX0, y_), 0.028, color=PALE))
# header of the panel: two small drugs and the buttons of a window
capsule(ax, PX0 + 0.2, PY1 - 0.18, np.deg2rad(-38), length=0.085, r=0.04, lw=0.6)
capsule(ax, PX0 + 0.42, PY1 - 0.18, np.deg2rad(-20), length=0.085, r=0.04, fill=(1.0, 0.72, 0.50), lw=0.6)
for k in range(3):
    ax.add_patch(plt.Circle((PX1 - 0.14 - 0.1 * k, PY1 - 0.18), 0.026, color=(*PALE, 0.55)))
ax.plot([PX0 + 0.1, PX1 - 0.1], [PY1 - 0.33, PY1 - 0.33], color=(*PALE, 0.35), lw=0.6)
# the matrix, the first drug is at the bottom left
for i in range(K):
    for j in range(K):
        if np.isnan(p_mat[i, j]):
            continue
        cx_, cy_ = PX0 + 0.12 + (j + 0.5) * pitch, PY0 + 0.12 + (i + 0.5) * pitch
        g_ = app_grey(p_mat[i, j])
        if i == j:
            ax.add_patch(plt.Circle((cx_, cy_), 0.42 * pitch, facecolor=(g_, g_, g_), edgecolor=(*PALE, 0.45), lw=0.5))
        else:
            h_ = 0.4 * pitch
            ax.add_patch(Polygon(round_rect(cx_ - h_, cx_ + h_, cy_ - h_, cy_ + h_, 0.3 * h_, n=5), closed=True,
                                 facecolor=(g_, g_, g_), edgecolor=(*PALE, 0.28), lw=0.4))
L_panel = grab(fig)


# ---------- composition ----------
h, w = L_edges.shape[:2]
yy, xx = np.mgrid[0:h, 0:w].astype(np.float32)
dist = np.sqrt(((xx / w - 0.5) / 0.75) ** 2 + ((yy / h - 0.5) / 0.75) ** 2)
img = np.clip(1 - dist, 0, 1)[..., None] ** 1.4 * np.array([0.035, 0.075, 0.16], np.float32) + np.array([0.012, 0.02, 0.045], np.float32)
del yy, xx, dist
img += glow(L_far, 0, (0.8, 0.012), (0.35, 0.08))
img += 0.75 * L_fill
img += glow(L_nuc, 0, (0.9, 0.012), (0.5, 0.12))
L_edges *= 1 - M_prot
L_sel *= 1 - M_prot
L_nodes *= 1 - M_prot
img += glow(L_edges, 1.1, (0.9, 0.02), (0.6, 0.09))
img += glow(L_sel, 1.0, (1.0, 0.015), (0.8, 0.07))
img += glow(L_nodes, 1.0, (1.3, 0.012), (1.2, 0.05))
img += glow(L_thread, 1.0, (0.8, 0.015), (0.5, 0.06))
img += glow(L_prot, 1.0, (0.9, 0.012), (0.7, 0.06))
img += glow(L_cell, 1.0, (0.9, 0.02), (0.8, 0.10))
img += glow(L_drug, 1.0, (0.6, 0.03), (0.5, 0.12))
# the panel is a pane of dark glass over the scene
img = img * (1 - 0.78 * M_panel) + 0.78 * M_panel * np.array([0.045, 0.085, 0.17], np.float32)
img += glow(L_panel, 1.0, (0.6, 0.012), (0.35, 0.06))
img = 1 - np.exp(-1.25 * img)
img = np.clip(img, 0, 1) ** (1 / 1.05)

im = Image.fromarray((img * 255 + 0.5).astype(np.uint8), "RGB")
im.save(out + ".png", dpi=(dpi, dpi))
# LZW with the horizontal predictor keeps the file below the 10 MB limit of PLOS
im.save(out + ".tif", dpi=(dpi, dpi), compression="tiff_lzw", tiffinfo={317: 2})
print(f"{out}.png and {out}.tif written, {N} proteins and {len(ei)} dependencies drawn, "
      f"{len(segs_s)} of them of the selected protein")
