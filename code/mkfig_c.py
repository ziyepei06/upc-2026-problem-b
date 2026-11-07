"""Figures 3 and 4: the answer in the goal plane, and its marginals."""
import json, os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
import figstyle as S
import bkick as B
import estimators as E
from mplhelp import hist2d_contour

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
ANS = json.load(open(os.path.join(RES, "answer.json")))
CH = np.load(os.path.join(RES, "chains", "base_full.npz"))

R_OBS = B.SCENARIO["r_obs"] if "r_obs" in B.SCENARIO else np.array(
    [[1.50, 0.45, 0.44], [2.96, 0.86, 0.74]])
T_OBS = np.array([0.05, 0.10])


def _samples():
    e1 = E.e1_predict(R_OBS, T_OBS, 0.010, n_mc=40000, seed=3)
    y1, z1 = e1["_y"], e1["_z"]
    y2, z2 = CH["y"], CH["z"]
    m = min(len(y1), len(y2))
    return (y1, z1), (y2[:m], z2[:m]), e1


def f03():
    (y1, z1), (y2, z2), e1 = _samples()
    fig, ax = plt.subplots(figsize=(S.TEXTW * 0.80, 2.85),
                           layout="constrained")
    b, h, r = 3.66, 2.44, 0.11
    ax.add_patch(Rectangle((-b, 0), 2 * b, h, facecolor="#f2f5fa",
                           edgecolor=PAL["ink"], lw=1.6, zorder=1))
    ax.add_patch(Rectangle((-b + r, r), 2 * (b - r), h - 2 * r,
                           facecolor="none", edgecolor=PAL["grey"],
                           lw=0.7, ls=(0, (4, 3)), zorder=2))
    hist2d_contour(ax, y1, z1, color=PAL["rust"], ylim=(-5, 5), zlim=(0, 3.6),
                   nbin=90, smooth=1.3)
    hist2d_contour(ax, y2, z2, color=PAL["blue"], ylim=(-5, 5), zlim=(0, 3.6),
                   nbin=90, smooth=1.3)
    pts = [("o", 2.4948, 1.9204, PAL["rust"], "constant $a$", "2.50, 1.92"),
           ("s", 2.764, 1.937, PAL["green"], "model average", "2.77, 1.93"),
           ("D", 3.029, 1.950, PAL["blue"], "physical (E2)", "3.03, 1.95")]
    for mk, y, z, c, lab, num in pts:
        ax.plot([y], [z], marker=mk, ms=4.6, mfc="white", mew=1.2, color=c,
                zorder=6)
    handles = [Line2D([], [], color=PAL["rust"], lw=1.2,
                      label="E1, 68 % contour"),
               Line2D([], [], color=PAL["blue"], lw=1.2,
                      label="E2, 68 % contour"),
               Line2D([], [], marker="o", ls="none", mfc="white", mew=1.2,
                      color=PAL["rust"],
                      label="constant $a$: (2.50, 1.92)"),
               Line2D([], [], marker="s", ls="none", mfc="white", mew=1.2,
                      color=PAL["green"],
                      label="model average: (2.76, 1.94)"),
               Line2D([], [], marker="D", ls="none", mfc="white", mew=1.2,
                      color=PAL["blue"],
                      label="physical (E2): (3.03, 1.95)")]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(0.005, 0.50),
              fontsize=6.9, handlelength=1.7, borderpad=0.3, labelspacing=0.45)
    ax.text(0.15, 3.02, "68, 87 and 95 % of the predictive mass",
            fontsize=6.6, color=PAL["grey"], ha="center")
    ax.set_xlim(-3.95, 3.95)
    ax.set_ylim(-0.30, 3.25)
    ax.set_aspect("equal")
    ax.set_xlabel("$y$ at the goal plane  (m)")
    ax.set_ylabel("$z$  (m)")
    ax.set_xticks([-3, -2, -1, 0, 1, 2, 3])
    ax.set_yticks([0, 1, 2])
    S.save(fig, "f03-goalplane")


def f04():
    (y1, z1), (y2, z2), e1 = _samples()
    rng = np.random.default_rng(11)
    m = min(len(y1), len(y2))
    pick = rng.random(m) < 0.5
    ym = np.where(pick, y1[:m], y2[:m])
    zm = np.where(pick, z1[:m], z2[:m])
    fig, axes = plt.subplots(1, 2, figsize=(S.TEXTW * 0.95, 2.5),
                             layout="constrained")
    for ax, dat, lab, cut, cutlab, xl, xr in (
            (axes[0], (y1, y2, ym), "$y$ at the goal plane  (m)", 3.66,
             "post", 0.2, 4.9),
            (axes[1], (z1, z2, zm), "$z$ at the goal plane  (m)", 2.44,
             "crossbar", 0.5, 3.1)):
        bins = np.linspace(0, 5 if ax is axes[0] else 3.4, 90)
        ax.hist(dat[0], bins=bins, density=True, histtype="step", lw=1.1,
                color=PAL["rust"], label="constant-$a$ (E1)")
        ax.hist(dat[1], bins=bins, density=True, histtype="step", lw=1.1,
                color=PAL["blue"], label="physical (E2)")
        ax.hist(dat[2], bins=bins, density=True, histtype="step", lw=1.0,
                ls=(0, (5, 2)), color=PAL["green"], label="average")
        ax.axvline(cut, color=PAL["ink"], lw=0.9, ls=":")
        ax.set_xlabel(lab)
        ax.set_xlim(xl, xr)
        ax.set_yticks([])
        ax.set_ylabel("predictive density")
    axes[0].set_ylim(0, 1.95)
    axes[1].set_ylim(0, 2.9)
    axes[0].annotate("post", (3.66, 0.95), xytext=(4.15, 1.55), fontsize=6.8,
                     arrowprops=dict(arrowstyle="-", lw=0.5,
                                     color=PAL["ink"]))
    axes[1].annotate("crossbar", (2.44, 1.5), xytext=(2.55, 2.35), fontsize=6.8,
                     arrowprops=dict(arrowstyle="-", lw=0.5,
                                     color=PAL["ink"]))
    axes[0].legend(fontsize=6.8, loc="upper left")
    S.panel_label(axes[0], "a")
    S.panel_label(axes[1], "b")
    S.save(fig, "f04-histograms")


PAL = S.PAL

if __name__ == "__main__":
    import sys
    for n in (sys.argv[1:] or ["f03", "f04"]):
        globals()[n]()
        print("ok", n)
