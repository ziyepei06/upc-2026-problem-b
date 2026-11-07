"""mkfig_a.py -- figures 1-4 (setup, trajectories, answer, histograms)."""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import bkick as B
import estimators as E
import figstyle as S
from mplhelp import hist2d_contour

PAL = S.PAL
R_OBS = np.array([[1.50, 0.45, 0.44], [2.96, 0.86, 0.74]])
T_OBS = (0.05, 0.10)


def chain():
    import os
    d = np.load(os.path.join(S.__file__.rsplit("/", 2)[0], "results", "chains",
                             "base_full.npz"))
    return d["chain"], d["t"], d["y"], d["z"]


# ---------------------------------------------------------------------------
def f01():
    fig, ax = plt.subplots(figsize=(S.TEXTW * 0.82, 2.62), layout="constrained")
    ax.plot([-0.3, 11], [0, 0], color=PAL["grey"], lw=0.8, zorder=1)
    ax.plot([11, 11], [0, B.GOAL_HEIGHT], color=PAL["ink"], lw=2.0, zorder=3)
    ax.plot([10.78, 11.22], [B.GOAL_HEIGHT, B.GOAL_HEIGHT], color=PAL["ink"], lw=1.4)
    ax.plot([0], [B.KICK_START[2]], marker="o", ms=5.5, color=PAL["ink"], zorder=5)
    ax.annotate("kick, $t=0$\n$(0,\\,0,\\,0.11)$", (0.04, 0.13), xytext=(0.55, 0.95),
                fontsize=8, arrowprops=dict(arrowstyle="-", lw=0.6,
                                            color=PAL["grey"], shrinkA=2, shrinkB=3))
    for i, (t, p, c) in enumerate(zip(T_OBS, R_OBS, (PAL["blue"], PAL["rust"]))):
        ax.plot([p[0]], [p[2]], marker="o", ms=5.8, mfc="white", mew=1.3,
                color=c, zorder=6)
        ax.text(p[0] + 0.22, p[2] - 0.34, f"$\\mathbf{{r}}_{i+1}$, $t={t:.2f}$ s",
                fontsize=8, color=c, zorder=7)
    v0, a = B.constant_accel_inverse(R_OBS, T_OBS)
    tt = np.linspace(0, 0.42, 200)
    ax.plot(v0[0] * tt + 0.5 * a[0] * tt**2,
            B.KICK_START[2] + v0[2] * tt + 0.5 * a[2] * tt**2,
            ls="-.", lw=1.0, color=PAL["green"], zorder=4,
            label="constant-$a$ extrapolation")
    ax.annotate("goal plane $x=11$ m", (11.0, 1.35), xytext=(7.4, 2.35),
                fontsize=8, color=PAL["ink"],
                arrowprops=dict(arrowstyle="-", lw=0.6, color=PAL["grey"],
                                shrinkA=2, shrinkB=3))
    ax.annotate("", xy=(11.45, 0), xytext=(11.45, 2.44),
                arrowprops=dict(arrowstyle="<->", lw=0.6, color=PAL["grey"]))
    ax.text(11.62, 1.22, "$h=2.44$ m", fontsize=8, rotation=90, va="center",
            ha="center")
    ax.set_xlim(-0.9, 12.9)
    ax.set_ylim(-0.35, 2.95)
    ax.set_xlabel("$x$  (m, towards the goal)")
    ax.set_ylabel("$z$  (m, height)")
    ax.legend(loc="upper left", fontsize=7.5)
    return S.save(fig, "f01-geometry")


# ---------------------------------------------------------------------------
def f02():
    ch, _, _, _ = chain()
    sub = ch[::40]
    times = list(np.arange(0.0, 0.505, 0.0025))
    kd, km = B.coeffs(np.exp(sub[:, 3]), 1.0, n=sub.shape[0])
    P = np.array(B.batch_positions(sub[:, 0:3], sub[:, 4:7], kd, km, times,
                                   dt=0.002, tau=6.0))
    q = np.percentile(P, [5, 16, 50, 84, 95], axis=0)
    v0, a = B.constant_accel_inverse(R_OBS, T_OBS)
    tt = np.linspace(0, 0.406, 400)

    fig, axes = plt.subplots(1, 2, figsize=(S.TEXTW, 2.95), layout="constrained")
    ax = axes[0]
    ax.fill_between(q[2, :, 0], q[0, :, 1], q[4, :, 1], color=PAL["blue"],
                    alpha=0.13, lw=0)
    ax.fill_between(q[2, :, 0], q[1, :, 1], q[3, :, 1], color=PAL["blue"],
                    alpha=0.33, lw=0)
    ax.plot(q[2, :, 0], q[2, :, 1], color=PAL["blue"], lw=1.0)
    ax.plot(v0[0] * tt + 0.5 * a[0] * tt**2, v0[1] * tt + 0.5 * a[1] * tt**2,
            ls="--", lw=1.2, color=PAL["rust"])
    for p in R_OBS:
        ax.errorbar([p[0]], [p[1]], xerr=0.01, yerr=0.01, fmt="o", ms=4.6,
                    color=PAL["ink"], mfc="none", mew=1.0, zorder=8)
    ax.axvline(11, color=PAL["ink"], lw=1.2)
    ax.set_xlim(-0.4, 11.5)
    ax.set_ylim(-0.2, 3.9)
    ax.set_xlabel("$x$  (m)")
    ax.set_ylabel("$y$  (m)")
    ax.set_title("horizontal plane")
    ax.plot([], [], color=PAL["blue"], alpha=0.33, lw=6, label="physical, 68 %")
    ax.plot([], [], color=PAL["blue"], alpha=0.13, lw=6, label="physical, 90 %")
    ax.plot([], [], color=PAL["rust"], ls="--", lw=1.2, label="constant $a$")
    ax.legend(loc="lower right", fontsize=7, handlelength=1.6)

    ax = axes[1]
    ax.add_patch(Rectangle((10.88, 0), 0.36, B.GOAL_HEIGHT, facecolor="#e3eaf4",
                           edgecolor=PAL["ink"], lw=0.8, zorder=1))
    ax.fill_between(q[2, :, 0], q[0, :, 2], q[4, :, 2], color=PAL["blue"],
                    alpha=0.13, lw=0)
    ax.fill_between(q[2, :, 0], q[1, :, 2], q[3, :, 2], color=PAL["blue"],
                    alpha=0.33, lw=0)
    ax.plot(q[2, :, 0], q[2, :, 2], color=PAL["blue"], lw=1.0)
    ax.plot(v0[0] * tt + 0.5 * a[0] * tt**2,
            B.KICK_START[2] + v0[2] * tt + 0.5 * a[2] * tt**2,
            ls="--", lw=1.2, color=PAL["rust"])
    for p in R_OBS:
        ax.errorbar([p[0]], [p[2]], xerr=0.01, yerr=0.01, fmt="o", ms=4.6,
                    color=PAL["ink"], mfc="none", mew=1.0, zorder=8)
    ax.text(11.35, 1.22, "goal\nmouth", fontsize=7, va="center", color=PAL["ink"])
    ax.set_xlim(-0.4, 11.9)
    ax.set_ylim(-0.1, 3.1)
    ax.set_xlabel("$x$  (m)")
    ax.set_ylabel("$z$  (m)")
    ax.set_title("vertical plane")


    S.panel_label(axes[0], "a", dx=0.03, dy=0.96)
    S.panel_label(ax, "b", dx=0.97, dy=0.96, loc="tr")
    return S.save(fig, "f02-trajectories")


# ---------------------------------------------------------------------------
def f03():
    ch, _, y2, z2 = chain()
    d1 = E.e1_predict(R_OBS, T_OBS, sigma=0.01, n_mc=40000)
    y1, z1 = d1["_y"], d1["_z"]
    fig, ax = plt.subplots(figsize=(5.6, 2.75), layout="constrained")
    ax.add_patch(Rectangle((-B.GOAL_HALF_WIDTH, 0), 2 * B.GOAL_HALF_WIDTH,
                           B.GOAL_HEIGHT, facecolor="#f4f6fa", edgecolor="none",
                           zorder=0))
    ax.add_patch(Rectangle((-B.GOAL_HALF_WIDTH, 0), 2 * B.GOAL_HALF_WIDTH,
                           B.GOAL_HEIGHT, facecolor="none", edgecolor=PAL["ink"],
                           lw=1.6, zorder=6))
    ax.axvline(0, color=PAL["light"], lw=0.6, zorder=0)
    rng = np.random.default_rng(1)
    s1 = rng.choice(y1.size, 5000, replace=False)
    s2 = rng.choice(y2.size, 5000, replace=False)
    ax.scatter(y1[s1], z1[s1], s=1.5, color=PAL["rust"], alpha=0.09, lw=0, zorder=1)
    ax.scatter(y2[s2], z2[s2], s=1.5, color=PAL["blue"], alpha=0.09, lw=0, zorder=2)
    for yy, zz, c in ((y1, z1, PAL["rust"]), (y2, z2, PAL["blue"])):
        hist2d_contour(ax, yy, zz, color=c, levels=(0.393, 0.865, 0.95), nbin=80)
    for lbl, yy, zz, c, tx in (
            ("constant $a$", 2.495, 1.920, PAL["rust"], (0.42, 2.86)),
            ("model average", 2.767, 1.933, PAL["green"], (0.42, 0.30)),
            ("physical prior", 3.012, 1.946, PAL["blue"], (3.75, 0.42))):
        ax.plot([yy], [zz], marker="+", ms=10, mew=1.5, color=c, zorder=9)
        ax.annotate(lbl, (yy, zz), xytext=tx, fontsize=8, color=c,
                    arrowprops=dict(arrowstyle="-", lw=0.7, color=c,
                                    shrinkA=2, shrinkB=3))
    ax.set_xlim(-4.3, 4.3)
    ax.set_ylim(-0.08, 3.10)
    ax.set_aspect("equal")
    ax.set_xlabel("$y$  (m)")
    ax.set_ylabel("$z$  (m)")
    ax.text(-4.05, 2.92, "68 / 90 / 95 % of the predictive mass", fontsize=7.5,
            color=PAL["grey"])
    return S.save(fig, "f03-goalplane")


# ---------------------------------------------------------------------------
def f04():
    d1 = E.e1_predict(R_OBS, T_OBS, sigma=0.01, n_mc=40000)
    y1, z1 = d1["_y"], d1["_z"]
    ch, t2, y2, z2 = chain()
    m = 20000
    rng = np.random.default_rng(7)
    pick1 = rng.random(m) < 0.5
    ym = np.where(pick1, rng.permutation(y1)[:m], rng.permutation(y2)[:m])
    zm = np.where(pick1, rng.permutation(z1)[:m], rng.permutation(z2)[:m])

    fig, axes = plt.subplots(1, 2, figsize=(S.TEXTW, 2.45), layout="constrained")
    ax = axes[0]
    bins = np.linspace(-0.2, 4.8, 60)
    for v, c, l, ls in ((y1, PAL["rust"], "constant $a$", "-"),
                        (y2, PAL["blue"], "physical prior", "-"),
                        (ym, PAL["green"], "model average", "--")):
        ax.hist(v, bins=bins, density=True, histtype="step", lw=1.3 if ls == "--"
                else 0.9, color=c, ls=ls, label=l)
    ax.axvline(B.GOAL_HALF_WIDTH, color=PAL["ink"], lw=1.0, ls=":")
    ax.text(B.GOAL_HALF_WIDTH + 0.06, 0.86, "post", fontsize=7.5, va="top")
    ax.set_xlabel("$y$ at the goal plane  (m)")
    ax.set_ylabel("density  (m$^{-1}$)")
    ax.legend(loc="upper left", fontsize=7)
    ax.set_title("lateral position")
    ax = axes[1]
    bins = np.linspace(0.4, 3.5, 60)
    for v, c, l, ls in ((z1, PAL["rust"], "constant $a$", "-"),
                        (z2, PAL["blue"], "physical prior", "-"),
                        (zm, PAL["green"], "model average", "--")):
        ax.hist(v, bins=bins, density=True, histtype="step", lw=1.3 if ls == "--"
                else 0.9, color=c, ls=ls, label=l)
    ax.axvline(B.GOAL_HEIGHT, color=PAL["ink"], lw=1.0, ls=":")
    ax.text(B.GOAL_HEIGHT - 0.08, 0.55, "crossbar", fontsize=7.5, ha="right")
    ax.set_xlabel("$z$ at the goal plane  (m)")
    ax.set_ylabel("density  (m$^{-1}$)")
    ax.set_title("height")
    S.panel_label(axes[0], "a", dx=0.03, dy=0.96)
    S.panel_label(ax, "b", dx=0.97, dy=0.96, loc="tr")
    return S.save(fig, "f04-histograms")


if __name__ == "__main__":
    import sys
    for n in (sys.argv[1:] or ["f01", "f02", "f03", "f04"]):
        globals()[n]()
