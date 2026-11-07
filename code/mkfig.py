"""mkfig.py -- every figure in the paper.  Usage: python3 mkfig.py f01 f02 ..."""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch, Polygon
from matplotlib.lines import Line2D

import bkick as B
import estimators as E
import figstyle as S

S.setup()
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
PAL = S.PAL
R_OBS = np.array([[1.50, 0.45, 0.44], [2.96, 0.86, 0.74]])
T_OBS = (0.05, 0.10)


def load(name):
    with open(os.path.join(RES, name)) as f:
        return json.load(f)


def chain():
    d = np.load(os.path.join(RES, "chains", "base_full.npz"))
    return d["chain"], d["t"], d["y"], d["z"]


def hist2d_contour(ax, y, z, levels=(0.28, 0.60, 0.90), color=None, nbin=60,
                   ylim=(-4.2, 4.2), zlim=(-0.2, 3.2), smooth=1.2):
    """Smoothed 2-D density contour enclosing the requested probability mass."""
    H, ye, ze = np.histogram2d(y, z, bins=(nbin, nbin),
                               range=((ylim[0], ylim[1]), (zlim[0], zlim[1])))
    # Gaussian smoothing by separable convolution
    k = int(max(3, round(smooth * 3))) | 1
    g = np.exp(-0.5 * (np.arange(k) - k // 2) ** 2 / smooth**2)
    g /= g.sum()
    Hp = np.apply_along_axis(lambda m: np.convolve(m, g, "same"), 0, H)
    Hp = np.apply_along_axis(lambda m: np.convolve(m, g, "same"), 1, Hp)
    tot = Hp.sum()
    if tot <= 0:
        return None
    # levels from cumulative mass, from the outside in
    flat = np.sort(Hp.ravel())[::-1]
    cum = np.cumsum(flat) / tot
    lv = []
    for l in levels:
        idx = np.searchsorted(cum, l)
        lv.append(flat[min(idx, flat.size - 1)])
    lv = sorted(set(lv))
    yc = 0.5 * (ye[:-1] + ye[1:])
    zc = 0.5 * (ze[:-1] + ze[1:])
    cs = ax.contour(yc, zc, Hp.T, levels=lv, colors=[color] * len(lv),
                    linewidths=(0.9, 0.7, 0.5)[:len(lv)])
    return cs


# ---------------------------------------------------------------------------
# f01 geometry of the problem
# ---------------------------------------------------------------------------
def f01():
    fig, axes = plt.subplots(1, 2, figsize=(S.TEXTW, 2.55),
                             gridspec_kw=dict(width_ratios=[1.55, 1.0], wspace=0.30))
    ax = axes[0]
    # side view x-z
    ax.add_patch(Rectangle((11, 0), 0.25, B.GOAL_HEIGHT, facecolor="#e8eef7",
                           edgecolor="none", zorder=0))
    ax.plot([11, 11], [0, B.GOAL_HEIGHT], color=PAL["ink"], lw=1.4, zorder=3)
    ax.plot([11, 11.25], [B.GOAL_HEIGHT, B.GOAL_HEIGHT], color=PAL["ink"], lw=1.4)
    ax.plot([0, 11], [0, 0], color=PAL["grey"], lw=0.8)
    ax.plot([0], [B.KICK_START[2]], marker="o", ms=4.5, color=PAL["ink"], zorder=5)
    ax.annotate("kick, $t=0$\n$(0,0,0.11)$", (0, 0.11), xytext=(0.5, 0.62),
                fontsize=7.5, arrowprops=dict(arrowstyle="-", lw=0.6,
                                              color=PAL["grey"]))
    for (t, p, c) in zip(T_OBS, R_OBS, (PAL["blue"], PAL["rust"])):
        ax.plot([p[0]], [p[2]], marker="o", ms=4.5, mfc="white", mew=1.0,
                color=c, zorder=6)
        ax.annotate(f"$t={t:.2f}$ s\n({p[0]:.2f}, {p[2]:.2f})", (p[0], p[2]),
                    xytext=(p[0] - 1.2, p[2] + 0.55), fontsize=7.5, color=c,
                    arrowprops=dict(arrowstyle="-", lw=0.6, color=c))
    v0, a = B.constant_accel_inverse(R_OBS, T_OBS)
    tt = np.linspace(0, 0.42, 200)
    xs = v0[0] * tt + 0.5 * a[0] * tt**2
    zs = B.KICK_START[2] + v0[2] * tt + 0.5 * a[2] * tt**2
    ax.plot(xs, zs, ls="-.", lw=1.0, color=PAL["green"], zorder=4)
    ax.set_xlim(-0.7, 12.6)
    ax.set_ylim(-0.15, 3.1)
    ax.set_xlabel("$x$  (m, towards the goal)")
    ax.set_ylabel("$z$  (m, height)")
    ax.set_title("side view")
    ax.annotate("goal plane\n$x=11$ m", (11, 1.2), xytext=(11.7, 1.55),
                fontsize=7.5, color=PAL["grey"],
                arrowprops=dict(arrowstyle="-", lw=0.6, color=PAL["grey"]))

    ax = axes[1]
    # goal face y-z
    ax.add_patch(Rectangle((-B.GOAL_HALF_WIDTH, 0), 2 * B.GOAL_HALF_WIDTH,
                           B.GOAL_HEIGHT, facecolor="#e8eef7", edgecolor=PAL["ink"],
                           lw=1.4, zorder=1))
    ax.axvline(0, color=PAL["light"], lw=0.6, zorder=0)
    for yy in np.arange(-3.5, 3.6, 0.5):
        ax.plot([yy, yy], [0, 0.05], color=PAL["grey"], lw=0.5)
    ax.set_xlim(-4.3, 4.3)
    ax.set_ylim(-0.12, 3.2)
    ax.set_aspect("equal")
    ax.set_xlabel("$y$  (m)")
    ax.set_ylabel("$z$  (m)")
    ax.set_title("goal face")
    ax.text(0, 1.22, "3.66 m", ha="center", fontsize=7.5, color=PAL["grey"])
    ax.annotate("", xy=(-3.66, 2.62), xytext=(3.66, 2.62),
                arrowprops=dict(arrowstyle="<->", lw=0.6, color=PAL["grey"]))
    ax.annotate("", xy=(4.05, 0), xytext=(4.05, 2.44),
                arrowprops=dict(arrowstyle="<->", lw=0.6, color=PAL["grey"]))
    ax.text(4.15, 1.2, "2.44 m", fontsize=7.5, rotation=90, va="center",
            color=PAL["grey"])
    ax.plot([0.12], [0.12], marker="o", ms=3.4, color=PAL["grey"])
    ax.text(0.3, 0.10, "$r=0.11$ m", fontsize=7, color=PAL["grey"])
    S.panel_label(axes[0], "a", dx=0.03, dy=0.96)
    S.panel_label(ax, "b", dx=0.97, dy=0.96, loc="tr")
    return S.save(fig, "f01-geometry")


# ---------------------------------------------------------------------------
# f02 reconstructed trajectories
# ---------------------------------------------------------------------------
def f02():
    ch, _, _, _ = chain()
    sub = ch[::800]
    times = list(np.arange(0.0, 0.46, 0.005))
    kd, km = B.coeffs(np.exp(sub[:, 3]), 1.0, n=sub.shape[0])
    P = np.array(B.batch_positions(sub[:, 0:3], sub[:, 4:7], kd, km, times,
                                   dt=0.002, tau=6.0)).transpose(1, 0, 2)
    q = np.percentile(P, [16, 84, 5, 95], axis=0)

    fig, axes = plt.subplots(1, 2, figsize=(S.TEXTW, 2.8))
    ax = axes[0]
    ax.fill_between(q[3, :, 0], q[3, :, 1], q[1, :, 1], color=PAL["blue"],
                    alpha=0.10, lw=0)
    ax.fill_between(q[0, :, 0], q[0, :, 1], q[1, :, 1], color=PAL["blue"],
                    alpha=0.30, lw=0)
    ax.plot(q[1, :, 0], q[1, :, 1], color=PAL["blue"], lw=0.8)
    ax.plot(q[0, :, 0], q[0, :, 1], color=PAL["blue"], lw=0.8)
    for s in P[::4]:
        ax.plot(s[:, 0], s[:, 1], color=PAL["blue"], lw=0.25, alpha=0.12)
    v0, a = B.constant_accel_inverse(R_OBS, T_OBS)
    tt = np.linspace(0, 0.406, 300)
    ax.plot(v0[0] * tt + 0.5 * a[0] * tt**2, v0[1] * tt + 0.5 * a[1] * tt**2,
            ls="--", lw=1.1, color=PAL["rust"])
    for i, (p, c) in enumerate(zip(R_OBS, (PAL["blue"], PAL["rust"]))):
        ax.errorbar([p[0]], [p[1]], yerr=0.01, xerr=0.01, fmt="o", ms=4.2,
                    color=PAL["ink"], mfc="none", mew=0.9, zorder=8)
    ax.axvline(11, color=PAL["ink"], lw=1.0)
    ax.set_xlim(-0.4, 11.6)
    ax.set_ylim(-0.2, 4.0)
    ax.set_xlabel("$x$  (m)")
    ax.set_ylabel("$y$  (m)")
    ax.set_title("horizontal plane")
    ax.text(0.2, 3.55, "dashed: constant-$a$ extrapolation", fontsize=7,
            color=PAL["rust"])
    ax.text(0.2, 3.25, "shaded: physical posterior", fontsize=7, color=PAL["blue"])

    ax = axes[1]
    ax.fill_between(q[3, :, 0], q[3, :, 2], q[1, :, 2], color=PAL["blue"],
                    alpha=0.10, lw=0)
    ax.fill_between(q[0, :, 0], q[0, :, 2], q[1, :, 2], color=PAL["blue"],
                    alpha=0.30, lw=0)
    ax.plot(q[1, :, 0], q[1, :, 2], color=PAL["blue"], lw=0.8)
    ax.plot(q[0, :, 0], q[0, :, 2], color=PAL["blue"], lw=0.8)
    for s in P[::4]:
        ax.plot(s[:, 0], s[:, 2], color=PAL["blue"], lw=0.25, alpha=0.12)
    zs = B.KICK_START[2] + v0[2] * tt + 0.5 * a[2] * tt**2
    ax.plot(v0[0] * tt + 0.5 * a[0] * tt**2, zs, ls="--", lw=1.1, color=PAL["rust"])
    for p in R_OBS:
        ax.errorbar([p[0]], [p[2]], yerr=0.01, xerr=0.01, fmt="o", ms=4.2,
                    color=PAL["ink"], mfc="none", mew=0.9, zorder=8)
    ax.add_patch(Rectangle((11, 0), 0.18, B.GOAL_HEIGHT, facecolor="#e8eef7",
                           edgecolor=PAL["ink"], lw=1.0, zorder=1))
    ax.set_xlim(-0.4, 12.6)
    ax.set_ylim(-0.1, 3.2)
    ax.set_xlabel("$x$  (m)")
    ax.set_ylabel("$z$  (m)")
    ax.set_title("vertical plane")
    S.panel_label(axes[0], "a", dx=0.03, dy=0.95)
    S.panel_label(ax, "b", dx=0.97, dy=0.95, loc="tr")
    return S.save(fig, "f02-trajectories")


# ---------------------------------------------------------------------------
# f03 the answer: predictive distribution in the goal plane
# ---------------------------------------------------------------------------
def f03():
    ch, t2, y2, z2 = chain()
    d1 = E.e1_predict(R_OBS, T_OBS, sigma=0.01, n_mc=30000)
    y1, z1 = d1["_y"], d1["_z"]
    fig, axes = plt.subplots(1, 2, figsize=(S.TEXTW, 3.15),
                             gridspec_kw=dict(width_ratios=[1.15, 1.0], wspace=0.26))
    ax = axes[0]
    ax.add_patch(Rectangle((-B.GOAL_HALF_WIDTH, 0), 2 * B.GOAL_HALF_WIDTH,
                           B.GOAL_HEIGHT, facecolor="#eef2f8", edgecolor="none",
                           zorder=0))
    ax.add_patch(Rectangle((-B.GOAL_HALF_WIDTH, 0), 2 * B.GOAL_HALF_WIDTH,
                           B.GOAL_HEIGHT, facecolor="none", edgecolor=PAL["ink"],
                           lw=1.4, zorder=6))
    rng = np.random.default_rng(1)
    s1 = rng.choice(y1.size, 4000, replace=False)
    s2 = rng.choice(y2.size, 4000, replace=False)
    ax.scatter(y1[s1], z1[s1], s=1.6, color=PAL["rust"], alpha=0.10, lw=0, zorder=1)
    ax.scatter(y2[s2], z2[s2], s=1.6, color=PAL["blue"], alpha=0.10, lw=0, zorder=2)
    hist2d_contour(ax, y1, z1, color=PAL["rust"], levels=(0.39, 0.86, 0.95))
    hist2d_contour(ax, y2, z2, color=PAL["blue"], levels=(0.39, 0.86, 0.95))
    pts = [("constant $a$", 2.495, 1.920, PAL["rust"], (1.7, 2.62)),
           ("physical", 3.012, 1.946, PAL["blue"], (3.35, 1.35))]
    for lbl, yy, zz, c, tx in pts:
        ax.plot([yy], [zz], marker="+", ms=9, mew=1.3, color=c, zorder=9)
        ax.annotate(lbl, (yy, zz), xytext=tx, fontsize=7.5, color=c,
                    arrowprops=dict(arrowstyle="-", lw=0.6, color=c))
    ax.set_xlim(-4.2, 4.2)
    ax.set_ylim(-0.05, 3.05)
    ax.set_aspect("equal")
    ax.set_xlabel("$y$  (m)")
    ax.set_ylabel("$z$  (m)")
    ax.set_title("predictive distributions")

    ax = axes[1]
    edges = np.linspace(0, 4.6, 61)
    cy, _ = np.histogram(np.abs(y2), bins=edges, density=True)
    cz, _ = np.histogram(z2, bins=np.linspace(0.4, 3.4, 61), density=True)
    cy1, _ = np.histogram(np.abs(y1), bins=edges, density=True)
    ax.step(edges[:-1], cy1, where="post", color=PAL["rust"], lw=1.0,
            label="constant $a$")
    ax.step(edges[:-1], cy, where="post", color=PAL["blue"], lw=1.0,
            label="physical prior")
    ax.axvline(B.GOAL_HALF_WIDTH, color=PAL["ink"], lw=1.0, ls=":")
    ax.annotate("post", (B.GOAL_HALF_WIDTH, 0.62), xytext=(2.6, 0.78),
                fontsize=7.5, arrowprops=dict(arrowstyle="->", lw=0.6))
    ax.set_xlabel("$|y|$ at the goal plane  (m)")
    ax.set_ylabel("probability density  (m$^{-1}$)")
    ax.set_xlim(0, 4.6)
    ax.legend(loc="upper right")
    ax.set_title("distance from the goal axis")
    S.panel_label(axes[0], "a", dx=0.03, dy=0.96)
    S.panel_label(ax, "b", dx=0.97, dy=0.96, loc="tr")
    return S.save(fig, "f03-goalplane")


# ---------------------------------------------------------------------------
# f04 predictive histograms of y and z with the interval bands
# ---------------------------------------------------------------------------
def f04():
    ans = load("answer.json")
    fig, axes = plt.subplots(1, 2, figsize=(S.TEXTW, 2.5))
    for ax, key, lab in ((axes[0], "y", "$y$ at the goal plane  (m)"),
                         (axes[1], "z", "$z$ at the goal plane  (m)")):
        for name, c, lbl in (("e1", PAL["rust"], "constant $a$"),
                             ("e2", PAL["blue"], "physical"),
                             ("mix", PAL["green"], "model average")):
            s = ans[name]["summary"]
            mu, sd = s["mean"][key], s["std"][key]
            x = np.linspace(mu - 3.6 * sd, mu + 3.6 * sd, 400)
            pdf = np.exp(-0.5 * ((x - mu) / sd) ** 2) / (np.sqrt(2 * np.pi) * sd)
            ax.plot(x, pdf, color=c, lw=1.2 if name != "mix" else 1.5,
                    ls="-" if name != "mix" else "--", label=lbl)
            ax.axvline(mu, color=c, lw=0.6, ls=":")
        ax.set_xlabel(lab)
        ax.set_ylabel("density  (m$^{-1}$)")
        ax.set_yticks([])
    axes[0].legend(loc="upper left")
    axes[0].axvline(B.GOAL_HALF_WIDTH, color=PAL["ink"], lw=1.0)
    axes[0].text(B.GOAL_HALF_WIDTH - 0.15, axes[0].get_ylim()[1] * 0.55, "post",
                 rotation=90, ha="right", va="center", fontsize=7.5)
    S.panel_label(axes[0], "a", dx=0.03, dy=0.96)
    S.panel_label(axes[1], "b", dx=0.97, dy=0.96, loc="tr")
    return S.save(fig, "f04-histograms")




# ---------------------------------------------------------------------------
# f05 identifiability: singular spectrum and the compensation ridge
# ---------------------------------------------------------------------------
def f05():
    import an_budget as AB
    b1 = AB.e1_parts(0.01)
    b2 = AB.e2_parts(0.01, 0.19, 0.15, 35.0)
    J1, e1e = b1[0], None
    # rebuild J1 with the same function, simpler: recompute here
    sig = 0.01
    J1 = np.zeros((6, 6))
    for i, t in enumerate(np.array([0.05, 0.10])):
        J1[3 * i:3 * i + 3, :3] = np.eye(3) * t
        J1[3 * i:3 * i + 3, 3:] = np.eye(3) * 0.5 * t * t
    scales1 = np.array([1.0] * 3 + [25.0] * 3)      # 1 m/s, 25 m/s^2-ish
    A1 = (J1 * scales1) .T @ (J1 * scales1) / sig**2
    s1 = np.sqrt(np.maximum(np.sort(np.linalg.eigvalsh(A1))[::-1], 1e-30))

    th = np.array(b2["theta"])
    P = E.PhysicalPosterior(R_OBS, T_OBS, sigma=0.01)
    J2 = np.zeros((6, 7))
    for k in range(7):
        h = 1e-6 * max(1.0, abs(th[k]))
        tp = th.copy(); tp[k] += h
        tm = th.copy(); tm[k] -= h
        J2[:, k] = (np.asarray(P.positions(tp[None, :])).reshape(-1)
                    - np.asarray(P.positions(tm[None, :])).reshape(-1)) / (2 * h)
    scales2 = np.array([1.0, 1.0, 1.0, 0.15, 35.0, 35.0, 35.0])
    A2 = (J2 * scales2).T @ (J2 * scales2) / sig**2
    s2 = np.sqrt(np.maximum(np.sort(np.linalg.eigvalsh(A2))[::-1], 1e-30))

    fig, axes = plt.subplots(1, 2, figsize=(S.TEXTW, 2.75),
                             gridspec_kw=dict(wspace=0.32, width_ratios=[1.0, 1.12]))
    ax = axes[0]
    ax.semilogy(np.arange(1, s2.size + 1), s2, "o-", ms=3.6, lw=0.9,
                color=PAL["blue"], label=r"7-parameter physical model")
    ax.semilogy(np.arange(1, s1.size + 1), s1, "s--", ms=3.4, lw=0.9,
                color=PAL["rust"], label=r"6-parameter constant-$a$ model")
    ax.set_xlabel("direction index")
    ax.set_ylabel(r"$\sigma_i/\sigma_{\rm data}$  (prior-sd units)")
    ax.set_title("scaled spectrum of the data Jacobian")
    ax.legend(loc="lower left")
    ax.annotate("", xy=(6.6, s2[-1]), xytext=(6.6, s2[0]),
                arrowprops=dict(arrowstyle="<->", lw=0.6, color=PAL["grey"]))
    ax.text(6.75, np.sqrt(s2[0] * s2[-1]), "10$^{16}$", rotation=90, fontsize=7.5,
            va="center", color=PAL["grey"])

    ax = axes[1]
    th0 = th.copy()
    ny, nz = 81, 81
    v0y = th0[1] + np.linspace(-2.2, 2.2, ny)
    wz = th0[6] + np.linspace(-110, 110, nz)
    Y, W = np.meshgrid(v0y, wz, indexing="ij")
    chi2 = np.zeros_like(Y)
    ygo = np.zeros_like(Y)
    for i in range(ny):
        for j in range(nz):
            th_ij = th0.copy()
            th_ij[1] = Y[i, j]; th_ij[6] = W[i, j]
            r = np.asarray(P.positions(th_ij[None, :])).reshape(-1)
            r0 = np.asarray(P.positions(th0[None, :])).reshape(-1)
            chi2[i, j] = np.sum(((r - R_OBS.reshape(-1)) / sig) ** 2)
            ygo[i, j] = B.fast_goal_plane(th_ij[:3], th_ij[4:7],
                                          float(np.exp(th_ij[3])), 1.0)[1]
    cs = ax.contour(Y, W, chi2, levels=[1.0, 4.0, 9.0],
                    colors=[PAL["green"], PAL["grey"], PAL["light"]],
                    linewidths=(1.0, 0.7, 0.5))
    ax.clabel(cs, fmt=r"$\Delta\chi^2=%.0f$", fontsize=6.5)
    c2 = ax.contour(Y, W, ygo, levels=[2.6, 2.8, 3.0, 3.2, 3.4],
                    colors=PAL["blue"], linewidths=0.6, linestyles="--")
    ax.clabel(c2, fmt=r"$y=%.1f$", fontsize=6.5)
    ax.plot([th0[1]], [th0[6]], "+", ms=8, color=PAL["rust"], mew=1.3)
    ax.set_xlabel(r"$v_{0y}$  (m s$^{-1}$)")
    ax.set_ylabel(r"$\omega_z$  (rad s$^{-1}$)")
    ax.set_title("data-equivalent ridge vs predicted entry")
    S.panel_label(axes[0], "a", dx=0.03, dy=0.96)
    S.panel_label(axes[1], "b", dx=0.97, dy=0.96, loc="tr")
    return S.save(fig, "f05-identifiability")


# ---------------------------------------------------------------------------
# f06 three members of the data-equivalent family
# ---------------------------------------------------------------------------
def f06():
    ch, _, _, _ = chain()
    # pick three samples spread along the sloppy direction v0y vs wz
    sc = ch[:, 1] / ch[:, 1].std() - ch[:, 6] / ch[:, 6].std()
    order = np.argsort(sc)
    picks = [order[int(0.10 * order.size)], order[int(0.50 * order.size)],
             order[int(0.90 * order.size)]]
    fig, axes = plt.subplots(1, 2, figsize=(S.TEXTW, 2.85),
                             gridspec_kw=dict(wspace=0.28))
    ax = axes[0]
    times = np.arange(0.0, 0.45, 0.0025)
    for i, (k, c) in enumerate(zip(picks, (PAL["rust"], PAL["green"], PAL["blue"]))):
        th = ch[k]
        tr = np.array(B.fast_positions(th[:3], th[4:7], float(np.exp(th[3])),
                                       1.0, list(times), dt=0.002, tau=6.0))
        yg = B.fast_goal_plane(th[:3], th[4:7], float(np.exp(th[3])), 1.0)[1]
        ax.plot(tr[:, 0], tr[:, 1], color=c, lw=1.2,
                label=f"$y={yg:.2f}$ m")
    for p in R_OBS:
        ax.errorbar([p[0]], [p[1]], xerr=0.01, yerr=0.01, fmt="o", ms=4.5,
                    color=PAL["ink"], mfc="none", mew=1.0, zorder=9)
    ax.axvline(11, color=PAL["ink"], lw=1.0)
    ax.set_xlim(-0.4, 11.6)
    ax.set_ylim(0.0, 4.2)
    ax.set_xlabel("$x$  (m)")
    ax.set_ylabel("$y$  (m)")
    ax.set_title("three fits, two fixes, three answers")
    ax.legend(loc="lower right", fontsize=7, title="entry $y$ (m)",
              title_fontsize=7)
    ax = axes[1]
    tp, yp, zp = chain(), None, None
    d = np.load(os.path.join(RES, "chains", "base_full.npz"))
    yy = d["y"]
    ax.hist(yy, bins=70, color=PAL["blue"], alpha=0.30, edgecolor="none",
            density=True)
    for k, c in zip(picks, (PAL["rust"], PAL["green"], PAL["blue"])):
        th = ch[k]
        res = B.fast_goal_plane(th[:3], th[4:7], float(np.exp(th[3])), 1.0)
        ax.axvline(res[1], color=c, lw=1.1)
    ax.set_xlabel("entry $y$  (m)")
    ax.set_ylabel("posterior predictive density")
    ax.set_yticks([])
    ax.set_title("the same fits land 0.5 m apart")
    S.panel_label(axes[0], "a", dx=0.03, dy=0.96)
    S.panel_label(axes[1], "b", dx=0.97, dy=0.96, loc="tr")
    return S.save(fig, "f06-family")


# ---------------------------------------------------------------------------
# f07 error budget
# ---------------------------------------------------------------------------
def f07():
    bud = load("budget.json")
    lab = {"obs0": r"$x_1$", "obs1": r"$y_1$", "obs2": r"$z_1$",
           "obs3": r"$x_2$", "obs4": r"$y_2$", "obs5": r"$z_2$",
           "data": "position noise (all)", "cd": r"$C_D$ prior",
           "spin_x": r"$\omega_x$ prior", "spin_y": r"$\omega_y$ prior",
           "spin_z": r"$\omega_z$ prior"}
    fig, axes = plt.subplots(1, 2, figsize=(S.TEXTW, 2.5),
                             gridspec_kw=dict(wspace=0.30))
    for ax, key, ttl in ((axes[0], "e1_sigma0.01", "constant-$a$ estimator"),
                         (axes[1], "e2_base", "physical estimator")):
        b = bud[key]
        share = b["share"]
        comps = [k for k in share if share[k][1] > 5e-3]
        vals = np.array([share[k][1] for k in comps])
        order = np.argsort(-vals)
        comps = [comps[i] for i in order]
        vals = vals[order]
        left = 0.0
        cols = [PAL["blue"], PAL["rust"], PAL["green"], PAL["purple"],
                PAL["gold"], PAL["grey"], PAL["light"]]
        for i, (k, v) in enumerate(zip(comps, vals)):
            ax.barh([0], [v], left=[left], height=0.55, color=cols[i % len(cols)],
                    edgecolor="white", lw=0.4, label=lab.get(k, k))
            if v > 0.06:
                ax.text(left + v / 2, 0, f"{v:.0%}", ha="center", va="center",
                        fontsize=7, color="white", fontweight="bold")
            left += v
        ax.set_xlim(0, 1.0)
        ax.set_yticks([])
        ax.set_xlabel("share of the variance of the entry-$y$ prediction")
        ax.set_title(ttl)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.34), ncol=2,
                  fontsize=6.8, handlelength=1.0, columnspacing=1.0)
    S.panel_label(axes[0], "a", dx=0.03, dy=0.96)
    S.panel_label(axes[1], "b", dx=0.97, dy=0.96, loc="tr")
    return S.save(fig, "f07-budget")


# ---------------------------------------------------------------------------
# f08 prior sensitivity
# ---------------------------------------------------------------------------
def f08():
    sp = load("sens_prior.json")
    base = sp["base"]
    sens = sp["sens"]
    rows = []
    for s in sens:
        rows.append(s)
    keys = [r.get("label", r.get("name", "?")) for r in rows]
    fig, axes = plt.subplots(1, 1, figsize=(S.TEXTW * 0.78, 0.34 * len(rows) + 1.3))
    ax = axes
    y0 = np.arange(len(rows))[::-1]
    for yy, r in zip(y0, rows):
        m, sd = r["mean"]["y"], r["std"]["y"]
        ax.errorbar([m], [yy], xerr=[2 * sd], fmt="o", ms=4.2, lw=1.2,
                    color=PAL["blue"], capsize=0)
    ax.axvline(base["mean"]["y"], color=PAL["ink"], lw=0.8, ls=":")
    ax.axvspan(base["mean"]["y"] - 2 * base["std"]["y"],
               base["mean"]["y"] + 2 * base["std"]["y"], color=PAL["light"],
               alpha=0.25, lw=0)
    ax.axvline(3.66, color=PAL["rust"], lw=1.0)
    ax.set_yticks(y0)
    ax.set_yticklabels(keys, fontsize=7.5)
    ax.set_xlabel(r"entry $y$  (m),  bars are $\pm2\sigma$")
    ax.set_xlim(2.2, 4.0)
    ax.set_ylim(-0.7, len(rows) - 0.3)
    ax.set_title("sensitivity to the prior widths")
    S.save(fig, "f08-sensitivity")


if __name__ == "__main__":
    todo = sys.argv[1:] or ["f01"]
    for name in todo:
        globals()[name]()
