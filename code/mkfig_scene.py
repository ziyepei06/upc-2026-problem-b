"""Axonometric illustration of the penalty: pitch, goal, flight path, entry cloud."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Circle
import figstyle as S
import fit6 as F
import paths

S.setup()
B, H, R = 3.66, 2.44, 0.110
KX, KZ = 0.52, 0.30          # axonometric coefficients of the x axis


def P(x, y, z):
    x = np.asarray(x, float); y = np.asarray(y, float); z = np.asarray(z, float)
    return np.stack([y + KX * x, z + KZ * x], axis=-1)


def seg(ax, p0, p1, **kw):
    q = P(*zip(p0, p1)) if False else P(*np.array([p0, p1]).T)
    ax.plot(q[..., 0], q[..., 1], **kw)


def traj(dt=0.0012, tmax=0.44):
    th, _, _, _ = F.fit_map()
    cd = th[3]; beta = F.beta_vec(th)
    p = F.R0.copy(); v = th[:3].copy()
    ps = [p.copy()]; t = 0.0
    while t < tmax:
        p, v = F.rk4_step(p, v, dt, cd, beta)
        t += dt
        ps.append(p.copy())
    ps = np.array(ps)
    return ps[ps[:, 0] <= 11.0]


def f18():
    fig, ax = plt.subplots(figsize=(S.TEXTW * 0.98, 3.05), layout="constrained")
    ax.set_axis_off()
    # pitch
    quad = np.array([[-4.5, -11.0, 0.0], [12.4, -11.0, 0.0],
                     [12.4, 11.0, 0.0], [-4.5, 11.0, 0.0]])
    ax.add_patch(Polygon(P(*quad.T), closed=True, fc="#edf4e6", ec="none", zorder=0))
    # markings
    for y in (-9.16, 9.16):
        seg(ax, (5.5, y, 0), (11, y, 0), color="#c6d6bd", lw=0.7, zorder=1)
    seg(ax, (5.5, -9.16, 0), (5.5, 9.16, 0), color="#c6d6bd", lw=0.7, zorder=1)
    # penalty arc
    th = np.linspace(-0.93, 0.93, 60)
    arcx = 11 - 9.15 * np.cos(th); arcy = 9.15 * np.sin(th)
    qq = P(arcx, arcy, 0.0 * th)
    ax.plot(qq[..., 0], qq[..., 1], color="#c6d6bd", lw=0.7, zorder=1)
    # goal line
    seg(ax, (11, -9.16, 0), (11, 9.16, 0), color="white", lw=2.6, zorder=5)
    # net (back panel at x = 11 + 1.3)
    D = 1.3
    for y in np.linspace(-B, B, 13):
        seg(ax, (11 + D, y, 0), (11 + D, y, H), color=S.PAL["light"], lw=0.35, zorder=3)
    for z in np.linspace(0.2, H, 6):
        seg(ax, (11 + D, -B, z), (11 + D, B, z), color=S.PAL["light"], lw=0.35, zorder=3)
    for y in (-B, B):
        for z in np.linspace(0.2, H, 6):
            seg(ax, (11, y, z), (11 + D, y, z), color=S.PAL["light"], lw=0.3, zorder=3)
        quadn = np.array([[11, y, 0], [11 + D, y, 0], [11 + D, y, H], [11, y, H]])
        ax.add_patch(Polygon(P(*quadn.T), closed=True, fc="#f5f7f9", ec="none",
                             alpha=0.6, zorder=2))
    # goal frame
    for y in (-B, B):
        seg(ax, (11, y, 0), (11, y, H), color=S.PAL["ink"], lw=3.0,
            solid_capstyle="round", zorder=7)
    seg(ax, (11, -B, H), (11, B, H), color=S.PAL["ink"], lw=3.0, zorder=7)
    # flight path and ground shadow
    ps = traj()
    q = P(ps[:, 0], ps[:, 1], ps[:, 2])
    sh = P(ps[:, 0], ps[:, 1], np.zeros(len(ps)))
    ax.plot(sh[:, 0], sh[:, 1], color=S.PAL["blue"], lw=0.8, alpha=0.20,
            ls=(0, (4, 4)), zorder=2)
    ax.plot(q[:, 0], q[:, 1], color=S.PAL["blue"], lw=2.0, zorder=8,
            solid_capstyle="round")
    ax.add_patch(Circle(q[0], 0.055, fc="white", ec=S.PAL["ink"], lw=0.8, zorder=9))
    # entry cloud on the goal plane
    nom = np.load(paths.npz("geometry_nominal.npz"))
    yz = nom["center_yz"]; w = nom["w"]
    ok = np.isfinite(yz).all(axis=1); yz = yz[ok]; w = w[ok] / w[ok].sum()
    idx = np.random.default_rng(11).choice(len(yz), 2500, replace=False, p=w / w.sum())
    qc = P(np.full(len(idx), 10.985), yz[idx, 0], yz[idx, 1])
    ax.scatter(qc[:, 0], qc[:, 1], s=2.2, color=S.PAL["rust"], alpha=0.18, lw=0,
               zorder=9)
    ax.add_patch(Circle(P(11, 2.785, 1.945), 0.075, fc="none", ec=S.PAL["rust"],
                        lw=1.1, zorder=10))
    # dimensions
    ax.text(*(P(11.0, 0.0, 0.0) - np.array([0, 0.95])), "$7.32$ m",
            fontsize=7.0, ha="center", color=S.PAL["grey"])
    ax.text(*(P(11.0, B + 0.30, H) + np.array([0.12, 0.06])), "$2.44$ m",
            fontsize=7.0, color=S.PAL["grey"], va="center")
    ax.text(*(P(5.5, -6.6, 0.0)), "$11$ m", fontsize=7.0, ha="center",
            color=S.PAL["grey"])
    # annotations
    ax.text(*(P(4.4, -3.4, 2.55)), "$0.40$ s of flight", fontsize=7.2,
            color=S.PAL["blue"], ha="center")
    ax.text(*(P(0, -1.5, 0.0) + np.array([0, -0.22])), "penalty mark", fontsize=7.0,
            color=S.PAL["grey"], ha="center")
    ax.text(*(P(11.0, 0.0, H + 0.55)), "goal plane $x=11$ m", fontsize=7.0,
            color=S.PAL["grey"], ha="center")
    ax.set_aspect("equal")
    xs = np.concatenate([ax.get_xlim()])
    ax.set_xlim(-12.6, 17.4)
    ax.set_ylim(-2.4, 9.2)
    return S.save(fig, "f18-scene")


if __name__ == "__main__":
    f18()
    print("ok f18", flush=True)
