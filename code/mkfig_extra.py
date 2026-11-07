"""Flow chart and scene figures for the main line."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
import figstyle as S
import fit6 as F
import paths

S.setup()
B, H, R = 3.66, 2.44, 0.110


def _box(ax, x, y, w, h, text, fc, ec, fs=7.0, tc=None, lw=0.8):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.018",
                                fc=fc, ec=ec, lw=lw, zorder=2))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
            color=tc or S.PAL["ink"], zorder=3, linespacing=1.35)


def _arrow(ax, p0, p1, style="-|>", color=None, lw=0.8, ls="-", rad=0.0):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle=style, mutation_scale=8,
                                 color=color or S.PAL["grey"], lw=lw, ls=ls,
                                 connectionstyle=f"arc3,rad={rad}", zorder=1))


def f17():
    fig, ax = plt.subplots(figsize=(S.TEXTW * 0.95, 3.45), layout="constrained")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    fs = 6.6
    xs, w = 0.035, 0.405
    main = [
        (0.870, "Two fixes $\\mathbf{r}_1,\\mathbf{r}_2$ at $t=0.05,\\ 0.10$ s\n"
                "assumed $\\sigma_{\\rm pos}=5$ mm (rounding: 2.9 mm)"),
        (0.680, "Six-parameter force law\n"
                "$a=-F C_D|\\mathbf{u}|\\mathbf{u}+|\\mathbf{u}|"
                "(\\beta\\times\\mathbf{u})-g\\hat{\\mathbf{z}}$"),
        (0.490, "Penalised Gauss-Newton fit\n"
                "(damped, numerical Jacobian, 12$\\times$6)"),
        (0.300, "Laplace covariance $\\times 1.15$\n"
                "Student-$t$ proposal, exact importance weights"),
        (0.110, "RK4 + cubic-Hermite crossing at $x=11$ m\n"
                "131 072 draws, ESS 109 585"),
    ]
    for y, txt in main:
        _box(ax, xs, y, w, 0.105, txt, "#eef3fb", S.PAL["blue"], fs=fs)
    for i in range(len(main) - 1):
        _arrow(ax, (xs + w / 2, main[i][0] - 0.002),
               (xs + w / 2, main[i + 1][0] + 0.107))
    xr, wr = 0.575, 0.395
    side = [
        (0.752, 0.108, "model ladder M0--M6\neach step adds one force term",
         "#f5f5f5", S.PAL["grey"], 0.900),
        (0.604, 0.108, "calibration: 1000 synthetic truths\n"
                       "$\\lambda_y=0.96$, $\\lambda_z=4.94$",
         "#f5f5f5", S.PAL["grey"], 0.720),
        (0.470, 0.096, "validation: 500 (M3) + 50 (M6) truths\ncoverage 0.92--0.94 ($y$), 0.93--0.94 ($z$)",
         "#f5f5f5", S.PAL["grey"], 0.352),
    ]
    for y, h, txt, fc, ec, yy in side:
        _box(ax, xr, y, wr, h, txt, fc, ec, fs=fs)
        _arrow(ax, (xr - 0.006, yy), (xs + w + 0.008, yy), color=ec,
               ls=(0, (2, 2)))
    _box(ax, xr, 0.292, wr, 0.138,
         "Whole-ball clearance test\n"
         "$|y|<3.66-0.11$, $0.11<z<2.33$\n"
         "$y=2.785$ m, $z=1.945$ m\n$P(\\mathrm{goal})=0.971$",
         "#eaf4ea", S.PAL["green"], fs=fs)
    _arrow(ax, (xs + w + 0.008, 0.150), (xr - 0.006, 0.330),
           color=S.PAL["green"])
    _box(ax, xr, 0.028, wr, 0.182,
         "assumption scenarios\n"
         "(25 one-factor, 4 joint, 1 replication)\n"
         "envelope $y\\in[2.259,3.935]$ m\n"
         "$z\\in[1.169,2.698]$ m",
         "#fbf0ec", S.PAL["rust"], fs=fs)
    _arrow(ax, (xr + 0.20, 0.286), (xr + 0.20, 0.212), color=S.PAL["rust"])
    return S.save(fig, "f17-workflow")


def traj(dt=0.0015, tmax=0.46):
    th, _, _, _ = F.fit_map()
    cd = th[3]; beta = F.beta_vec(th)
    p = F.R0.copy(); v = th[:3].copy()
    ps = [p.copy()]; ts = [0.0]
    t = 0.0
    while t < tmax:
        p, v = F.rk4_step(p, v, dt, cd, beta)
        t += dt
        ps.append(p.copy()); ts.append(t)
    return np.array(ps), np.array(ts)


def f18():
    ps, ts = traj()
    hit = int(np.argmax(ps[:, 0] >= 11.0))
    ps = ps[:hit + 1]
    fig = plt.figure(figsize=(S.TEXTW * 0.98, 3.55), layout="constrained")
    ax = fig.add_subplot(111, projection="3d")
    ax.set_axis_off()
    gx, gy = np.meshgrid([-1.5, 12.4], [-11.0, 11.0])
    ax.plot_surface(gx, gy, np.zeros_like(gx), color="#f0f5ec", alpha=0.85,
                    linewidth=0, shade=False, zorder=0)
    # penalty area (16.5 m x 40.3 m) and 6-yard box
    for y in (-20.16, 20.16):
        ax.plot([-5.5, 11], [y, y], [0, 0], color=S.PAL["light"], lw=0.6)
    ax.plot([11, 11], [-20.16, 20.16], [0, 0], color="white", lw=1.6)
    for y in (-9.16, 9.16):
        ax.plot([5.5, 11], [y, y], [0, 0], color=S.PAL["light"], lw=0.6)
    ax.plot([5.5, 5.5], [-9.16, 9.16], [0, 0], color=S.PAL["light"], lw=0.6)
    ax.plot([11, 11], [-3.66, 3.66], [0, 0], color=S.PAL["ink"], lw=2.4)
    # goal frame and net
    for y in (-B, B):
        ax.plot([11, 11], [y, y], [0, H], color=S.PAL["ink"], lw=2.4, zorder=6)
    ax.plot([11, 11], [-B, B], [H, H], color=S.PAL["ink"], lw=2.4, zorder=6)
    for z in np.linspace(0, H, 6):
        ax.plot([12.0, 12.0], [-B, B], [z, z], color=S.PAL["light"], lw=0.3)
    for y in np.linspace(-B, B, 13):
        ax.plot([12.0, 12.0], [y, y], [0, H], color=S.PAL["light"], lw=0.3)
    for z in np.linspace(0, H, 6):
        ax.plot([11, 12.0], [-B, B], [z, z], color=S.PAL["light"], lw=0.25)
    # trajectory
    ax.plot(ps[:, 0], ps[:, 1], ps[:, 2], color=S.PAL["blue"], lw=1.9, zorder=7)
    ax.plot([0], [0], [0.11], marker="o", ms=4.5, color=S.PAL["blue"], zorder=8)
    ax.plot([ps[-1, 0]], [ps[-1, 1]], [ps[-1, 2]], marker="X", ms=8,
            color=S.PAL["rust"], zorder=8)
    nom = np.load(paths.npz("geometry_nominal.npz"))
    yz = nom["center_yz"]
    idx = np.arange(0, len(yz), 60)
    ax.scatter(np.full(len(idx), 10.97), yz[idx, 0], yz[idx, 1], s=1.0,
               color=S.PAL["rust"], alpha=0.16, depthshade=False, zorder=7)
    ax.text(12.2, 1.2, 2.62, "$7.32\\times2.44$ m", fontsize=7.2,
            color=S.PAL["ink"])
    ax.text(5.6, 2.4, 1.05, "0.40 s of flight", fontsize=7.2, color=S.PAL["blue"])
    ax.text(0.1, -2.6, 0.16, "penalty mark", fontsize=7.2, color=S.PAL["grey"])
    ax.text(10.4, 6.6, 0.10, "$x=11$ m", fontsize=7.2, color=S.PAL["grey"])
    ax.text(12.2, 3.4, 1.15, "$2.44$ m", fontsize=7.0, color=S.PAL["grey"],
            rotation=90)
    ax.set_xlim(0, 13.0); ax.set_ylim(-7.0, 7.0); ax.set_zlim(0, 4.6)
    ax.view_init(elev=13, azim=-104)
    try:
        ax.set_box_aspect((13.0, 14.0, 4.6), zoom=1.05)
    except Exception:
        pass
    return S.save(fig, "f18-scene")


if __name__ == "__main__":
    import sys
    for n in (sys.argv[1:] or ["f17", "f18"]):
        globals()[n]()
        print("ok", n, flush=True)
