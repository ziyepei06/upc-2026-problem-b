"""f22: how precise is real ball tracking?  Residuals of real 3-D reconstructions."""
import csv, sys
import numpy as np
import matplotlib.pyplot as plt
import figstyle as S
import paths

CSV = paths.dataset("SoccerNet-v3D-parsed.csv")


def load():
    rows = list(csv.DictReader(open(CSV)))
    opt, rep = [], []
    for r in rows:
        for key, out in (("optimized_error", opt), ("rep_error", rep)):
            try:
                v = float(r[key])
            except (TypeError, ValueError):
                continue
            if np.isfinite(v) and v > 0:
                out.append(v)
    return np.array(opt), np.array(rep)


def f22():
    opt, rep = load()
    fig, ax = plt.subplots(figsize=(S.TEXTW * 0.72, 2.55), layout="constrained")
    bins = np.logspace(-3, 1, 70)
    ax.hist(opt * 1e3, bins=bins * 1e3, color=S.PAL["blue"], alpha=0.55,
            label=f"optimised 3-D residual ({opt.size} balls)")
    ax.hist(rep * 1e3, bins=bins * 1e3, color=S.PAL["rust"], alpha=0.35,
            label="raw reprojection residual")
    ax.set_xscale("log")
    ax.set_xlabel("residual of a real 3-D ball reconstruction  (mm)")
    ax.set_ylabel("count")
    ymax = ax.get_ylim()[1]
    for xx, lab, cc in ((5, "printed coordinates 2.9 mm", S.PAL["ink"]),
                        (15, "goal-line technology 15 mm", S.PAL["green"]),
                        (53, "real tracking 50--65 mm", S.PAL["blue"])):
        ax.axvline(xx, color=cc, lw=0.9, ls="--")
        ax.text(xx * 1.05, ymax * 0.96, lab, fontsize=6.2, color=cc,
                rotation=90, va="top", ha="left")
    ax.legend(fontsize=6.4, loc="upper right")
    return S.save(fig, "f22-realdata")


if __name__ == "__main__":
    S.setup()
    print("ok f22" if f22() else "fail")
