"""Generate figures/example.pdf from code/blowdown.py.

The figure must be reproducible from the data that produced the numbers in
the paper: this script imports the same model, sweeps the hole diameter,
and writes a vector PDF at the final print size. Run from the paper root:

    python3 scripts/make_example_fig.py
"""

import math
import os
import sys

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "code"))

from blowdown import blowdown  # noqa: E402

mpl.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.family": "serif",
    "font.serif": ["DejaVu Serif"],
    "font.size": 9,
    "axes.labelsize": 9,
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.linewidth": 0.8,
    "lines.linewidth": 1.4,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
})

OKABE = ["#0072B2", "#D55E00", "#009E73"]

fig, ax = plt.subplots(figsize=(4.5, 2.8))
for d_cm, colour, style in [(5, OKABE[0], "-"), (10, OKABE[1], "--"),
                            (20, OKABE[2], ":")]:
    t, p, _ = blowdown(d=d_cm / 100.0, dt=0.5)
    ax.plot([ti / 60.0 for ti in t], [pi / p[0] for pi in p], style,
            color=colour, label=r"$d$ = %g cm" % d_cm)
    ax.plot(t[::max(1, len(t) // 8)][1:] and
            [ti / 60.0 for ti in t[::max(1, len(t) // 8)]][1:],
            [pi / p[0] for pi in p[::max(1, len(t) // 8)]][1:],
            "o", ms=3.5, mfc="white", mec=colour, mew=0.9)

ax.axhline(0.3, color="0.45", lw=0.7, ls="-.")
ax.set_xlabel(r"$t$ (min)")
ax.set_ylabel(r"$p/p_0$")
ax.set_xlim(0, 20)
ax.set_ylim(0, 1.02)
ax.legend(frameon=False, loc="upper right")
fig.savefig(os.path.join(ROOT, "figures", "example.pdf"))
print("wrote figures/example.pdf")
