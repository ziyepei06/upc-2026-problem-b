"""figstyle.py -- one visual language for every figure in the paper."""
from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

FIGDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures")
os.makedirs(FIGDIR, exist_ok=True)

TEXTW = 6.5          # inches, matches the 1 in margins on letter paper

PAL = dict(
    ink="#1b1b1b",
    blue="#1A5FB4",
    rust="#A63A1E",
    green="#2E7D32",
    purple="#6A3D9A",
    grey="#6b6b6b",
    light="#b8b8b8",
    gold="#B8860B",
)


def setup():
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["STIX Two Text", "STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": 9,
        "axes.titlesize": 9.5,
        "axes.labelsize": 9,
        "legend.fontsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "axes.linewidth": 0.6,
        "axes.edgecolor": PAL["ink"],
        "axes.labelcolor": PAL["ink"],
        "text.color": PAL["ink"],
        "xtick.color": PAL["ink"],
        "ytick.color": PAL["ink"],
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.minor.width": 0.4,
        "ytick.minor.width": 0.4,
        "xtick.top": True,
        "ytick.right": True,
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        "axes.grid": False,
        "legend.frameon": False,
        "legend.handlelength": 1.6,
        "legend.borderaxespad": 0.2,
        "lines.linewidth": 1.1,
        "lines.solid_capstyle": "round",
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": None,
        
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def save(fig, name):
    eng = fig.get_layout_engine()
    if eng is None:
        fig.tight_layout(pad=0.25)
    else:
        eng.set(w_pad=0.035, h_pad=0.035, wspace=0.03, hspace=0.05)
    path = os.path.join(FIGDIR, name + ".pdf")
    fig.savefig(path)
    plt.close(fig)
    print("wrote", path)
    return path


def panel_label(ax, text, loc="tl", dx=0.02, dy=0.97):
    x, y = (0.02, 0.97) if loc == "tl" else (0.98, 0.97)
    ha = "left" if loc == "tl" else "right"
    ax.text(x + (dx - 0.02), y + (dy - 0.97), text, transform=ax.transAxes,
            ha=ha, va="top", fontsize=9, fontweight="bold")
