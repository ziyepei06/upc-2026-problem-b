"""mplhelp.py -- small plotting helpers shared by the figure scripts."""
from __future__ import annotations

import numpy as np


def hist2d_contour(ax, y, z, levels=(0.28, 0.60, 0.90), color=None, nbin=60,
                   ylim=(-4.2, 4.2), zlim=(-0.2, 3.2), smooth=1.2, w=None,
                   label=None):
    """Smoothed 2-D density contour enclosing the requested probability mass."""
    H, ye, ze = np.histogram2d(y, z, bins=(nbin, nbin), weights=w,
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


