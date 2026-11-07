"""an_data.py -- real tracking data: what position noise do we actually face?

Two open tracking datasets are used.  For every window of 9 consecutive
frames (0.36 s at 25 Hz) the ball position is fitted with a cubic in time; if
the motion is smooth (ball in free flight or rolling), the residual of that
fit is dominated by measurement noise, and sigma = sqrt(RSS/(n - p)) is an
unbiased estimate of the per-axis position error.
"""
from __future__ import annotations

import json
import os

import numpy as np

import paths

LIB = paths.DATASETS
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
os.makedirs(OUT, exist_ok=True)


def local_poly_sigma(t, x, win=9, deg=3, speed_min=6.0, speed_max=35.0):
    """Noise estimate from overlapping local polynomial fits."""
    t = np.asarray(t, float)
    x = np.asarray(x, float)
    if x.size < win + 2:
        return np.nan, 0, np.array([])
    dt = np.median(np.diff(t))
    v = np.gradient(x, dt)
    keep = []
    res_all = []
    for i in range(0, x.size - win, 2):
        sl = slice(i, i + win)
        if not (speed_min <= np.nanmax(np.abs(v[sl])) <= speed_max):
            continue
        tt = t[sl] - t[sl].mean()
        A = np.vander(tt, deg + 1)
        coef, *_ = np.linalg.lstsq(A, x[sl], rcond=None)
        r = x[sl] - A @ coef
        res_all.append(r)
        keep.append(i)
    if not res_all:
        return np.nan, 0, np.array([])
    R = np.vstack(res_all)
    rss = np.sum(R**2, axis=1)
    n_eff = win - (deg + 1)
    sig = np.sqrt(np.mean(rss) / n_eff)
    return float(sig), len(res_all), R.ravel()


def idsee():
    path = os.path.join(LIB, "IDSSE-J03WPY-ball.csv")
    import csv
    t, x, y = [], [], []
    with open(path) as f:
        rd = csv.reader(f)
        next(rd)
        for row in rd:
            try:
                t.append(float(row[1][11:19].replace(":", "")) if False else len(t) / 25.0)
                x.append(float(row[2]))
                y.append(float(row[3]))
            except (ValueError, IndexError):
                pass
    t = np.asarray(t)
    out = {}
    for name, arr in (("x", np.asarray(x)), ("y", np.asarray(y))):
        s, n, r = local_poly_sigma(t, arr)
        out[name] = dict(sigma=s, n_windows=n, rms=float(np.std(r)) if r.size else None)
    return out


def metrica():
    path = os.path.join(LIB, "Metrica-第1场-球轨迹.csv")
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    t = data[:, 2]
    bx = data[:, 3] * 105.0
    by = data[:, 4] * 68.0
    ok = np.isfinite(bx) & np.isfinite(by)
    t, bx, by = t[ok], bx[ok], by[ok]
    # split into contiguous segments (ball may leave the field of view)
    brk = np.where(np.diff(t) > 0.1)[0] + 1
    segs = np.split(np.arange(t.size), brk)
    print("Metrica segments:", len(segs), "longest", max(len(s) for s in segs))
    sig = {"x": [], "y": []}
    for s_ in segs:
        if len(s_) < 12:
            continue
        for name, arr in (("x", bx), ("y", by)):
            v, n, r = local_poly_sigma(t[s_], arr[s_])
            if np.isfinite(v):
                sig[name].append((v, n))
    out = {}
    for name in ("x", "y"):
        if sig[name]:
            w = np.array([n for _, n in sig[name]], float)
            v = np.array([a for a, _ in sig[name]], float)
            out[name] = dict(sigma=float(np.sum(v * w) / np.sum(w)),
                             n_windows=int(np.sum(w)), segments=len(sig[name]))
        else:
            out[name] = dict(sigma=float("nan"), n_windows=0)
    return out
    out = {}
    for name, arr in (("x", bx), ("y", by)):
        s, n, r = local_poly_sigma(t, arr)
        out[name] = dict(sigma=s, n_windows=n, rms=float(np.std(r)) if r.size else None)
    return out


if __name__ == "__main__":
    res = {"IDSSE": idsee(), "Metrica": metrica()}
    for k, v in res.items():
        print(k, {a: (round(b["sigma"], 4), b["n_windows"]) for a, b in v.items()})
    with open(os.path.join(OUT, "tracking_noise.json"), "w") as f:
        json.dump(res, f, indent=1)
