"""an_design.py -- design study with the linearised (Fisher) analysis.

Question: how accurately must the ball be tracked, and how many fixes are
needed, to pin the entry point down?  We linearise the physical model around
the MAP solution of the scenario and propagate the data covariance through the
posterior, which is exact to first order and cheap enough for large sweeps.
The linearisation is checked against the Monte-Carlo result in an_check.py.
"""
from __future__ import annotations

import json
import os

import numpy as np

import bkick as B
import estimators as E

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
R = np.array([[1.50, 0.45, 0.44], [2.96, 0.86, 0.74]])
T = (0.05, 0.10)


def build(theta, cl=1.0, tau=6.0, dt=0.002):
    cd = float(np.exp(theta[3]))
    return lambda t: np.asarray(B.fast_positions(theta[0:3], theta[4:7], cd, cl,
                                                [t], dt=dt, tau=tau)[0])


def jac_positions(theta, times, cl=1.0, tau=6.0, h=1e-5):
    f = build(theta, cl, tau)
    J = np.zeros((3 * len(times), 7))
    for k, t in enumerate(times):
        p0 = f(t)
        for i in range(7):
            tp = theta.copy()
            tp[i] += h * max(1.0, abs(theta[i]))
            J[3 * k:3 * k + 3, i] = (build(tp, cl, tau)(t) - p0) / (h * max(1.0, abs(theta[i])))
    return J


def grad_entry(theta, cl=1.0, tau=6.0, h=1e-5):
    cd = float(np.exp(theta[3]))
    g0 = np.array(B.fast_goal_plane(theta[0:3], theta[4:7], cd, cl, tau=tau))
    G = np.zeros((3, 7))
    for i in range(7):
        tp = theta.copy()
        tp[i] += h * max(1.0, abs(theta[i]))
        gi = np.array(B.fast_goal_plane(tp[0:3], tp[4:7], float(np.exp(tp[3])),
                                        cl, tau=tau))
        G[:, i] = (gi - g0) / (h * max(1.0, abs(theta[i])))
    return g0, G


PRIOR_SD = np.array([0, 0, 0, E.CD_LOG_SIGMA, E.OMEGA_SIGMA, E.OMEGA_SIGMA,
                     E.OMEGA_SIGMA])


def sigma_pred(theta, times, sigma_pos, prior_sd=PRIOR_SD, cl=1.0, tau=6.0):
    J = jac_positions(theta, times, cl, tau) / sigma_pos
    P = np.zeros((7, 7))
    for i in range(7):
        if prior_sd[i] > 0:
            P[i, i] = 1.0 / prior_sd[i]**2
    C = np.linalg.inv(J.T @ J + P)
    _, G = grad_entry(theta, cl, tau)
    Cyz = G @ C @ G.T
    return dict(sigma_t=float(np.sqrt(Cyz[0, 0])), sigma_y=float(np.sqrt(Cyz[1, 1])),
                sigma_z=float(np.sqrt(Cyz[2, 2])), cov=C.tolist())


def main():
    P = E.PhysicalPosterior(R, T)
    theta, _ = P.map()
    res = dict(theta=theta.tolist())

    # --- (a) number of fixes and position accuracy ------------------------
    grid = {}
    for n in (2, 3, 4, 5, 6, 8):
        times = tuple(np.round(0.05 * np.arange(1, n + 1), 4))
        for smm in (1, 2, 3, 5, 10, 15, 20, 30, 50):
            d = sigma_pred(theta, times, smm / 1000.0)
            grid[f"n{n}_s{smm}"] = d
    res["grid"] = dict(
        n=[2, 3, 4, 5, 6, 8],
        sigma_mm=[1, 2, 3, 5, 10, 15, 20, 30, 50],
        sigma_y=[[grid[f"n{n}_s{s}"]["sigma_y"] for s in (1, 2, 3, 5, 10, 15, 20, 30, 50)]
                 for n in (2, 3, 4, 5, 6, 8)],
        sigma_z=[[grid[f"n{n}_s{s}"]["sigma_z"] for s in (1, 2, 3, 5, 10, 15, 20, 30, 50)]
                 for n in (2, 3, 4, 5, 6, 8)],
    )
    for n in (2, 3, 4, 6):
        print(f"n={n}: sigma_y(mm) for sigma_pos=1,5,10,20,50 mm: " +
              ", ".join(f"{grid[f'n{n}_s{s}']['sigma_y']*1000:.0f}"
                        for s in (1, 5, 10, 20, 50)))

    # --- (b) knowing the spin would be equivalent to ... ------------------
    spin_known = {}
    for smm in (5, 10, 15, 20, 50):
        psd = PRIOR_SD.copy()
        psd[4:] = 0.6          # 0.1 rev/s
        d = sigma_pred(theta, (0.05, 0.10), smm / 1000.0, prior_sd=psd)
        spin_known[f"s{smm}"] = d
    res["spin_known"] = spin_known
    print("spin known (2 fixes):", {k: round(v["sigma_y"] * 1000) for k, v in spin_known.items()})

    # --- (c) where should the fixes be? -----------------------------------
    tgrid = np.round(np.arange(0.03, 0.31, 0.01), 3)
    SY = np.full((len(tgrid), len(tgrid)), np.nan)
    for i, t1 in enumerate(tgrid):
        for j, t2 in enumerate(tgrid):
            if t2 <= t1 + 0.015:
                continue
            SY[i, j] = sigma_pred(theta, (float(t1), float(t2)), 0.01)["sigma_y"]
    res["timing"] = dict(t=tgrid.tolist(), sigma_y=SY.tolist())
    ij = np.unravel_index(np.nanargmin(SY), SY.shape)
    print(f"best two-fix timing: t1={tgrid[ij[0]]:.2f} t2={tgrid[ij[1]]:.2f} "
          f"sigma_y={SY[ij]*1000:.0f} mm (at sigma_pos=10 mm)")
    print(f"problem timing 0.05/0.10: sigma_y={sigma_pred(theta,(0.05,0.10),0.01)['sigma_y']*1000:.0f} mm")

    # --- (d) base case, to compare with the Monte-Carlo -------------------
    res["base"] = sigma_pred(theta, (0.05, 0.10), 0.01)
    with open(os.path.join(OUT, "design.json"), "w") as f:
        json.dump(res, f, indent=1)
    print("wrote results/design.json")


if __name__ == "__main__":
    main()
