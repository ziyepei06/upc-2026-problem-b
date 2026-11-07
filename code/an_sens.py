"""an_sens.py -- parameter sweeps, identifiability, error budget."""
from __future__ import annotations

import json
import os

import numpy as np

import bkick as B
import estimators as E

R = np.array([[1.50, 0.45, 0.44], [2.96, 0.86, 0.74]])
T = (0.05, 0.10)
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
os.makedirs(OUT, exist_ok=True)


def bayes_point(sigma=0.01, n=8000, burn=1500, thin=2, seed=1, **kw):
    P = E.PhysicalPosterior(R, T, sigma=sigma, **kw)
    out = P.sample(n=n, burn=burn, thin=thin, seed=seed)
    t, y, z = P.predict(out["chain"])
    s = E.summarize(t, y, z)
    s["acc"] = out["acc"]
    s["theta_mean"] = out["chain"].mean(0).tolist()
    s["theta_sd"] = out["chain"].std(0).tolist()
    s["theta_map"] = out["theta_map"].tolist()
    s["omega_abs_mean"] = float(np.mean(np.linalg.norm(out["chain"][:, 4:7], axis=1)))
    return s, out


def main():
    res = {}

    # ---- 1. E1 noise sweep ------------------------------------------------
    e1 = []
    for sg in (0.0029, 0.005, 0.010, 0.015, 0.020, 0.030):
        d = E.e1_predict(R, T, sigma=sg, n_mc=20000)
        e1.append(dict(sigma=sg, t=d["t"], y=d["y"], z=d["z"],
                       std=d["std"], pct=d["pct"], p_inside=d["p_inside"],
                       p_center_inside=d["p_center_inside"], a=d["a"], v0=d["v0"]))
        print(f"E1 sigma={sg:.4f}  y={d['y']:.3f}+-{d['std']['y']:.3f} "
              f"z={d['z']:.3f}+-{d['std']['z']:.3f} P_in={d['p_inside']:.3f}")
    res["e1_sweep"] = e1

    # ---- 2. E2 noise sweep ------------------------------------------------
    e2 = []
    for sg in (0.005, 0.010, 0.015, 0.020):
        s, _ = bayes_point(sigma=sg)
        e2.append(dict(sigma=sg, **{k: s[k] for k in
                                    ("mean", "std", "pct", "p_inside",
                                     "p_center_inside", "acc", "theta_mean",
                                     "theta_sd", "omega_abs_mean")}))
        print(f"E2 sigma={sg:.3f}  y={s['mean']['y']:.3f}+-{s['std']['y']:.3f} "
              f"z={s['mean']['z']:.3f}+-{s['std']['z']:.3f} P_in={s['p_inside']:.3f}")
    res["e2_sweep"] = e2

    # ---- 3. prior sensitivity --------------------------------------------
    base, _ = bayes_point()
    res["e2_base"] = base
    print("base:", np.round(base["mean"]["y"], 3), np.round(base["std"]["y"], 3))
    sens = []
    for key, vals in (("cd_log_sigma", (0.15, 0.30, 0.50)),
                      ("cd_median", (0.17, 0.19, 0.22)),
                      ("omega_sigma", (17.0, 35.0, 70.0)),
                      ("cl", (0.7, 1.0, 1.3)),
                      ("tau", (2.0, 6.0, 12.0)),
                      ("sigma", (0.005, 0.010, 0.015))):
        for v in vals:
            if key == "cd_log_sigma":
                kw = {"cd": (E.CD_MEDIAN, v)}
            elif key == "cd_median":
                kw = {"cd": (v, E.CD_LOG_SIGMA)}
            else:
                kw = {key: v}
            s, _ = bayes_point(**kw)
            sens.append(dict(parkey=key, value=v, mean=s["mean"], std=s["std"],
                             p_inside=s["p_inside"],
                             omega_abs_mean=s["omega_abs_mean"],
                             cd=float(np.exp(s["theta_mean"][3]))))
            print(f"  {key:13s}={v:<6} y={s['mean']['y']:.3f}+-{s['std']['y']:.3f} "
                  f"z={s['mean']['z']:.3f} P_in={s['p_inside']:.3f}")
    res["prior_sens"] = sens

    # ---- 4. Fisher information / sloppy modes (E1 parameterisation) ------
    P = E.PhysicalPosterior(R, T)
    th, _ = P.map()
    Jm = np.zeros((6, 7))
    for i in range(7):
        h = 1e-6 * max(1.0, abs(th[i]))
        tp = th.copy()
        tp[i] += h
        Jm[:, i] = (P.positions(tp)[0] - P.positions(th)[0]).ravel() / h
    F = Jm.T @ Jm / 0.01**2
    w, V = np.linalg.eigh(F)
    order = np.argsort(w)[::-1]
    w, V = w[order], V[:, order]
    # semi-log sloppy spectrum; normalised eigenvector composition
    sloppy = dict(eigenvalues=w.tolist(),
                  sigma_eig=(1.0 / np.sqrt(np.maximum(w, 1e-30))).tolist(),
                  vectors=[[float(V[i, j]) for i in range(7)] for j in range(7)],
                  labels=["v0x", "v0y", "v0z", "lnCd", "wx", "wy", "wz"])
    # E1: (v0,a) Fisher from the data jacobian
    G = E.e1_data_jacobian(T)
    F1 = np.linalg.inv(G @ (0.01**2 * np.eye(6)) @ G.T)
    w1, V1 = np.linalg.eigh(F1)
    o1 = np.argsort(w1)[::-1]
    sloppy1 = dict(eigenvalues=w1[o1].tolist(),
                   sigma_eig=(1.0 / np.sqrt(np.maximum(w1[o1], 1e-30))).tolist(),
                   vectors=[[float(V1[i, j]) for i in range(6)] for j in o1],
                   labels=["v0x", "v0y", "v0z", "ax", "ay", "az"])
    res["fim"] = dict(full=sloppy, e1=sloppy1)
    print("FIM(7) eigen sigma:", np.round(sloppy["sigma_eig"], 3))
    print("FIM(6) eigen sigma:", np.round(sloppy1["sigma_eig"], 3))

    # ---- 5. error budget via the delta method on E1 ----------------------
    v0, a = B.constant_accel_inverse(R, T)
    base_t, base_y, base_z = E.quadratic_plane(v0[None, :], a[None, :])
    base = np.array([base_t[0], base_y[0], base_z[0]])
    sig_par = np.sqrt(np.diag(E.e1_data_jacobian(T) @ (0.01**2 * np.eye(6))
                              @ E.e1_data_jacobian(T).T))
    names = ["v0x", "v0y", "v0z", "ax", "ay", "az"]
    contrib = {}
    for target, idx in (("y", 1), ("z", 2), ("t", 0)):
        g = np.zeros(6)
        for i in range(6):
            p = np.concatenate([v0, a])
            h = 1e-6 * max(1.0, abs(p[i]))
            p[i] += h
            tq, yq, zq = E.quadratic_plane(p[None, :3], p[None, 3:])
            vec = np.array([tq[0], yq[0], zq[0]])
            g[i] = (vec[idx] - base[idx]) / h
        c = np.abs(g) * sig_par
        contrib[target] = dict(total=float(np.sqrt(np.sum(c**2))),
                               terms={n: float(x) for n, x in zip(names, c)},
                               grad={n: float(x) for n, x in zip(names, g)})
    res["budget"] = contrib
    print("budget y:", {k: round(v, 3) for k, v in contrib["y"]["terms"].items()},
          "total", round(contrib["y"]["total"], 3))

    with open(os.path.join(OUT, "sens.json"), "w") as f:
        json.dump(res, f, indent=1)
    print("wrote results/sens.json")


if __name__ == "__main__":
    main()
