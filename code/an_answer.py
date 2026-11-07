"""an_answer.py -- the headline answer: model-averaged predictive distribution.

The two estimators E1 (model-free constant acceleration) and E2 (physics
constrained, Bayesian) are both defensible and disagree by 0.5 m in y.  Since
six numbers cannot discriminate between them, we average their predictive
distributions with equal weights and quote the resulting (wider) uncertainty.
The difference between the two means is reported separately as the model
(misspecification) uncertainty.
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


def frac_region(y, z, kind):
    ay, az = np.abs(y), z
    if kind == "clean":
        return float(np.mean((ay < B.GOAL_HALF_WIDTH - B.BALL_RADIUS) &
                             (az > B.BALL_RADIUS) &
                             (az < B.GOAL_HEIGHT - B.BALL_RADIUS)))
    if kind == "center":
        return float(np.mean((ay < B.GOAL_HALF_WIDTH) & (az > 0) &
                             (az < B.GOAL_HEIGHT)))
    if kind == "wide":
        return float(np.mean(ay >= B.GOAL_HALF_WIDTH))
    if kind == "high":
        return float(np.mean(az >= B.GOAL_HEIGHT))
    if kind == "low":
        return float(np.mean(az <= 0))
    raise ValueError(kind)


def main():
    # --- E1 samples ------------------------------------------------------
    d1 = E.e1_predict(R, T, sigma=0.01, n_mc=40000)
    y1, z1, t1 = d1["_y"], d1["_z"], d1["_t"]

    # --- E2 posterior predictive -----------------------------------------
    P = E.PhysicalPosterior(R, T, sigma=0.01)
    out = P.sample(n=40000, burn=8000, thin=2, seed=20261107)
    t2, y2, z2 = P.predict(out["chain"])
    ok = np.isfinite(t2) & np.isfinite(y2) & np.isfinite(z2)
    t2, y2, z2 = t2[ok], y2[ok], z2[ok]
    chain = out["chain"]

    # diagnostics
    ch = chain
    ess = {}
    for i, nm in enumerate(["v0x", "v0y", "v0z", "lnCd", "wx", "wy", "wz"]):
        x = ch[:, i]
        x = x - x.mean()
        ac = np.correlate(x, x, "full")[x.size - 1:]
        ac /= ac[0]
        k = np.argmax(ac < 0.05)
        ess[nm] = float(x.size / max(k, 1))

    np.savez(os.path.join(OUT, "chains", "base_full.npz"),
             chain=chain, t=t2, y=y2, z=z2)
    dg = dict(acc=out["acc"], ess=ess, n_kept=int(chain.shape[0]),
              theta_mean=chain.mean(0).tolist(),
              theta_sd=chain.std(0).tolist(),
              cd_mean=float(np.exp(chain[:, 3]).mean()),
              cd_sd=float(np.exp(chain[:, 3]).std()),
              omega_abs_mean=float(np.linalg.norm(chain[:, 4:7], axis=1).mean()),
              omega_abs_sd=float(np.linalg.norm(chain[:, 4:7], axis=1).std()))

    # --- mixture ---------------------------------------------------------
    # equal weights; each component keeps its own spread
    m = min(y1.size, y2.size)
    rng = np.random.default_rng(0)
    pick1 = rng.random(m) < 0.5
    ym = np.where(pick1, rng.permutation(y1)[:m], rng.permutation(y2)[:m])
    zm = np.where(pick1, rng.permutation(z1)[:m], rng.permutation(z2)[:m])
    tm = np.where(pick1, rng.permutation(t1)[:m], rng.permutation(t2)[:m])

    def summ(y, z, t):
        return dict(
            mean=dict(t=float(np.mean(t)), y=float(np.mean(y)), z=float(np.mean(z))),
            std=dict(t=float(np.std(t, ddof=1)), y=float(np.std(y, ddof=1)),
                     z=float(np.std(z, ddof=1))),
            median=dict(t=float(np.median(t)), y=float(np.median(y)),
                        z=float(np.median(z))),
            pct={k: [float(np.percentile(v, p)) for p in (5, 16, 50, 84, 95)]
                 for k, v in (("t", t), ("y", y), ("z", z))},
            p_clean=frac_region(y, z, "clean"), p_center=frac_region(y, z, "center"),
            p_wide=frac_region(y, z, "wide"), p_high=frac_region(y, z, "high"),
            p_low=frac_region(y, z, "low"))

    res = dict(
        e1=dict(summary=summ(y1, z1, t1), point=dict(t=d1["t"], y=d1["y"], z=d1["z"]),
                std_analytic=d1["std_analytic"], v0=d1["v0"], a=d1["a"]),
        e2=dict(summary=summ(y2, z2, t2), diagnostics=dg),
        mix=dict(summary=summ(ym, zm, tm)),
        models=dict(
            ballistic_fit=dict(t=0.3706, y=3.217, z=1.957),
            ballistic_from_010=dict(t=0.3753, y=3.118, z=2.020),
            drag_only=dict(t=0.3881, y=3.217, z=1.919),
            full_physics=dict(t=0.3864, y=3.012, z=1.946),
            const_accel=dict(t=0.4050, y=2.495, z=1.920)),
    )
    res["spread"] = dict(
        y_min=2.495, y_max=3.217, y_range=0.722,
        z_min=1.919, z_max=2.020, z_range=0.101)

    with open(os.path.join(OUT, "answer.json"), "w") as f:
        json.dump(res, f, indent=1)

    for k in ("e1", "e2", "mix"):
        s = res[k]["summary"]
        print(f"{k}: y={s['mean']['y']:.3f}+-{s['std']['y']:.3f} "
              f"z={s['mean']['z']:.3f}+-{s['std']['z']:.3f} "
              f"t={s['mean']['t']:.4f}  P_clean={s['p_clean']:.3f} "
              f"P_wide={s['p_wide']:.3f} P_high={s['p_high']:.3f} P_low={s['p_low']:.3f}")
        print(f"   y68=[{s['pct']['y'][1]:.2f},{s['pct']['y'][3]:.2f}] "
              f"y95=[{s['pct']['y'][0]:.2f},{s['pct']['y'][4]:.2f}] "
              f"z68=[{s['pct']['z'][1]:.2f},{s['pct']['z'][3]:.2f}] "
              f"z95=[{s['pct']['z'][0]:.2f},{s['pct']['z'][4]:.2f}]")
    print("diag:", {k: (round(v, 2) if isinstance(v, float) else v) for k, v in dg.items() if k != "theta_mean" and k != "theta_sd"})


if __name__ == "__main__":
    main()
