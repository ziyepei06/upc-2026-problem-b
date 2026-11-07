"""an_valid.py -- synthetic (Monte-Carlo) validation of the two estimators.

Two truth models are used:

  physics truth   : Cd ~ lognormal(0.19, 0.15), spin magnitude ~ U(12.6, 63) rad/s
                    with an isotropic axis, kick parameters drawn so that the
                    ball enters the goal target drawn from the empirical
                    penalty placement distribution.
  kinematic truth : r(t) = r0 + v0 t + a t^2/2 exactly (the model that reproduces
                    the scenario fixes), i.e. a deliberate misspecification.

For each trial we simulate two noisy fixes (sigma = 0.01 m), run E1 and E2 and
compare the predictive distributions with the true entry point.  We report
bias, RMSE and the coverage of the 68/90 % intervals, which is what decides
whether the quoted uncertainty can be trusted.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

import bkick as B
import estimators as E

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
R_SCEN = np.array([[1.50, 0.45, 0.44], [2.96, 0.86, 0.74]])
T_OBS = (0.05, 0.10)
SIGMA = 0.010


# ---------------------------------------------------------------------------
def shoot(speed, cd, omega, y_t, z_t, cl=1.0, tau=6.0, dt=0.002, iters=12):
    """Initial velocity that sends the ball through (11, y_t, z_t)."""
    t_guess = 11.0 / max(speed * 0.92, 1.0)
    vy = y_t / t_guess
    vz = (z_t - B.KICK_START[2]) / t_guess + 0.5 * B.G * t_guess
    for _ in range(iters):
        vx = np.sqrt(max(speed**2 - vy**2 - vz**2, 1.0))
        v0 = np.array([vx, vy, vz])
        f = np.array(B.fast_goal_plane(v0, omega, cd, cl, dt=dt, tau=tau)[1:])
        res = f - np.array([y_t, z_t])
        if np.linalg.norm(res) < 1e-5:
            break
        Jm = np.zeros((2, 2))
        for j, h in enumerate((1e-4, 1e-4)):
            vp = v0.copy()
            vp[1 + j] += h
            vp[0] = np.sqrt(max(speed**2 - vp[1]**2 - vp[2]**2, 1.0))
            fp = np.array(B.fast_goal_plane(vp, omega, cd, cl, dt=dt, tau=tau)[1:])
            Jm[:, j] = (fp - f) / h
        try:
            step = np.linalg.solve(Jm, -res)
        except np.linalg.LinAlgError:
            break
        vy += np.clip(step[0], -2, 2)
        vz += np.clip(step[1], -2, 2)
    vx = np.sqrt(max(speed**2 - vy**2 - vz**2, 1.0))
    return np.array([vx, vy, vz])


def draw_truth(rng):
    speed = float(np.clip(rng.normal(27.5, 2.5), 21.0, 33.0))
    # empirical penalty placement: 31.3 % low, 36.0 % middle, 32.7 % high
    u = rng.random()
    if u < 0.313:
        z_t = rng.uniform(0.15, 0.813)
    elif u < 0.673:
        z_t = rng.uniform(0.813, 1.628)
    else:
        z_t = rng.uniform(1.628, 2.30)
    y_t = rng.choice([-1, 1]) * rng.uniform(0.3, 3.3)
    cd = float(rng.lognormal(np.log(0.19), 0.15))
    mag = rng.uniform(12.6, 63.0)
    axis = rng.normal(size=3)
    axis /= np.linalg.norm(axis)
    omega = axis * mag
    v0 = shoot(speed, cd, omega, y_t, z_t)
    return dict(speed=speed, cd=cd, omega=omega, v0=v0, y_t=y_t, z_t=z_t)


def truth_positions(tr, times=(0.05, 0.10)):
    return np.array(B.fast_positions(tr["v0"], tr["omega"], tr["cd"], 1.0,
                                     list(times), tau=6.0))


def run(n_trials=150, truth="physics", seed=7, sigma=SIGMA, n_mc=6000):
    rng = np.random.default_rng(seed)
    rec = []
    for i in range(n_trials):
        if truth == "physics":
            tr = draw_truth(rng)
            v0, a = tr["v0"], None
        else:
            # kinematic truth: uniform acceleration, drawn to resemble the scenario
            tr = dict(v0=np.array([30.4, 9.4, 6.9]) + rng.normal(0, 0.4, 3),
                      a=np.array([-16.0, -16.0, -12.0]) + rng.normal(0, 1.0, 3))
            tr["r"] = np.array([B.KICK_START + tr["v0"] * t + 0.5 * tr["a"] * t**2
                                for t in T_OBS])
        if truth == "physics":
            r_true = truth_positions(tr)
        else:
            r_true = tr["r"]
        r_obs = r_true + rng.normal(0, sigma, r_true.shape)

        d1 = E.e1_predict(r_obs, T_OBS, sigma=sigma, n_mc=n_mc)
        d2 = E.e2_predict_quick(r_obs, T_OBS, sigma=sigma, seed=i, n=3000,
                                burn=800, thin=1)
        # truth entry point
        if truth == "physics":
            t_t, y_t, z_t = B.fast_goal_plane(tr["v0"], tr["omega"], tr["cd"], 1.0)
        else:
            tq, yq, zq = E.quadratic_plane(tr["v0"][None, :], tr["a"][None, :])
            t_t, y_t, z_t = float(tq[0]), float(yq[0]), float(zq[0])
        rec.append(dict(i=i, y_true=float(y_t), z_true=float(z_t),
                        t_true=float(t_t),
                        e1_y=d1["mean"]["y"], e1_z=d1["mean"]["z"],
                        e1_sy=d1["std"]["y"], e1_sz=d1["std"]["z"],
                        e1_pin=d1["p_inside"], e1_pct=d1["pct"],
                        e2_y=d2["mean"]["y"], e2_z=d2["mean"]["z"],
                        e2_sy=d2["std"]["y"], e2_sz=d2["std"]["z"],
                        e2_pin=d2["p_inside"], e2_pct=d2["pct"]))
        if (i + 1) % 25 == 0:
            print(f"  {truth} trial {i+1}/{n_trials}", flush=True)
    return rec


def score(rec, tag):
    yt = np.array([r["y_true"] for r in rec])
    zt = np.array([r["z_true"] for r in rec])
    out = dict(tag=tag, n=len(rec))
    for est in ("e1", "e2"):
        y = np.array([r[f"{est}_y"] for r in rec])
        z = np.array([r[f"{est}_z"] for r in rec])
        sy = np.array([r[f"{est}_sy"] for r in rec])
        sz = np.array([r[f"{est}_sz"] for r in rec])
        pct = [r[f"{est}_pct"] for r in rec]
        cov = {}
        for lvl, (lo, hi) in ((68, (1, 3)), (90, (0, 4))):
            cov[f"y{lvl}"] = float(np.mean([(p["y"][lo] <= r["y_true"] <= p["y"][hi])
                                            for p, r in zip(pct, rec)]))
            cov[f"z{lvl}"] = float(np.mean([(p["z"][lo] <= r["z_true"] <= p["z"][hi])
                                            for p, r in zip(pct, rec)]))
        out[est] = dict(
            bias_y=float(np.mean(y - yt)), bias_z=float(np.mean(z - zt)),
            rmse_y=float(np.sqrt(np.mean((y - yt)**2))),
            rmse_z=float(np.sqrt(np.mean((z - zt)**2))),
            mean_sy=float(sy.mean()), mean_sz=float(sz.mean()),
            cov=cov,
            p_in_pred=float(np.mean([r[f"{est}_pin"] for r in rec])),
            p_in_true=float(np.mean((np.abs(yt) < B.GOAL_HALF_WIDTH - B.BALL_RADIUS) &
                                    (zt > B.BALL_RADIUS) &
                                    (zt < B.GOAL_HEIGHT - B.BALL_RADIUS))),
        )
    return out


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 120
    res = {}
    for truth in ("physics", "kinematic"):
        print(f"== {truth} truth ==", flush=True)
        rec = run(n_trials=n, truth=truth, seed=11 if truth == "physics" else 12)
        res[truth] = dict(trials=rec, score=score(rec, truth))
        s = res[truth]["score"]
        for est in ("e1", "e2"):
            e = s[est]
            print(f"  {est}: bias_y={e['bias_y']:+.3f} rmse_y={e['rmse_y']:.3f} "
                  f"mean_sy={e['mean_sy']:.3f} cov68y={e['cov']['y68']:.2f} "
                  f"cov90y={e['cov']['y90']:.2f} cov68z={e['cov']['z68']:.2f} "
                  f"Pin_pred={e['p_in_pred']:.2f} (true {e['p_in_true']:.2f})")
    with open(os.path.join(OUT, "valid.json"), "w") as f:
        json.dump(res, f)
    print("wrote results/valid.json")
