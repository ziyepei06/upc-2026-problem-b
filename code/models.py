"""models.py -- alternative flight models compared in Section (model comparison).

M0  gravity only                                 (3 parameters, v0)
M1  gravity + quadratic drag                     (4 parameters, v0, Cd)
M2  gravity + drag + Magnus, prior-regularised   (7 parameters) == estimator E2
M3  constant acceleration, exact interpolation   (6 parameters) == estimator E1

Every fit returns chi^2 over the six measured coordinates (two fixes) so the
models can be compared on the same footing, plus the goal-plane prediction.
"""

from __future__ import annotations

import numpy as np

import bkick as B
import estimators as E

R_OBS = np.array([[1.50, 0.45, 0.44], [2.96, 0.86, 0.74]])
T_OBS = (0.05, 0.10)


# ---------------------------------------------------------------------------
# M0: gravity only
# ---------------------------------------------------------------------------
def fit_gravity_only(r_obs=R_OBS, t_obs=T_OBS):
    t1, t2 = t_obs
    r1, r2 = np.asarray(r_obs, float)
    zhat = np.array([0.0, 0.0, 1.0])
    # weighted least squares in v0: v0 t_i = r_i - r0 + g t_i^2/2 zhat
    num = t1 * (r1 - B.KICK_START + 0.5 * B.G * t1**2 * zhat) \
        + t2 * (r2 - B.KICK_START + 0.5 * B.G * t2**2 * zhat)
    v0 = num / (t1**2 + t2**2)
    pred = np.array([B.KICK_START, B.KICK_START]) + np.outer(t_obs, v0) \
        - 0.5 * B.G * np.array(t_obs)[:, None]**2 * zhat
    chi2 = float(np.sum((pred - np.asarray(r_obs, float))**2))
    t, y, z, _ = B.Flight(cd=0.0).goal_plane(B.KICK_START, v0, np.zeros(3))
    return dict(name="M0 gravity only", n_par=3, v0=v0.tolist(), chi2=chi2 / E.SIGMA_POS**2,
                t=t, y=float(y), z=float(z), pred=pred)


# ---------------------------------------------------------------------------
# M1: gravity + drag (no spin)
# ---------------------------------------------------------------------------
def fit_drag(r_obs=R_OBS, t_obs=T_OBS, sigma=E.SIGMA_POS, cd0=0.19,
             cd_log_sigma=1e6, iters=200):
    r_obs = np.asarray(r_obs, float)
    theta = np.array([30.4, 9.4, 6.9, np.log(cd0)])
    ln_cd_mean = np.log(cd0)

    def resid(th):
        pred = np.array(B.fast_positions(th[0:3], np.zeros(3), float(np.exp(th[3])),
                                         1.0, t_obs))
        a = ((pred - r_obs) / sigma).ravel()
        b = (th[3] - ln_cd_mean) / cd_log_sigma
        return np.concatenate([a, [b]])

    mu = 1e-3
    R = resid(theta)
    cost = 0.5 * R @ R
    for _ in range(iters):
        J = np.zeros((R.size, 4))
        for i in range(4):
            h = 1e-6 * max(1.0, abs(theta[i]))
            tp = theta.copy()
            tp[i] += h
            J[:, i] = (resid(tp) - R) / h
        A = J.T @ J
        g = J.T @ R
        step = np.linalg.solve(A + mu * np.diag(np.diag(A) + 1e-12), -g)
        trial = theta + step
        Rt = resid(trial)
        ct = 0.5 * Rt @ Rt
        if ct < cost:
            theta, R, cost = trial, Rt, ct
            mu = max(mu * 0.3, 1e-12)
        else:
            mu *= 10.0
            if mu > 1e10:
                break
    pred = np.array(B.fast_positions(theta[0:3], np.zeros(3), float(np.exp(theta[3])),
                                     1.0, t_obs))
    chi2 = float(np.sum(((pred - r_obs) / sigma)**2))
    t, y, z = B.fast_goal_plane(theta[0:3], np.zeros(3), float(np.exp(theta[3])), 1.0)
    return dict(name="M1 gravity+drag", n_par=4, v0=theta[0:3].tolist(),
                cd=float(np.exp(theta[3])), chi2=chi2, t=t, y=float(y), z=float(z),
                pred=pred)


# ---------------------------------------------------------------------------
# M2: full physics, prior-regularised (E2)
# ---------------------------------------------------------------------------
def fit_full(r_obs=R_OBS, t_obs=T_OBS, sigma=E.SIGMA_POS, omega_sigma=E.OMEGA_SIGMA,
             cd=(E.CD_MEDIAN, E.CD_LOG_SIGMA)):
    P = E.PhysicalPosterior(r_obs, t_obs, sigma=sigma, cd=cd, omega_sigma=omega_sigma)
    theta, lp = P.map()
    cov = P.laplace_cov(theta)
    pred = P.positions(theta)[0]
    chi2 = float(np.sum(((pred - r_obs) / sigma)**2))
    t, y, z = B.fast_goal_plane(theta[0:3], theta[4:7], float(np.exp(theta[3])),
                                P.cl, tau=P.tau)
    return dict(name="M2 drag+Magnus", n_par=7, theta=theta.tolist(),
                v0=theta[0:3].tolist(), cd=float(np.exp(theta[3])),
                omega=theta[4:7].tolist(), chi2=chi2, t=t, y=float(y), z=float(z),
                pred=pred, cov=cov.tolist(), P=P)


# ---------------------------------------------------------------------------
# M3: constant acceleration (E1)
# ---------------------------------------------------------------------------
def fit_const_accel(r_obs=R_OBS, t_obs=T_OBS):
    r_obs = np.asarray(r_obs, float)
    v0, a = B.constant_accel_inverse(r_obs, t_obs)
    pred = np.array([B.KICK_START + v0 * t + 0.5 * a * t**2 for t in t_obs])
    t, y, z = E.quadratic_plane(v0[None, :], a[None, :])
    return dict(name="M3 constant accel", n_par=6, v0=v0.tolist(), a=a.tolist(),
                chi2=float(np.sum(((pred - r_obs) / E.SIGMA_POS)**2)),
                t=float(t[0]), y=float(y[0]), z=float(z[0]), pred=pred)
