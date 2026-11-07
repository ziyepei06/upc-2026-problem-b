"""estimators.py -- the two estimators used in the paper.

E1  kinematically exact (model-free) estimator
    r(t) = r0 + v0 t + a t^2/2 is inverted algebraically from the origin and
    the two fixes; the goal-plane crossing follows from the quadratic.  Its
    only assumption is that the acceleration is constant over the 0.10 s
    spanned by the data.

E2  physics-constrained Bayesian estimator
    Unknowns theta = (v0x, v0y, v0z, ln Cd, omega_x, omega_y, omega_z).
    Priors from wind-tunnel and tracking literature.  Posterior sampled with
    random-walk Metropolis started at the Levenberg-Marquardt MAP; every
    sample is propagated to x = 11 m, which yields the predictive
    distribution of the entry point directly.

Parameter defaults follow the evidence brief (notes/evidence-brief.md):
Cd = 0.19 +- 0.03 at 25-35 m/s (Carre 2005; Goff & Carre 2012),
spin decay time constant 6 s (Tsukada 2008: 5-15 % per 18 m),
Magnus slope CL = cL * S with cL = 1.0 (Asai 2020: CL = 0.21-0.29 at
S = 0.18-0.34), typical spin 2-10 rev/s.
"""

from __future__ import annotations

import numpy as np

import bkick as B

SIGMA_POS = 0.010          # m, 1-sigma position error of a fix (FIFA GLT 1.5 cm)
CD_MEDIAN = 0.19           # drag coefficient at 25-35 m/s
CD_LOG_SIGMA = 0.15        # log-normal shape -> ~0.16-0.22 at 1 sigma
OMEGA_SIGMA = 35.0         # rad/s per component (~5.6 rev/s)
CL_SLOPE = 1.0             # Magnus slope, CL = cL * S
SPIN_TAU = 6.0             # s
DT_FIT = 0.002


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _inside(y, z):
    """Ball entirely inside the frame (centre at least one radius from a post)."""
    return (np.abs(y) < B.GOAL_HALF_WIDTH - B.BALL_RADIUS) & \
           (z > B.BALL_RADIUS) & (z < B.GOAL_HEIGHT - B.BALL_RADIUS)


def _inside_center(y, z):
    return (np.abs(y) < B.GOAL_HALF_WIDTH) & (z > 0.0) & (z < B.GOAL_HEIGHT)


def quadratic_plane(v0, a, x_target=B.GOAL_X, r0z=B.KICK_START[2]):
    """Vectorised crossing of x = x_target for r(t) = r0 + v0 t + a t^2/2."""
    A = 0.5 * a[:, 0]
    Bq = v0[:, 0]
    C = -x_target
    disc = np.maximum(Bq**2 - 4.0 * A * C, 0.0)
    sq = np.sqrt(disc)
    t1 = (-Bq + sq) / (2.0 * A)
    t2 = (-Bq - sq) / (2.0 * A)
    with np.errstate(invalid="ignore"):
        t = np.where((t1 > 0) & ((t1 < t2) | (t2 <= 0)), t1, t2)
    y = v0[:, 1] * t + 0.5 * a[:, 1] * t**2
    z = r0z + v0[:, 2] * t + 0.5 * a[:, 2] * t**2
    return t, y, z


# ---------------------------------------------------------------------------
# E1
# ---------------------------------------------------------------------------


def e1_predict(r_obs, t_obs=(0.05, 0.10), sigma=SIGMA_POS, n_mc=20000, seed=0):
    r_obs = np.asarray(r_obs, float)
    v0, a = B.constant_accel_inverse(r_obs, t_obs)
    t0, y0, z0 = quadratic_plane(v0[None, :], a[None, :])
    t0, y0, z0 = float(t0[0]), float(y0[0]), float(z0[0])

    # --- Monte-Carlo over the two fixes ---
    rng = np.random.default_rng(seed)
    r1 = r_obs[0] + rng.normal(0.0, sigma, (n_mc, 3))
    r2 = r_obs[1] + rng.normal(0.0, sigma, (n_mc, 3))
    t1, t2 = t_obs
    a_s = 2.0 * ((r1 - B.KICK_START) * t2 - (r2 - B.KICK_START) * t1) / (t1 * t2 * (t1 - t2))
    v_s = (r1 - B.KICK_START) / t1 - 0.5 * a_s * t1
    t_s, y_s, z_s = quadratic_plane(v_s, a_s)

    # --- delta method ---
    J = np.zeros((3, 6))
    base = np.array([t0, y0, z0])
    for i in range(6):
        p = np.concatenate([v0, a])
        h = 1e-6 * max(1.0, abs(p[i]))
        p[i] += h
        tq, yq, zq = quadratic_plane(p[None, :3], p[None, 3:])
        J[:, i] = (np.array([tq[0], yq[0], zq[0]]) - base) / h
    G = e1_data_jacobian(t_obs)
    cov_pred = J @ G @ (sigma**2 * np.eye(6)) @ G.T @ J.T

    return dict(
        v0=v0.tolist(), a=a.tolist(), t=t0, y=y0, z=z0,
        mean=dict(t=float(t_s.mean()), y=float(y_s.mean()), z=float(z_s.mean())),
        std=dict(t=float(t_s.std(ddof=1)), y=float(y_s.std(ddof=1)),
                 z=float(z_s.std(ddof=1))),
        std_analytic=dict(t=float(np.sqrt(abs(cov_pred[0, 0]))),
                          y=float(np.sqrt(abs(cov_pred[1, 1]))),
                          z=float(np.sqrt(abs(cov_pred[2, 2])))),
        cov_yz=[[float(cov_pred[1, 1]), float(cov_pred[1, 2])],
                [float(cov_pred[2, 1]), float(cov_pred[2, 2])]],
        pct=dict(t=[float(np.percentile(t_s, q)) for q in (5, 16, 50, 84, 95)],
                 y=[float(np.percentile(y_s, q)) for q in (5, 16, 50, 84, 95)],
                 z=[float(np.percentile(z_s, q)) for q in (5, 16, 50, 84, 95)]),
        p_inside=float(np.mean(_inside(y_s, z_s))),
        p_center_inside=float(np.mean(_inside_center(y_s, z_s))),
        n_mc=n_mc, sigma=sigma,
        _t=t_s, _y=y_s, _z=z_s,
    )


def e1_data_jacobian(t_obs):
    """d(v0, a)/d(r1, r2), 6x6."""
    t1, t2 = t_obs
    ia = 2.0 / (t1 * t2 * (t1 - t2))
    Ja = ia * np.hstack([t2 * np.eye(3), -t1 * np.eye(3)])
    J = np.zeros((6, 6))
    J[3:, :] = Ja
    J[:3, :] = np.hstack([np.eye(3) / t1, np.zeros((3, 3))]) - 0.5 * t1 * Ja
    return J


# ---------------------------------------------------------------------------
# E2
# ---------------------------------------------------------------------------


class PhysicalPosterior:
    """Posterior of theta = (v0, ln Cd, omega) given two noisy fixes."""

    def __init__(self, r_obs, t_obs=(0.05, 0.10), sigma=SIGMA_POS,
                 cd=(CD_MEDIAN, CD_LOG_SIGMA), omega_sigma=OMEGA_SIGMA,
                 cl=CL_SLOPE, tau=SPIN_TAU, dt=DT_FIT):
        self.r_obs = np.asarray(r_obs, float)
        self.t_obs = tuple(t_obs)
        self.sigma = sigma
        self.cd_median, self.cd_log_sigma = cd
        self.omega_sigma = omega_sigma
        self.cl = cl
        self.tau = tau
        self.dt = dt
        self.ln_cd_mean = np.log(cd[0])

    # -- forward model -----------------------------------------------------
    def positions(self, theta):
        theta = np.atleast_2d(theta)
        if theta.shape[0] == 1:
            # fast path: pure-Python scalar RK4 (~0.1 ms vs ~19 ms for the
            # vectorised version, which matters for MCMC)
            return np.array([B.fast_positions(
                theta[0, 0:3], theta[0, 4:7], float(np.exp(theta[0, 3])),
                self.cl, self.t_obs, dt=self.dt, tau=self.tau)])
        cd = np.exp(theta[:, 3])
        kd, km = B.coeffs(cd, self.cl, n=theta.shape[0])
        return B.batch_positions(theta[:, 0:3], theta[:, 4:7], kd, km,
                                 self.t_obs, dt=self.dt, tau=self.tau)

    def log_prior(self, theta):
        theta = np.atleast_2d(theta)
        lp = -0.5 * ((theta[:, 3] - self.ln_cd_mean) / self.cd_log_sigma) ** 2
        lp = lp - 0.5 * np.sum(theta[:, 4:7] ** 2, axis=1) / self.omega_sigma**2
        return lp

    def log_like(self, theta):
        theta = np.atleast_2d(theta)
        pred = self.positions(theta)
        res = (pred - self.r_obs[None, :, :]).reshape(theta.shape[0], -1)
        return -0.5 * np.sum(res**2, axis=1) / self.sigma**2

    def log_post(self, theta):
        return self.log_like(theta) + self.log_prior(theta)

    def residuals(self, theta):
        return (self.positions(theta)[0] - self.r_obs).ravel()

    # -- MAP ---------------------------------------------------------------
    def prior_residuals(self, theta):
        theta = np.atleast_2d(theta)
        rp = np.empty((theta.shape[0], 4))
        rp[:, 0] = (theta[:, 3] - self.ln_cd_mean) / self.cd_log_sigma
        rp[:, 1:4] = theta[:, 4:7] / self.omega_sigma
        return rp

    def all_residuals(self, theta):
        """Joint residual vector (model misfit + prior pulls), length 10."""
        theta = np.atleast_2d(theta)
        pred = self.positions(theta)
        rm = ((pred - self.r_obs[None, :, :]).reshape(theta.shape[0], -1)
              / self.sigma)
        return np.hstack([rm, self.prior_residuals(theta)])

    @staticmethod
    def default_start(r_obs, t_obs=(0.05, 0.10)):
        v0, _ = B.constant_accel_inverse(r_obs, t_obs)
        return np.array([v0[0], v0[1], v0[2], np.log(CD_MEDIAN), 0.0, 0.0, -20.0])

    def map(self, theta0=None, iters=80, tol=1e-10):
        """Levenberg-Marquardt on the joint residual vector."""
        if theta0 is None:
            theta0 = self.default_start(self.r_obs, self.t_obs)
        theta = np.asarray(theta0, float).copy()
        n = 7
        R = self.all_residuals(theta)[0]
        cost = 0.5 * R @ R
        mu = 1e-3
        for _ in range(iters):
            J = np.zeros((R.size, n))
            for i in range(n):
                h = 1e-6 * max(1.0, abs(theta[i]))
                tp = theta.copy()
                tp[i] += h
                J[:, i] = (self.all_residuals(tp)[0] - R) / h
            A = J.T @ J
            g = J.T @ R
            try:
                step = np.linalg.solve(A + mu * np.diag(np.diag(A) + 1e-12), -g)
            except np.linalg.LinAlgError:
                break
            trial = theta + step
            Rt = self.all_residuals(trial)[0]
            ct = 0.5 * Rt @ Rt
            if ct < cost:
                done = (cost - ct) < tol * max(1.0, cost)
                theta, R, cost = trial, Rt, ct
                mu = max(mu * 0.3, 1e-12)
                if done:
                    break
            else:
                mu *= 10.0
                if mu > 1e10:
                    break
        self._J = J
        return theta, -cost

    def laplace_cov(self, theta):
        """Gauss-Newton (Laplace) covariance from the last Jacobian."""
        return np.linalg.inv(self._J.T @ self._J + 1e-12 * np.eye(7))

    # -- MCMC --------------------------------------------------------------
    def sample(self, n=30000, burn=5000, thin=3, seed=1, cov=None, theta0=None):
        theta_map, _ = self.map(theta0)
        if cov is None:
            cov = self.laplace_cov(theta_map)
        rng = np.random.default_rng(seed)
        prop_cov = (2.4**2 / 7.0) * cov
        keep = n // thin
        chain = np.empty((keep, 7))
        logp = np.empty(keep)
        cur = theta_map.copy()
        cur_lp = float(self.log_post(cur[None, :])[0])
        acc = 0
        k = 0
        for it in range(n + burn):
            cand = rng.multivariate_normal(cur, prop_cov)
            cand_lp = float(self.log_post(cand[None, :])[0])
            if np.log(rng.random()) < cand_lp - cur_lp:
                cur, cur_lp = cand, cand_lp
                acc += 1
            if it >= burn and (it - burn) % thin == 0 and k < keep:
                chain[k] = cur
                logp[k] = cur_lp
                k += 1
        return dict(chain=chain[:k], logp=logp[:k], acc=acc / (n + burn),
                    theta_map=theta_map, cov=cov)

    # -- prediction --------------------------------------------------------
    def predict(self, chain, dt=0.002):
        cd = np.exp(chain[:, 3])
        kd, km = B.coeffs(cd, self.cl, n=chain.shape[0])
        return B.batch_goal_plane(chain[:, 0:3], chain[:, 4:7], kd, km,
                                  dt=dt, tau=self.tau)


def summarize(t, y, z, w=None):
    ok = np.isfinite(t) & np.isfinite(y) & np.isfinite(z)
    t, y, z = t[ok], y[ok], z[ok]
    w = None if w is None else (w[ok] / w[ok].sum())

    def mean(v):
        return float(np.average(v, weights=w)) if w is not None else float(v.mean())

    def std(v):
        m = mean(v)
        return (float(np.sqrt(np.average((v - m) ** 2, weights=w)))
                if w is not None else float(v.std(ddof=1)))

    def q(v, p):
        if w is None:
            return float(np.percentile(v, p))
        o = np.argsort(v)
        return float(np.interp(p / 100.0, np.cumsum(w[o]), v[o]))

    def pin(fn):
        v = fn(y, z)
        return float(np.average(v, weights=w)) if w is not None else float(v.mean())

    return dict(
        mean=dict(t=mean(t), y=mean(y), z=mean(z)),
        std=dict(t=std(t), y=std(y), z=std(z)),
        median=dict(t=q(t, 50), y=q(y, 50), z=q(z, 50)),
        pct=dict(t=[q(t, p) for p in (5, 16, 50, 84, 95)],
                 y=[q(y, p) for p in (5, 16, 50, 84, 95)],
                 z=[q(z, p) for p in (5, 16, 50, 84, 95)]),
        p_inside=pin(_inside), p_center_inside=pin(_inside_center),
        frac_never=float(1.0 - ok.mean()), n=int(ok.sum()),
    )


def e2_predict_quick(r_obs, t_obs=(0.05, 0.10), sigma=SIGMA_POS, seed=0,
                     n=3000, burn=800, thin=1, **kw):
    """MAP + short MCMC + predictive summary (used by the validation runs)."""
    P = PhysicalPosterior(r_obs, t_obs, sigma=sigma, **kw)
    out = P.sample(n=n, burn=burn, thin=thin, seed=seed)
    t, y, z = P.predict(out["chain"])
    s = summarize(t, y, z)
    s["acc"] = out["acc"]
    s["theta_map"] = out["theta_map"].tolist()
    s["chain"] = out["chain"]
    return s
