"""bkick -- forward model for a football (soccer ball) in flight.

University Physics Competition 2026, Problem B.

The ball is treated as a rigid sphere of radius R and mass m moving through
still, dry air of density rho.  Two aerodynamic forces act in addition to
gravity:

  drag     F_D = -0.5 * rho * A * Cd * |v| * v
  Magnus   F_M =  0.5 * rho * A * cL * R * (omega x v)

The Magnus force is written in the linear form F_M = k (omega x v) with
k = 0.5*rho*A*cL*R/m; this is the small-spin-rate limit of
F_M = 0.5*rho*A*CL(S)*|v|^2 * (omega_hat x v_hat) with CL(S) = cL*S and
S = |omega| R / |v| the spin parameter.  It is smooth at omega = 0, which
matters because the estimator has to evaluate the model at zero spin.

Units: SI throughout (m, s, kg, rad).
"""

from __future__ import annotations

import math

import numpy as np

# ----------------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------------

G = 9.81                      # m s^-2, standard gravity (assumption A1)
RHO = 1.20                    # kg m^-3, dry air at 20 C, 1013 hPa (A2)
BALL_MASS = 0.430             # kg, mean of the IFAB band 0.410-0.450 kg
BALL_RADIUS = 0.110           # m, mean circumference 0.69 m -> r = 0.110 m
BALL_AREA = np.pi * BALL_RADIUS**2

GOAL_X = 11.0                 # m, goal plane
GOAL_HALF_WIDTH = 3.66        # m, -3.66 < y < 3.66
GOAL_HEIGHT = 2.44            # m, 0 < z < 2.44
POST_RADIUS = 0.06            # m, IFAB goalpost radius (either 0.06 or 0.06+)
KICK_START = np.array([0.0, 0.0, 0.11])   # m, ball centre at t = 0

# Measurement epochs of the problem
T_OBS = (0.05, 0.10)

# Scenario coordinates (Problem B statement)
SCENARIO = {
    0.05: np.array([1.50, 0.45, 0.44]),
    0.10: np.array([2.96, 0.86, 0.74]),
}


# ----------------------------------------------------------------------------
# Flight model
# ----------------------------------------------------------------------------


class Flight:
    """Drag + Magnus flight model with optional spin decay.

    Parameters
    ----------
    cd : float
        Drag coefficient (dimensionless).
    cl : float
        Magnus slope, CL = cl * S with S = |omega| R / |v|.
    tau : float or None
        Spin-decay time constant in seconds (omega ~ exp(-t/tau)).
        None means constant spin.
    """

    def __init__(self, cd=0.21, cl=1.0, tau=None):
        self.cd = float(cd)
        self.cl = float(cl)
        self.tau = None if tau is None else float(tau)
        self.k_d = 0.5 * RHO * BALL_AREA * self.cd / BALL_MASS      # 1/m
        self.k_m = 0.5 * RHO * BALL_AREA * self.cl * BALL_RADIUS / BALL_MASS

    # -- state layout: [r(3), v(3), omega(3), t] ----------------------------
    def accel(self, v, omega, t=0.0):
        """Acceleration (3,) for velocity v (3,) and spin omega (3,)."""
        speed = np.linalg.norm(v)
        a = -self.k_d * speed * v
        if omega is not None:
            a = a + self.k_m * np.cross(omega, v)
        a[2] -= G
        return a

    def deriv(self, state, t):
        v = state[3:6]
        omega = state[6:9]
        out = np.empty(9)
        out[0:3] = v
        out[3:6] = self.accel(v, omega, t)
        if self.tau is None:
            out[6:9] = 0.0
        else:
            out[6:9] = -omega / self.tau
        return out

    def rk4(self, state, t0, dt):
        """One classical RK4 step of length dt."""
        f = self.deriv
        k1 = f(state, t0)
        k2 = f(state + 0.5 * dt * k1, t0 + 0.5 * dt)
        k3 = f(state + 0.5 * dt * k2, t0 + 0.5 * dt)
        k4 = f(state + dt * k3, t0 + dt)
        return state + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

    def propagate(self, r0, v0, omega0, t_end, dt=0.001):
        """Integrate from t=0 to t_end; returns dense list (t, state)."""
        state = np.concatenate([np.asarray(r0, float), np.asarray(v0, float),
                                np.asarray(omega0, float)])
        n = int(np.ceil(t_end / dt))
        dt = t_end / n
        traj = [(0.0, state.copy())]
        t = 0.0
        for _ in range(n):
            state = self.rk4(state, t, dt)
            t += dt
            traj.append((t, state.copy()))
        return traj

    def position_at(self, r0, v0, omega0, t, dt=0.001):
        """Position at a single time t (<= ~1 s)."""
        traj = self.propagate(r0, v0, omega0, t, dt=dt)
        return traj[-1][1][0:3]

    def goal_plane(self, r0, v0, omega0, dt=0.0005, t_max=2.0):
        """Where and when the ball crosses x = GOAL_X.

        Returns (t_cross, y_cross, z_cross, trajectory).  y/z are None if the
        ball never reaches the plane within t_max.
        """
        r0 = np.asarray(r0, float)
        state = np.concatenate([r0, np.asarray(v0, float),
                                np.asarray(omega0, float)])
        t = 0.0
        traj = [(0.0, state.copy())]
        x_prev = state[0]
        while t < t_max:
            state_new = self.rk4(state, t, dt)
            if state_new[0] >= GOAL_X:
                # linear interpolation of the state where x = GOAL_X
                lam = (GOAL_X - x_prev) / (state_new[0] - x_prev)
                s = state + lam * (state_new - state)
                tc = t + lam * dt
                traj.append((tc, s.copy()))
                return tc, s[1], s[2], traj
            state, x_prev, t = state_new, state_new[0], t + dt
            traj.append((t, state.copy()))
        return None, None, None, traj


# ----------------------------------------------------------------------------
# Analytic limits (used for code verification, see validate.py)
# ----------------------------------------------------------------------------


def ballistic_plane(v0, r0=KICK_START, g=G):
    """Exact gravity-only crossing of x = GOAL_X.

    Solves 0.5 g t^2 - v0z t - (r0z - z_target) = 0 in the vertical, with
    linear motion in x and y.  Returns (t, y, z).
    """
    v0 = np.asarray(v0, float)
    t = (GOAL_X - r0[0]) / v0[0]
    return t, r0[1] + v0[1] * t, r0[2] + v0[2] * t - 0.5 * g * t * t


def constant_accel_inverse(r_obs, t_obs, r0=KICK_START):
    """Exact inversion of r(t) = r0 + v0 t + 0.5 a t^2 from two fixes.

    Returns (v0, a).  Requires equally spaced epochs t_obs = (t1, t2) with
    t2 = 2 t1 (true for 0.05/0.10 s).  Derivation in the appendix.
    """
    r1, r2 = np.asarray(r_obs[0], float), np.asarray(r_obs[1], float)
    r0 = np.asarray(r0, float)
    t1, t2 = float(t_obs[0]), float(t_obs[1])
    # r(t1)-r0 = v0 t1 + a t1^2/2 ; r(t2)-r0 = v0 t2 + a t2^2/2
    # eliminate v0 by t2*(first) - t1*(second):
    a = 2.0 * ((r1 - r0) * t2 - (r2 - r0) * t1) / (t1 * t2 * (t1 - t2))
    v0 = (r1 - r0) / t1 - 0.5 * a * t1
    return v0, a


def constant_accel_plane(v0, a, r0=KICK_START):
    """Crossing of x = GOAL_X under constant acceleration a (both roots)."""
    A, B, C = 0.5 * a[0], v0[0], r0[0] - GOAL_X
    roots = np.roots([A, B, C])
    out = []
    for rt in roots:
        if abs(rt.imag) < 1e-12 and rt.real > 0:
            t = float(rt.real)
            r = r0 + v0 * t + 0.5 * a * t * t
            out.append((t, r[1], r[2]))
    out.sort(key=lambda x: x[0])
    return out


# ----------------------------------------------------------------------------
# Vectorised integrator with per-sample coefficients
# ----------------------------------------------------------------------------
#
# The Monte-Carlo propagation needs N trajectories whose drag and Magnus
# coefficients differ from sample to sample (they are the sampled parameters).
# kd and km are therefore arrays of shape (N, 1).  The state is (N, 9).


def _accel_batch(v, omega, t, kd, km, tau):
    speed = np.linalg.norm(v, axis=1)[:, None]
    a = -kd * speed * v
    a = a + km * np.cross(omega, v)
    a[:, 2] -= G
    return a


def _deriv_batch(state, t, kd, km, tau):
    out = np.empty_like(state)
    out[:, 0:3] = state[:, 3:6]
    out[:, 3:6] = _accel_batch(state[:, 3:6], state[:, 6:9], t, kd, km, tau)
    if tau is None:
        out[:, 6:9] = 0.0
    else:
        out[:, 6:9] = -state[:, 6:9] / tau
    return out


def _rk4_batch(state, t, dt, kd, km, tau):
    k1 = _deriv_batch(state, t, kd, km, tau)
    k2 = _deriv_batch(state + 0.5 * dt * k1, t + 0.5 * dt, kd, km, tau)
    k3 = _deriv_batch(state + 0.5 * dt * k2, t + 0.5 * dt, kd, km, tau)
    k4 = _deriv_batch(state + dt * k3, t + dt, kd, km, tau)
    return state + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)


def batch_positions(v0, omega0, kd, km, times, r0=KICK_START, dt=0.002, tau=None):
    """Positions (N, len(times), 3) at the requested times."""
    v0 = np.atleast_2d(v0)
    N = v0.shape[0]
    state = np.empty((N, 9))
    state[:, 0:3] = np.broadcast_to(np.asarray(r0, float), (N, 3))
    state[:, 3:6] = v0
    state[:, 6:9] = np.broadcast_to(np.atleast_2d(omega0), (N, 3))
    out = np.empty((N, len(times), 3))
    t = 0.0
    for k, tk in enumerate(times):
        n = int(round((tk - t) / dt))
        for _ in range(n):
            state = _rk4_batch(state, t, dt, kd, km, tau)
            t += dt
        out[:, k, :] = state[:, 0:3]
    return out


def batch_goal_plane(v0, omega0, kd, km, r0=KICK_START, dt=0.002, t_max=1.5, tau=None):
    """Goal-plane crossing for N trajectories with individual coefficients.

    Returns (t, y, z), each of length N; NaN where no crossing happened.
    """
    v0 = np.atleast_2d(v0)
    N = v0.shape[0]
    state = np.empty((N, 9))
    state[:, 0:3] = np.broadcast_to(np.asarray(r0, float), (N, 3))
    state[:, 3:6] = v0
    state[:, 6:9] = np.broadcast_to(np.atleast_2d(omega0), (N, 3))
    out = np.full((N, 3), np.nan)
    done = np.zeros(N, bool)
    t = 0.0
    while t < t_max and not done.all():
        nxt = _rk4_batch(state, t, dt, kd, km, tau)
        hit = (~done) & (nxt[:, 0] >= GOAL_X)
        if hit.any():
            lam = (GOAL_X - state[hit, 0]) / (nxt[hit, 0] - state[hit, 0])
            s = state[hit] + lam[:, None] * (nxt[hit] - state[hit])
            out[hit, 0] = t + lam * dt
            out[hit, 1] = s[:, 1]
            out[hit, 2] = s[:, 2]
            done |= hit
        state = nxt
        t += dt
    return out[:, 0], out[:, 1], out[:, 2]


def coeffs(cd, cl=1.0, n=None):
    """Return (kd, km) arrays of shape (n,1) for arrays of cd (and cl)."""
    cd = np.atleast_1d(np.asarray(cd, float)).reshape(-1, 1)
    cl = np.atleast_1d(np.asarray(cl, float)).reshape(-1, 1)
    if n is not None:
        cd = np.broadcast_to(cd, (n, 1))
        cl = np.broadcast_to(cl, (n, 1))
    kd = 0.5 * RHO * BALL_AREA * cd / BALL_MASS
    km = 0.5 * RHO * BALL_AREA * cl * BALL_RADIUS / BALL_MASS
    return kd, km


# ----------------------------------------------------------------------------
# Fast scalar integrator (pure Python floats)
# ----------------------------------------------------------------------------
#
# The Bayesian sampler evaluates the forward model ~10^5 times, so the numpy
# version (which pays ~10 us of array overhead per call) is too slow.  These
# routines keep the same RK4 scheme but store the nine state variables as
# plain floats.  Agreement with the numpy integrator is checked in
# validate.py (test T1).


def _accel_f(vx, vy, vz, wx, wy, wz, kd, km):
    s = math.sqrt(vx * vx + vy * vy + vz * vz)
    ax = -kd * s * vx + km * (wy * vz - wz * vy)
    ay = -kd * s * vy + km * (wz * vx - wx * vz)
    az = -kd * s * vz + km * (wx * vy - wy * vx) - G
    return ax, ay, az


def fast_positions(v0, omega, cd, cl, times, dt=0.002, tau=None,
                   r0=KICK_START):
    """Positions at the requested (sorted) times, pure-Python RK4."""
    kd = 0.5 * RHO * BALL_AREA * cd / BALL_MASS
    km = 0.5 * RHO * BALL_AREA * cl * BALL_RADIUS / BALL_MASS
    x, y, z = float(r0[0]), float(r0[1]), float(r0[2])
    vx, vy, vz = float(v0[0]), float(v0[1]), float(v0[2])
    wx, wy, wz = float(omega[0]), float(omega[1]), float(omega[2])
    out = []
    t = 0.0
    for tk in times:
        n = int(round((tk - t) / dt))
        h = dt
        for _ in range(n):
            # spin decay: omega(t) = omega * exp(-t/tau); the three RK4
            # stages sit at t, t+h/2 (twice) and t+h
            if tau is None:
                f1 = f2 = f4 = 1.0
            else:
                f1 = math.exp(-t / tau)
                f2 = math.exp(-(t + 0.5 * h) / tau)
                f4 = math.exp(-(t + h) / tau)
            # stage 1
            ax1, ay1, az1 = _accel_f(vx, vy, vz, wx * f1, wy * f1, wz * f1, kd, km)
            # stage 2
            hvx, hvy, hvz = vx + 0.5 * h * ax1, vy + 0.5 * h * ay1, vz + 0.5 * h * az1
            ax2, ay2, az2 = _accel_f(hvx, hvy, hvz, wx * f2, wy * f2, wz * f2, kd, km)
            # stage 3
            hvx, hvy, hvz = vx + 0.5 * h * ax2, vy + 0.5 * h * ay2, vz + 0.5 * h * az2
            ax3, ay3, az3 = _accel_f(hvx, hvy, hvz, wx * f2, wy * f2, wz * f2, kd, km)
            # stage 4
            hvx, hvy, hvz = vx + h * ax3, vy + h * ay3, vz + h * az3
            ax4, ay4, az4 = _accel_f(hvx, hvy, hvz, wx * f4, wy * f4, wz * f4, kd, km)
            x += h * (vx + h * (ax1 + ax2 + ax3) / 6.0)
            y += h * (vy + h * (ay1 + ay2 + ay3) / 6.0)
            z += h * (vz + h * (az1 + az2 + az3) / 6.0)
            vx += h * (ax1 + 2.0 * ax2 + 2.0 * ax3 + ax4) / 6.0
            vy += h * (ay1 + 2.0 * ay2 + 2.0 * ay3 + ay4) / 6.0
            vz += h * (az1 + 2.0 * az2 + 2.0 * az3 + az4) / 6.0
            t += h
        out.append(np.array([x, y, z]))
    return out


def fast_goal_plane(v0, omega, cd, cl, dt=0.002, t_max=1.5, tau=None,
                    r0=KICK_START, x_target=GOAL_X):
    """Goal-plane crossing (t, y, z) for one trajectory, pure-Python RK4."""
    kd = 0.5 * RHO * BALL_AREA * cd / BALL_MASS
    km = 0.5 * RHO * BALL_AREA * cl * BALL_RADIUS / BALL_MASS
    x, y, z = float(r0[0]), float(r0[1]), float(r0[2])
    vx, vy, vz = float(v0[0]), float(v0[1]), float(v0[2])
    wx, wy, wz = float(omega[0]), float(omega[1]), float(omega[2])
    t = 0.0
    h = dt
    while t < t_max:
        if tau is None:
            f1 = f2 = f4 = 1.0
        else:
            f1 = math.exp(-t / tau)
            f2 = math.exp(-(t + 0.5 * h) / tau)
            f4 = math.exp(-(t + h) / tau)
        ax1, ay1, az1 = _accel_f(vx, vy, vz, wx * f1, wy * f1, wz * f1, kd, km)
        hvx, hvy, hvz = vx + 0.5 * h * ax1, vy + 0.5 * h * ay1, vz + 0.5 * h * az1
        ax2, ay2, az2 = _accel_f(hvx, hvy, hvz, wx * f2, wy * f2, wz * f2, kd, km)
        hvx, hvy, hvz = vx + 0.5 * h * ax2, vy + 0.5 * h * ay2, vz + 0.5 * h * az2
        ax3, ay3, az3 = _accel_f(hvx, hvy, hvz, wx * f2, wy * f2, wz * f2, kd, km)
        hvx, hvy, hvz = vx + h * ax3, vy + h * ay3, vz + h * az3
        ax4, ay4, az4 = _accel_f(hvx, hvy, hvz, wx * f4, wy * f4, wz * f4, kd, km)
        x_new = x + h * (vx + h * (ax1 + ax2 + ax3) / 6.0)
        y_new = y + h * (vy + h * (ay1 + ay2 + ay3) / 6.0)
        z_new = z + h * (vz + h * (az1 + az2 + az3) / 6.0)
        vx_new = vx + h * (ax1 + 2.0 * ax2 + 2.0 * ax3 + ax4) / 6.0
        vy_new = vy + h * (ay1 + 2.0 * ay2 + 2.0 * ay3 + ay4) / 6.0
        vz_new = vz + h * (az1 + 2.0 * az2 + 2.0 * az3 + az4) / 6.0
        if x_new >= x_target:
            lam = (x_target - x) / (x_new - x)
            return (t + lam * h, y + lam * (y_new - y), z + lam * (z_new - z))
        x, y, z, vx, vy, vz = x_new, y_new, z_new, vx_new, vy_new, vz_new
        t += h
    return (np.nan, np.nan, np.nan)
