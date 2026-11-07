"""Exploration: kinematic inversion vs drag+Magnus MAP fit for the scenario."""
import numpy as np
import bkick as B

r0 = B.KICK_START
r1 = B.SCENARIO[0.05]
r2 = B.SCENARIO[0.10]
t1, t2 = B.T_OBS

v0k, ak = B.constant_accel_inverse((r1, r2), (t1, t2))
print("constant-accel inversion: v0 =", v0k, " a =", ak)
print("ballistic check: ", B.constant_accel_plane(v0k, ak))
print("ballistic (gravity only) crossing:", B.ballistic_plane(v0k))

# ---- drag + Magnus model, fit to the three fixes -------------------------
CL = 1.0
flight = B.Flight(cd=0.21, cl=CL)


def model_positions(theta, dt=0.001):
    v0 = theta[0:3]
    cd = theta[3]
    om = theta[4:7]
    fl = B.Flight(cd=cd, cl=CL)
    traj = fl.propagate(r0, v0, om, t2, dt=dt)
    # pick positions at t1 and t2 (dt divides both epochs: 0.001)
    n1 = int(round(t1 / dt))
    n2 = int(round(t2 / dt))
    return np.array([traj[n1][1][0:3], traj[n2][1][0:3]])


def residuals(theta, sigma=0.01):
    pred = model_positions(theta)
    obs = np.vstack([r1, r2])
    return (pred - obs).ravel() / sigma


theta0 = np.concatenate([v0k, [0.21], [0.0, 0.0, -60.0]])

# quick check: what acceleration does the model produce at t=0.075?
fl = B.Flight(cd=0.21, cl=CL)
traj = fl.propagate(r0, v0k, theta0[4:7], 0.10, dt=0.0005)
print("model positions 0.05/0.10:", traj[100][1][0:3], traj[200][1][0:3])
print("target:", r1, r2)

# ---- Levenberg-Marquardt with priors -------------------------------------
priors = {
    "cd": (0.21, 0.06),
    "omega": (0.0, 30.0),   # per component, rad/s
}


def obj(theta, sigma=0.01):
    res = residuals(theta, sigma)
    val = 0.5 * res @ res
    cd = theta[3]
    val += 0.5 * ((cd - priors["cd"][0]) / priors["cd"][1]) ** 2
    for j in range(3):
        val += 0.5 * ((theta[4 + j] - 0.0) / priors["omega"][1]) ** 2
    return val


def num_jac(f, theta, h=1e-5):
    f0 = f(theta)
    J = np.zeros((f0.size, theta.size))
    for i in range(theta.size):
        th = theta.copy()
        th[i] += h * max(1.0, abs(theta[i]))
        J[:, i] = (f(th) - f0) / (h * max(1.0, abs(theta[i])))
    return J


def lm(theta, iters=40, mu=1e-3):
    for it in range(iters):
        r = residuals(theta)
        J = num_jac(lambda th: residuals(th), theta)
        # prior block
        rp = np.array([(theta[3] - priors["cd"][0]) / priors["cd"][1]] +
                      [(theta[4 + j]) / priors["omega"][1] for j in range(3)])
        Jp = np.zeros((4, theta.size))
        Jp[0, 3] = 1.0 / priors["cd"][1]
        for j in range(3):
            Jp[1 + j, 4 + j] = 1.0 / priors["omega"][1]
        r = np.concatenate([r, rp])
        J = np.vstack([J, Jp])
        A = J.T @ J
        g = J.T @ r
        step = np.linalg.solve(A + mu * np.diag(np.diag(A) + 1e-9), -g)
        trial = theta + step
        if obj(trial) < obj(theta):
            theta = trial
            mu = max(mu * 0.3, 1e-12)
        else:
            mu *= 10.0
    return theta


theta = lm(theta0.copy())
print("\nMAP theta: v0=", theta[0:3], " Cd=", theta[3], " omega=", theta[4:7],
      " |omega| =", np.linalg.norm(theta[4:7]), "rev/s", np.linalg.norm(theta[4:7]) / (2 * np.pi))
print("model positions:", model_positions(theta))
print("obs:", r1, r2)

fl = B.Flight(cd=theta[3], cl=CL)
tc, yc, zc, traj = fl.goal_plane(r0, theta[0:3], theta[4:7])
print("predicted goal crossing: t=%.4f s  y=%.3f  z=%.3f" % (tc, yc, zc))
