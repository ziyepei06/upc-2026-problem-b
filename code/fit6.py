"""Six-parameter transverse-force fit, our own least-squares implementation.

Model (GPT's M6_perp parameterisation, re-implemented independently):
  u = v (no wind),  a = -FACTOR*C_D*|u|*u + |u|*(beta x u) - g*zhat
  beta = b1*e1 + b2*e2, basis perpendicular to the launch velocity.
Fit: penalised (regularised) least squares on the six observed coordinates,
sigma_pos = 5 mm, solved by damped Gauss-Newton with our own integrator.
Posterior: Laplace covariance x 1.15, Student-t df=8 proposal, exact
importance weights, propagation to x = 11 m by RK4 + cubic Hermite crossing.
"""
import numpy as np, json, sys

G, R, MASS, RHO = 9.81, 0.110, 0.430, 1.20
A = np.pi*R*R
FACTOR = RHO*A/(2*MASS)          # 0.053042 1/m
R0 = np.array([0.0, 0.0, 0.11])
OBS = np.array([[1.50, 0.45, 0.44], [2.96, 0.86, 0.74]])
TOBS = np.array([0.05, 0.10])
GOAL_X, HALF_B, H = 11.0, 3.66, 2.44
SIGMA = 0.005
MU = np.array([30.0, 0.0, 7.0, 0.25, 0.0, 0.0])
SD = np.array([15.0, 15.0, 10.0, 0.15, 0.008, 0.008])
LOW = np.array([5.0, -30.0, 0.0, 0.0, -0.08, -0.08])
HIGH = np.array([60.0, 30.0, 30.0, 1.0, 0.08, 0.08])


def basis(v0):
    e0 = v0/np.linalg.norm(v0)
    e1 = np.cross(np.array([0.0, 0.0, 1.0]), e0)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(e0, e1)
    return e0, e1, e2


def beta_vec(th):
    e0, e1, e2 = basis(th[:3])
    return th[4]*e1 + th[5]*e2


def accel(p, v, cd, beta):
    sp = np.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
    a = -FACTOR*cd*sp*v + sp*np.cross(beta, v)
    a[2] -= G
    return a


def rk4_step(p, v, h, cd, beta):
    k1 = accel(p, v, cd, beta)
    q = v + 0.5*h*k1
    k2 = accel(p + 0.5*h*v, q, cd, beta)
    q2 = v + 0.5*h*k2
    k3 = accel(p + 0.5*h*q, q2, cd, beta)
    q3 = v + h*k3
    k4 = accel(p + h*q2, q3, cd, beta)
    return (p + h/6.0*(v + 2*q + 2*q2 + q3),
            v + h/6.0*(k1 + 2*k2 + 2*k3 + k4))


def positions(th, times=(0.05, 0.10), dt=0.005):
    cd = th[3]; beta = beta_vec(th)
    p = R0.copy(); v = th[:3].copy()
    n = int(round(times[-1]/dt)); out = []
    want = {int(round(t/dt)): i for i, t in enumerate(times)}
    for i in range(1, n+1):
        p, v = rk4_step(p, v, dt, cd, beta)
        if i in want:
            out.append((want[i], p.copy()))
    res = [None]*len(times)
    for i, q in out:
        res[i] = q
    return np.array(res)


def residuals(th):
    r = np.concatenate([(positions(th) - OBS).ravel()/SIGMA,
                        (th - MU)/SD])
    return r


def fit_map(th0=None, iters=200):
    th = (MU + np.array([0.4, 9.3, -0.1, 0.10, 0, 0])) if th0 is None else np.array(th0, float)
    th = np.clip(th, LOW+1e-9, HIGH-1e-9)
    lam = 1e-3
    r = residuals(th); cost = float(r@r)
    for _ in range(iters):
        J = np.zeros((len(r), 6))
        for j in range(6):
            h = 1e-6*max(1.0, abs(th[j]))
            tp = th.copy(); tp[j] += h
            J[:, j] = (residuals(tp) - r)/h
        A_ = J.T@J; g = J.T@r
        try:
            step = np.linalg.solve(A_ + lam*np.diag(np.diag(A_) + 1e-12), -g)
        except np.linalg.LinAlgError:
            break
        tn = th + step
        tn = np.clip(tn, LOW+1e-9, HIGH-1e-9)
        rn = residuals(tn); cn = float(rn@rn)
        if cn < cost:
            th, r, cost = tn, rn, cn
            lam = max(lam*0.3, 1e-12)
        else:
            lam *= 10
            if lam > 1e8:
                break
    return th, cost, J, r


def cross_plane(th, dt=0.0025):
    cd = th[3]; beta = beta_vec(th)
    p = R0.copy(); v = th[:3].copy()
    t = 0.0
    while p[0] < GOAL_X and t < 1.5:
        pn, vn = rk4_step(p, v, dt, cd, beta)
        if pn[0] >= GOAL_X:
            # cubic Hermite in x
            f0, f1 = p[0]-GOAL_X, pn[0]-GOAL_X
            s = f0/(f0-f1)
            h00 = 2*s**3-3*s**2+1; h10 = s**3-2*s**2+s
            h01 = -2*s**3+3*s**2; h11 = s**3-s**2
            def her(c0, c1, v0_, v1_):
                return h00*c0 + h10*dt*v0_ + h01*c1 + h11*dt*v1_
            y = her(p[1], pn[1], v[1], vn[1]); z = her(p[2], pn[2], v[2], vn[2])
            return t + s*dt, y, z
        p, v, t = pn, vn, t+dt
    return np.nan, np.nan, np.nan


def clear_full(y, z):
    # disc criterion: every point of the sphere crosses the aperture
    if not np.isfinite(y):
        return False
    if abs(y) > HALF_B - R:
        return False
    if z > H - R or z < R:
        return False
    return True


def t_logpdf(z, df):
    from math import lgamma, pi, log
    d = len(z)
    return (lgamma((df+d)/2)-lgamma(df/2)-0.5*d*np.log(df*pi)
            - 0.5*np.linalg.slogdet(np.eye(d))[1]
            -(df+d)/2*np.log1p((z@z)/df))


def run(n=131072, seed=950001, df=8):
    th, cost, J, r = fit_map()
    cov = np.linalg.inv(J.T@J)*1.15
    print('MAP  v0=(%.4f, %.4f, %.4f) C_D=%.4f b=(%.5f, %.5f)' % (*th[:3], th[3], th[4], th[5]))
    print('rms data residual = %.4f mm   chi2 = %.2f' % (1000*np.sqrt(np.mean(r[:6]**2))*SIGMA, r@r))
    t_map, y_map, z_map = cross_plane(th)
    print('MAP entry  t=%.4f  y=%.4f  z=%.4f' % (t_map, y_map, z_map))
    rng = np.random.default_rng(seed)
    L = np.linalg.cholesky(cov)
    Z = rng.standard_normal((n, 6))
    chi = np.sqrt(rng.chisquare(df, n)/df)
    ths = th + (Z@L.T)/chi[:, None]
    ok = np.all((ths > LOW) & (ths < HIGH), axis=1)
    # posterior density on the unconstrained scale (flat prior inside bounds)
    rr = np.array([residuals(x) for x in ths[ok]])
    logpost = -0.5*np.sum(rr*rr, axis=1)
    Zok = (ths[ok]-th)@np.linalg.inv(L).T*chi[ok][:, None]
    logq = np.array([t_logpdf(z_, df) for z_ in Zok]) - np.linalg.slogdet(cov)[1]/2
    lw = logpost - logq
    w = np.exp(lw - lw.max()); w /= w.sum()
    ess = 1.0/np.sum(w*w)
    print('samples kept %d  ESS %.1f' % (ok.sum(), ess))
    res = np.array([cross_plane(x) for x in ths[ok]])
    good = np.isfinite(res[:, 0])
    ww = w[good]/w[good].sum(); r2 = res[good]
    def wq(v, q):
        o = np.argsort(v); vv = v[o]; wc = ww[o]
        c = (np.cumsum(wc)-0.5*wc)/wc.sum()
        return np.interp(q, c, vv)
    out = {}
    for i, name in enumerate(['t', 'y', 'z']):
        v = r2[:, i]
        out[name] = dict(mean=float(np.sum(ww*v)), sd=float(np.sqrt(np.sum(ww*(v-np.sum(ww*v))**2))),
                         q90=[float(wq(v, 0.05)), float(wq(v, 0.95))],
                         q95=[float(wq(v, 0.025)), float(wq(v, 0.975))])
    clr = np.array([clear_full(y_, z_) for _, y_, z_ in r2])
    p = float(np.sum(ww*clr))
    # MCSE of the weighted probability
    mcse = float(np.sqrt(np.sum(ww**2*(clr-p)**2)))
    out['p_clear'] = p; out['p_clear_mcse'] = mcse
    out['n'] = int(good.sum()); out['ess'] = float(ess)
    out['map'] = dict(v0=list(map(float, th[:3])), cd=float(th[3]), b=[float(th[4]), float(th[5])],
                      rms_mm=float(1000*np.sqrt(np.mean(r[:6]**2))*SIGMA), t=t_map, y=y_map, z=z_map,
                      chi2=float(r@r))
    json.dump(out, open('../results/fit6.json', 'w'), indent=1)
    for k in ['t', 'y', 'z']:
        d = out[k]
        print('%s: %.4f +- %.4f  90%% [%.3f, %.3f]  95%% [%.3f, %.3f]' %
              (k, d['mean'], d['sd'], d['q90'][0], d['q90'][1], d['q95'][0], d['q95'][1]))
    print('P(full clearance) = %.5f +- %.5f' % (p, mcse))
    return out


if __name__ == '__main__':
    run(n=int(sys.argv[1]) if len(sys.argv) > 1 else 20000)
