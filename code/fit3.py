"""M3 main line: drag plus a constant lateral acceleration (five parameters).

Our own least-squares fit (damped Gauss-Newton) and our own sampler
(random-walk Metropolis with a Laplace proposal), pure-Python RK4.
Writes ../results/fit3.json.
"""
import json, math, numpy as np

G, RAD, XG = 9.81, 0.110, 11.0
FACT = 1.20*math.pi*0.110**2/(2*0.430)
MU = np.array([30., 0., 7., 0.25, 0.0])
SD = np.array([15., 15., 10., 0.15, 20.0])
LO = np.array([5., -30., 0., 0., -60.])
HI = np.array([60., 30., 30., 1., 60.])
SIG = 0.005
TIMES = (0.05, 0.10)


def advance(th, targets, dt=0.002):
    vx, vy, vz = th[0], th[1], th[2]
    kd = FACT*th[3]; al = th[4]
    px, py, pz = 0.0, 0.0, 0.11
    out = []; t = 0.0
    def acc(x, y, z):
        s = math.sqrt(x*x + y*y + z*z)
        return (-kd*s*x, -kd*s*y + al, -kd*s*z - G)
    for tgt in targets:
        while t < tgt - 1e-12:
            h = min(dt, tgt - t)
            a1 = acc(vx, vy, vz)
            a2 = acc(vx + .5*h*a1[0], vy + .5*h*a1[1], vz + .5*h*a1[2])
            a3 = acc(vx + .5*h*a2[0], vy + .5*h*a2[1], vz + .5*h*a2[2])
            a4 = acc(vx + h*a3[0], vy + h*a3[1], vz + h*a3[2])
            px += h/6*(vx + 2*(vx + .5*h*a1[0]) + 2*(vx + .5*h*a2[0]) + vx + h*a3[0])
            py += h/6*(vy + 2*(vy + .5*h*a1[1]) + 2*(vy + .5*h*a2[1]) + vy + h*a3[1])
            pz += h/6*(vz + 2*(vz + .5*h*a1[2]) + 2*(vz + .5*h*a2[2]) + vz + h*a3[2])
            vx += h/6*(a1[0] + 2*a2[0] + 2*a3[0] + a4[0])
            vy += h/6*(a1[1] + 2*a2[1] + 2*a3[1] + a4[1])
            vz += h/6*(a1[2] + 2*a2[2] + 2*a3[2] + a4[2])
            t += h
        out.append((px, py, pz))
    while px < XG and t < 1.0:
        h = dt
        a1 = acc(vx, vy, vz)
        a2 = acc(vx + .5*h*a1[0], vy + .5*h*a1[1], vz + .5*h*a1[2])
        a3 = acc(vx + .5*h*a2[0], vy + .5*h*a2[1], vz + .5*h*a2[2])
        a4 = acc(vx + h*a3[0], vy + h*a3[1], vz + h*a3[2])
        npx = px + h/6*(vx + 2*(vx + .5*h*a1[0]) + 2*(vx + .5*h*a2[0]) + vx + h*a3[0])
        npy = py + h/6*(vy + 2*(vy + .5*h*a1[1]) + 2*(vy + .5*h*a2[1]) + vy + h*a3[1])
        npz = pz + h/6*(vz + 2*(vz + .5*h*a1[2]) + 2*(vz + .5*h*a2[2]) + vz + h*a3[2])
        if npx >= XG:
            fr = (XG - px)/(npx - px)
            return out, (t + fr*h, py + fr*(npy - py), pz + fr*(npz - pz))
        px, py, pz = npx, npy, npz
        vx += h/6*(a1[0] + 2*a2[0] + 2*a3[0] + a4[0])
        vy += h/6*(a1[1] + 2*a2[1] + 2*a3[1] + a4[1])
        vz += h/6*(a1[2] + 2*a2[2] + 2*a3[2] + a4[2])
        t += h
    return out, (float('nan'),)*3


OBS = None


def set_obs(r1=(1.50, 0.45, 0.44), r2=(2.96, 0.86, 0.74)):
    global OBS
    OBS = np.array([r1, r2])


set_obs()


def resid(th):
    pos, _ = advance(th, TIMES)
    d = (OBS - np.array(pos)).ravel()/SIG
    pr = (th - MU)/SD
    return np.concatenate([d, pr])


def fit_map(iters=60):
    th = np.clip(np.array([30., 5., 7., 0.25, -5.]), LO, HI)
    r = resid(th); cost = r @ r; lam = 1e-2
    for _ in range(iters):
        J = np.zeros((r.size, 5))
        for i in range(5):
            h = 1e-5*max(1., abs(th[i]))
            tp = th.copy(); tp[i] += h; tm = th.copy(); tm[i] -= h
            J[:, i] = (resid(tp) - resid(tm))/(2*h)
        A = J.T @ J; g = J.T @ r; ok = False
        for _ in range(10):
            try:
                step = np.linalg.solve(A + lam*np.diag(np.diag(A) + 1e-12), -g)
            except np.linalg.LinAlgError:
                lam *= 10; continue
            tn = np.clip(th + step, LO, HI)
            rn = resid(tn); cn = rn @ rn
            if cn < cost:
                th, r, cost = tn, rn, cn; lam = max(lam*0.3, 1e-9); ok = True; break
            lam *= 10
        if not ok:
            break
    # Laplace covariance
    J = np.zeros((r.size, 5))
    for i in range(5):
        h = 1e-5*max(1., abs(th[i]))
        tp = th.copy(); tp[i] += h; tm = th.copy(); tm[i] -= h
        J[:, i] = (resid(tp) - resid(tm))/(2*h)
    cov = np.linalg.inv(J.T @ J)
    return th, cost, cov


def sample(th0, cov, n=12000, burn=2000, thin=2):
    rng = np.random.default_rng(5)
    L = np.linalg.cholesky(cov*2.4**2/5)
    def lnp(th):
        if np.any(th < LO) or np.any(th > HI):
            return -np.inf
        r = resid(th)
        return -0.5*(r @ r)
    th = th0.copy(); lp = lnp(th); keep = []; acc = 0
    for i in range(n):
        cand = th + L @ rng.normal(0, 1, 5)
        lq = lnp(cand)
        if math.log(rng.random() + 1e-300) < lq - lp:
            th, lp = cand, lq; acc += 1
        if i >= burn and i % thin == 0:
            keep.append(th.copy())
    return np.array(keep), acc/n


if __name__ == '__main__':
    th, cost, cov = fit_map()
    pos, tc = advance(th, TIMES)
    print('MAP', np.round(th, 5), 'chi2 %.2f' % cost)
    print('MAP entry t %.4f y %.4f z %.4f' % tc)
    _, tc = advance(th, TIMES)
    chain, acc = sample(th, cov)
    Y, Z, T = [], [], []
    for t in chain:
        _, c = advance(t, TIMES, dt=0.004)
        if np.isfinite(c[1]):
            T.append(c[0]); Y.append(c[1]); Z.append(c[2])
    T, Y, Z = map(np.array, (T, Y, Z))
    inb = (np.abs(Y) <= 3.55) & (Z >= 0.11) & (Z <= 2.33)
    out = dict(map=[float(x) for x in th], chi2=float(cost), acc=float(acc),
               n=int(len(Y)), t=float(T.mean()), y=float(Y.mean()), z=float(Z.mean()),
               sd_y=float(Y.std()), sd_z=float(Z.std()), sd_t=float(T.std()),
               q90_y=[float(np.quantile(Y, .05)), float(np.quantile(Y, .95))],
               q90_z=[float(np.quantile(Z, .05)), float(np.quantile(Z, .95))],
               q95_y=[float(np.quantile(Y, .025)), float(np.quantile(Y, .975))],
               q95_z=[float(np.quantile(Z, .025)), float(np.quantile(Z, .975))],
               p_clean=float(inb.mean()),
               p_wide=float((np.abs(Y) > 3.66).mean()),
               p_over=float(((Z > 2.44) & (np.abs(Y) <= 3.66)).mean()),
               p_clips=float(1 - inb.mean() - (np.abs(Y) > 3.66).mean()
                             - ((Z > 2.44) & (np.abs(Y) <= 3.66)).mean()))
    print(json.dumps({k: (np.round(v, 4).tolist() if isinstance(v, list) else
                          (round(v, 4) if isinstance(v, float) else v))
                      for k, v in out.items()}, indent=1))
    json.dump(out, open('../results/fit3.json', 'w'), indent=1)
