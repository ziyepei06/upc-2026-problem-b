"""Prior-sensitivity tornado and noise sweep, recomputed with the corrected
spin-decay model. Deltas are taken against a baseline run with the SAME
sampler settings, so the figure can be anchored at the paper's E2 answer."""
import json
import numpy as np
import estimators as E

R = np.array([[1.50, 0.45, 0.44], [2.96, 0.86, 0.74]])
T = np.array([0.05, 0.10])
PAPER_SPEC = dict(n=40000, burn=8000, thin=2, seed=1)   # same as predict.py
COARSE = dict(n=12000, burn=4000, thin=1, seed=3)

def run(sigma=0.010, spec=None, **kw):
    spec = spec or COARSE
    P = E.PhysicalPosterior(R, T, sigma=sigma, **kw)
    out = P.sample(**spec)
    t, y, z = P.predict(out["chain"])
    clean = (np.abs(y) <= 3.55) & (z >= 0.11) & (z <= 2.33)
    return dict(mean=float(y.mean()), std=float(y.std()),
                z=float(z.mean()), sz=float(z.std()), p=float(clean.mean()))

base_exact = run(spec=PAPER_SPEC)
print("base_exact", {k: round(v, 4) for k, v in base_exact.items()}, flush=True)
base_coarse = run()
print("base_coarse", {k: round(v, 4) for k, v in base_coarse.items()}, flush=True)

rows = []
for key, vals in (("omega_sigma", (17., 70.)), ("cl", (0.7, 1.3)),
                  ("cd_median", (0.17, 0.22)), ("cd_log_sigma", (0.30, 0.50)),
                  ("tau", (2., 12.))):
    for v in vals:
        if key == "cd_median":
            kw = dict(cd=(v, 0.15))
        elif key == "cd_log_sigma":
            kw = dict(cd=(0.19, v))
        else:
            kw = {key: v}
        r = run(**kw)
        rows.append(dict(parkey=key, value=v, mean=r["mean"], std=r["std"], p=r["p"],
                         dy=r["mean"] - base_coarse["mean"],
                         dstd=r["std"] - base_coarse["std"]))
        print(key, v, "dy=%+.3f dsd=%+.3f p=%.3f" % (rows[-1]["dy"], rows[-1]["dstd"], r["p"]), flush=True)
json.dump(dict(base_exact=base_exact, base_coarse=base_coarse, rows=rows),
          open("../results/sens_prior.json", "w"), indent=1)

sweep = []
for sg in (0.005, 0.010, 0.015, 0.020):
    r = run(sigma=sg, spec=PAPER_SPEC)
    sweep.append(dict(sigma=sg, **r))
    print("sweep", {k: round(x, 4) for k, x in r.items()}, flush=True)
json.dump(sweep, open("../results/e2_sweep.json", "w"), indent=1)
print("done")
