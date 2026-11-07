"""an_budget.py -- where the uncertainty in the entry point comes from.

Linear (delta-method) variance decomposition.  For a Gaussian posterior
covariance C = (A + P)^{-1}, with A = J^T J / sigma^2 the data information and
P the prior information, the identity C = C A C + C P C makes the predictive
variance exactly additive:
    Var[y] = sum_data g' C A C g  +  sum_prior g' C P_k C g.
"""
from __future__ import annotations

import json
import os

import numpy as np

import bkick as B
import estimators as E

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
R = np.array([[1.50, 0.45, 0.44], [2.96, 0.86, 0.74]])
T = np.array([0.05, 0.10])


def e1_parts(sigma):
    J = np.zeros((6, 6))
    for i, t in enumerate(T):
        J[3 * i:3 * i + 3, :3] = np.eye(3) * t
        J[3 * i:3 * i + 3, 3:] = np.eye(3) * 0.5 * t * t
    A = J.T @ J / sigma**2
    C = np.linalg.inv(A)
    v0, a = B.constant_accel_inverse(R, tuple(T))
    th = np.concatenate([v0, a])
    # prediction gradient by finite differences
    G = np.zeros((3, 6))
    for k in range(6):
        h = 1e-6 * max(1.0, abs(th[k]))
        tp = th.copy(); tp[k] += h
        tm = th.copy(); tm[k] -= h
        pp = E.quadratic_plane(tp[None, :3], tp[None, 3:])
        pm = E.quadratic_plane(tm[None, :3], tm[None, 3:])
        for j in range(3):
            G[j, k] = (pp[j][0] - pm[j][0]) / (2 * h)
    var = np.einsum("ij,jk,ik->i", G, C, G)
    e1eig = np.linalg.eigvalsh(A)
    # per-measurement contribution: A = sum_k j_k j_k^T / sigma^2
    comp = {}
    for k in range(6):
        jk = J[k]
        comp[f"obs{k}"] = np.array([(float((C @ jk) @ G[i]) / sigma) ** 2
                                    for i in range(3)])
    return var, comp, C, e1eig


def e2_parts(sigma, cd_med, cd_sig, om_sig, cl=1.0, tau=6.0):
    P = E.PhysicalPosterior(R, T, sigma=sigma, cd=(cd_med, cd_sig),
                            omega_sigma=om_sig, cl=cl, tau=tau)
    th0 = P.default_start(R)
    th, _ = P.map(th0)
    # numeric data Jacobian
    J = np.zeros((6, 7))
    for k in range(7):
        h = 1e-6 * max(1.0, abs(th[k]))
        tp = th.copy(); tp[k] += h
        tm = th.copy(); tm[k] -= h
        rp = np.asarray(P.positions(tp[None, :])).reshape(-1)
        rm = np.asarray(P.positions(tm[None, :])).reshape(-1)
        J[:, k] = (rp - rm) / (2 * h)
    A = J.T @ J / sigma**2
    # prediction gradient
    G = np.zeros((3, 7))
    for k in range(7):
        h = 1e-6 * max(1.0, abs(th[k]))
        tp = th.copy(); tp[k] += h
        tm = th.copy(); tm[k] -= h
        fp = np.array(B.fast_goal_plane(tp[:3], tp[4:7], float(np.exp(tp[3])), cl,
                                        tau=tau))
        fm = np.array(B.fast_goal_plane(tm[:3], tm[4:7], float(np.exp(tm[3])), cl,
                                        tau=tau))
        G[:, k] = (fp - fm) / (2 * h)
    blocks = {"cd": [3], "spin_x": [4], "spin_y": [5], "spin_z": [6]}
    prior_sd = {"cd": cd_sig, "spin_x": om_sig, "spin_y": om_sig, "spin_z": om_sig}
    Pmat = np.zeros((7, 7))
    for k, idx in blocks.items():
        for i in idx:
            Pmat[i, i] = 1.0 / prior_sd[k] ** 2
    C = np.linalg.inv(A + Pmat)
    var = np.einsum("ij,jk,ik->i", G, C, G)
    comp = {"data": np.einsum("ij,jk,ik->i", G, C @ A @ C, G)}
    for k in range(6):
        jk = J[k]
        comp[f"obs{k}"] = np.array([(float((C @ jk) @ G[i]) / sigma) ** 2
                                    for i in range(3)])
    for k in blocks.items():
        Pk = np.zeros((7, 7))
        for i in k[1]:
            Pk[i, i] = 1.0 / prior_sd[k[0]] ** 2
        comp[k[0]] = np.einsum("ij,jk,ik->i", G, C @ Pk @ C, G)
    var_data_share = comp["data"] / np.maximum(var, 1e-30)
    return dict(theta=th.tolist(), var=var, comp=comp, C=C, G=G,
                var_data_share=var_data_share, sigma_eig=np.linalg.eigvalsh(A),
                post_eig=np.linalg.eigvalsh(A + Pmat))


def main():
    out = {}
    for sig in (0.0029, 0.005, 0.01, 0.02):
        var, comp, C, e1e = e1_parts(sig)
        out[f"e1_sigma{sig}"] = dict(
            sd=[float(np.sqrt(max(v, 0))) for v in var],
            share={k: [float(x / max(v, 1e-30)) for x, v in zip(c, var)]
                   for k, c in comp.items()},
            total=[float(v) for v in var],
            fim_eig=[float(x) for x in e1e])
    b = e2_parts(0.01, 0.19, 0.15, 35.0)
    out["e2_base"] = dict(sd=[float(np.sqrt(max(v, 0))) for v in b["var"]],
                          share={k: [float(x / max(v, 1e-30))
                                     for x, v in zip(c, b["var"])]
                                 for k, c in b["comp"].items()},
                          var_data_share=[float(x) for x in b["var_data_share"]],
                          fim_eig=[float(x) for x in b["sigma_eig"]],
                          post_eig=[float(x) for x in b["post_eig"]])
    b2 = e2_parts(0.01, 0.19, 0.15, 0.6)
    out["e2_spinknown"] = dict(sd=[float(np.sqrt(max(v, 0))) for v in b2["var"]],
                               share={k: [float(x / max(v, 1e-30))
                                          for x, v in zip(c, b2["var"])]
                                      for k, c in b2["comp"].items()},
                               fim_eig=[float(x) for x in b2["sigma_eig"]],
                               post_eig=[float(x) for x in b2["post_eig"]])
    with open(os.path.join(OUT, "budget.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({k: {kk: (np.round(vv, 3).tolist() if isinstance(vv, list) else vv)
                          for kk, vv in v.items() if kk in ("sd", "share", "var_data_share")}
                      for k, v in out.items()}, indent=1)[:3000])


if __name__ == "__main__":
    main()
