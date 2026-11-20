"""PEAK_pos: czułość masy na zniekształcenia toru (MC sygnał, WP70):
 sagitta (weak mode, antysymetryczna w ładunku): q/pT → q/pT + δ, δ w TeV⁻¹  → dm/dδ dla D (q=+1); m(D)−m(D̄) = 2·dm/dδ·δ
 skala pędu (symetryczna): p → (1+α)p → dm/dα [MeV/‰]. Wynik: out/sagitta.json"""
import json
import numpy as np, pandas as pd
from common import load, PTB, OUT
MK, MPI = 0.493677, 0.13957039


def mass(trk, dlt=0.0, alpha=0.0):
    E = px = py = pz = 0
    for pt, eta, phi, q, m in trk:
        ptn = pt / (1 + q * dlt * pt) * (1 + alpha)
        a, b, c = ptn * np.cos(phi), ptn * np.sin(phi), ptn * np.sinh(eta)
        E = E + np.sqrt(a * a + b * b + c * c + m * m); px = px + a; py = py + b; pz = pz + c
    return np.sqrt(E * E - px * px - py * py - pz * pz)


res = {}
for sp in ("D0", "Dplus"):
    S = load(sp, "mc")
    if sp == "D0":
        S = S[S.hypo_swap == 0]
    nd = 2 if sp == "D0" else 3
    q = S.q.to_numpy(float)
    # ładunki córek: K = −q, π = +q
    trk = [(S[f"d_pt{i}"].to_numpy(float), S[f"d_eta{i}"].to_numpy(float), S[f"d_phi{i}"].to_numpy(float),
            (-q if i == 0 else q), MK if i == 0 else MPI) for i in range(nd)]
    m0 = mass(trk)
    d = 1e-4  # GeV⁻¹ = 0,1 TeV⁻¹
    dmdd = (mass(trk, dlt=d) - mass(trk, dlt=-d)) / (2 * d) * q    # dla D (q=+1) znak naturalny; antycząstka symetrycznie
    dmda = (mass(trk, alpha=1e-3) - mass(trk, alpha=-1e-3)) / 2e-3
    r = {}
    for name, sel in [(f"pT{lo}-{hi}", (S.pt >= lo) & (S.pt < hi)) for lo, hi in PTB] + [("pT3-30", S.pt > 0),
                      ("abseta<1", np.abs(S.eta) < 1), ("abseta>1", np.abs(S.eta) >= 1)]:
        s = np.asarray(sel)
        r[name] = dict(dm_ddelta_MeV_per_TeVinv=float(np.mean(dmdd[s]) * 1e3 * 1e-3),   # MeV na 1 TeV⁻¹
                       dm_dalpha_MeV_per_permille=float(np.mean(dmda[s]) * 1e3 * 1e-3))
    res[sp] = r
    print(sp, json.dumps(r, indent=0))
json.dump(res, open(f"{OUT}/sagitta.json", "w"), indent=1)
