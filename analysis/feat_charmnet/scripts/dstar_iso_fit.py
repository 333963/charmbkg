"""FEAT_charmnet: izolacja prawdziwego sygnału D⁰ w danych z tagu D* — fit Δm w klasach ptcone40.
Tło Δm: A·(Δm − m_π)^p (Poisson NLL w bokach 140–142,5 i 150–170 MeV), sygnał = N(143,5–147,5) − ∫tło.
Preselekcja wierzchołkowa (podnosi S/B; cechy wierzchołka ≈ niezależne od izolacji — sprawdzane na fd):
Lxy > LXY, cosθ_P > CP, pT(D*) 3–30, |m(D⁰) − m0| < 40 MeV. Porównanie z fd truth 0/1 po tej samej preselekcji.
"""
import glob
import json
import os

import h5py
import numpy as np
from scipy.optimize import minimize

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
COLS = ["Dmeson_decayType", "Dmeson_pt", "Dmeson_delta_mass", "Dmeson_mass", "Dmeson_ptcone40", "Dmeson_ptcone100",
        "Dmeson_truth_category", "Dmeson_Lxy", "Dmeson_cosThetaP"]
MPI = 0.13957
EDGES = [0, 1, 2, 4, 8, 1e4]
LXY, CP = 0.5, 0.995


def read(sample):
    out = {c: [] for c in COLS}
    for p in sorted(glob.glob(os.path.join(BASE, "data", sample, "*_none.h5"))):
        a = h5py.File(p, "r")["events"].fields(COLS)[:]
        s = ((a["Dmeson_decayType"] == 1) & (a["Dmeson_pt"] >= 3) & (a["Dmeson_pt"] < 30)
             & (np.abs(a["Dmeson_mass"] - a["Dmeson_delta_mass"] - 1.86484) < 0.040)
             & (a["Dmeson_Lxy"] > LXY) & (a["Dmeson_cosThetaP"] > CP))
        for c in COLS:
            out[c].append(a[c][s])
    return {c: np.concatenate(v) for c, v in out.items()}


def yield_fit(dm):
    e = np.linspace(0.140, 0.170, 61)
    h, _ = np.histogram(dm, e)
    c = 0.5 * (e[1:] + e[:-1])
    w = e[1] - e[0]
    side = (c < 0.1425) | (c > 0.150)
    x = c - MPI

    def nll(q):
        mu = np.exp(q[0]) * x ** q[1] * w
        return float(np.sum(mu[side] - h[side] * np.log(mu[side] + 1e-300)))
    r = minimize(nll, [np.log(max(h[side].mean(), 1) / w / 0.1 ** 0.5), 0.5], method="Nelder-Mead")
    A, p = np.exp(r.x[0]), r.x[1]
    lo, hi = 0.1435, 0.1475
    nsr = int(((dm > lo) & (dm < hi)).sum())
    nb = A / (p + 1) * ((hi - MPI) ** (p + 1) - (lo - MPI) ** (p + 1))
    return nsr - nb, np.sqrt(nsr + nb * 0.05 ** 2 * 0 + nb), nsr, nb


def main():
    D, M = read("data_428770"), read("fd_nd")
    res = {"preselection": f"Lxy>{LXY} mm, cosThetaP>{CP}", "edges_ptcone40": EDGES}
    for lo, hi in ((3, 30), (3, 8), (8, 30)):
        r = {}
        a = (D["Dmeson_pt"] >= lo) & (D["Dmeson_pt"] < hi)
        ys = []
        for i in range(len(EDGES) - 1):
            k = a & (D["Dmeson_ptcone40"] >= EDGES[i]) & (D["Dmeson_ptcone40"] < EDGES[i + 1])
            ys.append(yield_fit(D["Dmeson_delta_mass"][k]))
        tot = yield_fit(D["Dmeson_delta_mass"][a])
        r["data_total"] = dict(S=tot[0], dS=tot[1], N_sr=tot[2], B_sr=tot[3])
        r["data_frac"] = [dict(S=y[0], dS=y[1], frac=y[0] / tot[0], dfrac=y[1] / tot[0]) for y in ys]
        b = (M["Dmeson_pt"] >= lo) & (M["Dmeson_pt"] < hi) & np.isin(M["Dmeson_truth_category"], [0, 1]) \
            & (M["Dmeson_delta_mass"] > 0.1435) & (M["Dmeson_delta_mass"] < 0.1475)
        hm, _ = np.histogram(M["Dmeson_ptcone40"][b], EDGES)
        r["fd_frac"] = (hm / hm.sum()).tolist()
        r["fd_N"] = int(b.sum())
        bb = (D["Dmeson_delta_mass"] > 0.150) & a
        hb, _ = np.histogram(D["Dmeson_ptcone40"][bb], EDGES)
        r["data_bkg_sideband_frac"] = (hb / max(hb.sum(), 1)).tolist()
        res[f"{lo}-{hi}"] = r
        print(lo, hi, json.dumps(r, default=float), flush=True)
    json.dump(res, open(os.path.join(BASE, "out", "dstar_iso_fit.json"), "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
