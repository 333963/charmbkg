"""FEAT_charmnet: przeniesienie zysku z MC na dane (modele z testu mc).
ε_S^dane ≈ ε_S MC przeważone w klasach ptcone40 stosunkiem f_dane/f_MC zmierzonym na D*-tag (dstar_iso_fit.json; 3–8 i 8–30 GeV;
dla 3–5 i 5–8 użyty 3–8). ε_B^dane = dalekie SB danych (z wetami) ocenione modelem MC (OOF po foldach zdarzeń).
Progi WP: z sygnału MC (jak w analizie). Z_transfer = ε_S^dane/√ε_B^dane. Użycie: python transfer.py SPECIES
"""
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from study import BASE, BINS, WPS, boot_weights, cfg, load  # noqa: E402


def main(sp):
    S, Bmc, Bda, dall, mb = load(sp)
    F = json.load(open(os.path.join(BASE, "out", "dstar_iso_fit.json")))
    edges = F["edges_ptcone40"]
    far = {"D0": 0.100, "Dplus": 0.110}[sp]
    m0 = cfg["species"][sp]["m0"]
    dsel = np.abs(dall["m"].to_numpy() - m0) > far  # dall już po wetach
    D = dall[dsel]
    odir = os.path.join(BASE, "out", f"{sp}_mc")
    ptS, wS = S["pt"].to_numpy(float), S["weight"].to_numpy(float)
    ib = np.clip(np.digitize(S["ptcone40"].to_numpy(float), edges) - 1, 0, len(edges) - 2)
    rw = np.ones(len(S))
    for lo, hi in BINS:
        key = "8-30" if lo >= 8 else "3-8"
        fd = np.array(F[key]["fd_frac"])
        fda = np.array([x["frac"] for x in F[key]["data_frac"]])
        fda = np.clip(fda, 0, None)
        fda = fda / fda.sum()
        r = np.where(fd > 0, fda / fd, 1.0)
        a = (ptS >= lo) & (ptS < hi)
        rw[a] = r[ib[a]]
    ptD = D["pt"].to_numpy(float)
    bw = boot_weights(D["evkey"].to_numpy(np.int64), 200, 21)
    rows = []
    for f in sorted(os.listdir(odir)):
        if not (f.startswith("scores_") and f.endswith(".npz")):
            continue
        n = f[7:-4]
        z = np.load(os.path.join(odir, f))
        sS, sD = z["S"], z["dall"][dsel]
        for lo, hi in BINS:
            a = (ptS >= lo) & (ptS < hi)
            b = (ptD >= lo) & (ptD < hi)
            s, w = sS[a], wS[a]
            o = np.argsort(-s)
            c = np.cumsum(w[o]) / w.sum()
            for wp in WPS:
                t = s[o][min(np.searchsorted(c, wp), len(o) - 1)]
                es_mc = w[s >= t].sum() / w.sum()
                es_da = (w * rw[a])[s >= t].sum() / (w * rw[a]).sum()
                pas = sD[b] >= t
                eb = pas.mean()
                ebr = np.array([(bw[r][b] * pas).sum() / bw[r][b].sum() for r in range(200)])
                rows.append(dict(variant=n, bin=f"{lo}-{hi}", wp=wp, epsS_mc=es_mc, epsS_data_iso_rw=es_da, epsB_dataSB=eb,
                                 Z_mcS=es_mc / np.sqrt(eb), Z_transfer=es_da / np.sqrt(eb),
                                 Z_transfer_err=float(np.std(es_da / np.sqrt(np.maximum(ebr, 1e-9))))))
    T = pd.DataFrame(rows)
    for col in ("Z_mcS", "Z_transfer"):
        T[col + "_ratio_vs_H"] = T.apply(lambda r: r[col] / T[(T.variant == "H") & (T.bin == r["bin"]) & (T.wp == r["wp"])][col].iloc[0], axis=1)
    T.to_csv(os.path.join(BASE, "out", f"transfer_{sp}.csv"), index=False)
    print(T[T.wp == 0.7].to_string(index=False, float_format=lambda x: f"{x:.3g}"))


if __name__ == "__main__":
    main(sys.argv[1])
