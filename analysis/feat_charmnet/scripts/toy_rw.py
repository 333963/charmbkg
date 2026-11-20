# niepewność ε_S^dane (przeważenie ptcone40 z D*-tag) — toy po frakcjach danych
import json, os, sys
import numpy as np
sys.path.insert(0, "/home/k/Documents/charm/charm_work/FEAT_charmnet/scripts")
from study import BASE, BINS, cfg, load
sp = sys.argv[1]
S, *_ = load(sp)
F = json.load(open(os.path.join(BASE, "out", "dstar_iso_fit.json")))
edges = F["edges_ptcone40"]
ptS, wS = S["pt"].to_numpy(float), S["weight"].to_numpy(float)
ib = np.clip(np.digitize(S["ptcone40"].to_numpy(float), edges) - 1, 0, len(edges) - 2)
rng = np.random.default_rng(5)
odir = os.path.join(BASE, "out", f"{sp}_mc")
for n in ["H", "H+iso2", "H+ALL"]:
    sS = np.load(os.path.join(odir, f"scores_{n}.npz"))["S"]
    for lo, hi in BINS:
        key = "8-30" if lo >= 8 else "3-8"
        fd = np.array(F[key]["fd_frac"])
        a = (ptS >= lo) & (ptS < hi)
        s, w = sS[a], wS[a]
        o = np.argsort(-s); c = np.cumsum(w[o]) / w.sum()
        t = s[o][min(np.searchsorted(c, 0.7), len(o) - 1)]
        out = []
        for k in range(2000):
            fda = np.array([rng.normal(x["frac"], x["dfrac"]) for x in F[key]["data_frac"]])
            fda = np.clip(fda, 1e-3, None); fda /= fda.sum()
            r = fda / fd
            ww = w * r[ib[a]]
            out.append(ww[s >= t].sum() / ww.sum())
        out = np.array(out)
        print(sp, n, f"{lo}-{hi}", f"epsS_data = {np.median(out):.3f} [{np.percentile(out,16):.3f},{np.percentile(out,84):.3f}]")
