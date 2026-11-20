"""FIT_scorebins — toye z LOSOWANĄ prawdziwą wydajnością (delta_k ~ Gauss(0, sig) na każdej granicy kategorii):
błąd całkowity (stat ⊕ SF) fitu z ułamkami ograniczonymi vs ustalonymi z MC. Użycie: fsb_toys_sf.py [NTOY=500] → out/toys_sf.json"""
import os, sys, json
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fsb_common import NCELL, PTB, merge, build
from fsb_run import load, configs, truth, OUT

CASES = [(3, 0.05), (5, 0.05), (3, 0.10), (5, 0.10), (8, 0.10), (5, 0.30), (8, 0.30)]


def run(a):
    sp, ib, ntoy = a
    inp = load(sp, ib)
    cfg, asim, V, F = configs(inp)
    rng = np.random.default_rng(555 + 7 * ib + (0 if sp == "D0" else 100))
    out = {}
    for K, sig in CASES:
        b = cfg[f"fix{K}"]["bounds"]
        ed = b + [NCELL]
        mc, fc = build(inp, b, "constr", sig=sig)
        mf, _ = build(inp, b, "fixed")
        E = np.cumsum(fc[::-1])[::-1]
        Bm = merge(inp["B"], b)
        r = {"c": [], "ce": [], "f": []}
        for t in range(ntoy):
            while True:
                d = np.r_[0.0, rng.normal(0, sig, K - 1)]
                Ee = E * (1 + d)
                phi = Ee - np.r_[Ee[1:], 0.0]
                if np.all(phi > 0):
                    break
            fp = inp["f"].copy()
            for c in range(K):
                fp[ed[c]:ed[c + 1]] *= phi[c] / fc[c]
            S, B = truth(inp, f=fp)
            n = merge(rng.poisson(S + B).astype(float), b)
            _, N, e, ok, _ = mc.fit(n, mc.theta0(inp["NS"], Bm))
            _, N2, e2, ok2, _ = mf.fit(n, mf.theta0(inp["NS"], Bm))
            r["c"].append(N); r["ce"].append(e); r["f"].append(N2)
        c, ce, f = (np.array(r[k]) for k in ("c", "ce", "f"))
        NS = inp["NS"]
        pull = (c - NS) / ce
        out[f"K{K}_s{int(sig*100):02d}"] = dict(constr_std=float(c.std(ddof=1) / NS), constr_hesse=float(ce.mean() / NS), bias=float((c.mean() - NS) / c.std(ddof=1)),
                                               pull_w=float(pull.std(ddof=1)), cov=float(np.mean(np.abs(pull) < 1)), fixed_std=float(f.std(ddof=1) / NS),
                                               fixed_bias=float((f.mean() - NS) / f.std(ddof=1)), ntoy=ntoy)
    return f"{sp}_{ib}", out


if __name__ == "__main__":
    ntoy = int(sys.argv[1]) if len(sys.argv) > 1 else 500
    with Pool(4) as p:
        res = dict(p.map(run, [(sp, ib, ntoy) for sp in ("D0", "Dplus") for ib in range(3)]))
    json.dump(res, open(f"{OUT}/toys_sf.json", "w"), indent=1)
    for k, v in res.items():
        for c, o in v.items():
            print(k, c, " ".join(f"{a}={100*b:.2f}%" if a.endswith(("std", "hesse", "cov")) else f"{a}={b:.2f}" for a, b in o.items()))
