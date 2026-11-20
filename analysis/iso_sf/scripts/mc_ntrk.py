"""ISO_sf: zależność ε_iso(MC) od krotności n_trk (nTracksPt500Eta25) i efekt przeważenia w1(n_trk) = p_dane/p_MC(ND)
(etap1_definicje_masy_v20260925_01/weights_ntrk.csv; mapa 2D w2(n_trk, lead_pt) nie jest dostępna lokalnie — w2d.npy tylko na lxplus).
Wariant „w1×w_lead”: dodatkowo przybliżenie lead_pt — brak, zob. RAPORT.
Wejście: out/root_mc_dstar.parquet (D* truth, v9), out/root_mc_d0incl.parquet (D0 truth inkluzywne).
Też: przeważenie MC do rozkładu n_trkid sygnału w danych (z sf_fits.json) — wariant data-driven. → out/mc_ntrk.json"""
import json, os
import numpy as np, pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
W1 = pd.read_csv("/home/k/Documents/charm/etap1_definicje_masy_v20260925_01/weights_ntrk.csv")
E1 = np.r_[W1.lo.values, W1.hi.values[-1]]
CUTS = [("ptcone40", c) for c in (0.5, 1, 2, 4)] + [("ptcone100", c) for c in (1, 2, 4, 8)] \
    + [("rel40", c) for c in (0.1, 0.25, 0.5, 1.0)] + [("rel100", c) for c in (0.25, 0.5, 1.0, 2.0)]
PTB = {"3-30": (3, 30), "3-5": (3, 5), "5-8": (5, 8), "8-30": (8, 30)}
NT = [0, 10, 20, 30, 40, 60, 1e4]


def w1(n):
    return W1.w.values[np.clip(np.digitize(n, E1) - 1, 0, len(W1) - 1)]


def eff(x, w):
    s = w.sum(); e = (x * w).sum() / s; neff = s ** 2 / (w ** 2).sum()
    return float(e), float(np.sqrt(e * (1 - e) / neff)), float(neff)


def main():
    out = {}
    for nm in ("dstar", "d0incl"):
        X = pd.read_parquet(os.path.join(BASE, "out", f"root_mc_{nm}.parquet"))
        if nm == "dstar":
            X = X[(np.abs(X.m - 1.86484) < 0.040) & (X.dm > 0.1435) & (X.dm < 0.1475)]
        X = X.assign(rel40=X.ptcone40 / X.pt, rel100=X.ptcone100 / X.pt, w1=w1(X.n_trk.values))
        r = {"N": int(len(X)), "mean_n_trk": float(X.n_trk.mean()), "mean_n_trk_w1": float(np.average(X.n_trk, weights=X.w1)),
             "mean_n_trkid": float(X.n_trkid.mean()), "mean_n_trkid_w1": float(np.average(X.n_trkid, weights=X.w1)),
             "corr_ntrk_ntrkid": float(np.corrcoef(X.n_trk, X.n_trkid)[0, 1]), "cuts": {}}
        for pb, (lo, hi) in PTB.items():
            x = X[(X.pt >= lo) & (X.pt < hi)]
            for v, c in CUTS:
                p = (x[v] < c).to_numpy(float)
                k = f"{pb}|{v}<{c}"
                r["cuts"][k] = {"unw": eff(p, np.ones(len(x))), "w1": eff(p, x.w1.to_numpy()),
                                "by_ntrk": [eff(p[m], np.ones(m.sum()))[0] if m.sum() > 20 else None
                                            for m in [((x.n_trk >= a) & (x.n_trk < b)).to_numpy() for a, b in zip(NT[:-1], NT[1:])]]}
        out[nm] = r
        print(nm, json.dumps({k: r[k] for k in r if k != "cuts"}))
        for k in ("3-30|ptcone40<1", "3-30|ptcone40<2", "3-30|ptcone100<4", "3-30|rel40<0.25"):
            print("  ", k, r["cuts"][k]["unw"][:2], r["cuts"][k]["w1"][:2], r["cuts"][k]["by_ntrk"])
    out["NT"] = NT
    json.dump(out, open(os.path.join(BASE, "out", "mc_ntrk.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
