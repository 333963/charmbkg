"""ISO_sf: SF_iso = ε_dane/ε_MC dla cięć izolacji D0 (wiersz Kπ) z tagu D*.
Dane: out/data_428770_tag.parquet (2 lokalne pliki h5 v1.2 z ptcone), ε z jednoczesnego fitu Δm pass/fail
(charmtag.dstar_refit.fit_passfail, kształt sygnału z MC truth). MC: out/fd_nd_tag.parquet (8 lokalnych plików), ε z prawdy
(truth ∈ {0,1}, Δm 143,5–147,5 MeV), te same cięcia. Przeważenie krotności: w(n_trkid) = f_dane^sygn/f_MC^sygn, f_dane z fitów
pass/fail „w binie / poza binem” (proxy n_trk: liczba różnych torów córek wszystkich kandydatów w zdarzeniu).
Użycie: python sf_fit.py [NPROC]  → out/sf_fits.json
"""
import json, os, sys
from multiprocessing import Pool
import numpy as np, pandas as pd
sys.path.insert(0, "/home/k/Documents/charm/charm_work/AA_merge/src")
from charmtag import dstar_refit as R

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# 16.11 ISO-FULL: ISO_IN (katalog z data_428770_tag.parquet, fd_nd_tag.parquet) / ISO_OUT (wyniki); domyślnie out/ jak dotąd
IN = os.environ.get("ISO_IN", os.path.join(BASE, "out")); OUT = os.environ.get("ISO_OUT", os.path.join(BASE, "out"))
M0, MW = 1.86484, 0.040
E = np.arange(R.LO, R.HI + 1e-6, 0.1)
PRESEL = {"A_score0.6": lambda X: X.score > 0.6, "B_none": lambda X: X.score > -1}
PTB = {"3-30": (3, 30), "3-5": (3, 5), "5-8": (5, 8), "8-30": (8, 30)}
CUTS = [("ptcone40", c) for c in (0.5, 1, 2, 4)] + [("ptcone100", c) for c in (1, 2, 4, 8)] \
    + [("rel40", c) for c in (0.1, 0.25, 0.5, 1.0)] + [("rel100", c) for c in (0.25, 0.5, 1.0, 2.0)]
NTB = [3, 15, 25, 35, 50, 1000]     # biny n_trkid (do wag)
NTB2 = {"lo": (3, 25), "hi": (25, 1000)}   # podział do SF(krotność)


def load():
    D = pd.read_parquet(os.path.join(IN, "data_428770_tag.parquet"))
    M = pd.read_parquet(os.path.join(IN, "fd_nd_tag.parquet"))
    for X in (D, M):
        X["rel40"], X["rel100"] = X.ptcone40 / X.pt, X.ptcone100 / X.pt
    D = D[np.abs(D.m - M0) < MW].reset_index(drop=True)
    M = M[(np.abs(M.m - M0) < MW) & M.truth.isin([0, 1]) & (M.dm > 0.1435) & (M.dm < 0.1475)].reset_index(drop=True)
    return D, M


def _fit(job):
    key, hp, hf, shape = job
    try:
        r = R.fit_passfail(hp, hf, E, shape)
        out = {k: (float(r[k]) if isinstance(r.get(k), (int, float, np.floating, np.integer, bool, np.bool_)) else None)
               for k in ("eps", "eps_err", "eps_lo", "eps_hi", "N_tot", "N_tot_err", "scale", "eps_expo", "dnll_expo")}
        out["determined"] = bool(r["determined"]); out["n_pass"], out["n_fail"] = int(hp.sum()), int(hf.sum())
    except Exception as ex:  # noqa
        out = {"error": str(ex)}
    return key, out


def hist(dm):
    return np.histogram(dm * 1e3, E)[0]


def main(nproc):
    D, M = load()
    shape = R.shape_from_mc(M.dm.values * 1e3)
    jobs, mc = [], {}
    w_ok = {}
    for ps, f in PRESEL.items():
        d0, m0 = D[f(D)], M[f(M)]
        for pb, (lo, hi) in PTB.items():
            if ps != "A_score0.6" and pb != "3-30":
                continue
            d, m = d0[(d0.pt >= lo) & (d0.pt < hi)], m0[(m0.pt >= lo) & (m0.pt < hi)]
            for v, c in CUTS:
                k = f"{ps}|{pb}|{v}<{c}"
                p = d[v] < c
                jobs.append((k, hist(d.dm[p]), hist(d.dm[~p]), shape))
                mc[k] = {"eps_mc": float((m[v] < c).mean()), "n_mc": int(len(m)),
                         "eps_mc_by_ntrkid": [float((m[v][(m.n_trkid >= a) & (m.n_trkid < b)] < c).mean())
                                              for a, b in zip(NTB[:-1], NTB[1:])]}
            if pb == "3-30":
                # rozkład krotności sygnału w danych: frakcje w binach n_trkid i n_cand
                for var in ("n_trkid",):
                    for a, b in zip(NTB[:-1], NTB[1:]):
                        p = (d[var] >= a) & (d[var] < b)
                        k = f"{ps}|{pb}|frac_{var}[{a},{b})"
                        jobs.append((k, hist(d.dm[p]), hist(d.dm[~p]), shape))
                        mc[k] = {"eps_mc": float(((m[var] >= a) & (m[var] < b)).mean()), "n_mc": int(len(m))}
                # SF w podziale na krotność
                for tag, (a, b) in NTB2.items():
                    dd, mm = d[(d.n_trkid >= a) & (d.n_trkid < b)], m[(m.n_trkid >= a) & (m.n_trkid < b)]
                    for v, c in [("ptcone40", 1), ("ptcone40", 2), ("ptcone100", 4), ("rel40", 0.25)]:
                        k = f"{ps}|{pb}|ntrkid_{tag}|{v}<{c}"
                        p = dd[v] < c
                        jobs.append((k, hist(dd.dm[p]), hist(dd.dm[~p]), shape))
                        mc[k] = {"eps_mc": float((mm[v] < c).mean()), "n_mc": int(len(mm))}
    print("fitów:", len(jobs), flush=True)
    res = {}
    with Pool(nproc) as pool:
        for k, o in pool.imap_unordered(_fit, jobs):
            o.update(mc[k]); res[k] = o
            print(k, {x: (round(y, 4) if isinstance(y, float) else y) for x, y in o.items() if x in ("eps", "eps_err", "N_tot", "eps_mc", "determined")}, flush=True)
    json.dump({"shape": shape, "res": res, "NTB": NTB}, open(os.path.join(OUT, "sf_fits.json"), "w"), indent=1)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)
