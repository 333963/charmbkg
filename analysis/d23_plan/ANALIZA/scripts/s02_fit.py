"""D23/ANALIZA krok 2+3: fity m(Kππ) w 2–3 GeV (po wyborze WP w s01_wp.py), bootstrap po zdarzeniach, skan cięcia.
Histogramy nominalne: pełna statystyka (hist_*.npz). Bootstrap: wiersze score >= 0,93 (kompletne) -> Poisson(1) per zdarzenie, sparowany między WP.
Zapis: out/fit.json.  Użycie: python s02_fit.py [--nrep 300] [--nproc 7]"""
import argparse
import json
import os
import sys
import time
import multiprocessing as mp

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
TOP = os.path.dirname(HERE)
CW = os.path.dirname(os.path.dirname(TOP))
IN = os.path.join(os.path.dirname(TOP), "runs", "lxplus")
sys.path.insert(0, os.path.join(CW, "DLZ_dplus", "SCAN_2026-11-17", "K1c", "scripts"))
sys.path.insert(0, os.path.join(CW, "DLZ_dplus", "DATA_2026-11-17", "scripts"))
import fastfit as FF  # noqa: E402
import fitlib as FL  # noqa: E402

RUNS = {"427929": 1, "428770": 0}
NB = FF.NB
G = {}


def job(a):
    run, wi = a
    wpd = G["wps"][wi]
    sh = wpd["shape"]
    h = G["nom"][run][wi]
    out = dict(run=run, wp=wpd["wp"], n=float(h.sum()), fit={}, boot={})
    prev = None
    for o in (2, 3, 4):
        for free in (False, True):
            r = FF.fit(h, sh, order=o, free_shape=free, start=prev if not free else None)
            if free:
                r = FF.fit(h, sh, order=o, free_shape=True, start=parfix)
            else:
                prev = r["par"]
                parfix = r["par"]
            # błąd dmu/sc z hesjanu nie jest zwracany przez fastfit -> z bootstrapu
            out["fit"][f"{o}|{'free' if free else 'fix'}"] = {k: r[k] for k in ("N", "N_err", "chi2", "ndf", "valid", "at_limit", "SB", "dmu", "sc", "nb")}
            if o == 3 and run in G["boot"] and G["boot"][run] is not None and wi in G["bwi"]:
                hb = G["boot"][run][:, G["bwi"].index(wi), :]
                x = np.array([FF.fit_fast(hb[i], sh, r["par"], order=o, free_shape=free) for i in range(hb.shape[0])])
                out["boot"]["free" if free else "fix"] = dict(N=x[:, 0].tolist(), err=x[:, 1].tolist(), dmu=x[:, 2].tolist(), sc=x[:, 3].tolist(),
                                                              ok=((x[:, 4] > 0) & np.isfinite(x[:, 0]) & (x[:, 1] > 0)).astype(int).tolist())
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nrep", type=int, default=300)
    ap.add_argument("--nproc", type=int, default=7)
    a = ap.parse_args()
    t0 = time.time()
    WP = json.load(open(os.path.join(TOP, "out", "wp.json")))
    wps = WP["wp"]["grid"]
    G["wps"] = wps
    bwi = [i for i, w in enumerate(wps) if w["thr"] >= 0.93 - 1e-9]   # WP z kompletnymi wierszami
    G["bwi"] = bwi
    G["nom"], G["boot"] = {}, {}
    out = dict(nrep=a.nrep, info={}, cells=[], mcshape_boot={})
    rng = np.random.default_rng(20261118)

    # błąd stat. kształtu MC (mu, rms_eff) — bootstrap po kandydatach sygnału MC
    mc = pd.read_parquet(os.path.join(IN, "dp_fd_nd.parquet"), columns=["m", "score", "is_sig", "dstar_veto", "file_idx", "entry"])
    sig = mc[(mc.is_sig == 1) & (mc.dstar_veto == 0)]
    mm, ss = sig.m.to_numpy(float), sig.score.to_numpy(float)
    for wi in bwi:
        ps = ss >= wps[wi]["thr"]
        x = mm[ps]
        v = []
        for _ in range(40):
            sh = FL.mc_shape(x, rng.poisson(1.0, len(x)).astype(float))
            v.append((sh["mu"], sh["rms_eff"]))
        v = np.array(v)
        out["mcshape_boot"][str(wps[wi]["wp"])] = dict(mu_err=float(v[:, 0].std()), rms_err=float(v[:, 1].std()))
    print("błąd kształtu MC:", out["mcshape_boot"], flush=True)

    for run, cls in RUNS.items():
        hd = np.load(os.path.join(IN, f"hist_data_{run}.npz"))
        Hc = hd["dp_h"][cls].sum(0)
        H5 = Hc.reshape(Hc.shape[0], 100, 5).sum(2).astype(float)
        G["nom"][run] = [H5[w["ie"]:].sum(0) for w in wps]
        d = pd.read_parquet(os.path.join(IN, f"dp_data_{run}.parquet"),
                            columns=["run", "event", "m", "score", "dstar_veto", "grl", "trig_HLT_mb_mbts_L1MBTS_1"], filters=[("score", ">=", 0.93)])
        sel = (d.dstar_veto == 0) & (d.score >= 0.93)
        if cls == 1:
            sel &= (d.grl == 1) & (d["trig_HLT_mb_mbts_L1MBTS_1"] == 1)
        d = d[sel]
        ev = d.run.to_numpy(np.int64) * 10_000_000_000 + d.event.to_numpy(np.int64)
        ue, einv, cnt = np.unique(ev, return_inverse=True, return_counts=True)
        mbin = np.clip(((d.m.to_numpy(float) - FF.LO) / 0.005).astype(int), 0, NB - 1)
        s = d.score.to_numpy(np.float32)
        se = hd["s_edges"]
        sidx = np.clip(np.searchsorted(se, s, side="right") - 1, 0, len(se) - 2)   # jak ibin() w ekstrakcji
        BW = rng.poisson(1.0, (a.nrep, len(ue))).astype(np.float32)
        boot = np.zeros((a.nrep, len(bwi), NB))
        eq = True
        for j, wi in enumerate(bwi):
            ps = sidx >= wps[wi]["ie"]
            hrow = np.bincount(mbin[ps], minlength=NB).astype(float)
            eq &= bool(np.array_equal(hrow, G["nom"][run][wi]))
            e_, m_ = einv[ps], mbin[ps]
            for r in range(a.nrep):
                boot[r, j] = np.bincount(m_, weights=BW[r, e_], minlength=NB)
        G["boot"][run] = boot
        out["info"][run] = dict(rows_hi=int(len(d)), n_events=int(len(ue)), cand_per_event_mean=float(cnt.mean()), cand_per_event_max=int(cnt.max()),
                                frac_cand_in_multi=float(cnt[cnt > 1].sum() / cnt.sum()), rows_hist_equal=eq)
        print(f"{run}: wiersze score>=0.93 {len(d)}, zdarzeń {len(ue)}, kand./zdarzenie {cnt.mean():.3f} (max {cnt.max()}), wiersze==hist pełny: {eq}  ({time.time()-t0:.0f} s)", flush=True)
        del BW
    jobs = [(run, wi) for run in RUNS for wi in range(len(wps))]
    with mp.get_context("fork").Pool(a.nproc) as p:
        out["cells"] = p.map(job, jobs, chunksize=1)
    out["hists"] = {run: [h.tolist() for h in G["nom"][run]] for run in RUNS}
    json.dump(out, open(os.path.join(TOP, "out", "fit.json"), "w"))
    print(f"gotowe {time.time()-t0:.0f} s")


if __name__ == "__main__":
    main()
