"""HIGHSCORE krok 3 — cechy H: sygnał danych (okno Δm − przeskalowane pasma boczne) vs MC, ranking rozbieżności, test domknięcia
przez przeważenie MC w 1–2 cechach (TEST, nie kalibracja). EKSPLORACYJNE, PO FAKCIE. python hs_feats.py [NBOOT] -> out/feats.json"""
import json, os, sys
import multiprocessing as mp
import numpy as np
import hs_common as C
R = C.R
NB = int(sys.argv[1]) if len(sys.argv) > 1 else 200
NQ = 10
WIN = 2.0
ALLF = C.FEATS + C.EXTRA
JOBS = [(428770, (8, 30)), (428770, (8, 12)), (428770, (12, 30)), (427929, (8, 30))]
REG = {"all": None, "geWP70": "0.70", "geWP50": "0.50", "gtWP25": "q25"}


def job(a):
    run, b = a
    d = C.load(f"data_{run}"); M = C.load("fd_nd")
    mu = C.MU["8-30"]
    thr6 = C.thr_for("8-30", R.WP_MAIN)                       # progi WP95..25 binu 8–30
    tcut = {"all": -1.0, "geWP70": thr6[2], "geWP50": thr6[3], "gtWP25": thr6[4]}
    sd0 = C.base_mask(d, b, isdata=True, run=run)
    sm0 = C.base_mask(M, b, isdata=False) & (np.abs(M["dm"] - mu) < WIN)
    wpt_all, _ = C.pt_weights(run)
    cen = R.CEN
    inwin = np.abs(d["dm"] - mu) < WIN
    insb = ((d["dm"] > 141.0) & (d["dm"] < 143.0)) | ((d["dm"] > 149.0) & (d["dm"] < 155.0))
    sbb = ((cen > 141.0) & (cen < 143.0)) | ((cen > 149.0) & (cen < 155.0)); wb = np.abs(cen - mu) < WIN
    msc = M["score"][sm0]; wpt = wpt_all[sm0]; wpt = wpt / wpt.mean()
    mtop = msc >= thr6[4]; m70 = msc >= thr6[2]; m50 = msc >= thr6[3]
    out = {"run": run, "bin": list(b), "regions": {}}
    rng = np.random.default_rng(900 + run)
    BW = [rng.poisson(1.0, d["nev"]).astype(np.float32) for _ in range(NB)]

    def mc_obs(w):
        return {"Rtop": (w * mtop).sum() / w.sum(), "rho70": (w * mtop).sum() / (w * m70).sum(), "rho50": (w * mtop).sum() / (w * m50).sum()}
    for mcw_name, w0 in (("MC", np.ones(len(msc))), ("MCptw", wpt)):
        for rname, tc in tcut.items():
            base_obs = mc_obs((msc >= tc).astype(float))        # odniesienie: MC nieważone, ograniczone do regionu
            sd = sd0 & (d["score"] >= tc)
            ib, ok = C.dm_bin(d["dm"][sd])
            h = np.bincount(ib[ok], minlength=R.NB).astype(float)
            p = C.fit_bkg(h, np.abs(cen - mu) > 3.0, 2 if h[np.abs(cen - mu) > 3.0].sum() >= R.NSB_FULL else 1)
            bs = p[0] * R.bkg_shape(*p[1:])
            k = bs[wb].sum() / bs[sbb].sum()
            iw = np.where(sd & inwin)[0]; isb = np.where(sd & insb)[0]
            Nsig = len(iw) - k * len(isb)
            mreg = msc >= tc
            rows = {}
            keep = {}
            for f in ALLF:
                xm = M[f][sm0][mreg].astype(float)
                ed = np.unique(np.quantile(xm, np.linspace(0, 1, NQ + 1)[1:-1]))
                nb = len(ed) + 1
                bm_all = np.searchsorted(ed, M[f][sm0].astype(float), side="right")
                pm = np.bincount(bm_all[mreg], weights=w0[mreg], minlength=nb); 
                neff = w0[mreg].sum() ** 2 / (w0[mreg] ** 2).sum()
                pm = pm / pm.sum()
                bw_ = np.searchsorted(ed, d[f][iw].astype(float), side="right"); bsb = np.searchsorted(ed, d[f][isb].astype(float), side="right")
                hd = np.bincount(bw_, minlength=nb) - k * np.bincount(bsb, minlength=nb)
                pd_ = hd / hd.sum()
                reps = np.empty((NB, nb))
                for i, w in enumerate(BW):
                    hh = np.bincount(bw_, weights=w[d["inv"][iw]], minlength=nb) - k * np.bincount(bsb, weights=w[d["inv"][isb]], minlength=nb)
                    reps[i] = hh / hh.sum()
                sdv = reps.std(0, ddof=1)
                sm_ = np.sqrt(pm * (1 - pm) / neff)
                chi2 = float((((pd_ - pm) ** 2) / (sdv ** 2 + sm_ ** 2)).sum())
                rank = (np.arange(nb) + 0.5) / nb
                zr = float((pd_ @ rank - pm @ rank) / np.sqrt((reps @ rank).std(ddof=1) ** 2 + (pm @ (rank - pm @ rank) ** 2) / neff))
                ks = np.abs(np.cumsum(pd_) - np.cumsum(pm)).max()
                ks_null = np.abs(np.cumsum(reps - reps.mean(0), axis=1)).max(1)             # rozrzut KS z samych fluktuacji
                rows[f] = {"chi2": chi2, "ndf": nb - 1, "z_rank": zr, "ks": float(ks), "ks_p": float((ks_null >= ks).mean()),
                           "pd": pd_.tolist(), "pd_err": sdv.tolist(), "pm": pm.tolist(), "edges": ed.tolist()}
                # test domknięcia: wagi 1D z tego regionu nałożone na MC (kandydaci regionu), obserwable ogona
                def closure(pdv):
                    wr = np.clip(pdv / np.clip(pm, 1e-9, None), 0, None)
                    ww = np.where(mreg, w0 * wr[bm_all], 0.0)
                    o = mc_obs(ww)
                    return [o[q] / base_obs[q] for q in ("Rtop", "rho70", "rho50")]
                cn = closure(pd_); cr = np.array([closure(reps[i]) for i in range(NB)])
                rows[f]["closure"] = dict(zip(("Rtop", "rho70", "rho50"), cn)); rows[f]["closure_err"] = dict(zip(("Rtop", "rho70", "rho50"), cr.std(0, ddof=1).tolist()))
                keep[f] = (bm_all, pm, pd_, reps)
            # raking 2 cech o najwyższym chi2/ndf (poza log_pt przy MC pT-ważonym) — wybór z nominału
            order = sorted(C.FEATS, key=lambda f: -rows[f]["chi2"] / max(rows[f]["ndf"], 1))
            pair = order[:2]

            def rake(pds):
                ww = np.where(mreg, w0, 0.0).astype(float)
                for _ in range(6):
                    for f, pdv in zip(pair, pds):
                        bm_all = keep[f][0]; nbf = len(pdv)
                        cur = np.bincount(bm_all, weights=ww, minlength=nbf); cur = cur / cur.sum()
                        ww = ww * np.clip(pdv / np.clip(cur, 1e-9, None), 0, None)[bm_all]
                o = mc_obs(ww)
                return [o[q] / base_obs[q] for q in ("Rtop", "rho70", "rho50")]
            rn = rake([keep[f][2] for f in pair]); rr = np.array([rake([keep[f][3][i] for f in pair]) for i in range(0, NB, 2)])
            base_w = mc_obs(np.where(mreg, w0, 0.0))
            out["regions"][f"{rname}|{mcw_name}"] = {"N_sig": float(Nsig), "k": float(k), "n_win": len(iw), "n_sb": len(isb), "feats": rows,
                                                     "pair": pair, "rake": dict(zip(("Rtop", "rho70", "rho50"), rn)),
                                                     "rake_err": dict(zip(("Rtop", "rho70", "rho50"), rr.std(0, ddof=1).tolist())),
                                                     "mc_w_only": {q: base_w[q] / base_obs[q] for q in base_w}}
    return out


if __name__ == "__main__":
    C.load("fd_nd"); [C.load(f"data_{r}") for r in R.RUNS]
    with mp.get_context("fork").Pool(4) as P:
        out = P.map(job, JOBS, chunksize=1)
    json.dump(out, open(os.path.join(C.OUT, "feats.json"), "w"))
    for r in out:
        for rg, v in r["regions"].items():
            print(f"\n== {r['run']} pT {r['bin']} region {rg}: N_sig={v['N_sig']:.0f} (okno {v['n_win']}, k={v['k']:.3f}); sama waga MC -> Rtop {v['mc_w_only']['Rtop']:.3f} rho70 {v['mc_w_only']['rho70']:.3f} rho50 {v['mc_w_only']['rho50']:.3f}")
            fs = sorted(v["feats"], key=lambda f: -v["feats"][f]["chi2"] / max(v["feats"][f]["ndf"], 1))
            for f in fs[:9]:
                x = v["feats"][f]
                print(f"   {f:20s} chi2/ndf={x['chi2']:6.1f}/{x['ndf']} z_rank={x['z_rank']:+5.1f} KS={x['ks']:.3f} (p={x['ks_p']:.3f})  domkn.: Rtop={x['closure']['Rtop']:.3f}±{x['closure_err']['Rtop']:.3f} rho70={x['closure']['rho70']:.3f}±{x['closure_err']['rho70']:.3f} rho50={x['closure']['rho50']:.3f}")
            print(f"   raking {v['pair']}: Rtop={v['rake']['Rtop']:.3f}±{v['rake_err']['Rtop']:.3f} rho70={v['rake']['rho70']:.3f}±{v['rake_err']['rho70']:.3f} rho50={v['rake']['rho50']:.3f}±{v['rake_err']['rho50']:.3f}")
