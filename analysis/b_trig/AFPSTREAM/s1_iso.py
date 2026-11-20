#!/usr/bin/env python3
"""Krok 1: N_sig D* i SF izolacji per klasa trygera (strumień AFP), bez wag i z przeważeniem pT × lead_pt × nSelBL (ISO_origin),
kontrola na tle kombinatorycznym (Δm > 150 MeV, bez score) dane/MC. Wynik: out/s1.json.
Estymatory: SBS (zawsze) + fit Δm pass/fail (gdy S_SBS ≥ 120 ≈ N_sig(fit) ≥ 150)."""
from multiprocessing import Pool
import numpy as np
import afp_common as A
from iso_common import load_pairs, sbs_w, eff_w, fit_pf, C, DR
import s2_universal_rw as U

MWIN, STHR = 0.040, 0.6
CUTS = [("ptcone40", 1.0), ("ptcone40", 2.0), ("ptcone40", 4.0), ("rel40", 0.25), ("rel40", 0.5), ("ptcone100", 5.0), ("ptcone100", 10.0)]
MAIN = ("ptcone40", 2.0)
SCH = U.SCHEMES["pT×lead_pt×nSelBL"]
PTB = {"3-30": (3, 30), "3-5": (3, 5), "5-8": (5, 8), "8-30": (8, 30)}
CL8 = A.DISJ8 + ["AFPonly", "AFPsptrk4_any", "AFPnoalg_any", "AFPstream", "J12only", "MBTSbased", "other", "j20!J12", "all"]
CL9 = ["MB427929", "MB9&!j20&!J12", "MB9&j20&!J12", "MB9&AFPnoalg", "MB9&AFPsptrk4", "MB9&AFPsptrk4&!j20&!J12", "MB9&!AFP"]


def tpd(D):
    return {c: D["tp_" + c] == 1 for c in A.ALLCH}


def samples():
    D8 = load_pairs("bt_data_428770"); c8 = A.classes8(tpd(D8)); c8["AFPnoalg_any"] = c8["_an"]
    D9 = load_pairs("bt_data_427929"); c9 = A.classes9(tpd(D9), D9["e_grl"] == 1)
    S = {"8:" + k: (D8, c8[k]) for k in CL8}
    S.update({"9:" + k: (D9, c9[k]) for k in CL9})
    return D8, D9, S


def _fit(job):
    key, dp, df_, shape, interval = job
    return key, fit_pf(dp, df_, shape, interval=interval)


def bkg_sel(D, cm, lo, hi):
    return cm & (np.abs(D["m"] - C.M0_D0) < MWIN) & (D["pt"] >= lo) & (D["pt"] < hi) & (D["dm"] > 0.150) & (D["dm"] < 0.170)


def main():
    M = load_pairs("bt_fd_nd")                       # sygnał MC
    B = load_pairs("bt_mb_nd")                       # tło MC (inkluzywne pary mb_nd)
    mf_all = U.mc_sel(M, MWIN, STHR)
    shape = DR.shape_from_mc(M["dm"][mf_all] * 1e3)
    D8, D9, S = samples()
    res = {"cuts": [f"{v}<{c:g}" for v, c in CUTS], "rows": {}, "mc": {}}
    jobs = []
    for pb, (lo, hi) in PTB.items():
        mm = U.mc_sel(M, MWIN, STHR, lo, hi)
        bm = bkg_sel(B, np.ones(len(B["m"]), bool), lo, hi)
        res["mc"][pb] = {"n_sig": int(mm.sum()), "n_bkg": int(bm.sum()),
                         "eps": {f"{v}<{c:g}": float((M[v][mm] < c).mean()) for v, c in CUTS},
                         "f_bkg": {f"{v}<{c:g}": float((B[v][bm] < c).mean()) for v, c in CUTS}}
        for name, (D, cm) in S.items():
            md = U.base(D, MWIN, STHR, lo, hi) & cm
            bb = bkg_sel(D, cm, lo, hi)
            row = {"n_pairs_cls": int(cm.sum()), "ncand": int(md.sum()), "n_bkg": int(bb.sum()), "cuts": {}}
            # N_sig: SBS bez podziału pass/fail
            if md.sum() > 20:
                w0 = sbs_w(D["dm"][md], D["pt"][md], np.ones(int(md.sum()), bool))
                row["S"] = float(w0.sum()); row["S_err"] = float(np.sqrt((w0 ** 2).sum()))
                row["nwin"] = int((w0 == 1).sum())
            else:
                row["S"], row["S_err"], row["nwin"] = 0.0, 0.0, 0
            for v, c in CUTS:
                ck = f"{v}<{c:g}"
                pas = D[v] < c; pasM = M[v] < c; pasB = B[v] < c
                r = {}
                if md.sum() > 20:
                    w = sbs_w(D["dm"][md], D["pt"][md], pas[md])
                    e, er, T, eT = eff_w(w, pas[md].astype(float))
                    r.update({"eps": e, "err": er, "S": T, "S_err": eT})
                    emc = res["mc"][pb]["eps"][ck]
                    r["mc_eps"] = emc
                    r["SF"] = e / emc; r["SF_err"] = float(np.hypot(er / emc, e / emc * np.sqrt((1 - emc) / emc / mm.sum())))
                    rw = U.sf_rw(D, md, pas, M, mm, pasM, SCH)
                    if rw:
                        r["rw"] = rw
                    if T >= 120 and (pb == "3-30" or (v, c) == MAIN):
                        jobs.append(((name, pb, ck), D["dm"][md & pas], D["dm"][md & ~pas], shape, (v, c) == MAIN and pb == "3-30"))
                # tło: frakcja izolowana dane, MC bez wag, MC w komórkach (pT×lead_pt×nSelBL) rozkładu danych
                nb = int(bb.sum())
                if nb >= 30:
                    f = float(pas[bb].mean()); fe = float(np.sqrt(max(f * (1 - f), 1e-9) / nb))
                    fmc = res["mc"][pb]["f_bkg"][ck]
                    eps_c, var_c, fb = U.mc_eps_for(D, bb, B, bm, pasB, SCH)
                    frw = float(eps_c.mean())
                    # błąd MC komórek
                    key = np.round(eps_c, 12); u, inv = np.unique(key, return_inverse=True)
                    wc = np.bincount(inv).astype(float); vc = np.zeros(len(u)); vc[inv] = var_c
                    frw_e = float(np.sqrt((wc ** 2 * vc).sum()) / nb)
                    r["bkg"] = {"n": nb, "f": f, "err": fe, "f_mc": fmc, "R": f / fmc, "R_err": fe / fmc,
                                "f_mc_rw": frw, "R_rw": f / frw, "R_rw_err": float(f / frw * np.hypot(fe / f if f > 0 else 0, frw_e / frw)),
                                "fallback": float(fb.mean())}
                row["cuts"][ck] = r
            res["rows"][f"{name}|{pb}"] = row
    print("fity:", len(jobs), flush=True)
    with Pool(4) as P:
        for (name, pb, ck), f in P.imap_unordered(_fit, jobs):
            r = res["rows"][f"{name}|{pb}"]["cuts"][ck]
            r["fit"] = f
            emc = r["mc_eps"]
            if np.isfinite(f["eps"]):
                r["SF_fit"] = f["eps"] / emc; r["SF_fit_err"] = f["err"] / emc
    A.jdump(res, f"{A.OUT}/s1.json")
    print("zapisano out/s1.json")


if __name__ == "__main__":
    main()
