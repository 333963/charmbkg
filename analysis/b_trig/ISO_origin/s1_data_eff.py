#!/usr/bin/env python3
"""ISO_origin krok 1: ε_iso w DANYCH (bez MC) per klasa triggera / run / μ / LB / npu / nSelBL / lead_pt / pT / z,
kryminalistyka bitów triggera (strumień physics_Main), pile-up na tle kombinatorycznym. Wynik: out/s1.json."""
import sys
from multiprocessing import Pool
from iso_common import *

CUT = ("ptcone40", 2.0)
MWIN, STHR = 0.040, 0.6
VARS = {"pt": [3, 4, 5, 6.5, 8, 12, 30], "pt_ds": [3, 4, 5, 6.5, 8, 12, 32], "e_nSelBL": [0, 25, 35, 45, 60, 80, 1e9],
        "e_ntrk500": [0, 30, 45, 60, 80, 110, 1e9], "e_lead_pt": [0, 3, 4, 5, 6.5, 8, 12, 20, 1e9], "z": [0, 0.5, 0.8, 1.0, 1.2, 1.5, 1e9],
        "e_npu": [-0.5, 0.5, 1e9], "lbpar": [-0.5, 0.5, 1.5]}
MU = {427929: [0, 0.004, 0.0055, 0.008, 0.03], 428770: [0, 0.03, 0.045, 0.055, 0.07, 0.3]}
LBQ = {427929: None, 428770: [4, 26, 245, 275, 300, 324]}


def base(D, mwin=MWIN, sthr=STHR):
    return (np.abs(D["m"] - C.M0_D0) < mwin) & (D["score"] > sthr) & (D["pt"] >= 3) & (D["pt"] < 30)


def _fit(job):
    key, dp, df_, shape, interval = job
    return key, fit_pf(dp, df_, shape, interval=interval)


def measure(D, m, pas):
    w = sbs_w(D["dm"][m], D["pt"][m], pas[m])
    e, er, S, eS = eff_w(w, pas[m])
    return {"ncand": int(m.sum()), "S": S, "S_err": eS, "eps": e, "err": er}


def samples():
    D9 = load_pairs("bt_data_427929"); c9 = classes(D9, 427929)
    D8 = load_pairs("bt_data_428770"); c8 = classes(D8, 428770)
    for D in (D9, D8):
        D["lbpar"] = (D["e_lb"].astype(int) % 2).astype(float)
    tp9 = lambda c: D9["tp_" + c] == 1  # noqa
    mb = c9["MB427929"]; j9 = tp9("HLT_noalg_L1J12"); j20_9 = tp9("HLT_j20_L1MBTS_2"); s4_9 = tp9("HLT_mb_sptrk_pt4_L1MBTS_2")
    S = {"9:MB": (D9, mb), "9:MB&J12": (D9, mb & j9), "9:MB&j20&!J12": (D9, mb & j20_9 & ~j9), "9:MB&!j20&!J12": (D9, mb & ~j20_9 & ~j9),
         "9:MB&sptrk4&!j20&!J12": (D9, mb & s4_9 & ~j20_9 & ~j9), "9:MB&sptrk4": (D9, mb & s4_9),
         "8:all": (D8, c8["all"]), "8:J12only": (D8, c8["J12only"]), "8:J12any": (D8, c8["J12any"]), "8:MBTSbased": (D8, c8["MBTSbased"]),
         "8:j20&!J12": (D8, c8["j20_notJ12"]), "8:sptrk4&!j20&!J12": (D8, c8["sptrk4_notJ12"] & ~c8["j20MBTS"]),
         "8:other": (D8, c8["other"]), "8:unbMBTS(noalg_mb/mbts)": (D8, c8["unbMBTS"]), "8:MBTSsoft_only": (D8, c8["MBTSsoft_only"])}
    return D9, D8, S


def main():
    os.makedirs(os.path.join(HERE, "out"), exist_ok=True)
    M = load_pairs("bt_fd_nd")
    mf0 = np.isin(M["truth"], [0, 1]) & (np.abs(M["m"] - C.M0_D0) < MWIN) & (M["dm"] > 0.1435) & (M["dm"] < 0.1475) & (M["score"] > STHR) & (M["pt"] < 30)
    shape = DR.shape_from_mc(M["dm"][mf0] * 1e3)
    M["lbpar"] = np.zeros(len(M["m"]))
    D9, D8, S = samples()
    res = {"cut": f"{CUT[0]}<{CUT[1]}", "mc_eps_incl": float((M[CUT[0]][mf0] < CUT[1]).mean()), "mc_n": int(mf0.sum()), "classes": {}, "diff": {}, "mc_diff": {}}
    jobs = []
    PTB = {"3-30": (3, 30), "3-5": (3, 5), "5-8": (5, 8), "8-30": (8, 30)}
    # ---- klasy × bin pT: SBS + fit (profil) gdzie S>=150; tło: pasmo Δm>150 MeV
    for name, (D, cm) in S.items():
        b = base(D); pas = D[CUT[0]] < CUT[1]
        for pb, (lo, hi) in PTB.items():
            m = b & cm & (D["pt"] >= lo) & (D["pt"] < hi)
            r = measure(D, m, pas) if m.sum() > 20 else {"ncand": int(m.sum()), "S": 0.0, "S_err": 0.0, "eps": np.nan, "err": np.nan}
            # tło kombinatoryczne (Δm > 150 MeV): z cięciem score i bez
            for tag, bb in (("bkg", m & (D["dm"] > 0.150)),
                            ("bkg_noscore", cm & (np.abs(D["m"] - C.M0_D0) < MWIN) & (D["pt"] >= lo) & (D["pt"] < hi) & (D["dm"] > 0.150))):
                n = int(bb.sum()); f = float(pas[bb].mean()) if n else np.nan
                r[tag] = {"n": n, "f": f, "err": float(np.sqrt(max(f * (1 - f), 1e-9) / max(n, 1))) if n else np.nan}
            mm = mf0 & (M["pt"] >= lo) & (M["pt"] < hi)
            r["mc_eps"] = float((M[CUT[0]][mm] < CUT[1]).mean())
            res["classes"][f"{name}|{pb}"] = r
            if r["S"] >= 150:
                jobs.append(((f"{name}|{pb}", "cls"), D["dm"][m & pas], D["dm"][m & ~pas], shape, True))
    # ---- różniczkowo dla głównych próbek
    main_s = ["9:MB", "8:J12only", "8:MBTSbased", "8:all"]
    for name in main_s:
        D, cm = S[name]; run = 427929 if name[0] == "9" else 428770
        b = base(D) & cm; pas = D[CUT[0]] < CUT[1]
        vv = dict(VARS); vv["e_actual_mu"] = MU[run]
        if run == 428770:
            vv["e_lb"] = LBQ[run]
        else:
            q = np.quantile(D["e_lb"][b], [0, 0.25, 0.5, 0.75, 1]); q[-1] += 1; vv["e_lb"] = list(q)
        for var, ed in vv.items():
            for lo, hi in zip(ed[:-1], ed[1:]):
                m = b & (D[var] >= lo) & (D[var] < hi)
                key = f"{name}|{var}|{lo:g}|{hi:g}"
                if m.sum() < 30:
                    continue
                r = measure(D, m, pas)
                bb = cm & (np.abs(D["m"] - C.M0_D0) < MWIN) & (D["pt"] >= 3) & (D["pt"] < 30) & (D["dm"] > 0.150) & (D[var] >= lo) & (D[var] < hi)
                r["bkg_noscore"] = {"n": int(bb.sum()), "f": float(pas[bb].mean()) if bb.sum() else np.nan}
                res["diff"][key] = r
                if r["S"] >= 150:
                    jobs.append(((key, "diff"), D["dm"][m & pas], D["dm"][m & ~pas], shape, False))
    for var, ed in VARS.items():
        for lo, hi in zip(ed[:-1], ed[1:]):
            mm = mf0 & (M[var] >= lo) & (M[var] < hi)
            if mm.sum() > 20:
                e = float((M[CUT[0]][mm] < CUT[1]).mean())
                res["mc_diff"][f"{var}|{lo:g}|{hi:g}"] = {"n": int(mm.sum()), "eps": e, "err": float(np.sqrt(e * (1 - e) / mm.sum()))}
    print("fity:", len(jobs), flush=True)
    with Pool(6) as P:
        for (key, kind), f in P.imap_unordered(_fit, jobs):
            (res["classes"] if kind == "cls" else res["diff"])[key]["fit"] = f
    # ---- pile-up na tle (bez cięcia score, pT 3–5; wysoka statystyka)
    pu = {}
    for run, D, cm in ((427929, D9, np.ones(len(D9["m"]), bool)), (428770, D8, S["8:J12only"][1]), ("428770j20", D8, S["8:j20&!J12"][1]),
                       ("427929MB", D9, S["9:MB"][1])):
        bb = cm & (np.abs(D["m"] - C.M0_D0) < MWIN) & (D["pt"] < 5) & (D["dm"] > 0.150); pas = D[CUT[0]] < CUT[1]
        for var, ed in (("e_actual_mu", [0, 0.004, 0.0055, 0.008, 0.03, 0.045, 0.055, 0.07, 0.3, 2, 15, 25, 40]), ("e_npu", [-0.5, 0.5, 1.5, 2.5, 100])):
            for lo, hi in zip(ed[:-1], ed[1:]):
                m = bb & (D[var] >= lo) & (D[var] < hi)
                if m.sum() < 50:
                    continue
                f = float(pas[m].mean())
                pu[f"{run}|{var}|{lo:g}|{hi:g}"] = {"n": int(m.sum()), "f": f, "err": float(np.sqrt(f * (1 - f) / m.sum())),
                                                   "mu_mean": float(D["e_actual_mu"][m].mean()), "nSelBL_med": float(np.median(D["e_nSelBL"][m])),
                                                   "grl_frac": float((D["e_grl"][m] == 1).mean())}
    res["pileup_bkg"] = pu
    # ---- kryminalistyka bitów: tabele zdarzeń
    fo = {}
    z8 = np.load(os.path.join(R, "trg_data_428770.npz")); z9 = np.load(os.path.join(R, "trg_data_427929.npz"))
    tp8 = {k[7:]: z8[k] == 1 for k in z8.files if k.startswith("tp_")}
    ok9 = z9["grl"] == 1
    tp9 = {k[7:]: (z9[k] == 1) & ok9 for k in z9.files if k.startswith("tp_")}
    names = [k for k in tp8 if tp8[k].sum() > 0]
    anyc = np.any([tp8[k] for k in names], axis=0)
    fo["n_events_428770"] = int(len(anyc)); fo["none_of_stored"] = int((~anyc).sum())
    hard9 = ["noalg_L1J12", "j20_L1MBTS_2", "mb_sptrk_pt4_L1MBTS_2", "mb_sptrk_pt6_L1MBTS_2", "mb_sptrk_pt8_L1MBTS_2"]
    fo["chains"] = {}
    for k in names:
        oth = np.any([tp8[j] for j in names if j != k], axis=0)
        al = tp8[k] & ~oth
        d = {"n8": int(tp8[k].sum()), "alone8": int(al.sum()), "nSelBL_mean8": float(z8["nSelBL"][tp8[k]].mean()), "lead_med8": float(np.median(z8["lead_pt"][tp8[k]])),
             "ps8": float(np.median(z8["ps_HLT_" + k][tp8[k]]))}
        if k in tp9 and tp9[k].sum():
            o9 = np.any([tp9[j] for j in hard9 if j != k], axis=0)
            d.update({"n9": int(tp9[k].sum()), "frac_wo_other_hard9": float((tp9[k] & ~o9).sum() / tp9[k].sum()),
                      "nSelBL_mean9": float(z9["nSelBL"][tp9[k]].mean()), "lead_med9": float(np.median(z9["lead_pt"][tp9[k]]))})
        fo["chains"][k] = d
    j8 = tp8["noalg_L1J12"]
    fo["P_in_J12_events_428770"] = {k: float(tp8[k][j8].mean()) for k in ("j20_L1MBTS_2", "mb_sptrk_pt4_L1MBTS_2", "noalg_mb_L1MBTS_2")}
    fo["mean_inv_ps_j20"] = float((1 / z8["ps_HLT_j20_L1MBTS_2"]).mean())
    # domknięcie emulacji: oryginał (niezależne prescale’y, sptrk4 pisze do strumienia) vs model strumienia Main
    import btrig_emul as EM
    ref = tp9["mb_mbts_L1MBTS_1"]
    i9, j9_ = EM._bin2(z9["nSelBL"], z9["lead_pt"]); i8, j8_ = EM._bin2(z8["nSelBL"], z8["lead_pt"])
    shp = (len(EM.NSEL_EDGES) - 1, len(EM.LEAD_EDGES) - 1)

    def cellmean(mask):
        n = np.zeros(shp); d = np.zeros(shp)
        np.add.at(d, (i9[ref], j9_[ref]), 1.0); np.add.at(n, (i9[ref & mask], j9_[ref & mask]), 1.0)
        return n / np.maximum(d, 1), d
    J, j20, s4 = tp9["noalg_L1J12"], tp9["j20_L1MBTS_2"], tp9["mb_sptrk_pt4_L1MBTS_2"]
    eJ, nref = cellmean(J); eM_, _ = cellmean(j20 | (s4 & J)); eJM, _ = cellmean(J & (j20 | s4))
    e20, _ = cellmean(j20); e4, _ = cellmean(s4)
    q = 1.0 / z8["ps_HLT_j20_L1MBTS_2"]
    cl8 = C.classify({k: z8[k] for k in z8.files if k.startswith("tp_")}, 428770)
    sel = cl8["J12only"] | cl8["MBTSbased"]
    pm_new = q * eM_[i8, j8_]; pj_new = eJ[i8, j8_] - q * eJM[i8, j8_]
    nm = 1 - (1 - q * e20[i8, j8_]) * (1 - q * e4[i8, j8_]); pm_old = nm; pj_old = eJ[i8, j8_] * (1 - nm)
    out = {}
    for tag, pm, pj in (("oryginal_niezalezne_PS", pm_old, pj_old), ("model_strumienia_Main", pm_new, pj_new)):
        pred = pm / np.clip(pm + pj, 1e-12, None); x2 = 0.0; nd = 0; mx = 0.0
        for a in range(shp[0]):
            for b_ in range(shp[1]):
                c = sel & (i8 == a) & (j8_ == b_)
                if c.sum() < 50 or nref[a, b_] < 200:
                    continue
                obs = cl8["MBTSbased"][c].mean(); pr = pred[c].mean()
                x2 += ((obs - pr) / (np.sqrt(obs * (1 - obs) / c.sum()) + 1e-6)) ** 2; nd += 1
                if obs > 0:
                    mx = max(mx, abs(pr / obs - 1))
        out[tag] = {"chi2": float(x2), "ndf": nd, "max_rel_dev": float(mx), "pred_frac_MBTS": float(pred[sel].mean()), "obs_frac_MBTS": float(cl8["MBTSbased"][sel].mean())}
    fo["closure"] = out
    res["forensic"] = fo
    json.dump(res, open(os.path.join(HERE, "out", "s1.json"), "w"), indent=1, default=lambda o: None if (isinstance(o, float) and not np.isfinite(o)) else float(o))
    print("zapisano out/s1.json")


if __name__ == "__main__":
    main()
