#!/usr/bin/env python3
"""SIG_constr B(0): czy TRKJET2 odpowiada już na pytanie embeddingu? Z tj2_*.npz (stożki prostopadłe / losowy wokół kandydatów D0←D*):
ε(Σ pT w stożku R=0,4 < 2 GeV) dla (a) własnego stożka D0 [sygnał, SBS], (b) stożka prostopadłego i losowego w TYCH SAMYCH zdarzeniach sygnałowych (SBS)
= „embedding w kierunku prostopadłym/losowym w zdarzenia z prawdziwym D0”, (c) stożka prostopadłego/losowego dla WSZYSTKICH kandydatów klasy (głównie tło
kombinatoryczne = zdarzenia-nośniki bez prawdziwego D0). Δ_corr = ε_perp − ε_własny = aktywność skorelowana z D0. MC fd_nd: to samo na prawdzie."""
from sc_common import *

R2 = os.path.join(CW, "B_trig", "TRKJET2", "runs", "lxplus")
COLS = ["m", "pt", "eta", "dm", "score", "truth", "truth_ds", "event", "file_idx", "run", "e_grl", "ptcone40", "pis_pt", "p_s2_p", "p_s2_m", "r_s2", "r_ok", "v_s40_p2", "v_s40_p2x",
        "p_jet_p", "p_jet_m", "nmatch", "nd_lo", "same_entry", "e_nSelBL", "e_lead_pt", "e_ntrk500", "e_actual_mu"]


def ld(name):
    z = np.load(os.path.join(R2, f"tj2_{name}.npz"))
    d = {k: z[k] for k in COLS if k in z.files}
    for k in z.files:
        if k.startswith("e_tp_"):
            d[k[2:]] = z[k]
    par = (d["event"] % 2 == 0)
    d["perp1"] = np.where(par, d["p_s2_p"], d["p_s2_m"]); d["perp_avg"] = 0.5 * (d["p_s2_p"] + d["p_s2_m"])
    d["own"] = d["v_s40_p2x"]                                  # własny stożek bez wolnego pionu (porównywalny ze stożkiem bez D*)
    d["ok"] = (d["nmatch"] + d["nd_lo"] == 2) & (d["same_entry"] == 1)
    return d


rows = []
M = ld("fd_nd")
ms = np.isin(M["truth"], [0, 1]) & np.isin(M["truth_ds"], [0, 1]) & (np.abs(M["m"] - C.M0_D0) < MWIN) & (M["dm"] > 0.1435) & (M["dm"] < 0.1475) & (M["score"] > 0.6) & (M["pt"] >= 3) & (M["pt"] < 30) & M["ok"]
for name, run in (("data_427929", 427929), ("data_428770", 428770)):
    D = ld(name); cls = data_classes(D, run)
    for cn in (["9:MB", "9:MB&j20&!J12", "9:MB&!j20&!J12"] if run == 427929 else ["8:J12only", "8:j20&!J12", "8:MBTSbased", "8:other"]):
        for b, (lo, hi) in [("3-30", (3.0, 30.0))] + list(PTB.items()):
            sel = cls[cn] & base_sel(D, 0.6, lo, hi) & D["ok"]
            r = {"cls": cn, "bin": b}
            for v in ("own", "ptcone40", "perp1", "r_s2"):
                g = sel & ((D["r_ok"] == 1) if v == "r_s2" else True)
                pas = D[v][g] < 2.0
                w = IC.sbs_w(D["dm"][g], D["pt"][g], pas)
                e, ee, T, _ = IC.eff_w(w, pas.astype(float))
                r[v], r[v + "_err"], r["Ns"] = e, ee, T
                w0 = IC.sbs_w(D["dm"][g], D["pt"][g], np.ones(g.sum(), bool))
                r["mean_" + v] = float((w0 * D[v][g]).sum() / w0.sum())
            # wszystkie kandydaci klasy bez cięcia score (nośniki bez prawdziwego D0): pasmo górne Δm
            gb = cls[cn] & base_sel(D, None, lo, hi) & side_hi(D) & D["ok"]
            for v in ("perp1", "r_s2", "own"):
                g = gb & ((D["r_ok"] == 1) if v == "r_s2" else True)
                r["bkg_" + v] = float(np.mean(D[v][g] < 2.0)); r["bkg_n"] = int(g.sum())
            # MC przeważone w pT
            gm = ms & (M["pt"] >= lo) & (M["pt"] < hi)
            from s3_variants import pt_weights
            ws = IC.sbs_w(D["dm"][sel], D["pt"][sel], np.ones(sel.sum(), bool))
            wpt = pt_weights(M["pt"][gm], D["pt"][sel], ws)
            for v in ("own", "perp1", "r_s2"):
                r["mc_" + v] = float(np.average(M[v][gm] < 2.0, weights=wpt))
            r["dcorr"] = r["perp1"] - r["own"]; r["dcorr_err"] = float(np.hypot(r["perp1_err"], r["own_err"]))
            r["dcorr_mc"] = r["mc_perp1"] - r["mc_own"]
            rows.append(r)
            print(f"{cn:16s} {b:5s} N_sig {r['Ns']:5.0f} | własny(bez π_s) {r['own']:.3f}±{r['own_err']:.3f} (MC {r['mc_own']:.3f}) | prostop. {r['perp1']:.3f}±{r['perp1_err']:.3f} (MC {r['mc_perp1']:.3f}) | "
                  f"losowy {r['r_s2']:.3f}±{r['r_s2_err']:.3f} (MC {r['mc_r_s2']:.3f}) | nośniki-tło: prostop. {r['bkg_perp1']:.3f} losowy {r['bkg_r_s2']:.3f} własny {r['bkg_own']:.3f} | Δ_corr {r['dcorr']:.3f}±{r['dcorr_err']:.3f} (MC {r['dcorr_mc']:.3f})", flush=True)
jdump({"rows": rows}, "s6_perp.json")
