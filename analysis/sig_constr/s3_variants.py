#!/usr/bin/env python3
"""SIG_constr A3: warianty klasyfikatora na DANYCH (OOF, 5 foldów po zdarzeniach), per klasa triggera i zadanie:
  H: obszar PRE (bez cięcia score), pT 5–30 GeV, 18 cech H;   I: obszar H06 (score_H > 0,6), pT 3–30 GeV, 18 cech + izolacja + score.
Modele:  MC    — sygnał MC (przeważony w pT do sygnału danych) vs pasma Δm danych
         RW    — jak MC, sygnał MC dodatkowo przeważony w(x) = p_dane,syg(x)/p_MC,syg(x) (klasyfikator z wagami SBS; „SALAD dla sygnału”)
         DATA  — od zera na danych: okno Δm (etykieta 1) vs pasma Δm (0) — CWoLa, bez ujemnych wag
         STACK — dostrojenie: CWoLa na cechach + margines modelu MC jako dodatkowa cecha
         H     — sam score H (odniesienie; w zadaniu I = „bez izolacji”)
Miara: ε_S(dane) z odejmowania pasm Δm przy progu dającym to samo ε_B (pasmo górne OOF) co model MC w WP95/80/50; Z_rel = stosunek Z = ε_S/√ε_B."""
import sys, time
from sc_common import *
import sc_ml as ML

TASKS = {"H": dict(sthr=None, lo=5.0, hi=30.0, feats=F18, bins={"5-8": (5, 8), "8-30": (8, 30)}),
         "I": dict(sthr=0.6, lo=3.0, hi=30.0, feats=F18 + ISO + ["score"], bins=PTB)}
PT_RW = np.array([3, 3.5, 4, 4.5, 5, 6, 7, 8, 10, 12, 15, 20, 30.0])
WPS = (0.95, 0.8, 0.5)
MODELS = ["MC", "RW", "DATA", "STACK", "H"]
NB = 100


def dmreg(D, sel):
    dm = D["dm"][sel] * 1e3
    return (dm > IC.SIG[0]) & (dm < IC.SIG[1]), (dm > 140.6) & (dm < 143.4), (dm > 148.5) & (dm < 165.0)


def pt_weights(ptM, ptD, wD):
    hd = np.histogram(ptD, PT_RW, weights=wD)[0]; hm = np.histogram(ptM, PT_RW)[0].astype(float)
    r = np.clip(hd, 0, None) / np.maximum(hm, 1)
    w = r[np.clip(np.digitize(ptM, PT_RW) - 1, 0, len(PT_RW) - 2)]
    return w * len(w) / max(w.sum(), 1e-9)


def r2(x, y):
    if len(x) < 10 or np.std(x) == 0 or np.std(y) == 0:
        return 0.0
    return float(np.corrcoef(x, y)[0, 1] ** 2)


def run(D, cm, cn, M, task):
    T = TASKS[task]; feats = T["feats"]
    sel = cm & base_sel(D, T["sthr"], T["lo"], T["hi"])
    X = np.column_stack([D[f][sel] for f in feats]).astype(np.float32)
    kf = D["kf"][sel]; pt = D["pt"][sel]; dm = D["dm"][sel]; ek = D["evkey"][sel]; mm = D["m"][sel]
    win, slo, shi = dmreg(D, sel); side = slo | shi
    wsb = sbs_weights(D, sel)                                    # wagi sygnału (bez podziału pass/fail) — do pT i do RW
    ms = mc_sig(M, T["sthr"], T["lo"], T["hi"])
    XM = np.column_stack([M[f][ms] for f in feats]).astype(np.float32); kM = M["kf"][ms]; ptM = M["pt"][ms]
    wpt = pt_weights(ptM, pt, wsb)
    n = len(pt)
    S = {m: np.full(n, np.nan) for m in MODELS}; SM = {m: np.full(len(ptM), np.nan) for m in ("MC", "RW", "H")}
    S["H"] = D["score"][sel].astype(float); SM["H"] = M["score"][ms].astype(float)
    WRW = np.ones(len(ptM))
    rng = np.random.default_rng(11)
    t0 = time.time()
    for k in range(5):
        tr, te = kf != k, kf == k; trM, teM = kM != k, kM == k
        ib = np.nonzero(tr & side)[0]
        if len(ib) > 40000:
            ib = rng.choice(ib, 40000, replace=False)
        Xb = X[ib]
        # MC
        Xt = np.vstack([XM[trM], Xb]); yt = np.r_[np.ones(trM.sum()), np.zeros(len(ib))]
        b_mc = ML.fit_xgb(Xt, yt, np.r_[wpt[trM], np.ones(len(ib))], seed=k)
        S["MC"][te] = ML.pred(b_mc, X[te]); SM["MC"][teM] = ML.pred(b_mc, XM[teM])
        # RW: stosunek gęstości sygnału dane/MC (wagi SBS: okno +1, pasma −r), zastosowany do MC
        g = tr & (wsb != 0)
        rat = ML.fit_ratio(X[g], wsb[g], XM[trM], wpt[trM], seed=k)
        w_tr = ML.ratio_pred(rat, XM[trM]); WRW[teM] = ML.ratio_pred(rat, XM[teM])
        b_rw = ML.fit_xgb(Xt, yt, np.r_[wpt[trM] * w_tr, np.ones(len(ib))], seed=k)
        S["RW"][te] = ML.pred(b_rw, X[te]); SM["RW"][teM] = ML.pred(b_rw, XM[teM])
        # DATA (CWoLa)
        iw = np.nonzero(tr & win)[0]
        Xc = np.vstack([X[iw], Xb]); yc = np.r_[np.ones(len(iw)), np.zeros(len(ib))]
        b_d = ML.fit_xgb(Xc, yc, np.ones(len(yc)), seed=k)
        S["DATA"][te] = ML.pred(b_d, X[te])
        # STACK
        mg_c = ML.pred(b_mc, Xc, margin=True); mg_te = ML.pred(b_mc, X[te], margin=True)
        b_s = ML.fit_xgb(np.column_stack([Xc, mg_c]), yc, np.ones(len(yc)), seed=k, par=dict(max_depth=3))
        S["STACK"][te] = ML.pred(b_s, np.column_stack([X[te], mg_te]))
        print(f"  [{cn} {task}] fold {k}: okno {len(iw)}, tło {len(ib)}, MC {int(trM.sum())}, drzew MC/RW/DATA/STACK "
              f"{b_mc.best_iteration+1}/{b_rw.best_iteration+1}/{b_d.best_iteration+1}/{b_s.best_iteration+1}, ratio it {rat[1]}, {time.time()-t0:.0f}s", flush=True)
    # progi per fold przy wspólnym ε_B (z modelu MC w WP)
    out = {"cls": cn, "task": task, "n": int(n), "n_win": int(win.sum()), "Ns_sbs": float(wsb.sum()), "neff_sbs": neff(wsb), "n_mc": int(len(ptM)),
           "rw_weights": {"mean": float(WRW.mean()), "sd": float(WRW.std()), "neff_frac": float(neff(WRW * wpt) / max(neff(wpt), 1e-9)), "q01": float(np.quantile(WRW, 0.01)), "q99": float(np.quantile(WRW, 0.99))},
           "wp": {}, "leak": {}}
    boots = list(boot_ev(ek, NB, seed=5))
    for wp in WPS:
        epsB = []
        for k in range(5):
            tM = np.quantile(SM["MC"][kM == k], 1 - wp)        # WP na sygnale MC (bez wag pT — definicja progu jak w H)
            epsB.append(np.mean(S["MC"][(kf == k) & shi] > tM))
        eb = float(np.mean(epsB))
        P = {}
        for m in MODELS:
            pas = np.zeros(n, bool)
            for k in range(5):
                f = kf == k
                pas[f] = S[m][f] > ML.thr_at_epsb(S[m][f & shi], eb)
            P[m] = pas
        res = {"epsB_target": eb}
        binsx = dict(T["bins"]); binsx["all"] = (T["lo"], T["hi"])
        for b, (lo, hi) in binsx.items():
            g = (pt >= lo) & (pt < hi)
            rb = {}
            W = {m: IC.sbs_w(dm[g], pt[g], P[m][g]) for m in MODELS}
            for m in MODELS:
                e, ee, T_, _ = IC.eff_w(W[m], P[m][g].astype(float))
                ebm = float(np.mean(P[m][g & shi])); elo = float(np.mean(P[m][g & slo])) if (g & slo).sum() else np.nan
                rb[m] = {"epsS": e, "epsS_err": ee, "Ns": T_, "epsB": ebm, "epsB_lo": elo, "Z": e / np.sqrt(max(ebm, 1e-9))}
                if m in SM:
                    gm = (ptM >= lo) & (ptM < hi)
                    pm = np.zeros(len(ptM), bool)
                    for k in range(5):
                        # ten sam próg co na danych w foldzie k
                        thr = ML.thr_at_epsb(S[m][(kf == k) & shi], eb)
                        pm[kM == k] = SM[m][kM == k] > thr
                    rb[m]["epsS_mc"] = float(np.average(pm[gm], weights=wpt[gm]))
                    rb[m]["epsS_mc_rw"] = float(np.average(pm[gm], weights=(wpt * WRW)[gm]))
            # bootstrap po zdarzeniach: Z_rel względem MC oraz względem H
            eg = [bw[g] for bw in boots]
            shg = shi[g]
            def zb(m, bw):
                num = (bw * W[m] * P[m][g]).sum() / (bw * W[m]).sum()
                ebb = (bw[shg] * P[m][g][shg]).sum() / bw[shg].sum()
                return num / np.sqrt(max(ebb, 1e-9)), num
            for ref in ("MC", "H"):
                for m in MODELS:
                    if m == ref:
                        continue
                    zr = np.array([zb(m, bw)[0] / zb(ref, bw)[0] for bw in eg])
                    zr = zr[np.isfinite(zr)]
                    rb[m][f"Zrel_vs_{ref}"] = rb[m]["Z"] / rb[ref]["Z"] if rb[ref]["Z"] > 0 else np.nan
                    rb[m][f"Zrel_vs_{ref}_err"] = float(0.5 * (np.quantile(zr, 0.84) - np.quantile(zr, 0.16))) if len(zr) > 20 else np.nan
            res[b] = rb
        out["wp"][str(wp)] = res
        r = res["all"]
        print(f"  [{cn} {task}] WP{wp}: ε_B {eb:.4f} | " + " ".join(
            f"{m}: εS {r[m]['epsS']:.3f}±{r[m]['epsS_err']:.3f}" + (f" Zrel {r[m]['Zrel_vs_MC']:.3f}±{r[m]['Zrel_vs_MC_err']:.3f}" if m != "MC" else f" (MC {r[m]['epsS_mc']:.3f})") for m in MODELS), flush=True)
    for m in MODELS:
        out["leak"][m] = {"R2_dm_bkg": r2(S[m][shi], dm[shi]), "R2_m_bkg": r2(S[m][shi], mm[shi]), "R2_m_side_all": r2(S[m][side], mm[side])}
    np.savez_compressed(os.path.join(OUT, f"s3_scores_{cn.replace(':', '').replace('&', '_').replace('!', 'n')}_{task}.npz"),
                        **{"s_" + m: S[m].astype(np.float32) for m in MODELS}, dm=dm, pt=pt, kf=kf, evkey=ek, wsb=wsb, ptcone40=D["ptcone40"][sel],
                        **{"mc_" + m: SM[m].astype(np.float32) for m in SM}, mc_pt=ptM, mc_wpt=wpt, mc_wrw=WRW, mc_kf=kM, mc_ptcone40=M["ptcone40"][ms])
    return out


if __name__ == "__main__":
    want = sys.argv[1:] or ["9:MB", "8:J12only"]
    M = load("fd_nd", tj=False)
    res = []
    for name, run_ in (("data_427929", 427929), ("data_428770", 428770)):
        D = None
        for cn in want:
            if cn[0] != {427929: '9', 428770: '8'}[run_]:
                continue
            if D is None:
                D = load(name, tj=False); cls = data_classes(D, run_)
            for task in ("I", "H"):
                res.append(run(D, cls[cn], cn, M, task))
                jdump(res, "s3_variants_" + "_".join(w.replace(":", "").replace("&", "_").replace("!", "n") for w in want) + ".json")
