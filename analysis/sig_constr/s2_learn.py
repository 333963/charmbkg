#!/usr/bin/env python3
"""SIG_constr A2: krzywa uczenia — AUC i Z_proxy klasyfikatora w funkcji liczby zdarzeń sygnału MC N (tło: pasma Δm danych 9:MB, stałe).
Zadania: H (obszar PRE, 18 cech H) oraz I (obszar H06, 18 cech + izolacja + score). Test: sygnał MC i tło z foldu 0 (po zdarzeniach)."""
import sys, time
from sc_common import *
import sc_ml as ML

M = load("fd_nd", tj=False); D = load("data_427929", tj=False); cm = data_classes(D, 427929)["9:MB"]
NS = [100, 300, 1000, 3000, 10000, None]
res = {}
for task, sthr, feats in (("H", None, F18), ("I", 0.6, F18 + ISO + ["score"])):
    ms = mc_sig(M, sthr); mb = cm & base_sel(D, sthr) & side_hi(D)
    Xs = np.column_stack([M[f][ms] for f in feats]); ks = M["kf"][ms]; pts = M["pt"][ms]
    Xb = np.column_stack([D[f][mb] for f in feats]); kb = D["kf"][mb]; ptb = D["pt"][mb]
    Xs_tr, Xs_te, Xb_tr, Xb_te = Xs[ks != 0], Xs[ks == 0], Xb[kb != 0], Xb[kb == 0]
    rng = np.random.default_rng(3)
    if len(Xb_tr) > 40000:
        Xb_tr = Xb_tr[rng.choice(len(Xb_tr), 40000, replace=False)]
    print(task, "sygnał MC trening/test", len(Xs_tr), len(Xs_te), "tło trening/test", len(Xb_tr), len(Xb_te), flush=True)
    rows = []; full = None
    for N in NS[::-1]:
        nseed = 1 if N is None else (2 if N >= 5000 else 4)
        for sd in range(nseed):
            r = np.random.default_rng(100 + sd)
            idx = np.arange(len(Xs_tr)) if N is None else r.choice(len(Xs_tr), min(N, len(Xs_tr)), replace=False)
            X = np.vstack([Xs_tr[idx], Xb_tr]); y = np.r_[np.ones(len(idx)), np.zeros(len(Xb_tr))]
            t0 = time.time()
            bst = ML.fit_xgb(X, y, np.ones(len(y)), seed=sd)
            ps, pb = ML.pred(bst, Xs_te), ML.pred(bst, Xb_te)
            row = {"N": len(idx), "seed": sd, "auc": ML.auc_w(ps, pb), "trees": int(bst.best_iteration + 1)}
            if N is None:
                full = {wp: float(np.mean(pb > np.quantile(ps, 1 - wp))) for wp in (0.95, 0.8, 0.5)}   # ε_B modelu pełnego przy WP
            for wp, eb in full.items():
                row[f"epsS@B{wp}"] = float(np.mean(ps > ML.thr_at_epsb(pb, eb)))
                row[f"Zrel{wp}"] = row[f"epsS@B{wp}"] / wp
            rows.append(row)
            print(f"  N={len(idx):6d} seed {sd} AUC {row['auc']:.4f} Zrel95/80/50 {row['Zrel0.95']:.3f} {row['Zrel0.8']:.3f} {row['Zrel0.5']:.3f} drzew {row['trees']} ({time.time()-t0:.0f}s)", flush=True)
    res[task] = {"rows": rows, "epsB_full": full, "n_sig_train": len(Xs_tr), "n_sig_test": len(Xs_te), "n_bkg_train": len(Xb_tr), "n_bkg_test": len(Xb_te), "feats": feats}
    # podsumowanie i N*(x%)
    summ = []
    for N in sorted({r["N"] for r in rows}):
        rr = [r for r in rows if r["N"] == N]
        summ.append({"N": N, **{k: float(np.mean([r[k] for r in rr])) for k in ("auc", "Zrel0.95", "Zrel0.8", "Zrel0.5")},
                     **{k + "_sd": float(np.std([r[k] for r in rr])) for k in ("auc", "Zrel0.8")}})
    res[task]["summary"] = summ
    for x in (0.02, 0.05, 0.10):
        ns = None
        for s in summ[::-1]:                      # najmniejsze N, od którego wszystkie większe N spełniają warunek
            if s["Zrel0.8"] >= 1 - x:
                ns = s["N"]
            else:
                break
        res[task][f"Nstar_{int(x*100)}"] = ns
    print(task, {k: v for k, v in res[task].items() if k.startswith("Nstar")}, flush=True)
jdump(res, "s2_learn.json")
