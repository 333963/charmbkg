"""DLZ_dplus/DATA: kontrole na MC (bez danych): (G3) R_sculpt na mb truth≥8 w oknie sygnału (liniowy fit ε_B(m) w |Δm|>70 MeV),
(G5) closure fitu: mb_nd (cała próbka, dowolny truth) po cięciu WP → fit masy → N(WP)/N(WP90) vs zliczenie prawdy (truth 0/1),
(G1) AUC trening vs OOF z meta.json. Wyjście: out/mc_checks.json.  Użycie: python mc_checks.py"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fitlib as FL  # noqa: E402
from analyze import BASE, BINS, WPS, TAGS, SEED2, thresholds, bn  # noqa: E402

b = np.load(os.path.join(BASE, "oof", "base.npz"))
res = {"sculpt_mc": {}, "closure": {}, "auc": {}, "zproxy_mcS_dataSB": {}}
m, pt, tr = b["mb_m"].astype(float), b["mb_pt"].astype(float), b["mb_truth"]
ibm = np.digitize(pt, [3, 5, 8, 30]) - 1
e20 = np.linspace(FL.LO, FL.HI, 26)
c20 = 0.5 * (e20[1:] + e20[:-1])
for t in TAGS + SEED2:
    f = os.path.join(BASE, "oof", t + ".npz")
    if not os.path.exists(f):
        continue
    z = np.load(f)
    thr = thresholds(z["S"], b["S_w"].astype(float), b["S_pt"].astype(float))
    mf = os.path.join(BASE, "models", t, "meta.json")
    if os.path.exists(mf):
        me = json.load(open(mf))
        res["auc"][t] = dict(train=me["auc_train"], oof=me["auc_oof"], diff=float(np.mean(me["auc_train"]) - np.mean(me["auc_oof"])), best_it=me["best_it"],
                             fold_violations=me["fold_violations"], importance=me["importance_fold0"])
    bk = (tr >= 8) & (tr < 100)
    for wp in (0.8, 0.7):
        ok = z["mb"] >= np.array([thr[(lo, hi, wp)] for lo, hi in BINS])[ibm]
        n0, _ = np.histogram(m[bk], e20)
        n1, _ = np.histogram(m[bk & ok], e20)
        e = n1 / np.maximum(n0, 1)
        err = np.sqrt(np.maximum(e * (1 - e), 1e-12) / np.maximum(n0, 1))
        sb = (np.abs(c20 - FL.M0) > 0.07) & (n0 > 20)
        A = np.c_[np.ones(sb.sum()), c20[sb] - FL.M0]
        W = 1 / err[sb] ** 2
        cov = np.linalg.inv(A.T @ (A * W[:, None]))
        beta = cov @ (A.T @ (W * e[sb]))
        reg = bk & (np.abs(m - FL.M0) < FL.WIN)
        obs = ok[reg].mean()
        x = np.c_[np.ones(reg.sum()), m[reg] - FL.M0]
        pred = float((x @ beta).mean())
        vp = float(np.mean(x, 0) @ cov @ np.mean(x, 0))
        res["sculpt_mc"][f"{t}|{wp}"] = dict(R=float(obs / pred), R_err=float(obs / pred * np.hypot(np.sqrt(obs * (1 - obs) / reg.sum()) / obs, np.sqrt(vp) / pred)),
                                             slope=float(beta[1] * 0.1 / beta[0]), slope_err=float(np.sqrt(cov[1, 1]) * 0.1 / beta[0]),
                                             chi2_lin=float(np.sum(W * (e[sb] - A @ beta) ** 2)), ndf_lin=int(sb.sum() - 2), eps=e.tolist(), eps_err=err.tolist())
    # closure na mb (pT-integralnie i per bin)
    sig = np.isin(tr, [0, 1])
    for b_ in BINS + [(3, 30)]:
        inb = (pt >= b_[0]) & (pt < b_[1])
        fits, tru = {}, {}
        for wp in WPS:
            ok = z["mb"] >= np.array([thr[(lo, hi, wp)] for lo, hi in BINS])[ibm]
            sS = (b["S_pt"] >= b_[0]) & (b["S_pt"] < b_[1]) & (z["S"] >= np.array([thr[(lo, hi, wp)] for lo, hi in BINS])[np.digitize(b["S_pt"], [3, 5, 8, 30]) - 1])
            sh = FL.mc_shape(b["S_m"][sS].astype(float), b["S_w"][sS].astype(float))
            h, _ = np.histogram(m[inb & ok], FL.E)
            fits[wp] = FL.fit(h, sh)
            tru[wp] = int((inb & ok & sig).sum())
        for wp in WPS:
            if wp >= 0.9:
                continue
            ef = fits[wp]["N"] / max(fits[0.9]["N"], 1e-9)
            ee = ef * np.hypot(fits[wp]["N_err"] / max(fits[wp]["N"], 1e-9), 0)  # błąd tylko licznika (podzbiór) — przybliżenie dwumianowe niżej
            et = tru[wp] / max(tru[0.9], 1)
            eb = np.sqrt(max(ef * (1 - ef), 0) / max(fits[0.9]["N"], 1)) if ef < 1 else ee
            res["closure"][f"{t}|{bn(b_)}|{wp}"] = dict(eps_fit=float(ef), eps_truth=float(et), ratio=float(ef / max(et, 1e-9)),
                                                        ratio_err=float(np.hypot(fits[wp]["N_err"] / max(fits[wp]["N"], 1e-9), fits[0.9]["N_err"] / max(fits[0.9]["N"], 1e-9)) * ef / max(et, 1e-9)),
                                                        ratio_err_binom=float(eb / max(et, 1e-9)), N_fit=fits[wp]["N"], N_fit_err=fits[wp]["N_err"], N_fit90=fits[0.9]["N"], N_truth=tru[wp], N_truth90=tru[0.9])
os.makedirs(os.path.join(BASE, "out"), exist_ok=True)
json.dump(res, open(os.path.join(BASE, "out", "mc_checks.json"), "w"), indent=1)
for k, v in res["sculpt_mc"].items():
    print("sculpt MC", k, round(v["R"], 3), round(v["R_err"], 3), round(v["slope"], 3))
for k, v in res["closure"].items():
    if "|3-30|" in k:
        print("closure", k, round(v["ratio"], 3), round(v["ratio_err"], 3))
for k, v in res["auc"].items():
    print("AUC", k, round(v["diff"], 4))
