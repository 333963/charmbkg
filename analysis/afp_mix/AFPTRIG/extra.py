#!/usr/bin/env python3
"""AFPTRIG/extra.py — krok 4 (czego metoda nie daje): (a) test lustrzany w strumieniu AFP (zależność prawdopodobieństwa protonu
przypadkowego od ξ_cal strony przeciwnej), (b) widma ξ protonów: „niemożliwe” w strumieniu AFP vs nie-AFP vs wszystkie tagi
(obciążenie trygerem zależne od ξ, różnica widma sygnał/przypadkowe), (c) topologia trygera (proton po stronie przeciwnej),
(d) nadwyżka danych nad szablonem tuż nad cięciem (przeciek sygnału przez rozdzielczość ξ_cal). → extra.json"""
import json, numpy as np, lib
from analysis import Run, NCLS, nominal_pool, POOL_THR
OUT = {}
EO = np.array([0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.2, 0.3, 0.5, 10.0])
EX = np.array([0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10001])
for run in lib.RUNS:
    R = Run(run); d = R.d; O = OUT[str(run)] = {}
    for stream in ("afp_any", "pt4_L1AFP", "nonAFP"):
        pool = nominal_pool(run, stream); Os = O[stream] = {}
        for s in lib.SIDES:
            o = "C" if s == "A" else "A"; r = Os[s] = {}
            # (a) lustrzany: zdarzenia z D, ξ_cal,s > 0,14, bez protonu po stronie o; dane = tag_s; pula = np_s == 0
            g = R.groups(s, stream, 8, "lb"); cell = g * NCLS + R.cls; nc = int(cell.max()) + 1
            base = R.has & (d[f"xical_{s}"] > POOL_THR) & (d[f"np_{o}"] == 0)
            dm = base & R.S[stream] & d[f"tag_{s}"]; pm = base & R.S[pool] & (d[f"np_{s}"] == 0)
            bo_d = np.digitize(d[f"xical_{o}"][dm], EO) - 1; bo_p = np.digitize(d[f"xical_{o}"][pm], EO) - 1
            nb = len(EO) - 1; norm_bins = np.arange(nb) >= 7   # ξ_cal,o > 0,14
            Hd = np.zeros((nc, nb)); Hp = np.zeros((nc, nb))
            np.add.at(Hd, (cell[dm], bo_d), 1); np.add.at(Hp, (cell[pm], bo_p), 1)
            Nn = Hd[:, norm_bins].sum(1); Pn = Hp[:, norm_bins].sum(1); ok = Pn > 0
            pred = (Hp[ok] * (Nn[ok] / Pn[ok])[:, None]).sum(0); obs = Hd[ok].sum(0)
            # błąd pred: Poisson N_norm per komórka (dominujący) ⊕ pula
            vpred = ((Hp[ok] / Pn[ok][:, None]) ** 2 * Nn[ok][:, None]).sum(0) + (Hp[ok] * ((Nn[ok] / Pn[ok]) ** 2)[:, None]).sum(0)
            r["mirror"] = {"edges": EO.tolist(), "obs": obs.tolist(), "pred": pred.tolist(), "pred_err": np.sqrt(vpred).tolist()}
            lo = ~norm_bins
            ol, pl = obs[lo].sum(), pred[lo].sum(); el = np.sqrt(vpred[lo].sum())   # (przybliżenie: bez korelacji między binami od N_norm → dodaj wprost)
            el = np.sqrt((((Hp[ok][:, lo].sum(1) / Pn[ok]) ** 2) * Nn[ok]).sum() + (Hp[ok][:, lo].sum(1) * (Nn[ok] / Pn[ok]) ** 2).sum())
            r["mirror_low"] = {"obs": float(ol), "pred": float(pl), "ratio": float(ol / pl) if pl > 0 else None, "ratio_err": float(np.hypot(np.sqrt(ol) / pl, ol * el / pl ** 2)) if pl > 0 else None}
            # (b) widma ξ
            sp = {}
            for nm, m in (("impossible_stream", R.S[stream] & d[f"tag_{s}"] & (d[f"xical_{s}"] > POOL_THR)), ("impossible_stream_hasD", R.S[stream] & R.has & d[f"tag_{s}"] & (d[f"xical_{s}"] > POOL_THR)),
                          ("impossible_nonAFP", d["non"] & d[f"tag_{s}"] & (d[f"xical_{s}"] > POOL_THR)), ("alltags_stream_hasD", R.S[stream] & R.has & d[f"tag_{s}"]),
                          ("SR_stream_hasD", R.S[stream] & R.has & d[f"tag_{s}"] & (R.delta[s] <= lib.DELTA)), ("CR_stream_hasD", R.S[stream] & R.has & d[f"tag_{s}"] & (R.delta[s] > lib.DELTA + lib.M_NOM))):
                x = d[f"xi1_{s}"][m]; sp[nm] = {"n": int(m.sum()), "h": np.histogram(x, EX)[0].tolist(), "mean": float(x.mean()) if len(x) else None, "mean_err": float(x.std() / np.sqrt(max(len(x), 1))) if len(x) else None}
            r["xi_spectra"] = {"edges": EX.tolist(), **sp}
            # (c) topologia: udział zdarzeń z protonem po stronie przeciwnej
            t = R.S[stream] & R.has & d[f"tag_{s}"]
            cr = t & (R.delta[s] > lib.DELTA + lib.M_NOM); sr = t & (R.delta[s] <= lib.DELTA)
            pp = R.S[pool] & R.has & (d[f"np_{s}"] == 0)
            r["opp_proton_frac"] = {"CR": float((d[f"np_{o}"][cr] > 0).mean()) if cr.any() else None, "SR": float((d[f"np_{o}"][sr] > 0).mean()) if sr.any() else None,
                                    "pool": float((d[f"np_{o}"][pp] > 0).mean()), "n_CR": int(cr.sum()), "n_CR_opp": int((cr & (d[f"np_{o}"] > 0)).sum())}
            # estymator z wykluczeniem zdarzeń z protonem po stronie przeciwnej (dane i pula)
            w = (R.has & (d[f"np_{o}"] == 0)).astype(float)
            e1 = R.est(s, stream, pool, w=w); e0 = R.est(s, stream, pool)
            r["est_no_opp_proton"] = {k: e1[k] for k in ("n_sr", "n_cr", "n_acc", "err", "f_acc")}; r["est_nominal"] = {k: e0[k] for k in ("n_sr", "n_cr", "n_acc", "err", "f_acc")}
            print(run, stream, s, "mirror low", r["mirror_low"], "opp", r["opp_proton_frac"], "f_noopp", e1["f_acc"], "f_nom", e0["f_acc"], flush=True)
        R.cache.clear()
json.dump(OUT, open("extra.json", "w"), indent=1)
