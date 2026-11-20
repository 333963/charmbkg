#!/usr/bin/env python3
"""Uzupełnienie v1: (a) skład klasy „other” w 428770 wg bitów trygera (bez zmiennych protonowych), (b) próbka J12only ∪ MBTSbased
(inkluzywna bez klasy „other”): ε_veto zdarzeń z D0, sygnału D* (SBS) i fit Poissona, (c) liczby do normalizacji. Wyjście: out/addendum.{md,json}."""
import sys, os, json
sys.dont_write_bytecode = True
import numpy as np
import v1_veto_lib as L   # te same funkcje co v1 (bez pętli głównej)
I, C = L.I, L.C
z = np.load(os.path.join(L.R, "bt_data_428770.npz"), allow_pickle=False)
A = L.tab(z, "all_", ("pt",)); Pp = L.tab(z, "", ("m", "pt", "score", "dm"))
_, fi = np.unique(A["uid"], return_index=True); E = {k: v[fi] for k, v in A.items()}
g = (E["grl"] == 1) & (E["actual_mu"] > 0); c = C.classify(E, 428770); cp = C.classify(Pp, 428770)
out = {}; md = []
o = c["other"] & g
afp = (E["tp_HLT_noalg_L1AFP_A_OR_C"] == 1) | (E["tp_HLT_mb_sptrk_pt4_L1AFP_A_OR_C"] == 1)
out["other_frac_AFPchains"] = float(afp[o].mean()); out["other_frac_of_D0events"] = float(o.sum() / g.sum())
md.append(f"Klasa „other” (428770): {o.sum()} zdarzeń z D0 = {100*o.sum()/g.sum():.1f} % próbki; {100*afp[o].mean():.1f} % z nich ma bit HLT_noalg_L1AFP_A_OR_C lub HLT_mb_sptrk_pt4_L1AFP_A_OR_C.")
u = (c["J12only"] | c["MBTSbased"]) & g
ev = L.ev_frac(E["pas"], u); f = L.mu_bins_fit(E["actual_mu"][u], E["pas"][u])
base = (np.abs(Pp["m"] - C.M0_D0) < L.MWIN) & (Pp["score"] > L.STHR) & (Pp["pt"] >= 3) & (Pp["pt"] < 30) & (Pp["grl"] == 1) & (Pp["actual_mu"] > 0) & (cp["J12only"] | cp["MBTSbased"])
pas = Pp["pas"].astype(float); w = np.zeros(len(pas)); w[base] = I.sbs_w(Pp["dm"][base], Pp["pt"][base], pas[base] > 0)
sg = L.boot_eps(Pp["uid"][base], w[base], pas[base])
cpu = float(np.exp(-f["eps_v"] * E["actual_mu"][u]).mean())
out["J12orMBTS"] = {"D0cand_events": ev, "signal": sg, "pois": f, "C_PU": cpu}
md.append(f"J12only ∪ MBTSbased (428770): zdarzenia z D0 {L.fmt(ev)} % (N = {ev[2]}), sygnał D* (SBS) {L.fmt(sg)} % (N_S = {sg[2]:.0f}); "
          f"fit: f = {100*f['f']:.2f} ± {100*f['f_err']:.2f} %, ε_v = {f['eps_v']:.3f} ± {f['eps_v_err']:.3f}, χ²/ndf = {f['chi2']:.1f}/{f['ndf']}; C_PU = {cpu:.4f}.")
# różnica sygnał − zdarzenia z D0, oba runy łącznie (z out/veto.json)
V = json.load(open(os.path.join(L.OUT, "veto.json")))
d = [(V[s][k]["signal"][0] - V[s][k]["D0cand_events"][0], V[s][k]["signal"][1]) for s, k in (("data_428770", "all"), ("data_427929", "all427929"))]
wm = sum(x / e ** 2 for x, e in d) / sum(1 / e ** 2 for x, e in d); we = 1 / np.sqrt(sum(1 / e ** 2 for x, e in d))
out["signal_minus_D0events_comb"] = (float(wm), float(we))
md.append(f"Sygnał D* − zdarzenia z D0 (dane, oba runy, średnia ważona): {100*wm:+.2f} ± {100*we:.2f} %; MC: fd_nd sygnał − fd_nd zd. z D0 = "
          f"{100*(V['fd_nd']['MC']['signal'][0]-V['fd_nd']['MC']['D0cand_events'][0]):+.2f} ± {100*V['fd_nd']['MC']['signal'][1]:.2f} %.")
# ε_v: zgodność runów
a, b = V["data_428770"]["J12only"]["pois_D0cand"], V["data_427929"]["all427929"]["pois_D0cand"]
out["eps_v_diff_runs"] = (a["eps_v"] - b["eps_v"], float(np.hypot(a["eps_v_err"], b["eps_v_err"])))
md.append(f"ε_v(428770 J12only) − ε_v(427929) = {a['eps_v']-b['eps_v']:.3f} ± {np.hypot(a['eps_v_err'], b['eps_v_err']):.3f} ({(a['eps_v']-b['eps_v'])/np.hypot(a['eps_v_err'], b['eps_v_err']):.1f}σ).")
mu8 = E["actual_mu"][u]
md.append(f"Wpływ rozrzutu ε_v (0,66–0,79) na C_PU w 428770: {100*(np.exp(-0.656*mu8).mean()-np.exp(-0.789*mu8).mean()):.2f} %.")
open(os.path.join(L.OUT, "addendum.md"), "w").write("\n\n".join(md) + "\n"); json.dump(out, open(os.path.join(L.OUT, "addendum.json"), "w"), indent=1, default=float)
print("\n".join(md))
