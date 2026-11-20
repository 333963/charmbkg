#!/usr/bin/env python3
"""Krok 0: liczności klas (zdarzenia z kandydatem ev_*, pary D0←D*, inkluzywne all_*), macierz nakładania bitów, preskale.
Żadnych wydajności. Wynik: out/s0.json."""
import numpy as np
import afp_common as A

res = {}
for run in (428770, 427929):
    z = np.load(f"{A.R}/bt_data_{run}.npz", allow_pickle=True)
    r = {}
    for pref, lab in (("ev_", "zdarzenia"), ("e_", "pary"), ("all_e_", "all")):
        tp = A.bits(z, pref)
        n = len(tp[A.J12])
        grl = z[{"ev_": "ev_grl", "e_": "e_grl", "all_e_": "all_e_grl"}[pref]] == 1
        d = {"n": n, "grl": int(grl.sum()), "bit": {c: int(tp[c].sum()) for c in A.ALLCH}}
        anyb = np.any([tp[c] for c in A.ALLCH], axis=0)
        d["zaden_bit"] = int((~anyb).sum())
        d["alone"] = {c: int((tp[c] & ~np.any([tp[k] for k in A.ALLCH if k != c], axis=0)).sum()) for c in A.ALLCH}
        if run == 428770:
            cl = A.classes8(tp)
            d["cls"] = {k: int(v.sum()) for k, v in cl.items() if not k.startswith("_")}
            d["suma_rozlacznych"] = int(sum(cl[k].sum() for k in A.DISJ8))
            d["rd0_w_none"] = int((cl["none"] & cl["_rd0"]).sum())
            # preskale w klasach (mediana, min, max wartości > 0) i ⟨1/PS⟩
            ps = {}
            for c in (A.AN, A.AS, A.J12, "HLT_j20_L1MBTS_2", "HLT_mb_sptrk_pt4_L1MBTS_2", A.RD0):
                p = z[pref + "ps_" + c]
                ok = p > 0
                u, cnt = np.unique(np.round(p[ok], 2), return_counts=True)
                ps[c] = {"wartosci": {str(a): int(b) for a, b in zip(u, cnt)}, "frac_ps_le0": float((~ok).mean())}
                pb = p[tp[c] & ok]
                ps[c]["w_zd_z_bitem"] = {str(a): int(b) for a, b in zip(*np.unique(np.round(pb, 2), return_counts=True))}
            d["ps"] = ps
            lb = z[{"ev_": "ev_lb", "e_": "e_lb", "all_e_": "all_e_lb"}[pref]]
            d["lb_frac_lt17"] = {k: float((lb[cl[k]] < 17).mean()) if cl[k].any() else None for k in A.DISJ8}
            d["P_an_given_as"] = float(cl["_an"][cl["_as"]].mean()); d["P_as_given_an"] = float(cl["_as"][cl["_an"]].mean())
            d["P_A_given_J"] = float(cl["_A"][cl["_J"]].mean()); d["P_J_given_A"] = float(cl["_J"][cl["_A"]].mean())
            d["P_M_given_A"] = float(cl["_M"][cl["_A"]].mean())
        else:
            cl = A.classes9(tp, grl)
            d["cls"] = {k: int(v.sum()) for k, v in cl.items() if not k.startswith("_")}
        r[lab] = d
    res[str(run)] = r
A.jdump(res, f"{A.OUT}/s0.json")
for run, r in res.items():
    for lab, d in r.items():
        print(run, lab, "n", d["n"], "zaden_bit", d["zaden_bit"])
        print("  bit", d["bit"]); print("  alone", d["alone"]); print("  cls", d["cls"])
        if "ps" in d:
            print("  suma", d["suma_rozlacznych"], "rd0_w_none", d["rd0_w_none"])
            for c, p in d["ps"].items():
                print("  ps", c, p["wartosci"], "| z bitem", p["w_zd_z_bitem"], "ps<=0", round(p["frac_ps_le0"], 4))
            print("  P(an|as) %.4f P(as|an) %.4f P(A|J) %.4f P(J|A) %.4f P(M|A) %.4f" % (d["P_an_given_as"], d["P_as_given_an"], d["P_A_given_J"], d["P_J_given_A"], d["P_M_given_A"]))
            print("  lb<17", d["lb_frac_lt17"])
