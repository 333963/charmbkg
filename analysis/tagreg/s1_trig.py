#!/usr/bin/env python3
"""s1 (część A): skład trygerów w regionie tagowanym (T1, T2), per klasa/bit/LB; preskale zapisane i efektywne; ułamki L_eff/L."""
import numpy as np, os, common as C
res = {}
def pct(a, b): return 100.0 * a / b if b else float("nan")
for r in C.RUNS:
    d = C.load(r); cl, b = C.classes(d); T0, T1, T2 = C.regions(d); R = {}
    grl = d["grl"]
    # wariant ze score: D0 z bt all_ (pT>=3, okno 1.75-2.0, Lxy>=0.25), score>0.6
    z = np.load(os.path.join(C.BT, f"bt_data_{r}.npz"), allow_pickle=False)
    ae = z["all_event"]; asc = z["all_score"]
    evs = np.sort(d["evno"]); o = np.argsort(d["evno"])
    def mark(e):
        p = np.searchsorted(evs, e); p = np.clip(p, 0, len(evs) - 1); ok = evs[p] == e
        m = np.zeros(len(evs), bool); m[o[p[ok]]] = True; return m, ok.mean()
    inall, cov = mark(np.unique(ae)); hs, _ = mark(np.unique(ae[asc > 0.6]))
    R["bt_all_cov_match"] = float(cov); R["bt_all_frac_of_hasD0_grl"] = float((inall & grl & d["hasD0"]).sum() / max((grl & d["hasD0"]).sum(), 1))
    d["hasD0s"] = hs
    sel = {"zdarzenia": np.ones(len(grl), bool), "z kand. D (slim)": d["hasD"], "z kand. D0 pT>=3": d["hasD0"], "z D0 score>0.6 (bt all_)": hs}
    R["tot"] = {}
    for reg, m in (("GRL", grl), ("T0", T0), ("T1", T1), ("T2", T2)):
        R["tot"][reg] = {k: int((m & v).sum()) for k, v in sel.items()}
        R["tot"][reg]["sb_cands"] = {c: int(d["sbc"][m, i].sum()) for i, c in enumerate(C.SBC)}
    R["disj"] = {}; R["bits"] = {}
    for reg, m in (("GRL", grl), ("T1", T1), ("T2", T2)):
        R["disj"][reg] = {k: {s: int((m & v & sv).sum()) for s, sv in sel.items()} | {"sb": {c: int(d["sbc"][m & v, i].sum()) for i, c in enumerate(C.SBC)}}
                          for k, v in cl.items()}
        R["bits"][reg] = {c: {s: int((m & d["T_" + c] & sv).sum()) for s, sv in sel.items()} for c in C.CH}
        R["bits"][reg]["L1TBP"] = {c: {s: int((m & d["TBP_" + c] & sv).sum()) for s, sv in sel.items()} for c in C.L1}
    # strony w T2
    R["sides_T2"] = {k: {"A only": int((T2 & d["tag_A"] & ~d["tag_C"] & v).sum()), "C only": int((T2 & d["tag_C"] & ~d["tag_A"] & v).sum()),
                         "A&C": int((T2 & d["tag_A"] & d["tag_C"] & v).sum())} for k, v in sel.items()}
    R["sides_T2_none"] = {"A only": int((T2 & cl["none"] & d["tag_A"] & ~d["tag_C"]).sum()), "C only": int((T2 & cl["none"] & d["tag_C"] & ~d["tag_A"]).sum()), "A&C": int((T2 & cl["none"] & d["tag_A"] & d["tag_C"]).sum())}
    # LB: 8 bloków równolicznych w T1
    e = np.unique(np.quantile(d["lb"][T1], np.linspace(0, 1, 9)[1:-1])); g = np.searchsorted(e, d["lb"], side="right")
    R["lb_edges"] = e.tolist(); R["lb"] = []
    for k in range(len(e) + 1):
        mk = g == k
        row = {"lb_lo": int(d["lb"][mk & grl].min()), "lb_hi": int(d["lb"][mk & grl].max()), "mu": float(d["mu"][mk & grl].mean()), "nGRL": int((mk & grl).sum()),
               "T1": int((mk & T1).sum()), "T2": int((mk & T2).sum()), "T2_D0": int((mk & T2 & d["hasD0"]).sum()),
               "cls_T2": {c: int((mk & T2 & v).sum()) for c, v in cl.items()}, "cls_T2_D0": {c: int((mk & T2 & v & d["hasD0"]).sum()) for c, v in cl.items()}}
        R["lb"].append(row)
    # preskale zapisane per LB
    cols, ps = C.trigps(r); R["ps_cols"] = cols
    lbs = ps[:, 0].astype(int); inGRL = np.isin(lbs, np.unique(d["lb"][grl]))
    # wagi LB: (a) liczba zdarzeń łańcucha bez preskali, (b) <mu>
    ref = "T_" + (C.J12 if r == 428770 else "HLT_mb_mbts_L1MBTS_1")
    wa = np.array([((d["lb"] == l) & grl & d[ref]).sum() for l in lbs], float)
    wb = np.array([d["mu"][(d["lb"] == l) & grl].mean() if ((d["lb"] == l) & grl).any() else 0.0 for l in lbs])
    R["lumi_frac"] = {}; R["ps_values"] = {}
    for j, c in enumerate(cols[1:-1], 1):
        p = ps[:, j]; ok = inGRL & (p > 0)
        inv = np.where(p > 0, 1.0 / np.where(p > 0, p, 1), 0.0)
        Fa = (wa * inv)[inGRL].sum() / wa[inGRL].sum(); Fb = (wb * inv)[inGRL].sum() / wb[inGRL].sum()
        u, cn = np.unique(np.round(p[inGRL], 2), return_counts=True)
        R["lumi_frac"][c] = {"F_a": float(Fa), "F_b": float(Fb)}; R["ps_values"][c] = {str(x): int(y) for x, y in zip(u, cn)}
    R["nLB_GRL"] = int(inGRL.sum()); R["lb_grl_range"] = [int(lbs[inGRL].min()), int(lbs[inGRL].max())]
    # preskale efektywne: P(X | zapisane przez łańcuch niezależny Y ∧ TBP L1(X))
    psmap = {c: dict(zip(lbs.tolist(), ps[:, j].tolist())) for j, c in enumerate(cols[1:-1], 1)}
    def pstored(c, m):
        v = np.array([psmap[c].get(int(l), 0.0) for l in np.unique(d["lb"][m])]); cnt = np.array([(m & (d["lb"] == l)).sum() for l in np.unique(d["lb"][m])])
        return float((cnt * v).sum() / cnt.sum()) if cnt.sum() else float("nan")   # średnia arytm. po zdarzeniach odniesienia
    tests = [("HLT_noalg_L1AFP_A_OR_C", "L1_AFP_A_OR_C", C.J12), ("HLT_noalg_L1AFP_A_OR_C", "L1_AFP_A_OR_C", "MBTSany"), ("HLT_noalg_L1AFP_A_OR_C", "L1_AFP_A_OR_C", C.AS),
             ("HLT_noalg_L1J12", "L1_J12", "AFPany"), ("HLT_noalg_L1J12", "L1_J12", "MBTSany"),
             ("HLT_noalg_mb_L1MBTS_2", "L1_MBTS_2", C.J12), ("HLT_noalg_mb_L1MBTS_2", "L1_MBTS_2", "AFPany"), ("HLT_noalg_mb_L1MBTS_1", "L1_MBTS_1", C.J12),
             ("HLT_noalg_mb_L1MBTS_1", "L1_MBTS_1", "AFPany"), ("HLT_mb_mbts_L1MBTS_1", "L1_MBTS_1", C.J12)]
    R["ps_eff"] = []
    for X, l1, Y in tests:
        ym = b["M"] if Y == "MBTSany" else b["A"] if Y == "AFPany" else d["T_" + Y]
        for lab, extra in (("wszystkie", np.ones(len(grl), bool)), ("T1", T1)):
            m = grl & ym & d["TBP_" + l1] & extra; nref = int(m.sum()); k = int((m & d["T_" + X]).sum())
            R["ps_eff"].append({"X": X, "Y": Y, "sel": lab, "n_ref": nref, "n_pass": k, "p_eff": nref / k if k else None,
                                "p_eff_err": (nref / k) / np.sqrt(k) if k else None, "p_stored_mean": pstored(X, m) if nref else None})
    # P(TBP L1 | łańcuch przeszedł) — kontrola sensu bitu TBP
    R["tbp_given_pass"] = {X: float((grl & d["T_" + X] & d["TBP_" + l1]).sum() / max((grl & d["T_" + X]).sum(), 1))
                           for X, l1 in ((C.AN, "L1_AFP_A_OR_C"), (C.AS, "L1_AFP_A_OR_C"), (C.J12, "L1_J12"), ("HLT_noalg_mb_L1MBTS_1", "L1_MBTS_1"), ("HLT_mb_mbts_L1MBTS_1", "L1_MBTS_1"))}
    # N·p w T2 per łańcuch (równoważnik bez preskali, preskale zapisane)
    R["Nxp_T2"] = {}
    for c in C.CH:
        m = T2 & d["T_" + c]
        pv = np.array([psmap[c].get(int(l), 0.0) for l in d["lb"][m]])
        R["Nxp_T2"][c] = {"N": int(m.sum()), "Nxp": float(pv.sum()), "N_D0": int((m & d["hasD0"]).sum()), "Nxp_D0": float(pv[d["hasD0"][m]].sum())}
    res[str(r)] = R
C.jdump(res, os.path.join(C.OUT, "s1.json"))
# krótki wydruk
for r, R in res.items():
    print("==", r, "bt all_ pokrycie hasD0:", round(R["bt_all_frac_of_hasD0_grl"], 3))
    for reg in ("GRL", "T0", "T1", "T2"):
        print(" ", reg, {k: v for k, v in R["tot"][reg].items() if k != "sb_cands"})
    for reg in ("T1", "T2"):
        tot = R["tot"][reg]
        print("  klasy", reg)
        for k, v in R["disj"][reg].items():
            print(f"    {k:10s} ev {v['zdarzenia']:8d} {pct(v['zdarzenia'], tot['zdarzenia']):6.2f}%  D {v['z kand. D (slim)']:7d} {pct(v['z kand. D (slim)'], tot['z kand. D (slim)']):6.2f}%  D0 {v['z kand. D0 pT>=3']:7d} {pct(v['z kand. D0 pT>=3'], tot['z kand. D0 pT>=3']):6.2f}%  D0s {v['z D0 score>0.6 (bt all_)']:6d}")
    print("  bity T2:", {c: (v["zdarzenia"], v["z kand. D0 pT>=3"]) for c, v in R["bits"]["T2"].items() if c != "L1TBP"})
    print("  L1TBP T2:", {c: (v["zdarzenia"], v["z kand. D0 pT>=3"]) for c, v in R["bits"]["T2"]["L1TBP"].items()})
    print("  strony T2:", R["sides_T2"], R["sides_T2_none"])
    for row in R["lb"]:
        t = row["T2"]; print(f"    LB {row['lb_lo']}-{row['lb_hi']} mu {row['mu']:.4f} GRL {row['nGRL']} T1 {row['T1']} T2 {t} T2/T1 {t/row['T1']:.3f} D0 {row['T2_D0']}",
                             {c: round(pct(v, t), 1) for c, v in row["cls_T2"].items()}, "| D0:", {c: round(pct(v, row['T2_D0']), 1) for c, v in row["cls_T2_D0"].items()})
    print("  preskale:", R["ps_values"]); print("  F:", R["lumi_frac"], R["nLB_GRL"], R["lb_grl_range"])
    for t in R["ps_eff"]:
        print("   ", t)
    print("  TBP|pass", R["tbp_given_pass"]); print("  Nxp", R["Nxp_T2"])
