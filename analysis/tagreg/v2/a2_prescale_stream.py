"""TAGREG v2 §4: preskale drzewo vs COOL, efektywne preskale L1/HLT z bitów, test H-strumień. Tylko liczniki zdarzeń i bitów."""
import json
import numpy as np, pandas as pd
from v2common import *
pd.set_option("display.width", 250, "display.max_rows", 500, "display.max_columns", 60)
res = {}
def ratio(n, k):
    """p = n/k, błąd dwumianowy (k z n)."""
    if k == 0: return (np.nan, np.nan)
    f = k / n; return (1 / f, np.sqrt(f * (1 - f) / n) / f ** 2)
def fmt(pe): return "-" if np.isnan(pe[0]) else f"{pe[0]:.4g} ± {pe[1]:.2g}"

# ---------- 4a. drzewo vs COOL
for run in RUNS:
    p = pd.read_csv(f"{D}/ps_perLB_{run}.csv"); h = p[p.kind == "HLT"]
    c = pd.read_csv(f"{D}/coolps_{run}.csv")
    lbs = np.array(sorted(h.lumiBlock.unique()))
    def cool(name, level):
        d = c[(c.name == name) & (c.level == level)]; out = np.full(len(lbs), np.nan); en = np.zeros(len(lbs), bool); key = np.zeros(len(lbs), int)
        for r in d.itertuples():
            m = (lbs >= r.lb_first) & (lbs <= r.lb_last); out[m] = r.prescale; en[m] = r.enabled == 1; key[m] = r.key
        return out, en, key
    tree = {ch: g.set_index("lumiBlock").ps_mode.reindex(lbs).to_numpy() for ch, g in h.groupby("chain")}
    rows = []
    for name in sorted(c[c.level == "HLT"].name.unique()):
        l1 = "L1_" + name.split("_L1")[1]
        ph, eh, kh = cool(name, "HLT"); pl, el, kl = cool(l1, "L1")
        prod = np.where(eh & el, ph * pl, 0.0); t = tree[name]
        ok = np.isclose(prod, t, rtol=0.01, atol=1e-6)
        rows.append({"lancuch": name, "LB": len(lbs), "LB_COOL_wlaczone": int((eh & el).sum()), "LB_zgodne_<1%": int(ok.sum()),
                     "max|stos-1|": float(np.nanmax(np.abs(prod[t > 0] / t[t > 0] - 1))) if (t > 0).any() else 0.0,
                     "drzewo": sorted(set(np.round(t, 2))), "COOL_L1": sorted(set(np.round(pl, 2))), "COOL_HLT": sorted(set(np.round(ph, 2))),
                     "drzewo_L1": sorted(set(np.round(tree[l1], 2))), "L1_zgodne": bool(np.allclose(np.where(el, pl, 0), tree[l1], rtol=0.01))})
    t = pd.DataFrame(rows); t.to_csv(f"{OUT}/ps_drzewo_vs_cool_{run}.csv", index=False)
    print(f"\n==== {run}: drzewo vs COOL (L1×HLT), {len(lbs)} LB GRL; klucze HLT w GRL {sorted(set(cool('HLT_noalg_L1J12','HLT')[2]))}, L1 {sorted(set(cool('L1_J12','L1')[2]))}")
    print(t.to_string(index=False))
    # czynnik HLT z samego drzewa = ps(łańcuch)/ps(pozycja L1) dla wszystkich 75
    rows = []
    for name in sorted(tree):
        if not name.startswith("HLT_"): continue
        l1 = "L1_" + name.split("_L1")[1]
        if l1 not in tree: l1n = None
        t_, l_ = tree[name], tree.get(l1)
        f = sorted(set(np.round(t_[l_ > 0] / l_[l_ > 0], 3))) if l_ is not None and (l_ > 0).any() else None
        rows.append({"lancuch": name, "ps_calk": sorted(set(np.round(t_, 3))), "L1": l1 if l_ is not None else "(brak w drzewie)", "ps_L1": sorted(set(np.round(l_, 3))) if l_ is not None else None,
                     "czynnik_HLT": f, "n_pass_GRL": int(h[h.chain == name].n_pass.sum())})
    t = pd.DataFrame(rows); t.to_csv(f"{OUT}/ps_drzewo_rozklad_{run}.csv", index=False)
    print(t[t.n_pass_GRL > 0].sort_values("n_pass_GRL", ascending=False).to_string(index=False))

# ---------- 4b–d. efektywne preskale i strumień (428770; 427929 skrótowo)
for run in RUNS:
    d = load(run); H, hlt, AV, TBP, l1 = d["H"], d["hlt"], d["AV"], d["TBP"], d["l1"]
    ix = {c: j for j, c in enumerate(hlt)}; il = {c: j for j, c in enumerate(l1)}
    lb = d["ev"].lumiBlock.to_numpy(); nsel = d["ev"].nSelBL.to_numpy(); n = len(lb)
    b = lambda c: H[:, ix[c]]
    print(f"\n==== {run}: T0 = {n}; TAV ⊂ TBP: naruszeń {int((AV & ~TBP).sum())}")
    nb = H.sum(1)
    # zbiór S
    only = {c: int((b(c) & (nb == 1)).sum()) for c in hlt}
    tot = {c: int(b(c).sum()) for c in hlt}
    S0abs = [c for c in hlt if only[c] > 20]                      # kryterium z KRYTERIA.md (próg bezwzględny)
    # PO FAKCIE: próg względny — N_tylko > 1 % zdarzeń łańcucha (i ≥ 5); bezwzględny wpuszczał sptrk_pt4_L1AFP (28 z 177 876 = 0,016 %,
    # poniżej częstości zdarzeń bez żadnego bitu HLT, 0,06 %)
    rel = lambda v, c: v >= 5 and v > 0.01 * tot[c]
    S0 = [c for c in hlt if rel(only[c], c)]
    inS0 = np.array([c in S0 for c in hlt])
    anyS0 = H[:, inS0].any(1)
    ext = {c: int((b(c) & ~anyS0).sum()) for c in hlt if c not in S0}
    S1 = [c for c, v in ext.items() if rel(v, c)]
    S = S0 + S1; inS = np.array([c in S for c in hlt]); anyS = H[:, inS].any(1)
    print("S0 wg progu bezwzględnego > 20:", [(c, only[c], tot[c]) for c in S0abs])
    print("S0 (N_tylko > 1 % N_T0):", [(c, only[c], tot[c]) for c in S0])
    print("S1 (bit bez żadnego łańcucha S0, > 1 %):", [(c, ext[c], tot[c]) for c in S1])
    print("poza S, z bitem bez S0 (>0):", [(c, v, tot[c]) for c, v in ext.items() if v > 0 and c not in S1])
    print("poza S, nigdy bez S0:", [(c, tot[c]) for c, v in ext.items() if v == 0 and tot[c] > 0])
    print(f"zdarzenia T0 z ≥1 bitem S: {anyS.sum()} / {n} = {100*anyS.mean():.4f} %; bez S: {int((~anyS).sum())}, w tym bez żadnego bitu HLT {int((nb==0).sum())}")
    pd.DataFrame([{"lancuch": c, "N_T0": tot[c], "N_tylko": only[c], "N_bez_S0": ext.get(c, only[c]), "w_S": c in S} for c in hlt]).to_csv(f"{OUT}/strumien_S_{run}.csv", index=False)
    res[run] = {"S0": S0, "S1": S1, "frac_anyS": float(anyS.mean()), "n_noS": int((~anyS).sum()), "n_noHLT": int((nb == 0).sum())}
    if run != 428770:
        continue
    blocks = {"LB4-16": lb <= 16, "LB17-20": (lb >= 17) & (lb <= 20), "LB>=21": lb >= 21, "wszystkie": np.ones(n, bool)}
    # efektywne preskale L1
    J = b("HLT_noalg_L1J12"); NA = b("HLT_noalg_L1AFP_A_OR_C"); AAC = b("HLT_noalg_L1AFP_A_AND_C")
    p = pd.read_csv(f"{D}/ps_perLB_{run}.csv"); h = p[p.kind == "HLT"]
    tps = {ch: dict(zip(g.lumiBlock, g.ps_mode)) for ch, g in h.groupby("chain")}
    rows = []
    for item, ref, rn in [("L1_AFP_A_OR_C", J, "noalg_L1J12"), ("L1_AFP_A_AND_C", J, "noalg_L1J12"), ("L1_AFP_A", J, "noalg_L1J12"), ("L1_AFP_C", J, "noalg_L1J12"),
                          ("L1_AFP_FSA_BGRP12", J, "noalg_L1J12"), ("L1_AFP_NSC_BGRP12", J, "noalg_L1J12"),
                          ("L1_J12", NA | AAC, "noalg_L1AFP(OR|AND)"), ("L1_J15", NA | AAC, "noalg_L1AFP(OR|AND)"), ("L1_J20", NA | AAC, "noalg_L1AFP(OR|AND)"),
                          ("L1_MBTS_2", NA | AAC, "noalg_L1AFP(OR|AND)"), ("L1_MBTS_2", J, "noalg_L1J12"), ("L1_MBTS_1", NA | AAC | J, "noalg AFP|J12"),
                          ("L1_EM7", NA | AAC, "noalg_L1AFP(OR|AND)"), ("L1_EM3", NA | AAC, "noalg_L1AFP(OR|AND)"), ("L1_MU3V", NA | AAC | J, "noalg AFP|J12"),
                          ("L1_RD0_FILLED", NA | AAC | J, "noalg AFP|J12")]:
        j = il[item]
        for bn, bm in blocks.items():
            m = ref & bm & TBP[:, j]; k = int((m & AV[:, j]).sum()); nn = int(m.sum())
            ps = sorted(set(round(tps[item][x], 2) for x in np.unique(lb[bm]))) if bm.any() else []
            rows.append({"L1": item, "odniesienie": rn, "blok": bn, "N_TBP": nn, "N_TAV": k, "p_eff": fmt(ratio(nn, k)), "ps_drzewo": ps})
    t = pd.DataFrame(rows); t.to_csv(f"{OUT}/peff_L1_{run}.csv", index=False)
    print("\n-- efektywne preskale L1 (T0):"); print(t.to_string(index=False))
    # efektywne preskale HLT noalg_L1AFP_A_OR_C wobec różnych odniesień
    jA = il["L1_AFP_A_OR_C"]; avA, tbA = AV[:, jA], TBP[:, jA]
    inS_other = inS.copy(); inS_other[ix["HLT_noalg_L1AFP_A_OR_C"]] = False
    anyS_other = H[:, inS_other].any(1)
    refs = ["HLT_noalg_L1J12", "HLT_j20_L1AFP_A_OR_C", "HLT_j20f_L1AFP_A_OR_C", "HLT_2j10_20detaAA_L1AFP_A_OR_C", "HLT_noalg_L1AFP_A_AND_C", "HLT_j20_L1MBTS_2",
            "HLT_j10f_L1MBTS_2", "HLT_noalg_L1EM7", "HLT_noalg_L1EM3", "HLT_noalg_L1MU3V", "HLT_noalg_L1AFP_FSA_BGRP12",
            "HLT_mb_sptrk_pt4_L1AFP_A_OR_C", "HLT_mb_sptrk_pt6_L1AFP_A_OR_C", "HLT_mb_sptrk_pt8_L1AFP_A_OR_C", "HLT_mb_sptrk_pt2_L1AFP_A_OR_C",
            "HLT_mb_sptrk_pt2_L1AFP_A_AND_C", "HLT_mb_sptrk_pt4_L1MBTS_2", "HLT_mb_sptrk_pt6_L1MBTS_2", "HLT_noalg_mb_L1MBTS_2", "HLT_mb_sptrk_pt2_L1RD0_FILLED"]
    rows = []
    for r in refs:
        R = b(r)
        row = {"odniesienie": r, "w_S": r in S}
        for bn, bm in (("LB>=17", lb >= 17), ("LB4-16", lb <= 16)):
            m = R & bm & avA
            row[f"n_ref {bn}"] = int(m.sum()); row[f"p_HLT {bn}"] = fmt(ratio(int(m.sum()), int((m & NA).sum())))
            io = inS_other.copy(); io[ix[r]] = False            # S bez noalg_L1AFP i bez samego odniesienia
            m2 = m & H[:, io].any(1)
            row[f"n_ref∧S {bn}"] = int(m2.sum()); row[f"p_HLT|S {bn}"] = fmt(ratio(int(m2.sum()), int((m2 & NA).sum())))
            row[f"P(S|ref) {bn}"] = round(float(m2.sum() / max(m.sum(), 1)), 4)
        m = R & (lb >= 17) & tbA
        row["p_calk(TBP) LB>=17"] = fmt(ratio(int(m.sum()), int((m & NA).sum())))
        rows.append(row)
    t = pd.DataFrame(rows); t.to_csv(f"{OUT}/peff_noalgAFP_{run}.csv", index=False)
    print("\n-- p_eff HLT_noalg_L1AFP_A_OR_C = N(ref ∧ TAV_AFP)/N(ref ∧ TAV_AFP ∧ noalg); '|S' = dodatkowo ≥1 łańcuch z S inny niż noalg_L1AFP (T0):")
    print(t.to_string(index=False))
    # zależność od nSelBL i kandydata: sptrk4 odniesienie, z i bez warunku S
    rows = []
    R = b("HLT_mb_sptrk_pt4_L1AFP_A_OR_C") & (lb >= 17) & avA
    cuts = {"nSelBL<20": nsel < 20, "20-40": (nsel >= 20) & (nsel < 40), ">=40": nsel >= 40, "z kand. D": d["cand"]["D"], "bez kand.": ~d["cand"]["D"],
            "z L1_J12": TBP[:, il["L1_J12"]], "bez L1_J12": ~TBP[:, il["L1_J12"]]}
    for cn, cm in cuts.items():
        m = R & cm; m2 = m & anyS_other
        rows.append({"przekrój": cn, "n_ref": int(m.sum()), "p_eff": fmt(ratio(int(m.sum()), int((m & NA).sum()))), "P(S|ref)": round(float(m2.sum() / max(m.sum(), 1)), 3),
                     "n_ref∧S": int(m2.sum()), "p_eff|S": fmt(ratio(int(m2.sum()), int((m2 & NA).sum())))})
    t = pd.DataFrame(rows); t.to_csv(f"{OUT}/peff_noalgAFP_przekroje_{run}.csv", index=False)
    print("\n-- odniesienie sptrk_pt4_L1AFP, LB ≥ 17, przekroje:"); print(t.to_string(index=False))
    # kto towarzyszy sptrk4 (które łańcuchy S)
    R = b("HLT_mb_sptrk_pt4_L1AFP_A_OR_C")
    comp = sorted([(c, int((R & b(c)).sum())) for c in S], key=lambda x: -x[1])[:8]
    print("sptrk_pt4_L1AFP: N =", int(R.sum()), "towarzyszące łańcuchy S:", comp, "| bez żadnego S:", int((R & ~anyS).sum()))
    # inne łańcuchy noalg: p_HLT wobec odniesień w S
    rows = []
    for X, item in [("HLT_noalg_L1AFP_A_AND_C", "L1_AFP_A_AND_C"), ("HLT_noalg_L1J12", "L1_J12"), ("HLT_noalg_L1EM7", "L1_EM7"), ("HLT_noalg_L1EM3", "L1_EM3"),
                    ("HLT_noalg_L1MU3V", "L1_MU3V"), ("HLT_noalg_L1AFP_FSA_BGRP12", "L1_AFP_FSA_BGRP12"), ("HLT_noalg_L1AFP_NSA_BGRP12", "L1_AFP_NSA_BGRP12"),
                    ("HLT_noalg_mb_L1MBTS_2", "L1_MBTS_2"), ("HLT_noalg_mb_L1MBTS_1", "L1_MBTS_1"), ("HLT_noalg_L1AFP_A_OR_C_J12", "L1_AFP_A_OR_C_J12")]:
        refm = (NA | b("HLT_j20_L1AFP_A_OR_C")) if "AFP_A_OR_C" != item[3:] and not X.endswith("L1J12") else NA
        if X == "HLT_noalg_L1J12": refm = NA | AAC
        m = refm & ~b(X) if False else refm
        mt = m & TBP[:, il[item]]; ma = m & AV[:, il[item]]
        rows.append({"X": X, "N(ref∧TBP)": int(mt.sum()), "N(ref∧TAV)": int(ma.sum()), "N(X)": int((ma & b(X)).sum()),
                     "p_calk_eff": fmt(ratio(int(mt.sum()), int((mt & b(X)).sum()))), "p_HLT_eff": fmt(ratio(int(ma.sum()), int((ma & b(X)).sum()))),
                     "ps_drzewo": sorted(set(round(v, 2) for v in tps[X].values())), "w_S": X in S})
    t = pd.DataFrame(rows); t.to_csv(f"{OUT}/peff_inne_{run}.csv", index=False)
    print("\n-- inne łańcuchy noalg (odniesienie: noalg_L1AFP_A_OR_C ∨ j20_L1AFP_A_OR_C; dla J12: noalg AFP):"); print(t.to_string(index=False))
    # czy j20_L1J12 ≡ noalg_L1J12, j20 na AFP: ułamek TAV_AFP z j20 w próbce nieobciążonej (noalg)
    m = NA
    print("\nw próbce noalg_L1AFP (nieobciążona wzgl. HLT): udział z bitem",
          {c: round(float((m & b(c)).sum() / m.sum()), 4) for c in ("HLT_j20_L1AFP_A_OR_C", "HLT_j20f_L1AFP_A_OR_C", "HLT_2j10_20detaAA_L1AFP_A_OR_C",
                                                                   "HLT_mb_sptrk_pt2_L1AFP_A_OR_C", "HLT_mb_sptrk_pt4_L1AFP_A_OR_C", "HLT_mb_sptrk_pt6_L1AFP_A_OR_C", "HLT_noalg_L1AFP_A_AND_C", "HLT_noalg_L1J12")})
    m = NA & b("HLT_mb_sptrk_pt4_L1AFP_A_OR_C")
    print("w próbce noalg ∧ sptrk4: udział z j20_L1AFP", round(float((m & b('HLT_j20_L1AFP_A_OR_C')).sum() / m.sum()), 4), " z dowolnym S≠noalg", round(float((m & anyS_other).sum() / m.sum()), 4), "N", int(m.sum()))
    res[run]["S"] = S
json.dump(res, open(f"{OUT}/a2.json", "w"), indent=1)
