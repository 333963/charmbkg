"""TAGREG v2 §1–3: kontrola, klasa none11, rodziny. Tylko liczniki zdarzeń i bitów."""
import json, sys
import numpy as np, pandas as pd
from v2common import *
pd.set_option("display.width", 250, "display.max_rows", 500, "display.max_columns", 60)
S = json.load(open(f"{D}/summary.json"))["samples"]
EXP = {428770: dict(T0=1338608, T1=1233602, T2=1134024, D0=76030, A=570220, C=321855, AC=241949,
                    cls=[544031, 94359, 8769, 10623, 3095, 8725, 4760, 459662], clsD0=[5563, 39454, 2492, 6257, 1098, 3784, 1051, 16331], clsD=[11073, 53435, 3355, 7851, 1474, 5131, 1609, 30192]),
       427929: dict(T0=707688, T1=703350, T2=701183, D0=10438, A=354983, C=338503, AC=7697,
                    cls=[92, 0, 0, 0, 475972, 0, 225108, 11], clsD0=[1, 0, 0, 0, 7138, 0, 3299, 0], clsD=None)}
CL = ["AFPnoalg", "AFPsptrk4", "AFPboth", "AFP&J12", "AFP&MBTS", "J12only", "MBTS!A", "none"]
res = {}
def pct(x, n): return f"{100 * x / n:.2f}" if n else "-"
for run in RUNS:
    d = load(run); H, hlt, reg, cand = d["H"], d["hlt"], d["reg"], d["cand"]
    ix = {c: j for j, c in enumerate(hlt)}
    T2 = reg["T2"]; s = S[f"data_{run}"]; e = EXP[run]
    print(f"\n======== {run} ========")
    got = dict(T0=len(T2), T1=int(reg["T1"].sum()), T2=int(T2.sum()), D0=int((T2 & cand["D0"]).sum()), A=int((T2 & d["a"] & ~d["c"]).sum()),
               C=int((T2 & d["c"] & ~d["a"]).sum()), AC=int((T2 & d["a"] & d["c"]).sum()))
    print("kontrola liczb:", {k: (v, e[k], "OK" if v == e[k] else "ROZNE") for k, v in got.items()})
    print("pliki", s["n_files"], "dup", s["n_dup_event_number"], "identyczne listy", s["chain_lists_identical"], "ps_missing", s["ps_missing"], "nHLT", len(hlt), "nL1", len(d["l1"]))
    print("T2 z D", int((T2 & cand["D"]).sum()), "D+", int((T2 & cand["Dp"]).sum()), "D*", int((T2 & cand["Dst"]).sum()))
    # klasy 11-bitowe
    b = lambda c: H[:, ix[c]]
    na, s4, J = b("HLT_noalg_L1AFP_A_OR_C"), b("HLT_mb_sptrk_pt4_L1AFP_A_OR_C"), b("HLT_noalg_L1J12")
    M = np.zeros(len(T2), bool)
    for c in TRIG11:
        if "MBTS" in c: M |= b(c)
    A = na | s4
    cls = np.select([A & M, A & J, na & s4, na, s4, M, J], [4, 3, 2, 0, 1, 6, 5], 7)   # jak ../common.py: MBTS przed J12
    t = pd.DataFrame({"klasa": CL, "zdarz": [int((T2 & (cls == k)).sum()) for k in range(8)], "exp": e["cls"],
                      "zD0": [int((T2 & cand["D0"] & (cls == k)).sum()) for k in range(8)], "expD0": e["clsD0"],
                      "zD": [int((T2 & cand["D"] & (cls == k)).sum()) for k in range(8)]})
    if e["clsD"]: t["expD"] = e["clsD"]
    print(t.to_string(index=False)); print("klasy zgodne:", (t.zdarz == t.exp).all(), "D0:", (t.zD0 == t.expD0).all())
    none11 = cls == 7
    assert (none11 == ~H[:, [ix[c] for c in TRIG11]].any(1)).all()
    noneH = ~H.any(1)
    res[run] = {"control": got, "cls11": t.to_dict("list")}
    # ---- §2 klasa none11
    sets = {"wszystkie": T2, "D": T2 & cand["D"], "D0": T2 & cand["D0"], "Dp": T2 & cand["Dp"], "Dst": T2 & cand["Dst"]}
    nb = H.sum(1)
    rows = []
    for j, c in enumerate(hlt):
        if c in TRIG11: continue
        r = {"lancuch": c}
        for k, m in sets.items():
            mm = m & none11
            r[f"lacz_{k}"] = int(H[mm, j].sum()); r[f"tylko_{k}"] = int((H[mm, j] & (nb[mm] == 1)).sum())
        rows.append(r)
    t = pd.DataFrame(rows).sort_values("lacz_wszystkie", ascending=False)
    t = t[t.lacz_wszystkie > 0]
    t.to_csv(f"{OUT}/none11_ranking_{run}.csv", index=False)
    print("\n-- none11: N =", {k: int((m & none11).sum()) for k, m in sets.items()}, " noneHLT:", {k: int((m & noneH).sum()) for k, m in sets.items()})
    print(t.head(22).to_string(index=False))
    # zachłanny
    for k in ("wszystkie", "D0", "Dp", "Dst"):
        left = sets[k] & none11; n0 = left.sum(); out = []
        while left.sum() and len(out) < 8:
            cnt = H[left].sum(0); j = int(cnt.argmax())
            if cnt[j] == 0: break
            out.append((hlt[j], int(cnt[j]), round(100 * cnt[j] / max(n0, 1), 2))); left &= ~H[:, j]
        print(f"zachłannie none11 [{k}] N={n0}:", out, "reszta", int(left.sum()))
        res[run][f"greedy_none11_{k}"] = out
    # liczba bitów HLT w none11
    print("liczba bitów HLT w none11∧T2:", np.bincount(np.minimum(nb[T2 & none11], 8)).tolist())
    # reszta noneHLT
    m = T2 & noneH
    print("noneHLT∧T2:", int(m.sum()), "L1 po preskali:", sorted([(d["l1"][j], int(v)) for j, v in enumerate(d["AV"][m].sum(0)) if v], key=lambda x: -x[1])[:8],
          "| TBP:", sorted([(d["l1"][j], int(v)) for j, v in enumerate(d["TBP"][m].sum(0)) if v], key=lambda x: -x[1])[:8],
          "| LB min/max", (int(d["ev"].lumiBlock[m].min()), int(d["ev"].lumiBlock[m].max())) if m.sum() else None)
    if m.sum(): print("  noneHLT per LB (top):", d["ev"].lumiBlock[m].value_counts().head(6).to_dict(), " AV żaden:", int((~d["AV"][m].any(1)).sum()))
    # ---- §3 rodziny
    fc = famclass(H, d["fam"])
    for R in ("T2", "T1"):
        m0 = reg[R]
        rows = []
        for f in range(8):
            r = {"rodzina": FAM[f]}
            for k in ("wszystkie", "D", "D0", "Dp", "Dst"):
                mm = m0 & (cand[k] if k != "wszystkie" else True)
                n = int((mm & (fc == f)).sum()); r[k] = n; r[k + "_%"] = pct(n, mm.sum())
                inc = int((mm & H[:, d["fam"] == f].any(1)).sum()) if f < 7 else n
                r[k + "_lacz"] = inc
            rows.append(r)
        t = pd.DataFrame(rows)
        t.to_csv(f"{OUT}/rodziny_{R}_{run}.csv", index=False)
        print(f"\n-- rodziny rozłączne {R} (i łączne _lacz):"); print(t.to_string(index=False))
        if R == "T2": res[run]["fam_T2"] = t.to_dict("list")
    # nakładanie rodzin (T2, z D0)
    for k in ("wszystkie", "D0"):
        mm = T2 & (cand[k] if k != "wszystkie" else True)
        F = np.stack([H[mm][:, d["fam"] == f].any(1) for f in range(7)], 1).astype(int)
        print(f"macierz nakładania rodzin T2 [{k}]:"); print(pd.DataFrame(F.T @ F, index=FAM[:7], columns=FAM[:7]).to_string())
    # skład rodzin: łańcuchy w rodzinie (łącznie) T2 / D0
    rows = [{"lancuch": c, "rodzina": FAM[d["fam"][j]], "T2": int(H[T2, j].sum()), "T2_D": int(H[T2 & cand["D"], j].sum()), "T2_D0": int(H[T2 & cand["D0"], j].sum()),
             "T2_Dp": int(H[T2 & cand["Dp"], j].sum()), "T2_Dst": int(H[T2 & cand["Dst"], j].sum()),
             "T0_tylko": int((H[:, j] & (nb == 1)).sum()), "T0": int(H[:, j].sum())} for j, c in enumerate(hlt)]
    t = pd.DataFrame(rows).sort_values(["rodzina", "T2"], ascending=[True, False]); t.to_csv(f"{OUT}/lancuchy_{run}.csv", index=False)
    print("\n-- łańcuchy (T2>0):"); print(t[t.T0 > 0].to_string(index=False))
    # LB
    lb = d["ev"].lumiBlock.to_numpy(); e8 = lbblocks(lb, reg["T1"])
    rows = []
    for i in range(len(e8) - 1):
        m = (lb >= e8[i]) & (lb < e8[i + 1])
        r = {"LB": f"{e8[i]}-{e8[i+1]-1}", "T1": int((m & reg['T1']).sum()), "T2": int((m & T2).sum())}
        for k in ("wszystkie", "D0"):
            mm = m & T2 & (cand[k] if k != "wszystkie" else True)
            for f in range(8): r[f"{k}:{FAM[f]}"] = pct((mm & (fc == f)).sum(), mm.sum())
        rows.append(r)
    t = pd.DataFrame(rows); t.to_csv(f"{OUT}/rodziny_LB_{run}.csv", index=False)
    print("\n-- rodziny per blok LB (T2, %):"); print(t.to_string(index=False))
json.dump(res, open(f"{OUT}/a1.json", "w"), indent=1, default=int)
