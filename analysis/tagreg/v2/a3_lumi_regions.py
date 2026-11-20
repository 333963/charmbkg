"""TAGREG v2 §5: względna świetlność efektywna per łańcuch, pokrycie jednym łańcuchem o znanej preskali, regiony rozłączne. Tylko liczniki."""
import json
import numpy as np, pandas as pd
from v2common import *
pd.set_option("display.width", 250, "display.max_rows", 500, "display.max_columns", 60)
A2 = json.load(open(f"{OUT}/a2.json"))
COOL11 = set(TRIG11)
res = {}
for run in RUNS:
    d = load(run); H, hlt, reg, cand = d["H"], d["hlt"], d["reg"], d["cand"]
    ix = {c: j for j, c in enumerate(hlt)}; b = lambda c: H[:, ix[c]]
    T2 = reg["T2"]
    p = pd.read_csv(f"{D}/ps_perLB_{run}.csv"); h = p[p.kind == "HLT"]
    lbs = pd.read_csv(f"{D}/lb_summary_{run}.csv"); lbs = lbs[lbs.n_grl > 0].set_index("lumiBlock")
    refc = "HLT_noalg_L1J12" if run == 428770 else "HLT_mb_mbts_L1MBTS_1"
    wa = h[h.chain == refc].set_index("lumiBlock").n_pass.reindex(lbs.index).fillna(0).to_numpy().astype(float); wa /= wa.sum()
    wb = lbs.mean_mu_grl.to_numpy(); wb = wb / wb.sum()
    print(f"\n======== {run}: LB GRL {len(lbs)}; wagi (a) n_pass {refc}, (b) <mu>; max|a/b-1| = {np.abs(wa/wb-1).max():.3f}")
    S = set(A2[str(run)].get("S", A2[str(run)]["S0"] + A2[str(run)]["S1"]))
    sets = {"wszystkie": T2, "D": T2 & cand["D"], "D0": T2 & cand["D0"], "Dp": T2 & cand["Dp"], "Dst": T2 & cand["Dst"]}
    rows = []; F = {}
    for c in hlt:
        g = h[h.chain == c].set_index("lumiBlock").reindex(lbs.index)
        ps = g.ps_mode.to_numpy(); inv = np.where(ps > 0, 1 / np.where(ps > 0, ps, 1), 0.0)
        Fa, Fb = float((wa * inv).sum()), float((wb * inv).sum()); F[c] = Fa
        vals = sorted(set(np.round(ps, 2)))
        known = bool((ps > 0).all() and (g.ps_nvals == 1).all() and c in S)
        r = {"lancuch": c, "rodzina": FAM[d["fam"][ix[c]]], "ps": "/".join(f"{v:g}" for v in vals), "F_a": Fa, "F_b": Fb, "w_S": c in S, "COOL": "zgodne" if c in COOL11 else "tylko drzewo",
             "znana_ps": known}
        for k, m in sets.items():
            r[f"N_{k}"] = int((m & b(c)).sum()); r[f"%_{k}"] = round(100 * r[f"N_{k}"] / m.sum(), 2)
        r["N*p_wszystkie"] = float((ps[np.searchsorted(lbs.index.to_numpy(), d["ev"].lumiBlock.to_numpy()[T2 & b(c)])]).sum())
        r["N*p_D0"] = float((ps[np.searchsorted(lbs.index.to_numpy(), d["ev"].lumiBlock.to_numpy()[T2 & cand["D0"] & b(c)])]).sum())
        rows.append(r)
    t = pd.DataFrame(rows).sort_values("N_D0", ascending=False); t.to_csv(f"{OUT}/lancuchy_F_pokrycie_{run}.csv", index=False)
    print(t[t.N_wszystkie > 0].head(40).to_string(index=False))
    K = [c for c in t.lancuch if t.set_index("lancuch").znana_ps[c]]
    Kpos = [ix[c] for c in K]
    anyK = H[:, Kpos].any(1)
    print("\nłańcuchy o znanej preskali (K):", len(K))
    print("zdarzenia T2 z ≥1 bitem K [%]:", {k: round(100 * (m & anyK).sum() / m.sum(), 3) for k, m in sets.items()})
    # pokrycie zachłanne (wg D0) wśród K, R1 ustalone
    order = ["HLT_noalg_L1AFP_A_OR_C"] if run == 428770 else []
    left = sets["D0"].copy()
    for c in order: left &= ~b(c)
    while len(order) < 5:
        cnt = {c: int((left & b(c)).sum()) for c in K if c not in order}
        c = max(cnt, key=cnt.get)
        if cnt[c] == 0: break
        order.append(c); left &= ~b(c)
    print("kolejność regionów (R1 z góry, dalej zachłannie wg D0):", order)
    rows = []; taken = np.zeros(len(T2), bool)
    for i, c in enumerate(order):
        m = b(c) & ~taken
        r = {"region": f"R{i+1}", "lancuch": c + ("" if i == 0 else " ∧ ¬(poprzednie)"), "ps": t.set_index("lancuch").ps[c], "F": round(F[c], 5)}
        for k, s in sets.items():
            r[f"N_{k}"] = int((s & m).sum()); r[f"%_{k}"] = round(100 * (s & m).sum() / s.sum(), 2); r[f"%lacz_{k}"] = round(100 * (s & b(c)).sum() / s.sum(), 2)
        rows.append(r); taken |= b(c)
    r = {"region": "poza", "lancuch": "żaden z powyższych", "ps": "-", "F": np.nan}
    for k, s in sets.items():
        r[f"N_{k}"] = int((s & ~taken).sum()); r[f"%_{k}"] = round(100 * (s & ~taken).sum() / s.sum(), 2); r[f"%lacz_{k}"] = np.nan
    rows.append(r)
    tr = pd.DataFrame(rows); tr.to_csv(f"{OUT}/regiony_{run}.csv", index=False)
    print(tr.to_string(index=False))
    # co zostaje poza pierwszymi 3 regionami: rodziny/łańcuchy
    tk3 = np.zeros(len(T2), bool)
    for c in order[:3]: tk3 |= b(c)
    for k in ("wszystkie", "D0"):
        m = sets[k] & ~tk3
        topc = sorted([(c, int((m & b(c)).sum())) for c in hlt], key=lambda x: -x[1])[:6]
        print(f"poza R1–R3 [{k}] N={int(m.sum())} ({100*m.sum()/sets[k].sum():.2f} %):", topc)
    # stabilność regionów w LB (8 bloków)
    lb = d["ev"].lumiBlock.to_numpy(); e8 = lbblocks(lb, reg["T1"]); rows = []
    for i in range(len(e8) - 1):
        mb = (lb >= e8[i]) & (lb < e8[i + 1]); r = {"LB": f"{e8[i]}-{e8[i+1]-1}"}
        tk = np.zeros(len(T2), bool)
        for j, c in enumerate(order[:3]):
            m = b(c) & ~tk; tk |= b(c)
            r[f"R{j+1} zd%"] = round(100 * (T2 & mb & m).sum() / (T2 & mb).sum(), 2); r[f"R{j+1} D0%"] = round(100 * (sets['D0'] & mb & m).sum() / (sets['D0'] & mb).sum(), 2)
        rows.append(r)
    print(pd.DataFrame(rows).to_string(index=False))
    res[run] = {"order": order, "regions": tr.to_dict("list"), "anyK": {k: float((m & anyK).sum() / m.sum()) for k, m in sets.items()}}
json.dump(res, open(f"{OUT}/a3.json", "w"), indent=1, default=float)
