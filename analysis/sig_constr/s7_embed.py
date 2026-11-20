#!/usr/bin/env python3
"""SIG_constr B: analiza embeddingu (wyjście etapu ORCH `sigemb`: runs/lxplus/emb_<próbka>.npz; test: --in runs/emb_test).
Per klasa triggera × bin pT: ε_iso^embed (Σ pT torów w R = 0,4 wokół wstawionego D0 < 2 GeV) dla wariantów
  rnd, perpJ, oraz „jak prawdziwe D0”: wiersze rnd ∪ nearJ przeważone komórkami do rozkładu zmierzonego dla sygnału D0←D* (TRKJET, wagi SBS):
    WO (bez koła)  — max pT dżetu kol. O (z torów poza stożkiem) w ΔR < 0,8  × ΔR do najbliższego dżetu O z pT > 3
    WB             — najbliższy dżet torowy pT > 2 GeV: ΔR × pT (dla prawdziwego D0: dżet z towarzyszy, bez córek)
    WA [koło]      — pT towarzyszy w dżecie zawierającym D0 (a_jD_pt − pT_D0) × bin pT — sufit, nie test
vs ε_iso^tag = ε(własny stożek bez wolnego pionu < 2 GeV) sygnału z tagu D* w tej samej klasie. Δ_corr = ε^embed − ε^tag.
Kontrola: nośniki mb_nd vs sygnał prawdziwy fd_nd. Wynik: out/s7_embed.json, out/s7_embed.md."""
import argparse
from sc_common import *

ap = argparse.ArgumentParser(); ap.add_argument("--in", dest="inp", default=os.path.join(HERE, "runs", "lxplus")); ap.add_argument("--tag", default="")
a = ap.parse_args()
INF = 1e9
TJD = os.path.join(CW, "B_trig", "TRKJET", "runs", "lxplus")
TC = ["m", "pt", "eta", "dm", "score", "truth", "truth_ds", "event", "e_grl", "e_actual_mu", "ptcone40", "v_s40_p2x", "v_s40_p2", "o_ptmax08", "o_drn3", "b_drn2", "b_ptn2", "a_jD_pt", "nmatch", "nd_lo", "same_entry"]
SCH = {"WO (bez koła)": [("o_ptmax08", [0, 2, 3, 5, 8, 12, INF]), ("o_drn3", [0, 0.6, 0.8, 1.2, INF])],
       "WB": [("b_ptn2", [0, 2, 3, 5, 8, 12, 20, INF]), ("b_drn2", [0, 0.2, 0.4, 0.6, 0.8, 1.2, INF])],
       "WA [koło]": [("comp", [0, 0.5, 1, 2, 4, 8, 16, INF]), ("pt", [3, 5, 8, 12, 30])]}
CUT = 2.0


def tj(name):
    z = np.load(os.path.join(TJD, f"tj_{name}.npz"))
    d = {k: z[k] for k in TC}
    for k in z.files:
        if k.startswith("e_tp_"):
            d[k[2:]] = z[k]
    d["comp"] = np.maximum(d["a_jD_pt"] - d["pt"], 0.0)
    d["ok"] = (d["nmatch"] + d["nd_lo"] == 2) & (d["same_entry"] == 1)
    return d


def emb(name):
    p = os.path.join(a.inp, f"emb_{name}.npz")
    if not os.path.exists(p):
        return None
    z = np.load(p)
    d = {k: z[k] for k in z.files if k != "meta"}
    for k in list(d):
        if k.startswith("e_tp_"):
            d[k[2:]] = d[k]
    d["m"] = np.full(len(d["var"]), C.M0_D0); d["comp"] = np.maximum(d["a_jD_pt"] - d["pt"], 0.0)
    d.setdefault("e_grl", np.ones(len(d["var"]), np.int8)); d.setdefault("e_actual_mu", np.zeros(len(d["var"])))
    return d


def cells(X, sch):
    idx = np.zeros(len(X[sch[0][0]]), np.int64); nc = 1
    for v, ed in sch:
        ed = np.asarray(ed, float)
        idx = idx * (len(ed) - 1) + np.clip(np.digitize(X[v], ed) - 1, 0, len(ed) - 2); nc *= len(ed) - 1
    return idx, nc


def emb_rw(E, ge, S, ws, sch, min_n=20):
    """ε^embed po przeważeniu komórkami do rozkładu sygnału (wagi ws); zwraca (ε, błąd, udział sygnału w komórkach bez nośników)."""
    ie, nc = cells({v: E[v][ge] for v, _ in sch}, sch); isg, _ = cells(S, sch)
    pas = (E["v_s40_p2"][ge] < CUT).astype(float)
    cnt = np.bincount(ie, minlength=nc).astype(float); sp = np.bincount(ie, weights=pas, minlength=nc)
    e_c = sp / np.maximum(cnt, 1); good = cnt >= min_n
    g = good[isg]
    if ws[g].sum() <= 0:
        return np.nan, np.nan, 1.0
    e = float((ws[g] * e_c[isg[g]]).sum() / ws[g].sum())
    v_sig = float((ws[g] ** 2 * (e_c[isg[g]] - e) ** 2).sum() / ws[g].sum() ** 2)
    f_c = np.bincount(isg[g], weights=ws[g], minlength=nc) / ws[g].sum()
    v_emb = float(np.sum(f_c ** 2 * e_c * (1 - e_c) / np.maximum(cnt, 1)))
    return e, float(np.sqrt(v_sig + v_emb)), float(1 - ws[g].sum() / ws.sum())


rows = []
JOBS = [("data_427929", 427929, ["9:MB", "9:MB&j20&!J12", "9:MB&!j20&!J12"]), ("data_428770", 428770, ["8:J12only", "8:j20&!J12", "8:MBTSbased", "8:other"])]
for name, run, cl in JOBS:
    E = emb(name)
    if E is None:
        print("brak", name); continue
    S = tj(name); cS = data_classes(S, run); cE = data_classes(E, run)
    for cn in cl:
        for b, (lo, hi) in [("3-30", (3.0, 30.0))] + list(PTB.items()):
            ss = cS[cn] & base_sel(S, 0.6, lo, hi) & S["ok"]
            if ss.sum() < 50:
                continue
            pas = S["v_s40_p2x"][ss] < CUT
            w = IC.sbs_w(S["dm"][ss], S["pt"][ss], pas); et, ete, Ns, _ = IC.eff_w(w, pas.astype(float))
            ws = IC.sbs_w(S["dm"][ss], S["pt"][ss], np.ones(ss.sum(), bool))
            Ssub = {k: S[k][ss] for k in ("o_ptmax08", "o_drn3", "b_drn2", "b_ptn2", "comp", "pt")}
            r = {"cls": cn, "bin": b, "Ns": Ns, "eps_tag": et, "eps_tag_err": ete, "n_carriers": int((cE[cn] & (E["var"] == 0)).sum())}
            for v, vn in ((0, "rnd"), (1, "perpJ")):
                g = cE[cn] & (E["var"] == v) & (E["ok"] == 1)
                p = E["v_s40_p2"][g] < CUT
                r[vn] = float(p.mean()) if g.sum() else np.nan; r[vn + "_err"] = float(np.sqrt(p.mean() * (1 - p.mean()) / max(g.sum(), 1))) if g.sum() else np.nan
            ge = cE[cn] & (E["ok"] == 1) & np.isin(E["var"], [0, 2])
            for sn, sch in SCH.items():
                g2 = ge & ((E["pt"] >= lo) & (E["pt"] < hi) if sn.startswith("WA") else True)
                r[sn], r[sn + "_err"], r[sn + "_miss"] = emb_rw(E, g2, Ssub, ws, sch) if g2.sum() else (np.nan, np.nan, 1.0)
            for k in ("rnd", "perpJ") + tuple(SCH):
                r["d_" + k] = r[k] - et; r["d_" + k + "_err"] = float(np.hypot(r[k + "_err"], ete))
            rows.append(r)
            print(f"{cn:16s} {b:5s} N_sig {Ns:5.0f} nośników {r['n_carriers']:7d} | ε_tag {et:.3f}±{ete:.3f} | rnd {r['rnd']:.3f} perpJ {r['perpJ']:.3f} | " +
                  " | ".join(f"{sn} {r[sn]:.3f}±{r[sn+'_err']:.3f} (brak {r[sn+'_miss']:.2f})" for sn in SCH), flush=True)
# kontrola MC: nośniki mb_nd, sygnał prawdziwy fd_nd
E = emb("mb_nd")
if E is not None:
    S = tj("fd_nd")
    for b, (lo, hi) in [("3-30", (3.0, 30.0))] + list(PTB.items()):
        ss = np.isin(S["truth"], [0, 1]) & np.isin(S["truth_ds"], [0, 1]) & (np.abs(S["m"] - C.M0_D0) < MWIN) & (S["dm"] > 0.1435) & (S["dm"] < 0.1475) & (S["score"] > 0.6) & (S["pt"] >= lo) & (S["pt"] < hi) & S["ok"]
        pas = S["v_s40_p2x"][ss] < CUT; et = float(pas.mean()); ete = float(np.sqrt(et * (1 - et) / ss.sum())); ws = np.ones(int(ss.sum()))
        Ssub = {k: S[k][ss] for k in ("o_ptmax08", "o_drn3", "b_drn2", "b_ptn2", "comp", "pt")}
        r = {"cls": "MC: mb_nd ← fd_nd", "bin": b, "Ns": int(ss.sum()), "eps_tag": et, "eps_tag_err": ete, "n_carriers": int((E["var"] == 0).sum())}
        for v, vn in ((0, "rnd"), (1, "perpJ")):
            g = (E["var"] == v) & (E["ok"] == 1); p = E["v_s40_p2"][g] < CUT
            r[vn] = float(p.mean()); r[vn + "_err"] = float(np.sqrt(p.mean() * (1 - p.mean()) / g.sum()))
        ge = (E["ok"] == 1) & np.isin(E["var"], [0, 2])
        for sn, sch in SCH.items():
            g2 = ge & ((E["pt"] >= lo) & (E["pt"] < hi) if sn.startswith("WA") else True)
            r[sn], r[sn + "_err"], r[sn + "_miss"] = emb_rw(E, g2, Ssub, ws, sch)
        for k in ("rnd", "perpJ") + tuple(SCH):
            r["d_" + k] = r[k] - et; r["d_" + k + "_err"] = float(np.hypot(r[k + "_err"], ete))
        rows.append(r)
        print(f"{r['cls']:16s} {b:5s} N {r['Ns']:6d} | ε_tag {et:.3f} | rnd {r['rnd']:.3f} perpJ {r['perpJ']:.3f} | " + " | ".join(f"{sn} {r[sn]:.3f}±{r[sn+'_err']:.3f}" for sn in SCH), flush=True)
jdump({"rows": rows, "input": a.inp, "cut": CUT}, f"s7_embed{a.tag}.json")
L = [f"# SIG_constr — embedding: ε_iso^embed vs ε_iso^tag (wejście `{a.inp}`)", "",
     "Stożek: Σ pT wszystkich torów drzewa (pT > 0,2 GeV) w ΔR < 0,4 < 2 GeV. ε_tag: własny stożek prawdziwego D⁰ bez wolnego pionu (SBS). Δ = ε^embed − ε^tag.", "",
     "| klasa | bin pT | N_sig | nośników | ε_tag | rnd | perpJ | WO (bez koła) | WB | WA [koło] | Δ(rnd) | Δ(WO) | Δ(WB) |", "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
f = lambda x, e=None: "–" if x is None or not np.isfinite(x) else (f"{x:.3f}".replace(".", ",") + ("" if e is None else " ± " + f"{e:.3f}".replace(".", ",")))  # noqa: E731
for r in rows:
    L.append(f"| {r['cls']} | {r['bin']} | {r['Ns']:.0f} | {r['n_carriers']} | {f(r['eps_tag'], r['eps_tag_err'])} | {f(r['rnd'], r['rnd_err'])} | {f(r['perpJ'], r['perpJ_err'])} | "
             f"{f(r['WO (bez koła)'], r['WO (bez koła)_err'])} | {f(r['WB'], r['WB_err'])} | {f(r['WA [koło]'], r['WA [koło]_err'])} | {f(r['d_rnd'], r['d_rnd_err'])} | "
             f"{f(r['d_WO (bez koła)'], r['d_WO (bez koła)_err'])} | {f(r['d_WB'], r['d_WB_err'])} |")
open(os.path.join(OUT, f"s7_embed{a.tag}.md"), "w").write("\n".join(L) + "\n")
