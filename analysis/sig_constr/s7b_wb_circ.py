"""SIG_constr — na ile wariant WB embeddingu jest kołowy (8:J12only, 9:MB; H06, 3–30 GeV).
Komórki WB = (pT, ΔR) najbliższego dżetu torowego pT > 2 GeV. Podział: IN = dżet w ΔR < 0,4 (jego tory leżą w stożku), OUT = reszta."""
import sys, json, numpy as np
sys.argv = [sys.argv[0]]
src = open("s7_embed.py").read().split("\nrows = []")[0]
exec(src)
out = {}
for name, run, cn in [("data_428770", 428770, "8:J12only"), ("data_427929", 427929, "9:MB")]:
    E = emb(name); S = tj(name); cS = data_classes(S, run); cE = data_classes(E, run)
    ss = cS[cn] & base_sel(S, 0.6, 3.0, 30.0) & S["ok"]
    pasS = S["v_s40_p2x"][ss] < CUT; dm = S["dm"][ss]; pt = S["pt"][ss]
    ws = IC.sbs_w(dm, pt, np.ones(ss.sum(), bool))
    ge = cE[cn] & (E["ok"] == 1) & np.isin(E["var"], [0, 2]); pasE = (E["v_s40_p2"][ge] < CUT).astype(float)
    sch = SCH["WB"]
    ie, nc = cells({v: E[v][ge] for v, _ in sch}, sch); isg, _ = cells({v: S[v][ss] for v, _ in sch}, sch)
    cnt = np.bincount(ie, minlength=nc); e_c = np.bincount(ie, weights=pasE, minlength=nc) / np.maximum(cnt, 1)
    f_c = np.bincount(isg, weights=ws, minlength=nc) / ws.sum()
    eb = float((f_c * e_c).sum())
    det = 1 - float((f_c * e_c * (1 - e_c)).sum()) / (eb * (1 - eb))        # ułamek wariancji pass/fail ustalony przez komórkę
    inS = (S["b_drn2"][ss] < 0.4) & (S["b_ptn2"][ss] >= 2); inE = (E["b_drn2"][ge] < 0.4) & (E["b_ptn2"][ge] >= 2)
    r = {"eps_WB": eb, "determinizm": det, "f_sig_IN": float(ws[inS].sum() / ws.sum()), "f_emb_rnd_IN": float(inE[E["var"][ge] == 0].mean()),
         "f_sig_komorki_eps<0.1_lub_>0.9": float(f_c[(e_c < 0.1) | (e_c > 0.9)].sum())}
    for lab, mS, mE in (("IN", inS, inE), ("OUT", ~inS, ~inE)):
        w = IC.sbs_w(dm[mS], pt[mS], pasS[mS]); et, ete, Ns, _ = IC.eff_w(w, pasS[mS].astype(float))
        # embedding przeważony komórkami w obrębie podzbioru
        g = mS; fc = np.bincount(isg[g], weights=ws[g], minlength=nc) / ws[g].sum()
        c2 = np.bincount(ie[mE], minlength=nc); e2 = np.bincount(ie[mE], weights=pasE[mE], minlength=nc) / np.maximum(c2, 1)
        ee = float((fc * e2).sum()); r[lab] = {"N_sig": float(Ns), "eps_tag": float(et), "eps_tag_err": float(ete), "eps_emb": ee, "delta": ee - float(et),
                                               "eps_emb_bez_wag": float(pasE[mE].mean())}
    out[cn] = r
    print(cn, json.dumps(r, indent=1, ensure_ascii=False))
json.dump(out, open("out/s7b_wb_circ.json", "w"), indent=1, ensure_ascii=False)
