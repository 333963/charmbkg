#!/usr/bin/env python3
"""SIG_constr A(iii): test domknięcia przeważenia SYGNAŁU MC wagami w(x) = p_dane,syg(x)/p_MC,syg(x) (klasyfikator z wagami SBS, 5 foldów po zdarzeniach).
  ISO (obszar H06, pT 3–30): cięcia ptcone40 < 1/2/4 GeV; klasy 9:MB, 8:J12only, 8:j20&!J12.
  H   (obszar PRE, pT 5–8 i 8–30): cięcia score > próg(WP95/80/50); klasy 9:MB, 8:all.
ε_MC,rw = średnia ważona (waga pT × w(x), w z modelu NIE uczonego na danym foldzie MC; dane foldu k nie wchodzą do modelu k).
Błąd wag: jackknife grupowy po foldach danych (5 modeli, każdy zastosowany do całego MC). ε_dane: fit Δm pass/fail (HESSE) gdy N_sig ≥ 150, inaczej SBS.
Przenoszenie między klasami: model z klasy A + wagi pT klasy B → przewidywane ε klasy B."""
import sys, time
from sc_common import *
import sc_ml as ML
from s3_variants import pt_weights

EVX = ["log_pt", "abs_eta", "e_nSelBL", "e_lead_pt", "e_ntrk500"]
SETS_ISO = {"pT": None,
            "EV (bez koła)": EVX + ["ann", "o_ptmax08", "h_s_away", "o_drn5"],
            "EV+F18 (bez koła)": EVX + ["ann", "o_ptmax08", "h_s_away", "o_drn5"] + [f for f in F18 if f not in EVX],
            "ISO (z kołem)": ["log_pt", "ptcone40", "ptcone100"]}
SETS_H = {"pT": None,
          "EV+ISO (bez koła)": EVX + ["ptcone40", "ptcone100"],
          "F18 (wejścia H; koło pośrednie)": F18,
          "score (z kołem)": ["log_pt", "score"]}
thr = json.load(open(os.path.join(CW, "B_trig", "runs", "lxplus", "out", "btrig_sf.json")))["thresholds"]


def eps_data(D, sel, pas, shape):
    w = IC.sbs_w(D["dm"][sel], D["pt"][sel], pas)
    e, ee, T, _ = IC.eff_w(w, pas.astype(float))
    r = {"sbs": e, "sbs_err": ee, "Ns": T}
    if T >= 150:
        f = IC.fit_pf(D["dm"][sel][pas], D["dm"][sel][~pas], shape)
        if f["valid"] and np.isfinite(f["err"]):
            r.update({"fit": f["eps"], "fit_err": f["err"], "N_fit": f["N"]})
    r["eps"], r["err"], r["method"] = (r["fit"], r["fit_err"], "fit") if "fit" in r else (e, ee, "SBS")
    return r


def weights_models(D, sel, M, ms, feats, wpt):
    """5 modeli stosunku gęstości (dane bez foldu k); zwraca wagi OOF dla MC oraz macierz [5, n_MC] wag każdego modelu na całym MC."""
    X = np.column_stack([D[f][sel] for f in feats]).astype(np.float32); XM = np.column_stack([M[f][ms] for f in feats]).astype(np.float32)
    kf = D["kf"][sel]; kM = M["kf"][ms]
    wsb = sbs_weights(D, sel)
    allw = np.ones((5, len(XM))); oof = np.ones(len(XM)); models = []
    for k in range(5):
        g = (kf != k) & (wsb != 0)
        mod = ML.fit_ratio(X[g], wsb[g], XM[kM != k], wpt[kM != k], seed=k)
        allw[k] = ML.ratio_pred(mod, XM); oof[kM == k] = allw[k][kM == k]; models.append(mod)
    return oof, allw, models


def wavg(p, w):
    e = np.average(p, weights=w)
    return float(e), float(np.sqrt(np.sum(w ** 2 * (p - e) ** 2)) / w.sum())


def study(tag, D, cm, cn, M, sthr, lo, hi, sets, cuts, shape, store):
    sel = cm & base_sel(D, sthr, lo, hi); ms = mc_sig(M, sthr, lo, hi)
    wsb = sbs_weights(D, sel)
    wpt = pt_weights(M["pt"][ms], D["pt"][sel], wsb)
    ed = {cl: eps_data(D, sel, cf(D, sel), shape) for cl, cf in cuts.items()}
    pm = {cl: cf(M, ms).astype(float) for cl, cf in cuts.items()}
    rows = []
    for sn, feats in sets.items():
        t0 = time.time()
        if feats is None:
            oof = np.ones(int(ms.sum())); allw = np.ones((5, int(ms.sum()))); models = None
        else:
            oof, allw, models = weights_models(D, sel, M, ms, feats, wpt)
        store[(tag, cn, sn, lo, hi)] = (models, feats)
        for cl in cuts:
            e_rw, s_mc = wavg(pm[cl], wpt * oof)
            ek = np.array([np.average(pm[cl], weights=wpt * allw[k]) for k in range(5)])
            s_jk = float(np.sqrt(4 / 5 * np.sum((ek - ek.mean()) ** 2))) if feats is not None else 0.0
            d = ed[cl]
            s = float(np.sqrt(d["err"] ** 2 + s_mc ** 2 + s_jk ** 2))
            rows.append({"test": tag, "cls": cn, "bin": f"{lo:g}-{hi:g}", "set": sn, "cut": cl, "Ns": d["Ns"], "neff_sbs": neff(wsb), "eps_data": d["eps"], "eps_data_err": d["err"], "method": d["method"],
                         "eps_sbs": d["sbs"], "eps_mc_raw": float(pm[cl].mean()), "eps_mc_rw": e_rw, "err_mcstat": s_mc, "err_w_jk": s_jk, "pull": (e_rw - d["eps"]) / s,
                         "ratio": d["eps"] / e_rw, "ratio_err": s / e_rw, "neff_frac": neff(wpt * oof) / max(neff(wpt), 1e-9), "w_sd": float(oof.std())})
            r = rows[-1]
            print(f"[{tag} {cn} {lo:g}-{hi:g}] {sn:32s} {cl:14s} N_sig {r['Ns']:5.0f} ε_dane {r['eps_data']:.3f}±{r['eps_data_err']:.3f} ({r['method']}) ε_MC {r['eps_mc_raw']:.3f} → rw {e_rw:.3f} ±{s_mc:.3f}(MC) ±{s_jk:.3f}(wagi) "
                  f"dane/MCrw {r['ratio']:.3f}±{r['ratio_err']:.3f} pull {r['pull']:+.1f}  n_eff/N {r['neff_frac']:.2f} ({time.time()-t0:.0f}s)", flush=True)
    return rows, (sel, ms, ed, pm)


def transfer(tag, store, keyA, D_B, cmB, cnB, M, sthr, lo, hi, cuts, shape):
    """Model wag z klasy A (keyA) zastosowany do MC przeważonego w pT do klasy B → przewidywane ε(B) vs zmierzone ε(B)."""
    models, feats = store[keyA]
    sel = cmB & base_sel(D_B, sthr, lo, hi); ms = mc_sig(M, sthr, lo, hi)
    wpt = pt_weights(M["pt"][ms], D_B["pt"][sel], sbs_weights(D_B, sel))
    XM = np.column_stack([M[f][ms] for f in feats]).astype(np.float32)
    W = np.array([ML.ratio_pred(m, XM) for m in models])
    rows = []
    for cl, cf in cuts.items():
        d = eps_data(D_B, sel, cf(D_B, sel), shape); p = cf(M, ms).astype(float)
        ek = np.array([np.average(p, weights=wpt * W[k]) for k in range(5)])
        e = float(ek.mean()); s = float(np.sqrt(d["err"] ** 2 + 4 / 5 * np.sum((ek - e) ** 2)))
        rows.append({"test": tag, "from": keyA[1], "to": cnB, "set": keyA[2], "cut": cl, "eps_data_B": d["eps"], "eps_data_B_err": d["err"], "eps_pred": e, "pull": (e - d["eps"]) / s, "ratio": d["eps"] / e, "ratio_err": s / e})
        r = rows[-1]
        print(f"[przenoszenie {tag}] {keyA[1]} → {cnB} {keyA[2]:22s} {cl:14s} ε_dane(B) {r['eps_data_B']:.3f}±{r['eps_data_B_err']:.3f} przewidziane {e:.3f}  dane/przew. {r['ratio']:.3f}±{r['ratio_err']:.3f} pull {r['pull']:+.1f}", flush=True)
    return rows


if __name__ == "__main__":
    M = load("fd_nd"); shapes = mc_shapes(M)
    D9 = load("data_427929"); c9 = data_classes(D9, 427929); D8 = load("data_428770"); c8 = data_classes(D8, 428770)
    store = {}; out = {"iso": [], "H": [], "transfer": []}
    cuts_iso = {f"ptcone40<{c:g}": (lambda X, s, c=c: X["ptcone40"][s] < c) for c in (1.0, 2.0, 4.0)}
    CL = [("9:MB", D9, c9), ("8:J12only", D8, c8), ("8:j20&!J12", D8, c8)]
    for cn, D, c in CL:
        r, _ = study("ISO", D, c[cn], cn, M, 0.6, 3.0, 30.0, SETS_ISO, cuts_iso, shapes["3-30"], store)
        out["iso"] += r; jdump(out, "s4_closure.json")
        for b, (lo, hi) in PTB.items():                     # per bin pT: tylko wagi „z kołem” i „EV” (błąd per bin dla K(b)/K(c2))
            sub = {k: SETS_ISO[k] for k in ("pT", "EV (bez koła)", "ISO (z kołem)")}
            r, _ = study("ISO", D, c[cn], cn, M, 0.6, lo, hi, sub, {"ptcone40<2": cuts_iso["ptcone40<2"]}, shapes[b], store)
            out["iso"] += r; jdump(out, "s4_closure.json")
    for (a, Da, ca) in CL:
        for (b, Db, cb) in CL:
            if a != b:
                for sn in ("EV (bez koła)", "ISO (z kołem)"):
                    out["transfer"] += transfer("ISO", store, ("ISO", a, sn, 3.0, 30.0), Db, cb[b], b, M, 0.6, 3.0, 30.0, {"ptcone40<2": cuts_iso["ptcone40<2"]}, shapes["3-30"])
    jdump(out, "s4_closure.json")
    for cn, D, c in (("9:MB", D9, c9), ("8:all", D8, c8)):
        for b, (lo, hi) in (("5-8", (5.0, 8.0)), ("8-30", (8.0, 30.0))):
            cuts_h = {f"WP{wp}": (lambda X, s, t=thr[b][wp]: X["score"][s] > t) for wp in ("0.95", "0.80", "0.50")}
            r, _ = study("H", D, c[cn], cn, M, None, lo, hi, SETS_H, cuts_h, shapes[b], store)
            out["H"] += r; jdump(out, "s4_closure.json")
