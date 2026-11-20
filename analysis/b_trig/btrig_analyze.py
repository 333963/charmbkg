#!/usr/bin/env python3
"""B_trig — analiza odporności wyników B na trigger i niezależny run (lokalnie, po ekstrakcji).

Podkomendy:
  sf      SF klasyfikatora H z tagu D* per klasa triggera × bin pT × WP (fit Δm pass/fail jak EV / charmbkg.sf_rel),
          ε_MC per klasa z emulacją triggera (btrig_emul), test χ² zgodności klas rozłącznych + pulle vs „all” i vs nominał.
  epsb    ε_B pod pikiem per klasa: SB symetryczne i fit exp(pol2) (charmbkg.val_epsb), MC tła (truth ≥ 8) z emulacją,
          R_c = ε_B,dane / ε_B,MC; decyzja „reweighting MC→dane per trigger” (bramka G7).
  iso     SF izolacji (ptcone40, ptcone40/pT) per klasa (jak ISO_sf; preselekcja score_H > 0,6).
  effmap  mapy ε_k(nSelBL, lead_pt) z 427929 (odniesienie MB427929) + test domknięcia ułamków klas w 428770 (G5).

Wejścia danych: --data PLIK@RUN (powtarzalne). PLIK = npz z btrig_extract (kolumny e_tp_*) albo ev_data.npz z EV
+ --event-table trg_<run>.npz (trg_root_extract; złączenie po (run, event)). --toy-trig SEED: syntetyczne bity (testy).
MC: --mc fd npz (pary z prawdą; e_nSelBL/e_lead_pt → emulacja, inaczej wagi 1 i flaga), --mc-bkg mb npz (epsb).
Wyjście: --out KATALOG → btrig_<cmd>.json + btrig_<cmd>.md.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import btrig_common as C  # noqa: E402
import btrig_emul as EM  # noqa: E402

C.add_paths()


# ----------------------------------------------------------------------------------------------------------- wejście
def load(path):
    z = np.load(path, allow_pickle=False)
    d = {k: z[k] for k in z.files if k != "meta"}
    d["meta"] = json.loads(str(z["meta"])) if "meta" in z.files else {}
    return d


def _lead_gev(v):
    v = np.asarray(v, float)
    return v / 1000.0 if len(v) and np.nanmedian(v) > 200 else v


def toy_bits(run, event, seed, n_act=None):
    """Syntetyczne bity (TYLKO testy): J12 częściej w aktywnych zdarzeniach (n_act), MBTS-chains z prescale."""
    a = np.asarray(n_act, float) if n_act is not None else np.zeros(len(event))
    r = np.random.default_rng(seed).random((4, len(event)))
    pj = 0.45 + 0.4 * np.tanh(a / 30.0)
    return {"tp_HLT_noalg_L1J12": (r[0] < pj).astype(np.int8),
            "tp_HLT_j20_L1MBTS_2": (r[1] < 0.18).astype(np.int8),
            "tp_HLT_mb_sptrk_pt4_L1MBTS_2": (r[2] < 0.13).astype(np.int8),
            "tp_HLT_mb_mbts_L1MBTS_1": (r[3] < 0.001).astype(np.int8)}


def bits_for(D, prefix, run_key, ev_key, args, run):
    """Bity trigger (+grl) dla wierszy: e_tp_* z npz → tabela zdarzeń → toy. prefix '' (pary) lub 'all_'."""
    keys = [k for k in D if k.startswith(prefix + "e_tp_") or k == prefix + "e_grl"]
    if keys:
        return {k[len(prefix) + 2:]: D[k] for k in keys}, "npz"
    if args.event_table:
        for p in args.event_table:
            T = load(p)
            if not len(T.get("run", [])) or int(np.median(T["run"])) != int(run):
                continue
            cols = [k for k in T if k.startswith("tp_") or k == "grl"]
            J = C.join_event_table(D[run_key], D[ev_key], T, cols)
            nf = float((~J["found"]).mean()) if len(J["found"]) else 0.0
            out = {k: np.nan_to_num(J[k]).astype(np.int8) for k in cols}
            out["_notfound"] = nf
            return out, f"event_table:{os.path.basename(p)}"
    if args.toy_trig is not None:
        return toy_bits(D[run_key], D[ev_key], args.toy_trig, D.get(prefix + "n_act")), "TOY"
    return None, "brak"


def data_sets(args):
    """[(run, D, classes(dict masek dla par), src)]"""
    out = []
    for spec in args.data:
        path, run = spec.rsplit("@", 1)
        D = load(path)
        bits, src = bits_for(D, "", "run", "event", args, run)
        if bits is None:
            raise SystemExit(f"STOP: {path}: brak bitów triggera (e_tp_* / --event-table / --toy-trig)")
        nf = bits.pop("_notfound", 0.0)
        cls = C.classify(bits, int(run))
        out.append({"run": int(run), "D": D, "cls": cls, "src": src, "notfound": nf, "path": path})
    return out


def mc_weights(M, cls, args, prefix=""):
    """Wagi emulacji dla MC (metoda A, mapa z --emul-map). Brak zmiennych/ mapy → 1 i flaga."""
    n = len(M[prefix + "pt"])
    ns, ld = M.get(prefix + "e_nSelBL"), M.get(prefix + "e_lead_pt")
    if ns is None or ld is None or not args.emul_map:
        return np.ones(n), "none"
    em = EM.load_map(args.emul_map)
    return EM.class_prob(em, ns, _lead_gev(ld), cls), "effmap"


def thresholds(M):
    sc, pt, tr, m = M["all_score"], M["all_pt"], M["all_truth"], M["all_m"]
    sig = np.isin(tr, [0, 1]) & (np.abs(m - C.M0_D0) < C.NOM_WIN)
    thr = {}
    for b, (lo, hi) in C.BINS.items():
        s = sc[sig & (pt >= lo) & (pt < hi)]
        thr[b] = {wp: float(np.quantile(s, 1 - float(wp))) if len(s) else np.nan for wp in C.WPS}
    return thr


def _fit_job(job):
    key, dp, df_, shape, interval = job
    return key, C.fit_eps(dp, df_, shape, interval=interval)


# ----------------------------------------------------------------------------------------------------------- sf
def cmd_sf(args):
    from charmtag import dstar_refit as DR
    M = load(args.mc)
    thr = thresholds(M)
    f_sig = np.isin(M["truth"], [0, 1]) & np.isin(M["truth_ds"], [0, 1]) & (np.abs(M["m"] - C.M0_D0) < C.NOM_WIN)
    shapes = {}
    for b, (lo, hi) in C.BINS.items():
        bf = f_sig & (M["pt"] >= lo) & (M["pt"] < hi)
        shapes[b] = DR.shape_from_mc(M["dm"][bf] * 1e3)
    wps = args.wps.split(",")
    sets = data_sets(args)
    jobs, eps_mc, meta_cls = [], {}, {}
    for S in sets:
        D, run = S["D"], S["run"]
        nom = np.abs(D["m"] - C.M0_D0) < C.NOM_WIN
        for c, mask in S["cls"].items():
            if c.startswith("_"):
                continue
            w, wsrc = mc_weights(M, c, args)
            meta_cls[f"{run}:{c}"] = {"n_pairs": int(mask.sum()), "frac": float(mask.mean()), "mc_weights": wsrc}
            for b, (lo, hi) in C.BINS.items():
                md = mask & nom & (D["pt"] >= lo) & (D["pt"] < hi)
                mf = f_sig & (M["pt"] >= lo) & (M["pt"] < hi)
                for wp in wps:
                    t = thr[b][wp]
                    pas = D["score"] >= t
                    jobs.append(((run, c, b, wp), D["dm"][md & pas], D["dm"][md & ~pas], shapes[b], not args.fast))
                    ww = w[mf]
                    pm = (M["score"][mf] >= t)
                    e = float((ww * pm).sum() / ww.sum()) if ww.sum() > 0 else np.nan
                    e0 = float(pm.mean()) if mf.any() else np.nan
                    eps_mc[(run, c, b, wp)] = (e, e0, EM.n_eff(ww) if len(ww) else 0.0)
    with Pool(args.nproc) as P:
        fits = dict(P.map(_fit_job, jobs))
    rows = []
    for (run, c, b, wp), f in fits.items():
        e, e0, ne = eps_mc[(run, c, b, wp)]
        err = C.eps_sym_err(f)
        sf = f["eps"] / e if e and np.isfinite(f["eps"]) else np.nan
        mc_err = np.sqrt(e * (1 - e) / ne) / e if ne > 0 and 0 < e < 1 else 0.0
        rows.append({"run": run, "cls": c, "bin": b, "wp": wp, "eps_data": f["eps"], "eps_data_err": err,
                     "determined": f["determined"], "N_sig": f["N"], "N_sig_err": f["N_err"], "eps_mc": e, "eps_mc_unw": e0,
                     "mc_neff": ne, "SF": sf, "SF_err": float(np.hypot(err / e, sf * mc_err)) if e else np.nan,
                     "n_pass": f["n_pass"], "n_fail": f["n_fail"]})
    tests = sf_tests(rows, wps)
    res = {"inputs": {"data": [{k: S[k] for k in ("run", "src", "notfound", "path")} for S in sets], "mc": args.mc,
                      "emul_map": args.emul_map}, "thresholds": thr, "shapes": shapes, "classes": meta_cls, "rows": rows,
           "tests": tests}
    write(args, "sf", res, md_sf(res))


def _get(rows, run, c, b, wp):
    for r in rows:
        if r["run"] == run and r["cls"] == c and r["bin"] == b and r["wp"] == wp:
            return r
    return None


def _ok(r):
    return r is not None and r["determined"] and np.isfinite(r["SF"]) and np.isfinite(r["SF_err"]) and r["SF_err"] > 0


def sf_tests(rows, wps):
    """χ² zgodności klas rozłącznych (428770: J12only, MBTSbased, other; + MB427929), pulle podzbiór vs all, run vs run,
    vs nominał B (SF95)."""
    T = []
    runs = sorted({r["run"] for r in rows})
    for b in C.BINS:
        for wp in wps:
            disj = [_get(rows, 428770, c, b, wp) for c in C.DISJOINT_428770] + [_get(rows, 427929, "MB427929", b, wp)]
            disj = [r for r in disj if _ok(r)]
            t = {"bin": b, "wp": wp, "chi2_disjoint": C.chi2_consistency([r["SF"] for r in disj], [r["SF_err"] for r in disj]),
                 "members": [f'{r["run"]}:{r["cls"]}' for r in disj], "pulls": {}}
            for run in runs:
                ra = _get(rows, run, "all" if run != 427929 else "all427929", b, wp)
                for r in rows:
                    if r["run"] == run and r["bin"] == b and r["wp"] == wp and _ok(r) and _ok(ra) and r is not ra:
                        t["pulls"][f"{run}:{r['cls']}_vs_all"] = C.pull_subset(r["SF"], r["SF_err"], ra["SF"], ra["SF_err"])
            a8, a9 = _get(rows, 428770, "all", b, wp), _get(rows, 427929, "MB427929", b, wp)
            if _ok(a8) and _ok(a9):
                t["pulls"]["MB427929_vs_428770all"] = C.pull_indep(a9["SF"], a9["SF_err"], a8["SF"], a8["SF_err"])
            s8, s9 = _get(rows, 428770, "sptrk4", b, wp), _get(rows, 427929, "sptrk4_427929", b, wp)
            if _ok(s8) and _ok(s9):
                t["pulls"]["sptrk4_run_vs_run"] = C.pull_indep(s9["SF"], s9["SF_err"], s8["SF"], s8["SF_err"])
            nv, ne = C.NOMINAL_SF95.get(b, (np.nan, np.nan))
            if wp == "0.95" and np.isfinite(nv):
                if _ok(a9):
                    t["pulls"]["MB427929_vs_nominal"] = C.pull_indep(a9["SF"], a9["SF_err"], nv, ne)
                    t["MB427929_vs_1"] = (a9["SF"] - 1.0) / a9["SF_err"]
                if _ok(a8):
                    t["closure_all_vs_nominal"] = {"SF_all": a8["SF"], "nominal": nv, "diff": a8["SF"] - nv}
            T.append(t)
    return T


def md_sf(res):
    L = ["# B_trig — SF klasyfikatora H z tagu D* per klasa triggera", "",
         f"Wejścia: {[(d['run'], d['src'], round(d['notfound'], 5)) for d in res['inputs']['data']]}; MC {res['inputs']['mc']}; "
         f"emulacja {res['inputs']['emul_map'] or 'BRAK (wagi 1)'}.", "",
         "| run | klasa | bin | WP | N_sig | ε_dane | ε_MC (emul / bez) | SF | det |", "|---|---|---|---|---|---|---|---|---|"]
    for r in sorted(res["rows"], key=lambda r: (r["bin"], r["wp"], r["run"], r["cls"])):
        L.append(f"| {r['run']} | {r['cls']} | {r['bin']} | {r['wp']} | {r['N_sig']:.0f}±{r['N_sig_err']:.0f} | "
                 f"{r['eps_data']:.3f}±{r['eps_data_err']:.3f} | {r['eps_mc']:.3f} / {r['eps_mc_unw']:.3f} | "
                 f"**{r['SF']:.2f}±{r['SF_err']:.2f}** | {'tak' if r['determined'] else 'NIE'} |")
    L += ["", "## Testy zgodności", "", "| bin | WP | χ²/ndf (p) klasy rozłączne | pulle |", "|---|---|---|---|"]
    for t in res["tests"]:
        c = t["chi2_disjoint"]
        pl = ", ".join(f"{k} {v:+.1f}" for k, v in t["pulls"].items())
        L.append(f"| {t['bin']} | {t['wp']} | {c['chi2']:.1f}/{c['ndf']} ({c['p']:.2f}) | {pl} |")
    return "\n".join(L) + "\n"


# ----------------------------------------------------------------------------------------------------------- epsb
TOL_EPSB = {"0.80": 0.05, "0.70": 0.10, "0.50": 0.25}     # systematyka kształtu SALAD (§3a) — próg istotności różnic


def cmd_epsb(args):
    from charmbkg import val_epsb as V
    M = load(args.mc)
    thr = thresholds(M)
    Mb = load(args.mc_bkg) if args.mc_bkg else None
    reg = V.default_regions("D0")
    salad = json.load(open(args.salad_json)) if args.salad_json else {}
    wps = args.wps.split(",")
    rows = []
    for spec in args.data:
        path, run = spec.rsplit("@", 1)
        D = load(path)
        if "all_run" not in D:
            raise SystemExit(f"STOP: {path} bez all_run/all_event (EV ev_data.npz) — epsb wymaga btrig_extract")
        bits, src = bits_for(D, "all_", "all_run", "all_event", args, run)
        bits.pop("_notfound", None)
        cls = C.classify(bits, int(run))
        ref_c = "all" if int(run) != 427929 else "all427929"
        shp = {}                                   # kształty exp(pol2) z klasy odniesienia per bin × WP
        for b, (lo, hi) in C.BINS.items():
            s0 = cls[ref_c] & (D["all_pt"] >= lo) & (D["all_pt"] < hi)
            for wp in wps:
                shp[(b, wp)] = C.pol2_shapes(D["all_m"][s0], D["all_score"][s0] >= thr[b][wp], C.M0_D0, reg["sr_half"],
                                             reg["sb_min"], reg["lo"], reg["hi"])
        for c, mask in cls.items():
            if c.startswith("_"):
                continue
            for b, (lo, hi) in C.BINS.items():
                sel = mask & (D["all_pt"] >= lo) & (D["all_pt"] < hi)
                m = D["all_m"][sel]
                for wp in wps:
                    passed = D["all_score"][sel] >= thr[b][wp]
                    # liczba binów adaptacyjna: fit log N na binach z małą liczbą wpisów jest obciążony (test toy 17.11:
                    # 100 binów przy ~300 pass w SB → +50 % ε_B w małych klasach); ≥ ~30 wpisów pass / bin
                    nb = int(np.clip(passed.sum() / 30, 10, 100))
                    pol2 = V.eps_b_sr_pol2(m, passed, C.M0_D0, reg["sr_half"], reg["lo"], reg["hi"], nbins=nb)
                    sym = V.eps_b_sb_sym(m, passed, C.M0_D0, reg["sb_min"], reg["sb_max"])
                    ml = C.eps_b_pol2_ml(m, passed, C.M0_D0, reg["sr_half"], reg["sb_min"], reg["lo"], reg["hi"], shapes=shp[(b, wp)])
                    mlf = C.eps_b_pol2_ml(m, passed, C.M0_D0, reg["sr_half"], reg["sb_min"], reg["lo"], reg["hi"])
                    r = {"run": int(run), "cls": c, "bin": b, "wp": wp, "src": src, "n": int(sel.sum()), "nbins": nb,
                         "pol2_lsq": pol2["eps"],                 # val_epsb (LSQ log N) — tylko porównanie
                         "pol2_free": mlf["eps"],                 # kształt wolny per klasa (niestabilny w małych klasach)
                         "pol2": ml["eps"], "pol2_stat_rel": ml["err"] / ml["eps"] if ml["eps"] else np.nan, "sym": sym["eps"],
                         "sym_over_pol2": sym["eps"] / ml["eps"] if ml["eps"] else np.nan,
                         "salad": salad.get(b, {}).get(wp, np.nan)}
                    if Mb is not None:
                        mt = (Mb["all_truth"] >= 8) & (Mb["all_pt"] >= lo) & (Mb["all_pt"] < hi) & \
                             (np.abs(Mb["all_m"] - C.M0_D0) < reg["sr_half"])
                        w, wsrc = mc_weights({k[4:]: v for k, v in Mb.items() if k.startswith("all_")} | {"pt": Mb["all_pt"]}, c, args)
                        ww = w[mt]
                        r["mc_bkg"] = float((ww * (Mb["all_score"][mt] >= thr[b][wp])).sum() / max(ww.sum(), 1e-12))
                        r["mc_bkg_w"] = wsrc
                        r["R_data_mc"] = r["pol2"] / r["mc_bkg"] if r["mc_bkg"] > 0 else np.nan
                    if np.isfinite(r["salad"]):
                        r["pol2_over_salad"] = r["pol2"] / r["salad"] if r["salad"] else np.nan
                    rows.append(r)
    # decyzja G7: rozrzut pol2/ref między klasami rozłącznymi vs tolerancja kształtu SALAD
    dec = []
    for b in C.BINS:
        for wp in wps:
            key = "R_data_mc" if args.mc_bkg else ("pol2_over_salad" if salad else "pol2")
            R = [r for r in rows if r["bin"] == b and r["wp"] == wp and r["cls"] in C.DISJOINT_428770 + ["MB427929"]
                 and np.isfinite(r.get(key, np.nan))]
            ra = [r for r in rows if r["bin"] == b and r["wp"] == wp and r["cls"] == "all"]
            if len(R) < 2 or not ra:
                continue
            ref = ra[0].get(key)
            dev = [abs(r[key] / ref - 1) for r in R]
            chi = C.chi2_consistency([r[key] for r in R], [abs(r[key]) * r["pol2_stat_rel"] for r in R])
            tol = TOL_EPSB.get(wp, 0.05)
            need = (chi["p"] < 0.01) and (max(dev) > tol)
            dec.append({"bin": b, "wp": wp, "quantity": key, "max_rel_dev": max(dev), "tol": tol, "chi2": chi,
                        "per_trigger_reweighting_needed": bool(need)})
    write(args, "epsb", {"rows": rows, "decision": dec, "thresholds": thr}, md_epsb(rows, dec))


def md_epsb(rows, dec):
    L = ["# B_trig — ε_B pod pikiem per klasa triggera", "", "| run | klasa | bin | WP | n | pol2 | sym/pol2 | MC tła | R=dane/MC | pol2/SALAD |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    f = lambda x: f"{x:.3f}" if isinstance(x, float) and np.isfinite(x) else "–"  # noqa: E731
    for r in rows:
        L.append(f"| {r['run']} | {r['cls']} | {r['bin']} | {r['wp']} | {r['n']} | {f(r['pol2'])}±{f(r['pol2'] * r['pol2_stat_rel'])} | "
                 f"{f(r['sym_over_pol2'])} | {f(r.get('mc_bkg', np.nan))} | {f(r.get('R_data_mc', np.nan))} | {f(r.get('pol2_over_salad', np.nan))} |")
    L += ["", "## Decyzja G7 (reweighting per trigger?)", "", "| bin | WP | wielkość | max |Δ|/all | tol | χ² p | per trigger? |", "|---|---|---|---|---|---|---|"]
    for d in dec:
        L.append(f"| {d['bin']} | {d['wp']} | {d['quantity']} | {d['max_rel_dev']:.3f} | {d['tol']} | {d['chi2']['p']:.3f} | "
                 f"{'TAK' if d['per_trigger_reweighting_needed'] else 'nie'} |")
    return "\n".join(L) + "\n"


# ----------------------------------------------------------------------------------------------------------- iso
ISO_CUTS = [("ptcone40", 1.0), ("ptcone40", 2.0), ("ptcone40", 4.0), ("rel40", 0.25), ("rel40", 0.5)]
ISO_PTB = {"3-30": (3, 30), **C.BINS}


def cmd_iso(args):
    from charmtag import dstar_refit as DR
    M = load(args.mc)
    if "ptcone40" not in M:
        msg = ("# B_trig iso — POMINIĘTE\n\nBrak kolumny ptcone40 w MC/danych (h5 bez Dmeson_ptcone40). WYMAGANIE dla konwersji h5_v20: "
               "Dmeson_ptcone{40,60,80,100} (GeV) na wierszu kandydata (PROCEDURA §5).\n")
        write(args, "iso", {"skipped": "brak ptcone40"}, msg)
        return
    for X in (M,):
        X["rel40"] = X["ptcone40"] / X["pt"]
    mf0 = np.isin(M["truth"], [0, 1]) & (np.abs(M["m"] - C.M0_D0) < 0.040) & (M["dm"] > 0.1435) & (M["dm"] < 0.1475) & (M["score"] > 0.6)
    shape = DR.shape_from_mc(M["dm"][mf0 & (M["pt"] >= 3)] * 1e3)
    jobs, emc = [], {}
    for S in data_sets(args):
        D, run = S["D"], S["run"]
        if "ptcone40" not in D:
            continue
        D["rel40"] = D["ptcone40"] / D["pt"]
        base = (np.abs(D["m"] - C.M0_D0) < 0.040) & (D["score"] > 0.6)
        for c, mask in S["cls"].items():
            if c.startswith("_"):
                continue
            w, wsrc = mc_weights(M, c, args)
            for pb, (lo, hi) in ISO_PTB.items():
                md = mask & base & (D["pt"] >= lo) & (D["pt"] < hi)
                mf = mf0 & (M["pt"] >= lo) & (M["pt"] < hi)
                for v, cut in ISO_CUTS:
                    pas = D[v] < cut
                    jobs.append(((run, c, pb, f"{v}<{cut}"), D["dm"][md & pas], D["dm"][md & ~pas], shape, not args.fast))
                    ww = w[mf]
                    emc[(run, c, pb, f"{v}<{cut}")] = (float((ww * (M[v][mf] < cut)).sum() / max(ww.sum(), 1e-12)), wsrc)
    with Pool(args.nproc) as P:
        fits = dict(P.map(_fit_job, jobs))
    rows = []
    for k, f in fits.items():
        e, wsrc = emc[k]
        err = C.eps_sym_err(f)
        rows.append({"run": k[0], "cls": k[1], "bin": k[2], "cut": k[3], "eps_data": f["eps"], "eps_data_err": err, "eps_mc": e,
                     "mc_w": wsrc, "SF": f["eps"] / e if e else np.nan, "SF_err": err / e if e else np.nan,
                     "determined": f["determined"], "N_sig": f["N"]})
    L = ["# B_trig — SF izolacji per klasa triggera (preselekcja score_H > 0,6, |m−m0| < 40 MeV)", "",
         "| run | klasa | pT | cięcie | N_sig | ε_dane | ε_MC | SF | det |", "|---|---|---|---|---|---|---|---|---|"]
    for r in sorted(rows, key=lambda r: (r["cut"], r["bin"], r["run"], r["cls"])):
        L.append(f"| {r['run']} | {r['cls']} | {r['bin']} | {r['cut']} | {r['N_sig']:.0f} | {r['eps_data']:.3f}±{r['eps_data_err']:.3f} | "
                 f"{r['eps_mc']:.3f} ({r['mc_w']}) | {r['SF']:.2f}±{r['SF_err']:.2f} | {'tak' if r['determined'] else 'NIE'} |")
    write(args, "iso", {"rows": rows}, "\n".join(L) + "\n")


# ----------------------------------------------------------------------------------------------------------- effmap
def cmd_effmap(args):
    """Mapy ε_k z 427929 (tabela zdarzeń: trg_root_extract lub ev_* z btrig_extract) + domknięcie w 428770."""
    def table(p):
        T = load(p)
        if "ev_run" in T:
            T = {k[3:]: v for k, v in T.items() if k.startswith("ev_")}
        return T
    T9 = table(args.ref_table)
    ref = T9["tp_HLT_mb_mbts_L1MBTS_1"].astype(bool) & (T9["grl"].astype(bool) if "grl" in T9 else True)
    fired = {k: T9[f"tp_{ch}"].astype(bool) for k, ch in EM.CHAINS_EMUL.items() if f"tp_{ch}" in T9}
    em = EM.effmap(T9["nSelBL"], _lead_gev(T9["lead_pt"]), ref, fired)
    os.makedirs(args.out, exist_ok=True)
    EM.save_map(em, os.path.join(args.out, "effmap_427929.json"), {"ref": "MB427929", "n_ref": int(ref.sum())})
    res = {"n_ref": int(ref.sum()), "eps_global": {k: v["eps_global"] for k, v in em.items()}}
    if args.closure_table:
        T8 = table(args.closure_table)
        ns, ld = T8["nSelBL"], _lead_gev(T8["lead_pt"])
        cls = C.classify({k: v for k, v in T8.items() if k.startswith("tp_")}, 428770)
        ps = {}
        for k, ch in EM.CHAINS_EMUL.items():
            if f"ps_{ch}" in T8:
                ps[k] = float(np.median(T8[f"ps_{ch}"]))
        pj, pm = EM.class_prob(em, ns, ld, "J12only", ps), EM.class_prob(em, ns, ld, "MBTSbased", ps)
        pred = pm / np.clip(pj + pm, 1e-12, None)
        sel = cls["J12only"] | cls["MBTSbased"]
        i, j = EM._bin2(ns, ld)
        cells = []
        for a in range(len(EM.NSEL_EDGES) - 1):
            for bb in range(len(EM.LEAD_EDGES) - 1):
                c = sel & (i == a) & (j == bb)
                if c.sum() < 50:
                    continue
                obs = float(cls["MBTSbased"][c].mean())
                pr = float(pred[c].mean())
                er = np.sqrt(obs * (1 - obs) / c.sum()) + 1e-6
                cells.append({"nsel_bin": a, "lead_bin": bb, "n": int(c.sum()), "obs": obs, "pred": pr, "pull": (obs - pr) / er})
        pulls = np.array([c["pull"] for c in cells])
        res["closure_428770"] = {"cells": cells, "chi2": float((pulls ** 2).sum()), "ndf": len(cells),
                                 "max_rel_dev": float(max([abs(c["pred"] / c["obs"] - 1) for c in cells if c["obs"] > 0] or [np.nan]))}
    L = ["# B_trig — mapy emulacji triggera (427929, odniesienie MB427929)", "", f"n_ref = {res['n_ref']}; ε_k globalne: {res['eps_global']}"]
    if "closure_428770" in res:
        cl = res["closure_428770"]
        L.append(f"\nDomknięcie 428770 (ułamek MBTSbased wśród J12only∪MBTSbased per komórka x): χ²/ndf = {cl['chi2']:.1f}/{cl['ndf']}, "
                 f"max |pred/obs − 1| = {cl['max_rel_dev']:.2f}")
    write(args, "effmap", res, "\n".join(L) + "\n")


# ----------------------------------------------------------------------------------------------------------- wspólne
def write(args, name, res, md):
    os.makedirs(args.out, exist_ok=True)

    def conv(o):
        if isinstance(o, dict):
            return {str(k): conv(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [conv(x) for x in o]
        if isinstance(o, np.ndarray):
            return conv(o.tolist())
        if isinstance(o, (np.floating, float)):
            return float(o) if np.isfinite(o) else None
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, np.bool_):
            return bool(o)
        return o
    json.dump(conv(res), open(os.path.join(args.out, f"btrig_{name}.json"), "w"), indent=1)
    open(os.path.join(args.out, f"btrig_{name}.md"), "w").write(md)
    print(f"zapisano {args.out}/btrig_{name}.{{json,md}}")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["sf", "epsb", "iso", "effmap"])
    ap.add_argument("--data", action="append", default=[], help="PLIK.npz@RUN (powtarzalne)")
    ap.add_argument("--event-table", action="append", default=[], help="trg_<run>.npz (trg_root_extract)")
    ap.add_argument("--toy-trig", type=int, default=None, help="TYLKO testy: syntetyczne bity")
    ap.add_argument("--mc", help="fd npz (pary z prawdą)")
    ap.add_argument("--mc-bkg", help="mb npz (all_* z truth ≥ 8) — epsb")
    ap.add_argument("--salad-json", help="{bin: {wp: ε_B SALAD}} — epsb")
    ap.add_argument("--emul-map", help="effmap_427929.json (btrig_analyze effmap)")
    ap.add_argument("--ref-table", help="effmap: tabela zdarzeń 427929")
    ap.add_argument("--closure-table", help="effmap: tabela zdarzeń 428770 (test domknięcia)")
    ap.add_argument("--wps", default="0.95,0.90,0.80,0.70,0.50")
    ap.add_argument("--nproc", type=int, default=4)
    ap.add_argument("--fast", action="store_true", help="bez skanu profilu (testy)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    {"sf": cmd_sf, "epsb": cmd_epsb, "iso": cmd_iso, "effmap": cmd_effmap}[args.cmd](args)


if __name__ == "__main__":
    sys.exit(main())
