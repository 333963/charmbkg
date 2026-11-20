#!/usr/bin/env python3
"""DNET — analiza porównawcza klasyfikatorów (XGB H vs sieci) na tych samych zdarzeniach. Projekt B, tylko pojedyncze D.

  analyze.py compute --data DIR --sc DIR[,DIR2] --out dnet_results.json [--nproc 8] [--nboot 200] [--fast]
      (węzeł `ana` na lxplus; DIR = dn_<próbka>.npz z dn_train merge; sc = katalogi z sc/<tag>.npz + meta/<tag>.json; etap 2 dokłada drugi)
  analyze.py report --res dnet_results.json[,dnet2_results.json] --report RAPORT.md
      (lokalnie: sekcja wyników między znacznikami <!-- WYNIKI:BEGIN/END --> + werdykt wg kryteriów zapisanych w RAPORT.md §4)

Metryki (per gatunek, bin pT 3–5/5–8/8–30, WP80/70/50; progi z sygnału MC OOF w oknie sygnału):
  Z_mc = ε_S/√ε_B (tło mb truth ≥ 8 w oknie), Z_tr = ε_S(MC)/√ε_B(dalekie SB danych 428770) [+ 427929], stosunki do XGB H z bootstrapem
  Poissona po zdarzeniach (te same repliki dla wszystkich modeli → błąd sparowany); R_sculpt (liniowy fit ε_B(m) w SB → okno; MC truth ≥ 8
  oraz pasmo danych [okno, dalekie SB)), nachylenie ε_B(m), corr(score, m | tło MC); SF = ε_S(dane)/ε_S(MC) z fitu Δm pass/fail (tag D*:
  428770 all/J12only, 427929 MB427929; tylko D0); C2ST_1D(score sygnału dane vs MC) = 0,5 + 0,5·max_WP|ε_dane − ε_MC|; rozrzut seedów.
Dla modeli z etapu 2 (score = NaN poza podzbiorem z torami) wszystko liczone na wierszach wspólnych, odniesienie przeliczane na tych samych."""
import argparse
import datetime
import glob
import json
import math
import os
import sys
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
for p in (HERE, os.path.join(os.path.dirname(HERE), "Q_lxplus", "dnet1")):
    if os.path.isfile(os.path.join(p, "dn_feat.py")):
        sys.path.insert(0, p); sys.path.insert(0, os.path.join(p, "src"))
        break
BINS = [(3, 5), (5, 8), (8, 30)]
WPS = [0.80, 0.70, 0.50]
SFWPS = [0.80, 0.70, 0.50]
NEED = {"m", "pt", "uid", "truth", "ok", "swap", "m_swap", "dm", "truth_ds", "e_grl"}


def bname(b):
    return f"{b[0]}-{b[1]}"


# ----------------------------------------------------------------------------------------------------------- narzędzia
def closure(m, ok, m0, lo, hi, sbm, region_lo, region_hi):
    """FEAT_charmnet/scripts/sculpt.py: linia ε_B(m) w SB (|Δm| > sbm) → przewidywanie w paśmie |Δm| ∈ [region_lo, region_hi); R = obs/pred."""
    if len(m) < 200:
        return None
    edges = np.linspace(lo, hi, 26)
    n0, _ = np.histogram(m, edges)
    n1, _ = np.histogram(m[ok], edges)
    e = n1 / np.maximum(n0, 1)
    err = np.sqrt(np.maximum(e * (1 - e), 1e-12) / np.maximum(n0, 1))
    c = 0.5 * (edges[1:] + edges[:-1])
    sb = (np.abs(c - m0) > sbm) & (n0 > 20) & (n1 > 0)
    d = np.abs(m - m0)
    reg = (d >= region_lo) & (d < region_hi)
    if sb.sum() < 3 or reg.sum() < 50:
        return None
    A = np.c_[np.ones(sb.sum()), c[sb] - m0]
    W = 1 / err[sb] ** 2
    cov = np.linalg.inv(A.T @ (A * W[:, None]))
    beta = cov @ (A.T @ (W * e[sb]))
    obs = ok[reg].mean()
    x = np.c_[np.ones(reg.sum()), m[reg] - m0].mean(0)
    pred = float(x @ beta)
    if pred <= 0 or obs <= 0:
        return None
    r = obs / pred
    return {"R": float(r), "R_err": float(r * math.hypot(math.sqrt(obs * (1 - obs) / reg.sum()) / obs, math.sqrt(float(x @ cov @ x)) / pred)),
            "slope": float(beta[1] * 0.1 / beta[0])}


def dcor2(x, y, nmax=3000, seed=1):
    n = len(x)
    if n < 10:
        return float("nan")
    if n > nmax:
        i = np.random.default_rng(seed).choice(n, nmax, replace=False)
        x, y = x[i], y[i]
    a = np.abs(x[:, None] - x[None, :]); b = np.abs(y[:, None] - y[None, :])
    A = a - a.mean(0) - a.mean(1)[:, None] + a.mean(); B = b - b.mean(0) - b.mean(1)[:, None] + b.mean()
    return float((A * B).mean() / math.sqrt(max((A * A).mean() * (B * B).mean(), 1e-30)))


def boot_w(uid, nb, seed):
    u, inv = np.unique(uid, return_inverse=True)
    return np.random.default_rng(seed).poisson(1.0, (nb, len(u))).astype(np.uint8)[:, inv]


def _fit_job(job):
    import btrig_common as BC
    key, dp, df_, shape, interval = job
    return key, BC.fit_eps(dp, df_, shape, interval=interval)


def classes(D, run):
    import dn_feat as DF
    n = len(D["p_m"])
    tp = lambda c: (D["p_e_tp_" + c] == 1) if "p_e_tp_" + c in D else np.zeros(n, bool)  # noqa: E731
    mbts = np.zeros(n, bool)
    for c in DF.MBTS_CHAINS:
        mbts |= tp(c)
    grl = (D["p_e_grl"] == 1) if "p_e_grl" in D else np.ones(n, bool)
    if run == 427929:
        return {"MB427929": tp("HLT_mb_mbts_L1MBTS_1") & grl}
    return {"all": np.ones(n, bool), "J12only": tp(DF.J12) & ~mbts}


# ----------------------------------------------------------------------------------------------------------- compute
def compute(a):
    import dn_feat as DF
    import dn_train as DT
    metas, scf = {}, {}
    for d in a.sc.split(","):
        for f in sorted(glob.glob(os.path.join(d, "meta", "*.json"))):
            t = os.path.basename(f)[:-5]
            if os.path.isfile(os.path.join(d, "sc", t + ".npz")):
                metas[t] = json.load(open(f)); scf[t] = os.path.join(d, "sc", t + ".npz")
    if a.keep:
        import re
        import dn_feat as _D
        base = set(_D.REF.values()) | {f"D0.r128.HK.mc.d{_D.LAM[0]}.s1", "D0.xgb.HX.mc.n.s1", "D0.xgb.H.mc.n.s1", f"D0.r128.HKia.mc.d{_D.LAM[0]}.s1"}
        metas = {t: m for t, m in metas.items() if t in base or re.search(a.keep, t)}
    print(f"konfiguracji: {len(metas)}", flush=True)
    res = {"generated": datetime.datetime.now().isoformat(timespec="seconds"), "nboot": a.nboot, "fast": bool(a.fast), "meta": {}, "species": {}}
    for t, m in metas.items():
        res["meta"][t] = {k: m.get(k) for k in ("sp", "model", "feats", "bkg", "dec", "seed", "n_feat", "nparams", "npar_unit", "sec", "threads",
                                                "device", "n_sig", "n_bkg", "onnx", "arch", "env")}
        res["meta"][t]["best"] = [f.get("best_epoch", f.get("best_it")) for f in m.get("fit", [])]
    fitjobs, fitctx = [], {}
    for sp in DF.SPECIES:
        tags = [t for t in metas if metas[t]["sp"] == sp]
        if not tags:
            continue
        ref = DF.REF[sp] if DF.REF[sp] in tags else sorted(tags)[0]
        scfg = DF.cfg["species"][sp]
        m0, win, (lo, hi), sbm, far = scfg["m0"], scfg["signal_window"], scfg["mass_range"], scfg["sideband_min"], DF.FAR[sp]
        keys = NEED | {"e_tp_" + c for c in DF.TRIG_CHAINS}
        S = DT.load_species(a.data, sp, keys)
        sig, bmc, bsb = DT.masks(S, sp)
        fd, mb = S["fd_nd"], S["mb_nd"]
        SC = {t: dict(np.load(scf[t])) for t in tags}
        # zespoły (E): średnia logitów drzewa i sieci nominalnej
        lam0 = DF.LAM[0]
        nn0 = {"D0": f"D0.r128.HK.mc.d{lam0}.s1", "Dplus": f"Dplus.r128.HDLZ.mc.d{lam0}.s1"}[sp]
        lg = lambda s: np.log(np.clip(s, 1e-7, 1 - 1e-7) / (1 - np.clip(s, 1e-7, 1 - 1e-7)))  # noqa: E731
        for x in ([ref] + ([f"{sp}.xgb.HX.mc.n.s1"] if sp == "D0" else [f"{sp}.xgb.HDLZ.mc.n.s1"])):
            if x in SC and nn0 in SC:
                t = f"{sp}.ens.{x.split('.', 1)[1]}+{nn0.split('.', 1)[1]}"
                SC[t] = {k: (1 / (1 + np.exp(-0.5 * (lg(SC[x][k]) + lg(SC[nn0][k]))))).astype(np.float32) for k in SC[x]}
                res["meta"][t] = {"sp": sp, "model": "ens", "feats": "+", "bkg": "+", "dec": "+", "seed": 1, "nparams": None, "sec": 0}
        tags = list(SC)
        # zbiory oceny
        sets = {"sig": ("fd_nd", sig & (np.abs(fd["a_m"] - m0) < win)), "bmc": ("mb_nd", bmc & (np.abs(mb["a_m"] - m0) < win)),
                "bmc_all": ("mb_nd", bmc)}
        for run in (428770, 427929):
            D = S.get(f"data_{run}")
            if D is None:
                continue
            okd = D["a_ok"] & ((np.abs(D["a_m_swap"] - m0) > DF.REFL_VETO) if sp == "D0" else True)
            sets[f"sb{run}"] = (f"data_{run}", okd & (np.abs(D["a_m"] - m0) > far))
            sets[f"dall{run}"] = (f"data_{run}", okd)
        idx = {k: np.flatnonzero(mk) for k, (s, mk) in sets.items()}
        col = lambda k, c: S[sets[k][0]]["a_" + c][idx[k]]  # noqa: E731
        pt = {k: col(k, "pt") for k in sets}
        mm = {k: col(k, "m").astype(np.float64) for k in sets}
        ib = {k: np.digitize(pt[k], [3, 5, 8, 30]) - 1 for k in sets}
        BW = {k: boot_w(col(k, "uid"), a.nboot, 11 + i) for i, k in enumerate(sets) if k in ("sig", "bmc") or k.startswith("sb")}
        def sco(t, k):
            key = sets[k][0] + "__a"
            if key not in SC[t]:                 # próbka bez wyniku w tej konfiguracji (np. etap 2 bez danych) → wiersze nieważne
                return np.full(len(idx[k]), np.nan, np.float32)
            return SC[t][key][idx[k]]

        def one(t, valid):
            """Metryki konfiguracji t na wierszach valid[k] (maski per zbiór) → dict + repliki Z."""
            o, reps, thr = {}, {}, {}
            for bi, b in enumerate(BINS):
                inb = {k: (ib[k] == bi) & valid[k] for k in sets}
                ss = sco(t, "sig")[inb["sig"]]
                if len(ss) < 50:
                    continue
                for wp in WPS + [0.95]:
                    th = float(np.quantile(ss, 1 - wp)); thr[(bi, wp)] = th
                    if wp == 0.95:
                        continue
                    r = {"thr": th, "n_sig": int(len(ss))}
                    es_r = (BW["sig"][:, inb["sig"]] @ (ss >= th).astype(np.float32)) / np.maximum(BW["sig"][:, inb["sig"]].sum(1), 1)
                    r["eps_s"] = float((ss >= th).mean())
                    for k in sets:
                        if not (k == "bmc" or k.startswith("sb")) or inb[k].sum() < 20:
                            continue
                        p = (sco(t, k)[inb[k]] >= th)
                        eb = float(p.mean())
                        eb_r = (BW[k][:, inb[k]] @ p.astype(np.float32)) / np.maximum(BW[k][:, inb[k]].sum(1), 1)
                        nm = {"bmc": "mc"}.get(k, k)
                        r["eps_b_" + nm], r["n_b_" + nm] = eb, int(inb[k].sum())
                        r["Z_" + nm] = r["eps_s"] / math.sqrt(eb) if eb > 0 else float("nan")
                        reps[(bi, wp, nm)] = es_r / np.sqrt(np.maximum(eb_r, 1e-9))
                        r["Z_" + nm + "_err"] = float(np.std(reps[(bi, wp, nm)]))
                    o[f"{bname(b)}|{wp}"] = r
                sb_ = sco(t, "bmc")[inb["bmc"]]
                if len(sb_) > 20:
                    o[f"{bname(b)}|auc"] = DF.auc(np.concatenate([ss, sb_]), np.concatenate([np.ones(len(ss), bool), np.zeros(len(sb_), bool)]))
            # rzeźbienie i korelacja (wszystkie biny pT razem + per bin dla WP70)
            for wp in WPS:
                for k, nm, rl, rh, fm in (("bmc_all", "mc_SR", 0.0, win, sbm), ("dall428770", "data_band", win, far, far),
                                           ("dall427929", "data427_band", win, far, far)):
                    if k not in sets:
                        continue
                    v = valid[k] & (ib[k] >= 0) & (ib[k] < 3)
                    th = np.array([thr.get((bi, wp), np.inf) for bi in range(3)])
                    ok = sco(t, k)[v] >= th[ib[k][v]]
                    c = closure(mm[k][v], ok, m0, lo, hi, fm, rl, rh)
                    if c:
                        o[f"sculpt|{nm}|{wp}"] = c
                    if wp == 0.7 and k == "bmc_all":
                        for bi, b in enumerate(BINS):
                            vb = v & (ib[k] == bi)
                            c = closure(mm[k][vb], sco(t, k)[vb] >= th[bi], m0, lo, hi, fm, rl, rh)
                            if c:
                                o[f"sculpt|{nm}|{wp}|{bname(b)}"] = c
            v = valid["bmc_all"] & (ib["bmc_all"] >= 0) & (ib["bmc_all"] < 3)
            s_, m_ = sco(t, "bmc_all")[v].astype(np.float64), mm["bmc_all"][v]
            o["corr_score_m"] = float(np.corrcoef(s_, m_)[0, 1]) if len(s_) > 10 else None
            o["dcor2_score_m"] = dcor2(s_, m_)
            return o, reps, thr
        allv = {k: np.ones(len(idx[k]), bool) for k in sets}
        R = {"ref": ref, "n": {k: int(len(idx[k])) for k in sets}, "cfg": {}}
        o_ref, reps_ref, thr_ref = one(ref, allv)
        THR = {}
        for t in tags:
            valid = {k: np.isfinite(sco(t, k)) for k in sets}
            part = any(not v.all() for v in valid.values())
            o, reps, thr = one(t, valid) if t != ref else (dict(o_ref), reps_ref, thr_ref)
            orf, rrf = (one(ref, valid)[:2] if part else (o_ref, reps_ref))
            for key, r in o.items():
                if "|" in key and key.split("|")[1] in ("0.8", "0.7", "0.5") and key in orf:
                    bi = [bname(b) for b in BINS].index(key.split("|")[0]); wp = float(key.split("|")[1])
                    for nm in [k[2:] for k in r if k.startswith("Z_") and not k.endswith("_err")]:
                        if (bi, wp, nm) in reps and (bi, wp, nm) in rrf and orf[key].get("Z_" + nm):
                            r["G_" + nm] = r["Z_" + nm] / orf[key]["Z_" + nm]
                            r["G_" + nm + "_err"] = float(np.std(reps[(bi, wp, nm)] / rrf[(bi, wp, nm)]))
            o["partial_rows"] = bool(part)
            o["frac_rows"] = {k: float(valid[k].mean()) if len(valid[k]) else None for k in ("sig", "bmc")}
            if part:
                o["ref_same_rows"] = {k: orf[k] for k in orf if k.startswith("sculpt|mc_SR|0.7") or k.endswith("|0.7")}
            R["cfg"][t] = o
            THR[t] = thr
            print(f"[{sp}] {t}: " + " ".join(f"{bname(b)} G_mc {o.get(bname(b) + '|0.7', {}).get('G_mc', float('nan')):.3f}" for b in BINS)
                  + f"  R {o.get('sculpt|mc_SR|0.7', {}).get('R', float('nan')):.3f}", flush=True)
        # SF z tagu D* (D0): fit Δm pass/fail
        if sp == "D0" and "p_m" in fd and not a.no_sf:
            import btrig_common as BC
            from charmtag import dstar_refit as DR
            fp = np.isin(fd["p_truth"], DF.SIGNAL_CATS) & np.isin(fd["p_truth_ds"], DF.SIGNAL_CATS) & (np.abs(fd["p_m"] - m0) < win)
            pb = lambda T: np.digitize(T["p_pt"], [3, 5, 8, 30]) - 1  # noqa: E731
            shapes = {bi: DR.shape_from_mc(fd["p_dm"][fp & (pb(fd) == bi)] * 1e3) for bi in range(3) if (fp & (pb(fd) == bi)).sum() > 200}
            sft = [t for t in tags if a.sf_all or res["meta"][t].get("seed", 1) == 1]
            R["sf_eps_mc"], R["sf_ncls"] = {}, {}
            for t in sft:
                sfd = SC[t]["fd_nd__p"]
                for bi in shapes:
                    for wp in SFWPS:
                        th = THR[t].get((bi, wp)) if THR[t] else None
                        if th is None:
                            continue
                        sel = fp & (pb(fd) == bi) & np.isfinite(sfd)
                        k, n = int((sfd[sel] >= th).sum()), int(sel.sum())
                        R["sf_eps_mc"][f"{t}|{bname(BINS[bi])}|{wp}"] = [k / max(n, 1), math.sqrt(max(k / max(n, 1) * (1 - k / max(n, 1)), 1e-12) / max(n, 1)), n]
                        for run in (428770, 427929):
                            D = S.get(f"data_{run}")
                            if D is None or "p_m" not in D:
                                continue
                            sd = SC[t][f"data_{run}__p"]
                            for c, mk in classes(D, run).items():
                                md = mk & (np.abs(D["p_m"] - m0) < win) & (pb(D) == bi) & np.isfinite(sd)
                                R["sf_ncls"][f"{run}|{c}|{bname(BINS[bi])}"] = int(md.sum())
                                fitjobs.append(((sp, t, run, c, bname(BINS[bi]), wp), D["p_dm"][md & (sd >= th)], D["p_dm"][md & (sd < th)], shapes[bi],
                                                not a.fast))
        res["species"][sp] = R
        del S, SC
    if fitjobs:
        print(f"fitów Δm pass/fail: {len(fitjobs)} (nproc {a.nproc})", flush=True)
        import btrig_common as BC
        with Pool(a.nproc) as P:
            fits = dict(P.map(_fit_job, fitjobs, chunksize=1))
        for (sp, t, run, c, b, wp), f in fits.items():
            R = res["species"][sp]
            em, eme, _ = R["sf_eps_mc"][f"{t}|{b}|{wp}"]
            err = BC.eps_sym_err(f)
            sf = f["eps"] / em if em and np.isfinite(f["eps"]) else float("nan")
            sfe = float(np.hypot(err / em, sf * eme / em)) if em and np.isfinite(sf) else float("nan")
            R.setdefault("sf", {})[f"{t}|{run}|{c}|{b}|{wp}"] = {"eps_data": f["eps"], "eps_data_err": err, "eps_mc": em, "SF": sf, "SF_err": sfe,
                                                              "N_sig": f["N"], "determined": bool(f["determined"]), "valid": bool(f["valid"])}
    json.dump(res, open(a.out, "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    print("zapisano", a.out, flush=True)


# ----------------------------------------------------------------------------------------------------------- report
def pm(x, e=None, nd=3):
    if x is None or not np.isfinite(x):
        return "–"
    return f"{x:.{nd}f}" + (f" ± {e:.{nd}f}" if e is not None and np.isfinite(e) else "")


def fam(t):
    return t.rsplit(".s", 1)[0]


def report(a):
    res = {"meta": {}, "species": {}}
    for f in a.res.split(","):
        r = json.load(open(f))
        res["meta"].update(r["meta"])
        for sp, R in r["species"].items():
            if sp not in res["species"]:
                res["species"][sp] = R
            else:
                for k in ("cfg", "sf", "sf_eps_mc"):
                    res["species"][sp].setdefault(k, {}).update(R.get(k, {}))
        res["generated"] = r["generated"]; res["fast"] = r.get("fast")
    L = [f"Wygenerowano {res['generated']} przez `DNET/analyze.py report` z `{a.res}`." + (" **Tryb --fast (błędy HESSE).**" if res.get("fast") else ""), ""]
    verdicts = []
    for sp, R in res["species"].items():
        ref, C, M = R["ref"], R["cfg"], res["meta"]
        L += [f"### {sp} — odniesienie `{ref}`; wierszy: " + ", ".join(f"{k} {v}" for k, v in R["n"].items()), ""]
        seeds = {}
        for t in C:
            seeds.setdefault(fam(t), []).append(t)
        L += ["#### Zproxy względem XGB H (G = Z/Z_H, bootstrap sparowany), WP70; R_sculpt i corr na tle MC", "",
              "| konfiguracja | cech | G_mc 3–5 | G_mc 5–8 | G_mc 8–30 | G_tr(SB 428770) 3–5 | 5–8 | 8–30 | R_sculpt MC (WP70) | R dane (pasmo) | nachylenie /100 MeV | corr(s, m) | par. | czas [s] |",
              "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
        for t in sorted(C, key=lambda x: (M.get(x, {}).get("model", "") != "xgb", x)):
            o = C[t]
            g = lambda nm: [pm(o.get(f"{bname(b)}|0.7", {}).get("G_" + nm), o.get(f"{bname(b)}|0.7", {}).get("G_" + nm + "_err")) for b in BINS]  # noqa: E731
            s1, s2 = o.get("sculpt|mc_SR|0.7", {}), o.get("sculpt|data_band|0.7", {})
            L.append(f"| `{t}`{' (podzbiór)' if o.get('partial_rows') else ''} | {M.get(t, {}).get('n_feat', '–')} | " + " | ".join(g("mc") + g("sb428770")) +
                     f" | {pm(s1.get('R'), s1.get('R_err'))} | {pm(s2.get('R'), s2.get('R_err'))} | {pm(s1.get('slope'), nd=2)} | {pm(o.get('corr_score_m'))} | "
                     f"{M.get(t, {}).get('nparams') or '–'} | {M.get(t, {}).get('sec', '–')} |")
        L += ["", "#### WP80 / WP50 (G_mc; G_tr)", "", "| konfiguracja | WP80 3–5 | 5–8 | 8–30 | WP50 3–5 | 5–8 | 8–30 |", "|---|---|---|---|---|---|---|"]
        for t in sorted(C):
            L.append(f"| `{t}` | " + " | ".join(f"{pm(C[t].get(f'{bname(b)}|{wp}', {}).get('G_mc'), nd=2)}; {pm(C[t].get(f'{bname(b)}|{wp}', {}).get('G_sb428770'), nd=2)}"
                                                 for wp in (0.8, 0.5) for b in BINS) + " |")
        L += ["", "#### Rozrzut seedów (odchylenie std G_mc WP70 między seedami)", "", "| rodzina | seedy | 3–5 | 5–8 | 8–30 |", "|---|---|---|---|---|"]
        sstd = {}
        for f_, ts in seeds.items():
            if len(ts) >= 2:
                v = [[C[t].get(f"{bname(b)}|0.7", {}).get("G_mc", np.nan) for t in ts] for b in BINS]
                sstd[f_] = [float(np.nanstd(x, ddof=1)) for x in v]
                L.append(f"| `{f_}` | {len(ts)} | " + " | ".join(f"{np.nanmean(x):.3f} ± {np.nanstd(x, ddof=1):.3f}" for x in v) + " |")
        if R.get("sf"):
            L += ["", "#### SF = ε_S(dane)/ε_S(MC) z tagu D* (fit Δm pass/fail) i C2ST_1D sygnału", "",
                  "| konfiguracja | run:klasa | bin | SF WP80 | SF WP70 | SF WP50 | C2ST_1D |", "|---|---|---|---|---|---|---|"]
            keys = sorted({k.rsplit("|", 1)[0] for k in R["sf"]})
            for k in keys:
                t, run, c, b = k.split("|")
                x = {wp: R["sf"].get(f"{k}|{wp}") for wp in SFWPS}
                de = [abs(v["eps_data"] - v["eps_mc"]) for v in x.values() if v and v["determined"] and np.isfinite(v["eps_data"])]   # tylko fity określone
                L.append(f"| `{t}` | {run}:{c} | {b} | " + " | ".join(pm(x[wp]["SF"], x[wp]["SF_err"], 2) + ("" if x[wp]["determined"] else " (n)") if x[wp] else "–"
                                                                    for wp in (0.8, 0.7, 0.5)) + f" | {pm(0.5 + 0.5 * max(de) if de else None)} |")
        # werdykt wg §4
        L += ["", "#### Werdykt wg kryteriów §4", ""]
        for t in sorted(C):
            if t == ref or M.get(t, {}).get("model") == "xgb" or M.get(t, {}).get("seed", 1) != 1:
                continue
            o, orf = C[t], C[ref]
            nbin = 0
            for b in BINS:
                r = o.get(f"{bname(b)}|0.7", {})
                gm, gt = r.get("G_mc"), r.get("G_sb428770")
                if gm and gt and gm >= 1.10 and gt >= 1.10 and gm - 1 > 2 * r.get("G_mc_err", 0):
                    nbin += 1
            k1 = nbin >= 2
            s, sr = o.get("sculpt|mc_SR|0.7"), orf.get("sculpt|mc_SR|0.7")
            k2 = bool(s and sr and abs(s["R"] - 1) <= abs(sr["R"] - 1) + 2 * math.hypot(s["R_err"], sr["R_err"]))
            k3 = None
            if R.get("sf"):
                ncmp, bad = 0, 0
                for c_ in ("428770|all", "428770|J12only", "427929|MB427929"):
                    for b in ("5-8", "8-30"):
                        x, y = R["sf"].get(f"{t}|{c_}|{b}|0.7"), R["sf"].get(f"{ref}|{c_}|{b}|0.7")
                        if x and y and x["determined"] and y["determined"] and np.isfinite(x["SF"]) and np.isfinite(y["SF"]):
                            ncmp += 1
                            bad += abs(x["SF"] - 1) > abs(y["SF"] - 1) + math.hypot(x["SF_err"], y["SF_err"])
                k3 = None if ncmp < 6 else (bad == 0)          # wszystkie 6 pomiarów (3 klasy × 2 biny) muszą być określone
                if bad:
                    k3 = False
            sd = sstd.get(fam(t))
            gains = [o.get(f"{bname(b)}|0.7", {}).get("G_mc", 1) - 1 for b in BINS]
            k4 = None if sd is None else bool(all(s_ < max(g_, 1e-9) / 3 for s_, g_ in zip(sd[:2], gains[:2])))
            ok = k1 and k2 and (k3 is not False) and (k4 is not False)
            verdicts.append((t, ok))
            L.append(f"- `{t}`: K1 zysk ≥ 10 % w ≥ 2 binach (MC i transfer): **{'tak' if k1 else 'nie'}** ({nbin}/3); K2 rzeźbienie nie gorsze niż H: "
                     f"**{'tak' if k2 else 'nie'}**; K3 SF: **{ {True: 'tak', False: 'nie', None: 'brak pomiaru'}[k3]}**; K4 seedy: "
                     f"**{ {True: 'tak', False: 'nie', None: 'brak (1 seed)'}[k4]}** → **{'PRZYJĘTA' if ok and k3 and k4 else ('warunkowo' if ok else 'odrzucona')}**")
        L.append("")
    txt = "\n".join(L)
    if a.report:
        old = open(a.report).read() if os.path.isfile(a.report) else ""
        b, e = "<!-- WYNIKI:BEGIN -->", "<!-- WYNIKI:END -->"
        if b in old and e in old:
            new = old[:old.index(b) + len(b)] + "\n" + txt + "\n" + old[old.index(e):]
        else:
            new = old + f"\n{b}\n{txt}\n{e}\n"
        open(a.report, "w").write(new)
        print("zapisano", a.report)
    else:
        print(txt)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["compute", "report"])
    ap.add_argument("--data"); ap.add_argument("--sc"); ap.add_argument("--out")
    ap.add_argument("--nproc", type=int, default=4); ap.add_argument("--nboot", type=int, default=200)
    ap.add_argument("--fast", action="store_true"); ap.add_argument("--sf-all", action="store_true"); ap.add_argument("--no-sf", action="store_true")
    ap.add_argument("--keep", default="", help="wyrażenie regularne: tylko te konfiguracje (+ odniesienia) — etap 2")
    ap.add_argument("--res"); ap.add_argument("--report")
    a = ap.parse_args(argv)
    return compute(a) if a.cmd == "compute" else report(a)


if __name__ == "__main__":
    sys.exit(main())
