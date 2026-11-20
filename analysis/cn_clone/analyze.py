#!/usr/bin/env python3
"""CN_clone — analiza lokalna klonu CharmNet D0 (wariant „iso” z ptcone40+ptcone100 vs „noiso”) → RAPORT.md + out/cn_results.json.
Wejście: cn_{data_428770,data_427929,fd_nd,mb_nd}.npz + cn_summary.json (etap ORCH cnclone → CN_clone/runs/lxplus).
  .venv_h5/bin/python charm_work/CN_clone/analyze.py [--in runs/lxplus] [--out out] [--report RAPORT.md] [--nproc 8] [--fast]
Narzędzia fitu: B_trig/btrig_common.py (fit Δm pass/fail = charmtag.dstar_refit z AA_merge; ε_B = ML exp(pol2) na SB).
Projekt B (charmbkg). TO KLON, nie CharmNet. Kryteria rozstrzygnięcia (KRYTERIA poniżej) zapisane 17.11.2026 PRZED danymi."""
import argparse
import datetime
import json
import math
import os
import sys
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CW = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(CW, "B_trig"))
sys.path.insert(0, os.path.join(CW, "Q_lxplus", "cnclone"))
import btrig_common as C  # noqa: E402
import cn_feat as CF  # noqa: E402

VARS = ["iso", "noiso"]
BINS = C.BINS                                   # 3-5 (informacyjnie), 5-8, 8-30
MAIN_BINS = ["5-8", "8-30"]
AFP_CHAINS = ["HLT_noalg_L1AFP_A_OR_C", "HLT_mb_sptrk_pt4_L1AFP_A_OR_C"]
CLS = {428770: ["all", "J12only", "MBTSbased", "j20", "AFPonly", "other", "AFPtag"], 427929: ["MB427929"]}
REG = {"sr_half": C.NOM_WIN, "sb_min": 0.075, "lo": 1.75, "hi": 2.00}

KRYTERIA = """\
Zapisane 17.11.2026, przed pobraniem danych pełnej statystyki (test mechaniki tylko na 1 pliku h5 na próbkę, bez mocy statystycznej).
Punkt odniesienia: WP80 (pomocniczo WP50, WP65 ≈ punkt pracy D0 v37 ε_S = 0,64), biny 5–8 i 8–30 GeV; bin 3–5 tylko informacyjnie.
Klasy rozstrzygające w 428770: `all` i `J12only`; kontrola: `MB427929` (run 427929, MBTS).

- **K0 (bramka poprawności klonu, MC).** AUC OOF w oknie piku: iso > noiso; zysk Z = ε_S/√ε_B z izolacji na MC w WP80 zgodny co do znaku
  i rzędu z FEAT_charmnet (+14–19 % dla XGB); |corr(score, m | tło komb.)| < 0,1; AUC(pełny) − AUC(OOF) < 0,01. Niespełnione → wynik
  tylko jakościowy.
- **K1 (koszt: błąd ε_S z MC).** R_SF = SF_iso/SF_noiso. „Izolacja psuje ε_S z MC”, jeśli w 428770 (`all` lub `J12only`) w co najmniej
  jednym z binów 5–8 / 8–30 R_SF < 0,90 z istotnością ≥ 2σ, a w drugim binie R_SF ≤ 1 (ten sam znak). „Brak kosztu”, jeśli w obu binach
  |R_SF − 1| < 0,05 i < 1σ. Inaczej „nierozstrzygnięte”. Błąd R_SF liczony jak dla pomiarów niezależnych (zachowawczo — te same zdarzenia).
  W `MB427929` oczekiwane R_SF ≈ 1 (kontrola; odchylenie ≥ 2σ = sygnał problemu metody albo efekt nie-trygerowy).
- **K2 (czy przewaga przeżywa w danych).** G = Z_iso/Z_noiso, Z = ε_S/√ε_B; G_MC z MC (sygnał fd, tło mb), G_dane z ε_S(dane) z fitu Δm
  i ε_B(dane) z pasm bocznych. „Przewaga przeżywa”, jeśli w 428770 `all` w obu binach G_dane > 1 z ≥ 2σ i G_dane − 1 ≥ 0,5·(G_MC − 1).
  „Zanika”, jeśli G_dane ≤ 1 albo G_dane − 1 < 1σ w obu binach. Pomiędzy: „częściowo”.
- **K3 (rekomendacja).** K1 = „psuje” → sieć z ptcone wymaga SF z danych per klasa trygera (albo wariant bez ptcone, jeśli K2 ≠ „przeżywa”).
  K1 = „brak kosztu” i K2 = „przeżywa” → izolacja w sieci do utrzymania. Reszta → potrzebna większa statystyka / score prawdziwego CharmNet.
- Zakres wniosku: klon (inne wagi niż v60), tylko D0, selekcja par nasza, brak wymogu tagu AFP (klasa `AFPtag` tylko informacyjnie)."""


# ------------------------------------------------------------------------------------------------------------------ wejście
def load(d, s):
    p = os.path.join(d, f"cn_{s}.npz")
    if not os.path.isfile(p):
        return None
    z = np.load(p, allow_pickle=False)
    return {k: z[k] for k in z.files if k != "meta"}


def classes(D, pre, run):
    n = len(D[pre + "m"])

    def tp(c):
        k = f"{pre}e_tp_{c}"
        return (D[k] == 1) if k in D else np.zeros(n, bool)
    mbts = np.zeros(n, bool)
    for c in C.MBTS_CHAINS:
        mbts |= tp(c)
    j12 = tp(C.J12)
    afp = np.zeros(n, bool)
    for c in AFP_CHAINS:
        afp |= tp(c)
    grl = (D[pre + "e_grl"] == 1) if pre + "e_grl" in D else np.ones(n, bool)
    if run == 427929:
        return {"MB427929": tp("HLT_mb_mbts_L1MBTS_1") & grl}
    tag = ((D[pre + "e_afpA"] == 1) | (D[pre + "e_afpC"] == 1)) if pre + "e_afpA" in D else np.zeros(n, bool)
    return {"all": np.ones(n, bool), "J12only": j12 & ~mbts, "MBTSbased": mbts, "j20": tp("HLT_j20_L1MBTS_2"),
            "AFPonly": afp & ~j12 & ~mbts, "other": ~j12 & ~mbts & ~afp, "AFPtag": tag}


def inbin(pt, b):
    lo, hi = BINS[b]
    return (pt >= lo) & (pt < hi)


def binom(k, n):
    if n <= 0:
        return float("nan"), float("nan")
    e = k / n
    return e, math.sqrt(max(e * (1 - e), 1e-12) / n)


def _fit_job(job):
    key, dp, df_, shape, interval = job
    return key, C.fit_eps(dp, df_, shape, interval=interval)


def rel(a, ea, b, eb):
    """a/b z błędem (niezależne)."""
    if not (np.isfinite(a) and np.isfinite(b)) or b == 0 or a == 0:
        return float("nan"), float("nan")
    r = a / b
    return r, abs(r) * math.hypot(ea / a if np.isfinite(ea) else 0, eb / b if np.isfinite(eb) else 0)


# ------------------------------------------------------------------------------------------------------------------ analiza
def analyse(a):
    C.add_paths()
    from charmtag import dstar_refit as DR
    summ = json.load(open(os.path.join(a.inp, "cn_summary.json")))
    FD, MB = load(a.inp, "fd_nd"), load(a.inp, "mb_nd")
    if FD is None or MB is None:
        sys.exit("brak cn_fd_nd.npz / cn_mb_nd.npz")
    DATA = {r: load(a.inp, f"data_{r}") for r in (428770, 427929)}
    DATA = {r: d for r, d in DATA.items() if d is not None}
    wps = a.wps.split(",")
    res = {"inputs": a.inp, "wps": wps, "summary": {k: summ[k] for k in ("backend", "substitute", "folds", "n_train", "n_train_sig", "n_train_bkg",
                                                                       "samples", "versions", "time_s", "hp")}}
    # --- MC: progi (sygnał inkluzywny fd, okno nominalne), ROC/AUC, ε_B^MC, Z
    fs = np.isin(FD["a_truth"], CF.SIGNAL_CATS) & (np.abs(FD["a_m"] - C.M0_D0) < C.NOM_WIN)
    mbk = ~np.isin(MB["a_truth"], CF.SIGNAL_CATS) & (np.abs(MB["a_m"] - C.M0_D0) < C.NOM_WIN)
    thr, mc = {}, {}
    for v in VARS:
        for b in BINS:
            s_sig = FD["a_s_" + v][fs & inbin(FD["a_pt"], b)]
            s_bkg = MB["a_s_" + v][mbk & inbin(MB["a_pt"], b)]
            sf_sig = FD["a_sfull_" + v][fs & inbin(FD["a_pt"], b)]
            sf_bkg = MB["a_sfull_" + v][mbk & inbin(MB["a_pt"], b)]
            auc = CF.auc(np.concatenate([s_sig, s_bkg]), np.concatenate([np.ones(len(s_sig), bool), np.zeros(len(s_bkg), bool)]))
            aucf = CF.auc(np.concatenate([sf_sig, sf_bkg]), np.concatenate([np.ones(len(sf_sig), bool), np.zeros(len(sf_bkg), bool)]))
            mc[(v, b)] = {"n_sig": int(len(s_sig)), "n_bkg": int(len(s_bkg)), "auc_oof": auc, "auc_full": aucf, "wp": {}}
            for wp in wps:
                t = float(np.quantile(s_sig, 1 - float(wp))) if len(s_sig) else float("nan")
                thr[(v, b, wp)] = t
                eb, eeb = binom(int((s_bkg >= t).sum()), len(s_bkg))
                es, ees = binom(int((s_sig >= t).sum()), len(s_sig))
                # efekt braku OOF: ten sam próg na score modelu pełnego (in-sample)
                esf = float((sf_sig >= t).mean()) if len(sf_sig) else float("nan")
                mc[(v, b)]["wp"][wp] = {"thr": t, "eps_s": es, "eps_s_err": ees, "eps_b": eb, "eps_b_err": eeb, "eps_s_fullmodel": esf,
                                        "Z": es / math.sqrt(eb) if eb and eb > 0 else float("nan")}
    # --- SF: fit Δm pass/fail na parach
    fp = np.isin(FD["p_truth"], CF.SIGNAL_CATS) & np.isin(FD["p_truth_ds"], CF.SIGNAL_CATS) & (np.abs(FD["p_m"] - C.M0_D0) < C.NOM_WIN)
    mp = np.isin(MB["p_truth"], CF.SIGNAL_CATS) & np.isin(MB["p_truth_ds"], CF.SIGNAL_CATS) & (np.abs(MB["p_m"] - C.M0_D0) < C.NOM_WIN)
    shapes = {b: DR.shape_from_mc(FD["p_dm"][fp & inbin(FD["p_pt"], b)] * 1e3) for b in BINS}
    eps_mc = {}
    for v in VARS:
        for b in BINS:
            for wp in wps:
                t = thr[(v, b, wp)]
                sel = fp & inbin(FD["p_pt"], b)
                e, ee = binom(int((FD["p_s_" + v][sel] >= t).sum()), int(sel.sum()))
                selm = mp & inbin(MB["p_pt"], b)
                em, eem = binom(int((MB["p_s_" + v][selm] >= t).sum()), int(selm.sum()))
                eps_mc[(v, b, wp)] = {"fd": e, "fd_err": ee, "mb": em, "mb_err": eem, "n_fd": int(sel.sum()), "n_mb": int(selm.sum())}
    jobs, ncls = [], {}
    for run, D in DATA.items():
        cl = classes(D, "p_", run)
        nom = np.abs(D["p_m"] - C.M0_D0) < C.NOM_WIN
        for c in CLS[run]:
            if a.classes and c not in a.classes.split(","):
                continue
            ncls[f"{run}:{c}"] = int(cl[c].sum())
            for b in BINS:
                md = cl[c] & nom & inbin(D["p_pt"], b)
                for v in VARS:
                    for wp in wps:
                        pas = D["p_s_" + v] >= thr[(v, b, wp)]
                        jobs.append(((run, c, b, v, wp), D["p_dm"][md & pas], D["p_dm"][md & ~pas], shapes[b], not a.fast))
    print(f"fitów Δm pass/fail: {len(jobs)} (nproc {a.nproc}, {'HESSE' if a.fast else 'profil'})", flush=True)
    with Pool(a.nproc) as P:
        fits = dict(P.map(_fit_job, jobs, chunksize=1))
    sf = {}
    for (run, c, b, v, wp), f in fits.items():
        em = eps_mc[(v, b, wp)]
        err = C.eps_sym_err(f)
        s = f["eps"] / em["fd"] if em["fd"] and np.isfinite(f["eps"]) else float("nan")
        se = float(np.hypot(err / em["fd"], s * em["fd_err"] / em["fd"])) if em["fd"] and np.isfinite(s) else float("nan")
        sf[(run, c, b, v, wp)] = {"eps_data": f["eps"], "eps_data_err": err, "determined": bool(f["determined"]), "valid": bool(f["valid"]),
                                  "N_sig": f["N"], "N_sig_err": f["N_err"], "n_pass": f["n_pass"], "n_fail": f["n_fail"],
                                  "eps_mc": em["fd"], "eps_mc_mb": em["mb"], "SF": s, "SF_err": se}
    # --- ε_B w danych (SB → SR, ML exp(pol2); kształt z klasy odniesienia runu)
    epsb = {}
    for run, D in DATA.items():
        cl = classes(D, "a_", run)
        ref = CLS[run][0]
        for b in BINS:
            for v in VARS:
                for wp in wps:
                    t = thr[(v, b, wp)]
                    s0 = cl[ref] & inbin(D["a_pt"], b)
                    try:
                        shp = C.pol2_shapes(D["a_m"][s0], D["a_s_" + v][s0] >= t, C.M0_D0, REG["sr_half"], REG["sb_min"], REG["lo"], REG["hi"])
                    except Exception:  # noqa: BLE001
                        shp = None
                    for c in CLS[run]:
                        sel = cl[c] & inbin(D["a_pt"], b)
                        try:
                            r = C.eps_b_pol2_ml(D["a_m"][sel], D["a_s_" + v][sel] >= t, C.M0_D0, REG["sr_half"], REG["sb_min"], REG["lo"], REG["hi"],
                                                shapes=shp if shp and shp.get("all") and shp.get("pass") else None)
                        except Exception:  # noqa: BLE001
                            r = {"eps": float("nan"), "err": float("nan"), "n_sb": 0, "n_sb_pass": 0}
                        epsb[(run, c, b, v, wp)] = r
    # --- zysk z izolacji
    gain = {}
    for b in BINS:
        for wp in wps:
            zi, zn = mc[("iso", b)]["wp"][wp], mc[("noiso", b)]["wp"][wp]
            g = zi["Z"] / zn["Z"] if zn["Z"] else float("nan")
            ge = g * 0.5 * math.hypot(zi["eps_b_err"] / zi["eps_b"] if zi["eps_b"] else 0, zn["eps_b_err"] / zn["eps_b"] if zn["eps_b"] else 0) \
                if np.isfinite(g) else float("nan")
            gain[("MC", "mc", b, wp)] = {"G": g, "G_err": ge, "rej_ratio": zn["eps_b"] / zi["eps_b"] if zi["eps_b"] else float("nan")}
    for (run, c, b, v, wp) in list(sf):
        if v != "iso":
            continue
        out = {}
        Z = {}
        for vv in VARS:
            s_, e_ = sf[(run, c, b, vv, wp)], epsb.get((run, c, b, vv, wp), {})
            es, ees, eb, eeb = s_["eps_data"], s_["eps_data_err"], e_.get("eps", float("nan")), e_.get("err", float("nan"))
            ok = np.isfinite(es) and np.isfinite(eb) and eb > 0 and es > 0
            z = es / math.sqrt(eb) if ok else float("nan")
            ze = z * math.hypot(ees / es if np.isfinite(ees) else 0, 0.5 * eeb / eb if np.isfinite(eeb) else 0) if ok else float("nan")
            Z[vv] = (z, ze)
            out[f"Z_{vv}"], out[f"Z_{vv}_err"] = z, ze
        out["G"], out["G_err"] = rel(Z["iso"][0], Z["iso"][1], Z["noiso"][0], Z["noiso"][1])
        out["R_SF"], out["R_SF_err"] = rel(sf[(run, c, b, "iso", wp)]["SF"], sf[(run, c, b, "iso", wp)]["SF_err"],
                                           sf[(run, c, b, "noiso", wp)]["SF"], sf[(run, c, b, "noiso", wp)]["SF_err"])
        # odrzucenie tła przy stałej ε_S(dane): krzywa noiso (ε_S, ln ε_B) po WP, interpolowana w ε_S(dane) klonu iso
        xs = [(sf[(run, c, b, "noiso", w)]["eps_data"], epsb.get((run, c, b, "noiso", w), {}).get("eps", float("nan"))) for w in wps]
        xs = sorted([(x, y) for x, y in xs if np.isfinite(x) and np.isfinite(y) and y > 0])
        es_i, eb_i = sf[(run, c, b, "iso", wp)]["eps_data"], epsb.get((run, c, b, "iso", wp), {}).get("eps", float("nan"))
        rr = float("nan")
        if len(xs) >= 2 and np.isfinite(es_i) and np.isfinite(eb_i) and eb_i > 0 and xs[0][0] <= es_i <= xs[-1][0]:
            rr = float(math.exp(np.interp(es_i, [x for x, _ in xs], [math.log(y) for _, y in xs])) / eb_i)
        out["rej_ratio_fixed_epsS"] = rr
        gain[(run, c, b, wp)] = out
    res.update({"thr": thr, "mc": mc, "eps_mc": eps_mc, "sf": sf, "epsb": epsb, "gain": gain, "n_class_pairs": ncls, "diag": summ["diag"],
                "feat_map": summ["feat_map"], "missing": summ["missing"], "fast": bool(a.fast)})
    res["verdict"] = verdict(res, a.ref_wp if a.ref_wp in wps else wps[min(1, len(wps) - 1)])
    return res


def verdict(res, wp):
    sf, gain, mc, diag = res["sf"], res["gain"], res["mc"], res["diag"]
    V = {"ref_wp": wp}
    # K0
    k0 = []
    for b in MAIN_BINS:
        k0.append(mc[("iso", b)]["auc_oof"] > mc[("noiso", b)]["auc_oof"])
        g = gain[("MC", "mc", b, wp)]["G"]
        k0.append(np.isfinite(g) and g > 1.0)
        k0 += [abs(mc[(v, b)]["auc_full"] - mc[(v, b)]["auc_oof"]) < 0.01 for v in VARS]
    k0 += [abs(diag[v].get("corr_score_mass_comb") or 0) < 0.1 for v in VARS]
    V["K0"] = "spełnione" if all(k0) else "NIESPEŁNIONE"
    # K1
    k1 = {}
    for c in ("all", "J12only", "MB427929"):
        run = 427929 if c == "MB427929" else 428770
        rs = [gain.get((run, c, b, wp), {}) for b in MAIN_BINS]
        if not all("R_SF" in r and np.isfinite(r["R_SF"]) and np.isfinite(r["R_SF_err"]) and r["R_SF_err"] > 0 for r in rs):
            k1[c] = "brak pomiaru"
            continue
        r, e = [x["R_SF"] for x in rs], [x["R_SF_err"] for x in rs]
        bad = any(r[i] < 0.90 and (1 - r[i]) / e[i] >= 2 and r[1 - i] <= 1 for i in (0, 1))
        none = all(abs(r[i] - 1) < 0.05 and abs(r[i] - 1) / e[i] < 1 for i in (0, 1))
        k1[c] = "psuje" if bad else "brak kosztu" if none else "nierozstrzygnięte"
    V["K1"] = k1
    V["K1_428770"] = "psuje" if "psuje" in (k1.get("all"), k1.get("J12only")) else \
        "brak kosztu" if k1.get("all") == "brak kosztu" and k1.get("J12only") == "brak kosztu" else \
        "brak pomiaru" if k1.get("all") == "brak pomiaru" else "nierozstrzygnięte"
    # K2
    gs = [(gain.get((428770, "all", b, wp), {}), gain[("MC", "mc", b, wp)]) for b in MAIN_BINS]
    if not all("G" in g and np.isfinite(g["G"]) and np.isfinite(g["G_err"]) and g["G_err"] > 0 for g, _ in gs):
        V["K2"] = "brak pomiaru"
    else:
        surv = all(g["G"] > 1 and (g["G"] - 1) / g["G_err"] >= 2 and (g["G"] - 1) >= 0.5 * (m["G"] - 1) for g, m in gs)
        gone = all(g["G"] <= 1 or (g["G"] - 1) / g["G_err"] < 1 for g, m in gs)
        V["K2"] = "przeżywa" if surv else "zanika" if gone else "częściowo"
    if V["K1_428770"] == "psuje":
        V["K3"] = "sieć z ptcone wymaga SF z danych per klasa trygera" + ("" if V["K2"] == "przeżywa" else " albo wariantu bez ptcone")
    elif V["K1_428770"] == "brak kosztu" and V["K2"] == "przeżywa":
        V["K3"] = "izolacja w sieci do utrzymania"
    else:
        V["K3"] = "nierozstrzygnięte — potrzebna większa statystyka albo score prawdziwego CharmNet"
    return V


# ------------------------------------------------------------------------------------------------------------------ raport
def pm(x, e, nd=2):
    if x is None or not np.isfinite(x):
        return "–"
    return f"{x:.{nd}f} ± {e:.{nd}f}" if e is not None and np.isfinite(e) else f"{x:.{nd}f}"


def report(res, a):
    wps, sf, epsb, gain, mc = res["wps"], res["sf"], res["epsb"], res["gain"], res["mc"]
    s = res["summary"]
    V = res["verdict"]
    L = []
    w = L.append
    w("# CN_clone — klon CharmNet D⁰ z izolacją (ptcone40 + ptcone100) i bez: SF z tagu D*, koszt i zysk z izolacji")
    w("")
    w(f"Wygenerowano {datetime.datetime.now():%Y-%m-%d %H:%M} przez `CN_clone/analyze.py`. Wejście: `{res['inputs']}`. Projekt B, pojedyncze D (nie ślepe)."
      + (" **Błędy fitów: HESSE (--fast).**" if res["fast"] else ""))
    if a.note:
        w("")
        w(f"**{a.note}**")
    w("")
    w("## Zastrzeżenia")
    w("")
    w("- **To klon, nie CharmNet.** Te same wejścia, architektura i funkcja straty co w `charmpy`, ale własne wagi; zamrożonego modelu D⁰ v60 nie mamy.")
    w("- **Tylko D⁰.** D⁺ i D*⁺ poza zakresem.")
    w("- **Selekcja nasza, nie z charmpy:** pary D⁰←D* dobierane kinematycznie (jak B_trig), okno |m − m₀| < 55 MeV, L_xy ≥ "
      f"{next(iter(s['samples'].values())).get('lxy_min')} mm; preselekcja treningu z `training_D0.yaml` (pT ≥ 3 GeV, |η| < 2,2, min pT córki > 0,8 GeV).")
    w("- **Brak tagu AFP.** Pomiar w zdarzeniach nietagowanych; klasy `AFPonly` (łańcuchy AFP) i `AFPtag` (proton w oknie ξ) tylko informacyjnie.")
    w("- Trening K-fold po zdarzeniach (w charmpy `n_folds: 1`); ε_S(MC) bez wag emulacji trygera; ε_S(MC) z fd_nd.")
    if s["substitute"]:
        w(f"- **ZAMIENNIK:** brak torch w jobie → backend `{s['backend']}` zamiast MLP.")
    w("")
    w("## Kryteria rozstrzygnięcia (zapisane przed danymi)")
    w("")
    w(KRYTERIA)
    w("")
    w("## Werdykt")
    w("")
    w(f"Punkt odniesienia WP{V['ref_wp']}.")
    w("")
    w(f"- K0 (poprawność klonu na MC): **{V['K0']}**")
    w("- K1 (R_SF = SF_iso/SF_noiso): " + ", ".join(f"`{c}` **{v}**" for c, v in V["K1"].items()) + f" → 428770: **{V['K1_428770']}**")
    w(f"- K2 (przewaga izolacji w danych 428770): **{V['K2']}**")
    w(f"- K3: **{V['K3']}**")
    w("")
    w("## Wejścia klonu i braki")
    w("")
    w("| wejście CharmNet D⁰ (`training_D0.yaml`) | źródło w h5 v2.0 | wariant |")
    w("|---|---|---|")
    for f, src in res["feat_map"]:
        w(f"| `Dmeson_{f}` | {src} | {'tylko iso' if f in CF.ISO else 'oba'} |")
    w("")
    w("Wszystkie 14 wejść D⁰ jest dostępnych. Braki względem przepisu:")
    w("")
    for m in res["missing"]:
        w(f"- {m}")
    w("")
    w("## Trening")
    w("")
    hp = s["hp"]
    w(f"- Backend `{s['backend']}`, MLP 3 × {hp['hidden']} ReLU + sigmoid, Adam lr {hp['lr']}, {hp['epochs']} epok, batch {hp['batch']}, λ(dCor) {hp['lam']}, "
      f"podział {hp['train_fraction']:.2f}/{1 - hp['train_fraction']:.2f}.")
    w(f"- Próbka: {s['n_train']} kandydatów (sygnał {s['n_train_sig']}, tło {s['n_train_bkg']}), {s['folds']} foldów po zdarzeniach + model pełny. Czas jobu {s['time_s']} s.")
    w("")
    w("| próbka | grupy | inkluzywne D⁰ | pary D⁰←D* | ułamek zdarzeń (inkl.) |")
    w("|---|---|---|---|---|")
    for k, v in s["samples"].items():
        w(f"| {k} | {v['n_groups']} | {v['n_a']} | {v['n_p']} | {v['all_frac']} |")
    w("")
    w("## MC: ROC/AUC (okno |m − m₀| < 55 MeV; sygnał fd_nd, tło mb_nd bez sygnału)")
    w("")
    w("| bin | wariant | N_S / N_B | AUC (OOF) | AUC (model pełny, in-sample) | corr(score, m \\| tło komb.) |")
    w("|---|---|---|---|---|---|")
    for b in BINS:
        for v in VARS:
            m = mc[(v, b)]
            cr = res["diag"][v].get("corr_score_mass_comb")
            w(f"| {b} | {v} | {m['n_sig']} / {m['n_bkg']} | {m['auc_oof']:.4f} | {m['auc_full']:.4f} | {cr:+.3f} |" if cr is not None else
              f"| {b} | {v} | {m['n_sig']} / {m['n_bkg']} | {m['auc_oof']:.4f} | {m['auc_full']:.4f} | – |")
    w("")
    w("### Punkty pracy z MC (próg na score OOF per bin pT) i odrzucenie tła na MC")
    w("")
    w("| bin | WP | próg iso / noiso | ε_B(MC) iso | ε_B(MC) noiso | odrzucenie iso/noiso | G_MC = Z_iso/Z_noiso | ε_S przy tym progu, model pełny (iso / noiso) |")
    w("|---|---|---|---|---|---|---|---|")
    for b in BINS:
        for wp in wps:
            i, n, g = mc[("iso", b)]["wp"][wp], mc[("noiso", b)]["wp"][wp], gain[("MC", "mc", b, wp)]
            w(f"| {b} | {wp} | {i['thr']:.3f} / {n['thr']:.3f} | {pm(i['eps_b'], i['eps_b_err'], 4)} | {pm(n['eps_b'], n['eps_b_err'], 4)} | "
              f"{g['rej_ratio']:.2f} | {pm(g['G'], g['G_err'])} | {i['eps_s_fullmodel']:.3f} / {n['eps_s_fullmodel']:.3f} |")
    w("")
    w("Ostatnia kolumna = efekt braku out-of-fold: ε_S sygnału MC, gdy score pochodzi z modelu trenowanego na tych samych zdarzeniach (nominał = WP).")
    w("")
    w("### Ważność wejść (permutacyjna, spadek AUC OOF w oknie piku)")
    w("")
    w("| wejście | iso: ΔAUC | noiso: ΔAUC |")
    w("|---|---|---|")
    pi, pn = res["diag"]["iso"]["perm_importance"], res["diag"]["noiso"]["perm_importance"]
    for f in sorted(pi, key=lambda k: -pi[k][0]):
        w(f"| {f} | {pi[f][0]:.4f} ± {pi[f][1]:.4f} | " + (f"{pn[f][0]:.4f} ± {pn[f][1]:.4f}" if f in pn else "–") + " |")
    w("")
    w("## SF = ε_S(dane)/ε_S(MC) z fitu Δm pass/fail (tag D*)")
    w("")
    w("ε_S(MC): pary z prawdą z fd_nd. „(n)” = fit nieokreślony. R_SF = SF_iso/SF_noiso (błąd jak dla niezależnych — zachowawczy).")
    for wp in wps:
        w("")
        w(f"### WP {wp}")
        w("")
        w("| run | klasa | bin | N_sig (iso) | ε_S MC iso / noiso | SF iso | SF noiso | R_SF |")
        w("|---|---|---|---|---|---|---|---|")
        for run in CLS:
            for c in CLS[run]:
                for b in BINS:
                    if (run, c, b, "iso", wp) not in sf:
                        continue
                    i, n, g = sf[(run, c, b, "iso", wp)], sf[(run, c, b, "noiso", wp)], gain[(run, c, b, wp)]
                    fl = lambda r: "" if r["determined"] else " (n)"  # noqa: E731
                    w(f"| {run} | {c} | {b}{' (inf.)' if b not in MAIN_BINS else ''} | {i['N_sig']:.0f} | {i['eps_mc']:.3f} / {n['eps_mc']:.3f} | "
                      f"{pm(i['SF'], i['SF_err'])}{fl(i)} | {pm(n['SF'], n['SF_err'])}{fl(n)} | {pm(g['R_SF'], g['R_SF_err'])} |")
    w("")
    w("## ε_B pod pikiem: dane (SB → okno, ML exp(pol2)) i MC (tło mb_nd w oknie)")
    w("")
    w("| run | klasa | bin | WP | ε_B dane iso | ε_B dane noiso | dane/MC iso | dane/MC noiso |")
    w("|---|---|---|---|---|---|---|---|")
    for run in CLS:
        for c in CLS[run]:
            for b in BINS:
                for wp in wps:
                    if (run, c, b, "iso", wp) not in epsb or (run, c, b, "iso", wp) not in sf:
                        continue
                    i, n = epsb[(run, c, b, "iso", wp)], epsb[(run, c, b, "noiso", wp)]
                    mi, mn = mc[("iso", b)]["wp"][wp], mc[("noiso", b)]["wp"][wp]
                    ri = rel(i["eps"], i["err"], mi["eps_b"], mi["eps_b_err"])
                    rn = rel(n["eps"], n["err"], mn["eps_b"], mn["eps_b_err"])
                    w(f"| {run} | {c} | {b} | {wp} | {pm(i['eps'], i['err'], 4)} | {pm(n['eps'], n['err'], 4)} | {pm(*ri)} | {pm(*rn)} |")
    w("")
    w("## Zysk z izolacji: MC vs dane")
    w("")
    w("Z = ε_S/√ε_B (proxy istotności S/√B przy stałych S₀, B₀). G = Z_iso/Z_noiso. „Odrzucenie przy stałej ε_S(dane)” = ε_B(noiso)/ε_B(iso), "
      "gdzie ε_B(noiso) interpolowane po WP w ε_S(dane) klonu iso (– = poza zakresem WP).")
    w("")
    w("| run | klasa | bin | WP | Z_iso (dane) | Z_noiso (dane) | G_dane | G_MC | odrzucenie przy stałej ε_S: dane | MC (stały WP) |")
    w("|---|---|---|---|---|---|---|---|---|---|")
    for run in CLS:
        for c in CLS[run]:
            for b in BINS:
                for wp in wps:
                    g = gain.get((run, c, b, wp))
                    if not g:
                        continue
                    gm = gain[("MC", "mc", b, wp)]
                    rr = g["rej_ratio_fixed_epsS"]
                    w(f"| {run} | {c} | {b} | {wp} | {pm(g['Z_iso'], g['Z_iso_err'])} | {pm(g['Z_noiso'], g['Z_noiso_err'])} | {pm(g['G'], g['G_err'])} | "
                      f"{pm(gm['G'], gm['G_err'])} | {f'{rr:.2f}' if np.isfinite(rr) else '–'} | {gm['rej_ratio']:.2f} |")
    w("")
    w("## Liczba par per klasa")
    w("")
    w(", ".join(f"{k}: {v}" for k, v in res["n_class_pairs"].items()))
    w("")
    w("## Pliki")
    w("")
    w("- Paczka: `charm_work/Q_lxplus/cnclone/` (`cn_feat.py` — wejścia i mapowanie, `cn_extract.py`, `cn_train.py`), etap `ORCH/stages/cnclone.sh`.")
    w("- Wyniki liczbowe: `CN_clone/out/cn_results.json`; modele: `runs/lxplus/cn_models.npz`; diagnostyka treningu: `runs/lxplus/cn_summary.json`.")
    return "\n".join(L) + "\n"


def jsonable(o):
    if isinstance(o, dict):
        return {("|".join(map(str, k)) if isinstance(k, tuple) else str(k)): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(x) for x in o]
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return o


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", default=os.path.join(HERE, "runs", "lxplus"))
    ap.add_argument("--out", default=os.path.join(HERE, "out"))
    ap.add_argument("--report", default=os.path.join(HERE, "RAPORT.md"))
    ap.add_argument("--wps", default="0.95,0.80,0.65,0.50")
    ap.add_argument("--ref-wp", default="0.80")
    ap.add_argument("--classes", default="", help="ogranicz klasy (przecinki), np. all,J12only,MB427929")
    ap.add_argument("--nproc", type=int, default=8)
    ap.add_argument("--fast", action="store_true", help="bez skanu profilu (błędy HESSE)")
    ap.add_argument("--note", default="", help="uwaga wstawiana na początku raportu (np. TEST)")
    a = ap.parse_args(argv)
    res = analyse(a)
    os.makedirs(a.out, exist_ok=True)
    json.dump(jsonable(res), open(os.path.join(a.out, "cn_results.json"), "w"), indent=1)
    open(a.report, "w").write(report(res, a))
    print(f"zapisano {a.report} i {a.out}/cn_results.json")
    print("werdykt:", json.dumps(res["verdict"], ensure_ascii=False))


if __name__ == "__main__":
    sys.exit(main())
