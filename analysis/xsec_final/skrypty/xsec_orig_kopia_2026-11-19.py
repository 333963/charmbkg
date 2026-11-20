#!/usr/bin/env python3
"""XSEC (19.11.2026, projekt B) — inkluzywny przekrój czynny D0 (->K pi) i D+ (->K pi pi), pojedyncze D, bez AFP.

  sigma(bin) = N / (C * SF * BR * L),   C = A x eps_reco x eps_sel = N_reco^prawda(mb_nd, WP70) / N_gen^prawda(mb_nd, mode, bin pT, |eta| < 2.2)

Praca lokalna, tylko odczyt wejsc (nic poza katalogiem XSEC_2026-11-19 nie jest zapisywane; uruchamiac z `python -B`).
Wejscia (sciezki w stalych nizej): uzyski i domkniecie fitu — WPSCAN wps_fit.json (drzewo H, WP70); histogramy — wps_hist.npz (tylko do
wariacji odbicia K<->pi); prawda generatora mb_nd — GEN_mcstat truthD_mb.parquet; SF z tagu D* — DNET dnet4_results.json; swietlnosc —
LUMI_2026-11-19/lumi.json (jesli istnieje; inaczej wynik = sigma x L). Podmiana fitu po FITBIAS: plik `fitbias_override.json` w tym katalogu.

  python -B xsec.py [--sf-scen 0.10] [--no-refl-fit] [--out out]
Wyjscie: out/xsec_results.json, out/tabele.md.
"""
import argparse
import json
import math
import os
import sys

import numpy as np

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
CW = os.path.dirname(HERE)
P = {
    "wps_fit": os.path.join(CW, "WPSCAN_2026-11-19", "wps_fit.json"),
    "wps_hist": os.path.join(CW, "WPSCAN_2026-11-19", "runs", "lxplus", "wps_hist.npz"),
    "wps_dir": os.path.join(CW, "WPSCAN_2026-11-19"),
    "truth": os.path.join(CW, "GEN_mcstat", "runs", "lxplus_v20", "truthD_mb.parquet"),
    "sf": os.path.join(CW, "DNET", "runs", "lxplus4", "dnet4_results.json"),
    "lumi": os.path.join(CW, "LUMI_2026-11-19", "lumi.json"),
    "fitbias_report": os.path.join(CW, "FITBIAS_2026-11-19", "RAPORT.md"),
    "override": os.path.join(HERE, "fitbias_override.json"),
}

# ------------------------------------------------------------------------------------------------ stale (ze zrodlem)
BINS = [(3, 5), (5, 8), (8, 30)]
ETA_FID = 2.2
SPEC = {"D0": {"pdg": 421, "tag": "D0.xgb.H.sb.n.s1"}, "Dplus": {"pdg": 411, "tag": "Dplus.xgb.H.sb.n.s1"}}
RUNS = {"427929": "427929|MB427929", "428770": "428770|all"}
# BR: PDG 2026 wg charm_work/BR_frag_2026-11-18/RAPORT.md §1 (nota §6.9). D0: K-pi+ 3,936 +- 0,030 % + DCS K+pi- 0,0148 % (D0 + anty-D0 bez tagu
# zapachu; prawda MC `mode` = decayType 10 = corki {K, pi} bez wzgledu na ladunek).
BR = {"D0": (0.03936 + 0.000148, 0.00030), "Dplus": (0.0931, 0.0014)}
# Pythia 8.317 A3 samodzielna, 13,6 TeV: sigma_ND [mb] (charm_work/SD_var_2026-11-19/RAPORT.md §3 pkt 3; DIFF_model/RAPORT.md §1). NIE z AMI.
SIGMA_ND_MB = 55.587
# udzial ND wsrod D w akceptancji (A3, suma gatunkow) — DIFF_model/RAPORT.md §2, biny 3-5 / 5-8 / 8-30
F_ND = {"3-5": 0.9712, "5-8": 0.9768, "8-30": 0.9797}
N_EV_MB = 19220986          # zdarzenia mb_nd w 78 plikach — GEN_mcstat/runs/lxplus_v20/summary.md
N_EV_427929_GRL = 23.87e6   # zdarzenia w drzewie 427929 w GRL (458 plikow) — B_trig/PROCEDURA.md §6 (wszystkie w drzewie: 28 962 571)
# odbicie K<->pi: R/S dane/MC = 0,96 +- 0,04 przy WP70, L_xy >= 0,25 mm (nota §5.4) -> wariacja normalizacji +-5 %
REFL_VAR = 0.05
# nieprompt: wariacja udzialu b->D o +-50 % (nota §6.9 proponuje +-30 %; samodzielna A3 i mb_nd roznia sie mocniej, zob. DEFINICJA §3)
FNP_VAR = 0.5
SF_SCEN = (0.05, 0.10, 0.15)
LUMI_RUN, LUMI_TRIG = "427929", "HLT_mb_mbts_L1MBTS_1"
L_VARIANT_427929 = 0.392    # nb^-1: gdyby plateau N_drzewo/TAV = 0,84 bylo ubytkiem plikow, a nie wymogiem wierzcholka (LUMI_2026-11-19/RAPORT.md §4)
LUMI_REL_CONS = 0.07        # wariant konserwatywny bledu kalibracji (rozrzut tagow 001-003; LUMI/RAPORT.md §3)
L_MBTS = (0.4016, 0.4191)   # nb^-1: swietlnosc z licznika MBTS_1 przy sigma_wid 72 / 69 mb (lumi.json _sprawdzian_MBTS) — osobna pozycja budzetu
# CMS 13 TeV, arXiv:2107.01476 tab. 5 (prompt, D + anty-D, |eta| < 2,1), sumy binow [ub] (stat, syst) — arytmetyka z przeglądu literatury
# --- XSEC2 (19.11 noc): rozstrzygniecia; zrodla: WYNIKI.md §8, xsec2/ (skrypty, lxplus, out, fonll)
L_COMPL_REL = 0.003         # kompletnosc plikow w GRL: plateau N_drzewo/TAV = selekcja wierzcholka (nPV = 1 w 100 % zdarzen; udzialy plikow per LB = k/6), rozrzut plateau
EVSEL_REL = 0.002           # wymog wierzcholka (>= 0,998 dla zdarzen z D) i L1_MBTS_1 (>= 0,999 z trygera losowego HLT_*_L1RD0_FILLED) — nie ma ich w C; przyjete 1, blad 0,2 %
# FONLL (formularz sieciowy, NNPDF30_nlo_as_0118, |eta| < 2,2, prompt c->D, FF = 1, tylko czastka; 13,6 TeV: v1.4.0, 13 TeV: v1.3.4) [ub]: (centr., min, max)
FONLL136 = {"D0": {"3-5": (508.4, 218.5, 915.5), "5-8": (181.2, 113.4, 282.2), "8-30": (64.61, 49.23, 87.66)},
            "Dplus": {"3-5": (513.6, 216.6, 930.9), "5-8": (186.0, 114.9, 292.0), "8-30": (67.28, 50.94, 91.92)}}
FONLL_E_RATIO = {"3-5": 1.0285, "5-8": 1.036, "8-30": 1.0455}    # sigma(13,6)/sigma(13) prompt D, |eta| < 2,2 (skale/masy/PDF znosza sie do < 1 %)
FF = {"D0": 0.609, "Dplus": 0.240}                                # f(c->D) uzyte przez CMS w porownaniu z FONLL (arXiv:2107.01476 Sec. 6); ALICE: 0,542 dla D0
R1_L = (0.46758, 0.0224)    # nb^-1: 428770 HLT_noalg_L1AFP_A_OR_C, nasza kopia plikow (lumi.json), blad wzgl.
XS2_FIT = os.path.join(HERE, "xsec2", "out", "xs2_fit.json")
CMS = {"D0": {"5-8": (432, 22, 39, "5–8"), "8-30": (107.8, 2.9, 10.3, "8–40")}, "Dplus": {"5-8": (158.5, 11, 22, "5–8"), "8-30": (50.0, 1.2, 6.3, "8–40")}}


def bname(b):
    return f"{b[0]}-{b[1]}"


def q(*x):
    return math.sqrt(sum(float(v) ** 2 for v in x))


# ------------------------------------------------------------------------------------------------ wejscia
def load_fits():
    f = json.load(open(P["wps_fit"]))
    out = {}
    for x in f["res"]:
        if abs(x["eps_nom"] - 0.70) > 1e-6:
            continue
        for sp, c in SPEC.items():
            if x["sp"] == sp and x["tag"] == c["tag"]:
                out[(sp, x["src"], x["bin"])] = x
    return out


def load_truth():
    import pandas as pd
    T = pd.read_parquet(P["truth"], columns=["pdg", "pt", "eta", "fromB", "mode", "matched", "n_legs", "w2"])
    T["apdg"] = T.pdg.abs()
    T["aeta"] = T.eta.abs()
    return T


def load_sf():
    s = json.load(open(P["sf"]))["species"]["D0"]["sf"]
    out = {}
    for run, key in RUNS.items():
        for b in BINS:
            r = s.get(f"{SPEC['D0']['tag']}|{key}|{bname(b)}|0.7")
            if r:
                out[(run, bname(b))] = r
    return out


def load_lumi():
    """lumi.json (LUMI_2026-11-19, schemat: run -> nazwa trygera -> {"L" [nb^-1], "blad_wzgl", "L_nasza_kopia_plikow", ...}).
    Tylko 427929 / HLT_mb_mbts_L1MBTS_1: nominal = L naszej kopii plikow; wariant = L_VARIANT_427929 (plateau N_drzewo/TAV jako ubytek plikow).
    428770: przekroju bezwzglednego nie liczymy (regiony trygera, LUMI/RAPORT.md §5). Zwraca {run: (L [nb^-1], blad wzgl.)}."""
    if not os.path.exists(P["lumi"]):
        return {}, "brak pliku lumi.json"
    try:
        J = json.load(open(P["lumi"]))
        d = J[LUMI_RUN][LUMI_TRIG]
        L = float(d["L_nasza_kopia_plikow"])
        rel = q(d.get("błąd_wzgl_kalibracja_ZAŁOŻENIE", 0.022), d.get("błąd_wzgl_live", 0.003), L_COMPL_REL)
    except Exception as e:  # noqa: BLE001
        return {}, f"lumi.json: nie rozpoznano schematu ({e!r}) — wynik jako sigma x L"
    note = (f"{LUMI_RUN} {LUMI_TRIG}: L(nasza kopia plikow) = {L:.4f} nb^-1 ± {100 * rel:.1f} % (kalibracja 2,2 % = ZALOZENIE; konserwatywnie {100 * LUMI_REL_CONS:.0f} %); "
            f"L pelne GRL = {d['L']:.4f}; wariant {L_VARIANT_427929} nb^-1 ZAMKNIETY (XSEC2: plateau 0,84 = wymog wierzcholka, nie ubytek plikow); tag {d.get('tag')}; {d.get('n_LB')} LB")
    return {LUMI_RUN: (L, rel)}, note


# ------------------------------------------------------------------------------------------------ wariacja odbicia (lokalny refit)
def refl_variation(fits):
    """Wzgledna zmiana N(D0) przy f_refl x (1 +- REFL_VAR): powtorzenie nominalnego fitu WPSCAN (ten sam kod, te same histogramy)."""
    sys.path.insert(0, P["wps_dir"])
    import wps_fit as W
    A = W.A
    W._init(P["wps_hist"], 0)
    Z, INFO = W.Z, W.INFO
    sp, c = "D0", A.SP["D0"]
    lo, hi = INFO["sp"][sp]["mr"]
    j = INFO["eps_desc"].index(0.7)
    cum = lambda h: h[..., j + 1:, :].sum(-2)  # noqa: E731
    out = {}
    for b in BINS:
        key = f"{sp}|{SPEC[sp]['tag']}|{bname(b)}"
        sw, rw = cum(Z[key + "|fd|sig_w"]), cum(Z[key + "|fd|refl_w"])
        sigT, refT, fr = W.tmpl(sw, lo, hi), W.tmpl(rw, lo, hi), float(rw.sum() / sw.sum())
        flo, fhi = c["fit"]
        edges = np.arange(flo, fhi + A.BW / 2, A.BW)
        for run, src in RUNS.items():
            H = W.rebin(cum(Z[f"{key}|{src}"])[0], lo, flo, fhi)
            n = {}
            for lab, k in (("nom", 1.0), ("up", 1 + REFL_VAR), ("dn", 1 - REFL_VAR)):
                f = W.one(H, sigT, refT, fr * k, edges, c["m0"], c["nbkg"], c["win"])
                n[lab] = f["N"] if f["valid"] else float("nan")
            ref = fits[(sp, src, bname(b))]["nom"]["N"]
            out[(run, bname(b))] = {"N_nom_refit": n["nom"], "N_wps": ref, "repro": n["nom"] / ref - 1,
                                    "rel": max(abs(n["up"] / n["nom"] - 1), abs(n["dn"] / n["nom"] - 1)), "up": n["up"] / n["nom"] - 1, "dn": n["dn"] / n["nom"] - 1}
    return out


# ------------------------------------------------------------------------------------------------ obliczenia
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--sf-scen", type=float, default=None, help="scenariusz bledu nieznanego SF (domyslnie 0,10 dla 427929, 0,15 dla 428770 — nota §6.8)")
    ap.add_argument("--no-refl-fit", action="store_true")
    ap.add_argument("--out", default=os.path.join(HERE, "out"))
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    fits, T, SFD = load_fits(), load_truth(), load_sf()
    lumi, lumi_note = load_lumi()
    override = json.load(open(P["override"])) if os.path.exists(P["override"]) else {}
    fitbias_ready = os.path.exists(P["fitbias_report"])
    refl = {} if a.no_refl_fit else refl_variation(fits)
    R = {"inputs": P, "lumi": {k: list(v) for k, v in lumi.items()}, "lumi_note": lumi_note, "fitbias_report_exists": fitbias_ready,
         "override_used": bool(override), "cells": {}, "mc": {}}

    # ---- prawda mb_nd: mianowniki, BR w MC, nieprompt, przewidywanie Pythii
    for sp, c in SPEC.items():
        s = T[T.apdg == c["pdg"]]
        fid = s[(s.aeta < ETA_FID)]
        allpt = fid[fid.pt >= 3]
        k = int((allpt["mode"] == 1).sum())
        R["mc"][sp] = {"BR_mc": k / len(allpt), "BR_mc_err": math.sqrt(k) / len(allpt), "BR_pdg": BR[sp][0], "bins": {}}
        for b in BINS:
            bn = bname(b)
            inb = fid[(fid.pt >= b[0]) & (fid.pt < b[1])]
            sig = inb[inb["mode"] == 1]
            n_all, n_sig = len(inb), len(sig)
            p_, np_ = sig[sig.fromB == 0], sig[sig.fromB == 1]
            r_np = float((np_.n_legs > 0).mean() / (p_.n_legs > 0).mean())      # proxy: nogi selekcji `direct` (charmtag), nie WP70
            f_np = len(np_) / n_sig
            sig25 = s[(s.aeta < 2.5) & (s.pt >= b[0]) & (s.pt < b[1]) & (s["mode"] == 1)]
            out_fid = s[(s.aeta >= ETA_FID) & (s.pt >= b[0]) & (s.pt < b[1]) & (s["mode"] == 1)]
            R["mc"][sp]["bins"][bn] = {
                "N_gen_sig": n_sig, "N_gen_all": n_all, "f_np": f_np, "f_np_all": float(inb.fromB.mean()), "r_np_proxy": r_np,
                "matched_frac": float(sig.matched.mean()), "matched_frac_w2": float(np.average(sig.matched, weights=sig.w2)),
                "frac_matched_outside_fid": float(out_fid.matched.sum() / max(sig.matched.sum() + out_fid.matched.sum(), 1)),
                "N_gen_sig_eta25": len(sig25),
                "sigma_pythia_ND_ub": n_all / N_EV_MB * SIGMA_ND_MB * 1e3, "sigma_pythia_ND_ub_err": math.sqrt(n_all) / N_EV_MB * SIGMA_ND_MB * 1e3,
                "sigma_pythia_inel_ub": n_all / N_EV_MB * SIGMA_ND_MB * 1e3 / F_ND[bn], "per_event_ND": n_all / N_EV_MB}

    # ---- komorki gatunek x run x bin
    for sp in SPEC:
        for run, src in RUNS.items():
            for b in BINS:
                bn = bname(b)
                x, mb = fits[(sp, src, bn)], fits[(sp, "mb|all", bn)]
                mc = R["mc"][sp]["bins"][bn]
                N, hesse = x["nom"]["N"], x["nom"]["err"]
                boot = float(np.std([v for v in x["boot_N"] if v is not None], ddof=1))
                stat = max(hesse, boot)
                var = {k: v["N"] / N - 1 for k, v in x["var"].items()}
                fit_model = max(abs(v) for v in var.values())
                fit_model_pm1 = max(abs(v) for k, v in var.items() if k.startswith("exp"))
                # obciazenie fitu: domkniecie na pseudo-danych mb_nd (jedna realizacja tla)
                b_fit = mb["nom"]["N"] / mb["N_truth"] - 1
                b_err = mb["nom"]["err"] / mb["N_truth"]
                fit_bias = q(b_fit, b_err)
                ov = override.get(f"{sp}|{run}|{bn}")
                if ov:  # podmiana po FITBIAS: {"N":..., "stat":..., "fit_model_rel":..., "fit_bias_rel":...}
                    N, stat = ov.get("N", N), ov.get("stat", stat)
                    fit_model, fit_bias = ov.get("fit_model_rel", fit_model), ov.get("fit_bias_rel", fit_bias)
                # wydajnosc calkowita z mb_nd
                C = mb["N_truth"] / mc["N_gen_sig"]
                C_err = math.sqrt(C * (1 - C) / mc["N_gen_sig"]) / C
                # SF
                sfm = SFD.get((run, bn)) if sp == "D0" else None
                if sfm and sfm.get("determined") and b[0] >= 5:
                    SF, SF_rel, SF_kind = sfm["SF"], sfm["SF_err"] / sfm["SF"], "zmierzony (tag D*, WP70)"
                else:
                    scen = a.sf_scen if a.sf_scen is not None else (0.10 if run == "427929" else 0.15)
                    SF, SF_rel, SF_kind = 1.0, scen, f"SCENARIUSZ ±{scen:.0%} (SF nieznany, przyjeto 1)"
                br, br_e = BR[sp]
                rf = refl.get((run, bn), {}).get("rel", float("nan")) if sp == "D0" else 0.0
                fnp = FNP_VAR * mc["f_np"] * abs(mc["r_np_proxy"] - 1) / (1 + mc["f_np"] * (mc["r_np_proxy"] - 1))
                Ncorr = N / (C * SF * br)                       # = sigma x L  [liczba mezonow D + anty-D w obszarze fiducjalnym]
                rel = {"stat": stat / N, "fit_model": fit_model, "fit_bias": fit_bias, "SF": SF_rel, "mc_stat": C_err, "BR": br_e / br,
                       "refl": 0.0 if (sp != "D0" or math.isnan(rf)) else rf, "nonprompt": fnp, "evsel": EVSEL_REL}
                syst = q(*[v for k, v in rel.items() if k != "stat"])
                cell = {"N": N, "stat_hesse": hesse, "stat_boot": boot, "var": var, "fit_model_pm1": fit_model_pm1, "closure_b": b_fit, "closure_b_err": b_err,
                        "N_reco_truth_mb": mb["N_truth"], "N_gen_sig_mb": mc["N_gen_sig"], "C": C, "eps_cls_fd": x["eps"], "SF": SF, "SF_kind": SF_kind,
                        "SF_measured_any": sfm, "BR": br, "rel": rel, "syst": syst, "tot": q(rel["stat"], syst), "sigmaL": Ncorr,
                        "dsigmaL_dpt": Ncorr / (b[1] - b[0]), "sigmaL_closure_corrected": Ncorr / (1 + b_fit),
                        "syst_closure_corrected": q(*[v for k, v in rel.items() if k not in ("stat", "fit_bias")], b_err / (1 + b_fit)),
                        "SB_win": N * 0.0 if x["nom"]["bkg_win"] == 0 else None, "bkg_win": x["nom"]["bkg_win"], "chi2": x["nom"]["chi2"], "ndf": x["nom"]["ndf"],
                        "refl_detail": refl.get((run, bn)) if sp == "D0" else None}
                for sc in SF_SCEN:
                    cell[f"tot_sf{int(sc * 100)}"] = q(rel["stat"], *[v for k, v in rel.items() if k not in ("stat", "SF")], sc if SF_kind.startswith("SCEN") else SF_rel)
                if run in lumi:
                    L, Lrel = lumi[run]
                    cell["sigma_nb"] = Ncorr / L
                    cell["dsigma_dpt_nb"] = Ncorr / L / (b[1] - b[0])
                    cell["lumi_rel"] = Lrel
                R["cells"][f"{sp}|{run}|{bn}"] = cell

    # ---- wielkosci niezalezne od swietlnosci
    def cell(sp, run, bn):
        return R["cells"][f"{sp}|{run}|{bn}"]
    indep = ("stat", "fit_model", "fit_bias", "mc_stat", "BR", "refl", "nonprompt")
    R["ratio_DpD0"] = {}
    for run in RUNS:
        for b in BINS:
            bn = bname(b)
            p, z = cell("Dplus", run, bn), cell("D0", run, bn)
            # (i) SF+/SF0 = 1: SF sie skraca; blad scenariuszowy stosunku osobno
            r1 = (p["N"] / (p["C"] * p["BR"])) / (z["N"] / (z["C"] * z["BR"]))
            e = {k: q(p["rel"][k], z["rel"][k]) for k in indep}
            syst = q(*[v for k, v in e.items() if k != "stat"])
            # (ii) D0 ze zmierzonym SF (>= 5 GeV), D+ z SF = 1 +- scenariusz
            r2 = r1 * z["SF"] if z["SF_kind"].startswith("zmierz") else None
            # (iii) poprawka domkniecia MC na oba uzyski
            r3 = r1 * (1 + z["closure_b"]) / (1 + p["closure_b"])
            e3 = dict(e); e3["fit_bias"] = q(p["closure_b_err"] / (1 + p["closure_b"]), z["closure_b_err"] / (1 + z["closure_b"]))
            mcp, mcz = R["mc"]["Dplus"]["bins"][bn], R["mc"]["D0"]["bins"][bn]
            pyth = mcp["N_gen_all"] / mcz["N_gen_all"]
            R["ratio_DpD0"][f"{run}|{bn}"] = {
                "R_sf_equal": r1, "rel": e, "syst_noSF": syst, "tot_noSF": q(e["stat"], syst),
                **{f"tot_sfratio{int(s * 100)}": q(e["stat"], syst, s) for s in SF_SCEN},
                "R_sf0_measured": r2, "sf0_rel": z["rel"]["SF"] if r2 else None,
                "R_closure_corrected": r3, "tot_closure_corrected_noSF": q(*e3.values()),
                "pythia_A3_mbnd": pyth, "pythia_err": pyth * q(1 / math.sqrt(mcp["N_gen_all"]), 1 / math.sqrt(mcz["N_gen_all"]))}
    R["pt_shape"] = {}
    for sp in SPEC:
        for run in RUNS:
            for (b1, b2) in ((BINS[1], BINS[0]), (BINS[2], BINS[1])):
                n1, n2 = cell(sp, run, bname(b1)), cell(sp, run, bname(b2))
                r = n1["sigmaL"] / n2["sigmaL"]
                e = {k: q(n1["rel"][k], n2["rel"][k]) for k in ("stat", "fit_model", "fit_bias", "mc_stat", "refl", "nonprompt")}
                sfm = q(n1["rel"]["SF"], n2["rel"]["SF"]) if (n1["SF_kind"].startswith("zmierz") and n2["SF_kind"].startswith("zmierz")) else None
                m1, m2 = R["mc"][sp]["bins"][bname(b1)], R["mc"][sp]["bins"][bname(b2)]
                R["pt_shape"][f"{sp}|{run}|{bname(b1)}/{bname(b2)}"] = {
                    "ratio": r, "rel": e, "tot_noSF": q(*e.values()), "sf_measured_rel": sfm,
                    "ratio_dsdpt": r * (b2[1] - b2[0]) / (b1[1] - b1[0]),
                    "pythia": m1["N_gen_all"] / m2["N_gen_all"], "pythia_err": m1["N_gen_all"] / m2["N_gen_all"] * q(1 / math.sqrt(m1["N_gen_all"]), 1 / math.sqrt(m2["N_gen_all"]))}
    R["run_ratio"] = {}
    for b in BINS:
        bn = bname(b)
        rr = {}
        for sp in SPEC:
            a1, a2 = cell(sp, "427929", bn), cell(sp, "428770", bn)
            r = a1["N"] / a2["N"]                                 # ta sama wydajnosc MC w obu runach -> skraca sie; SF osobno
            e = q(a1["rel"]["stat"], a2["rel"]["stat"], a1["rel"]["fit_model"], a2["rel"]["fit_model"])
            sfr = (a2["SF"] / a1["SF"]) if (a1["SF_kind"].startswith("zmierz") and a2["SF_kind"].startswith("zmierz")) else None
            rr[sp] = {"ratio": r, "rel": e, "ratio_with_measured_SF": r * sfr if sfr else None,
                      "sf_rel": q(a1["rel"]["SF"], a2["rel"]["SF"]) if sfr else None}
        dr = rr["Dplus"]["ratio"] / rr["D0"]["ratio"]
        rr["double_ratio_Dp_over_D0"] = {"ratio": dr, "rel": q(rr["Dplus"]["rel"], rr["D0"]["rel"])}
        R["run_ratio"][bn] = rr
    # uzysk na zdarzenie (diagnostyka niezalezna od L; mianownik = zdarzenia drzewa w GRL, zob. zastrzezenie w DEFINICJA §6)
    R["per_event_427929"] = {}
    for sp in SPEC:
        for b in BINS:
            bn = bname(b)
            c_, mc = cell(sp, "427929", bn), R["mc"][sp]["bins"][bn]
            y = c_["sigmaL"] / N_EV_427929_GRL
            R["per_event_427929"][f"{sp}|{bn}"] = {"data": y, "rel_tot": c_["tot"], "pythia_ND": mc["per_event_ND"],
                                                    "pythia_inel": mc["per_event_ND"] * SIGMA_ND_MB / 78.58 / F_ND[bn], "data_over_pythia_ND": y / mc["per_event_ND"]}
    R["r1"] = r1_block(R)
    json.dump(R, open(os.path.join(a.out, "xsec_results.json"), "w"), indent=1, default=float)
    write_tables(R, os.path.join(a.out, "tabele.md"))
    print(open(os.path.join(a.out, "tabele.md")).read())
    return 0


def r1_block(R):
    """XSEC2: uzyski w klasach trygera AFP (xsec2/out/xs2_fit.json) -> sigma x P(L1_AFP | D) i P(AFP | D)."""
    if not os.path.exists(XS2_FIT) or LUMI_RUN not in R["lumi"]:
        return None
    J = json.load(open(XS2_FIT))
    L427 = R["lumi"][LUMI_RUN][0]
    out, acc = {"cells": {}, "tot": {}}, {"R1 (428770)": [], "MBafp (427929)": []}
    for sp in SPEC:
        m = "M8" if sp == "D0" else "M1"
        for b in BINS:
            bn = bname(b)
            c = R["cells"][f"{sp}|427929|{bn}"]

            def ne(cl):
                x = J[f"{sp}|{bn}|{cl}"][m]
                floor = math.sqrt(J[f"{sp}|{bn}|{cl}"]["n_win"]) if cl in ("R1", "MBafp") else 0.0   # male liczby: HESSE/bootstrap (30 replik) zanizaja -> podloga sqrt(N w oknie)
                return x["N"], max(x["err"] or 0.0, x.get("boot_rms") or 0.0, floor)
            (n1, e1), (na, ea), (nm, em) = ne("R1"), ne("MBafp"), ne("MB427929")
            k = 1.0 / (c["C"] * c["BR"]) / 1e3            # SF = 1
            s427 = nm * k / L427
            x = {"N_R1": n1, "e_R1": e1, "N_afp": na, "e_afp": ea, "N_MB": nm, "e_MB": em, "model": m,
                 "sig_R1": n1 * k / R1_L[0], "esig_R1": e1 * k / R1_L[0], "sig_afp": na * k / L427, "esig_afp": ea * k / L427,
                 "P_R1": n1 / R1_L[0] / (nm / L427), "P_afp": na / nm}
            x["eP_R1"] = math.sqrt((e1 / R1_L[0] / (nm / L427)) ** 2 + (x["P_R1"] * em / nm) ** 2)
            x["eP_afp"] = math.sqrt((ea / nm) ** 2 + (x["P_afp"] * em / nm) ** 2)
            out["cells"][f"{sp}|{bn}"] = x
            acc["R1 (428770)"].append((x["P_R1"], x["eP_R1"])); acc["MBafp (427929)"].append((x["P_afp"], x["eP_afp"]))
    for n, v in acc.items():
        w = np.array([1 / e ** 2 for _, e in v]); p = np.array([p_ for p_, _ in v])
        mu = float((w * p).sum() / w.sum())
        out["tot"][n] = (mu, float(1 / math.sqrt(w.sum())), float((w * (p - mu) ** 2).sum()), len(v) - 1)
    return out


def pc(v, n=1):
    return "—" if v is None or (isinstance(v, float) and math.isnan(v)) else f"{100 * v:.{n}f}"


def write_tables(R, path):
    L = ["# XSEC — tabele (generowane przez xsec.py; nie edytowac recznie)", ""]
    L.append(f"Swietlnosc: {R['lumi_note']}" + ("" if R["lumi"] else " → wyniki jako σ × L (liczba mezonow D + anty-D w obszarze fiducjalnym), L nie jest zgadywane."))
    L.append(f"FITBIAS/RAPORT.md istnieje: {R['fitbias_report_exists']}; podmiana `fitbias_override.json` uzyta: {R['override_used']}.")
    L += ["", "## T0. Prawda MC mb_nd (Pythia 8 A3 ND, mc23 801451, 78 plikow, 19 220 986 zdarzen, bez wag): mianowniki i kontrola",
          "", "| gatunek | bin pT | N_gen (mod sygnalu, \\|η\\| < 2,2) | N_gen (wszystkie mody) | f_nieprompt | ε_np/ε_p (proxy) | ε presel. | ε presel. z wagami w2 (428770) | reco spoza \\|η\\| < 2,2 |",
          "|---|---|---|---|---|---|---|---|---|"]
    for sp, m in R["mc"].items():
        for bn, x in m["bins"].items():
            L.append(f"| {sp} | {bn} | {x['N_gen_sig']} | {x['N_gen_all']} | {x['f_np']:.3f} | {x['r_np_proxy']:.2f} | {x['matched_frac']:.3f} | {x['matched_frac_w2']:.3f} | {pc(x['frac_matched_outside_fid'])} % |")
    L += ["", "BR w MC (ulamek modu sygnalu wsrod prawdziwych D, pT > 3, |η| < 2,2) wobec PDG: " +
          "; ".join(f"{sp}: {100 * m['BR_mc']:.2f} ± {100 * m['BR_mc_err']:.2f} % (PDG {100 * m['BR_pdg']:.3f} %, MC/PDG = {m['BR_mc'] / m['BR_pdg']:.3f})" for sp, m in R["mc"].items())]
    L += ["", "## T1. Skorygowane uzyski σ × L i budzet bledu (WP70, drzewo H) — bledy wzgledne w %",
          "", "| gatunek | run | bin | N_fit | C = A·ε [%] | SF | σ × L [10³] | dσ/dpT × L [10³/GeV] | stat | fit: model tla | fit: obciazenie (domkn. MC) | SF | stat. MC | BR | odbicie | nieprompt | wybor zdarzen + tryger | syst | razem | swietlnosc |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for k, c in R["cells"].items():
        sp, run, bn = k.split("|")
        r = c["rel"]
        sfl = f"{c['SF']:.3f} (zm.)" if c["SF_kind"].startswith("zmierz") else "1 (SCEN.)"
        L.append(f"| {sp} | {run} | {bn} | {c['N']:.0f} | {100 * c['C']:.2f} | {sfl} | {c['sigmaL'] / 1e3:.1f} | {c['dsigmaL_dpt'] / 1e3:.2f} | {pc(r['stat'])} | {pc(r['fit_model'])} | "
                 f"{pc(r['fit_bias'])} ({100 * c['closure_b']:+.0f} ± {100 * c['closure_b_err']:.0f}) | {pc(r['SF'])}{'' if c['SF_kind'].startswith('zmierz') else '*'} | {pc(r['mc_stat'])} | {pc(r['BR'])} | "
                 f"{pc(r['refl'])} | {pc(r['nonprompt'])} | {pc(r['evsel'])} | {pc(c['syst'])} | {pc(c['tot'])} | {pc(c.get('lumi_rel')) if 'lumi_rel' in c else 'brak L'} |")
    L += ["", "\\* SF nieznany: przyjeto 1, blad = scenariusz (427929: 10 %, 428770: 15 %). Razem przy scenariuszach 5 / 10 / 15 %:", "",
          "| gatunek | run | bin | razem (5 %) | razem (10 %) | razem (15 %) | σ × L po poprawce domkniecia MC [10³] | syst po poprawce |", "|---|---|---|---|---|---|---|---|"]
    for k, c in R["cells"].items():
        sp, run, bn = k.split("|")
        L.append(f"| {sp} | {run} | {bn} | {pc(c['tot_sf5'])} | {pc(c['tot_sf10'])} | {pc(c['tot_sf15'])} | {c['sigmaL_closure_corrected'] / 1e3:.1f} | {pc(c['syst_closure_corrected'])} |")
    if R["lumi"]:
        L += ["", "## T1b. Przekroj bezwzgledny, run 427929 (D + anty-D, prompt + nieprompt, |η| < 2,2); L = nasza kopia plikow (ostateczna po XSEC2)",
              "", "| gatunek | bin | SF | σ [µb] | ± stat | ± syst (bez L) | ± L | dσ/dpT [µb/GeV] | razem [%] | razem przy L ± 7 % [%] | σ gdyby L z licznika MBTS (0,402–0,419; NIE wariant nominalny) [µb] | Pythia A3 [µb] | dane/Pythia | CMS 13 TeV prompt \\|η\\| < 2,1 [µb] | dane w definicji CMS, 13,6 TeV [µb] | to samo przeskalowane do 13 TeV (FONLL) [µb] | dane/CMS (po poprawce na energie) | FONLL 13,6 TeV prompt, 2·f·σ [µb] (pasmo) | dane prompt / FONLL centr. |",
              "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
        for k, c in R["cells"].items():
            if "sigma_nb" not in c:
                continue
            sp, run, bn = k.split("|")
            sg = c["sigma_nb"] / 1e3
            mc = R["mc"][sp]["bins"][bn]
            py = mc["sigma_pythia_inel_ub"]
            tot, totc = q(c["tot"], c["lumi_rel"]), q(c["tot"], LUMI_REL_CONS)
            Lnom = R["lumi"][run][0]
            cm = CMS[sp].get(bn)
            conv = sg * (2.1 / ETA_FID) * (1 - mc["f_np"])
            conv13 = conv / FONLL_E_RATIO[bn]
            cms = "—" if not cm else f"{cm[0]} ± {cm[1]} ± {cm[2]} ({cm[3]})"
            rc = "—" if not cm else f"{conv13 / cm[0]:.2f} ± {conv13 / cm[0] * q(tot, q(cm[1], cm[2]) / cm[0]):.2f}"
            fo = [2 * FF[sp] * v for v in FONLL136[sp][bn]]
            prm = sg * (1 - mc["f_np"])
            c["sigma_cmsdef_13TeV_ub"], c["fonll_prompt_ub"], c["sigma_prompt_ub"] = conv13, fo, prm
            sfl = "zmierzony" if c["SF_kind"].startswith("zmierz") else "**SCENARIUSZ**"
            L.append(f"| {sp} | {bn} | {sfl} | {sg:.0f} | {sg * c['rel']['stat']:.0f} | {sg * c['syst']:.0f} | {sg * c['lumi_rel']:.0f} | {sg / (float(bn.split('-')[1]) - float(bn.split('-')[0])):.1f} | "
                     f"{pc(tot)} | {pc(totc)} | {sg * Lnom / L_MBTS[1]:.0f}–{sg * Lnom / L_MBTS[0]:.0f} | {py:.0f} | {sg / py:.2f} ± {sg / py * tot:.2f} | {cms} | {conv:.0f} | {conv13:.0f} | {rc} | "
                     f"{fo[0]:.0f} ({fo[1]:.0f}–{fo[2]:.0f}) | {prm / fo[0]:.2f} ± {prm / fo[0] * tot:.2f} |")
        L += ["", "Dane w definicji CMS = σ × (2,1/2,2) × (1 − f_nieprompt z mb_nd); poprawka na energie: ÷ σ_FONLL(13,6)/σ_FONLL(13) = 1,029 / 1,036 / 1,046; bin 8–30 porownany z CMS 8–40 (24–40 wnosi 1,6 µb D⁰).",
              "FONLL: formularz sieciowy, NNPDF30_nlo_as_0118, |η| < 2,2, prompt, × 2 (czastka + antyczastka) × f(c→D⁰) = 0,609, f(c→D⁺) = 0,240 (wartosci uzyte przez CMS); pasmo = skale ⊕ masy ⊕ PDF. Dane prompt = σ × (1 − f_nieprompt z mb_nd).",
              "428770 „all”: przekroju bezwzglednego nie liczymy (regiony R1 0,469 / R2 27,65 / R3 29,41 nb⁻¹). Region R1: T6."]
        if R.get("r1"):
            L += ["", "## T6. Klasy trygera AFP (XSEC2): 428770 R1 = HLT_noalg_L1AFP_A_OR_C ∧ GRL (preskala 60/120, L = 0,4676 nb⁻¹) i 427929 MBafp = MB427929 ∧ L1_AFP_A_OR_C (preskala 1, L = 0,4648 nb⁻¹)",
                  "", "To NIE jest przekroj inkluzywny: σ_AFP = σ × P(L1_AFP_A_OR_C | zdarzenie z D w akceptancji). Ta sama C, SF = 1 (bez bledu SF), BR. Fit: D⁰ M8, D⁺ M1 (ksztalt z ograniczeniem Gaussa); blad = max(HESSE, bootstrap, √N_okno dla podklas).",
                  "", "| gatunek | bin | N R1 (428770) | N MBafp (427929) | N MB427929 (ten sam model) | σ_AFP R1 [µb] | σ_AFP MBafp [µb] | P(AFP \\| D) R1 = σ_R1/σ_427929 [%] | P(AFP \\| D) MBafp = N_MBafp/N_MB [%] |",
                  "|---|---|---|---|---|---|---|---|---|"]
            for k, x in R["r1"]["cells"].items():
                sp, bn = k.split("|")
                L.append(f"| {sp} | {bn} | {x['N_R1']:.0f} ± {x['e_R1']:.0f} | {x['N_afp']:.0f} ± {x['e_afp']:.0f} | {x['N_MB']:.0f} ± {x['e_MB']:.0f} | {x['sig_R1']:.1f} ± {x['esig_R1']:.1f} | {x['sig_afp']:.1f} ± {x['esig_afp']:.1f} | "
                         f"{100 * x['P_R1']:.2f} ± {100 * x['eP_R1']:.2f} | {100 * x['P_afp']:.2f} ± {100 * x['eP_afp']:.2f} |")
            t = R["r1"]["tot"]
            L += ["", "Srednie wazone po binach i gatunkach: " + "; ".join(f"{n}: {100 * v[0]:.2f} ± {100 * v[1]:.2f} % (χ²/ndf = {v[2]:.1f}/{v[3]})" for n, v in t.items()),
                  "Odniesienie (wszystkie zdarzenia z wierzcholkiem): 427929 P(L1_AFP | MBTS) = 2,40 % (0,53 % dla 20–40 torow, 0,13 % dla 40–80; zdarzenia z kandydatem D: 0,57 %); 428770 strumien J12: 0,58 %.",
                  "Poziom przypadkowy (AFP na przeciecie): 427929 ≈ 0,045 %, 428770 ≈ 0,19 %."]
    L += ["", "## T2. Stosunek D⁺/D⁰ (σ(D⁺ + D⁻)/σ(D⁰ + anty-D⁰), |η| < 2,2; swietlnosc sie skraca)",
          "", "| run | bin | R (SF⁺ = SF⁰) | stat | fit model | fit obciazenie | stat. MC | BR | odbicie | nieprompt | razem bez SF | razem, SF⁺/SF⁰ = 1 ± 5 / 10 / 15 % | R z SF⁰ zmierzonym (SF⁺ = 1) | R po poprawce domkniecia MC (± razem bez SF) | Pythia A3 (prawda mb_nd) |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for k, x in R["ratio_DpD0"].items():
        run, bn = k.split("|")
        e = x["rel"]
        r2 = "—" if x["R_sf0_measured"] is None else f"{x['R_sf0_measured']:.3f} (SF⁰ ± {pc(x['sf0_rel'])} %)"
        L.append(f"| {run} | {bn} | {x['R_sf_equal']:.3f} | {pc(e['stat'])} | {pc(e['fit_model'])} | {pc(e['fit_bias'])} | {pc(e['mc_stat'])} | {pc(e['BR'])} | {pc(e['refl'])} | {pc(e['nonprompt'])} | "
                 f"{pc(x['tot_noSF'])} (± {x['R_sf_equal'] * x['tot_noSF']:.3f}) | {pc(x['tot_sfratio5'])} / {pc(x['tot_sfratio10'])} / {pc(x['tot_sfratio15'])} | {r2} | "
                 f"{x['R_closure_corrected']:.3f} ± {x['R_closure_corrected'] * x['tot_closure_corrected_noSF']:.3f} | {x['pythia_A3_mbnd']:.3f} ± {x['pythia_err']:.3f} |")
    L += ["", "## T3. Ksztalt widma: stosunki σ × L sasiednich binow (SF: przyjety plaski tam, gdzie nieznany — blad SF nie wliczony)",
          "", "| gatunek | run | biny | stosunek σ | stosunek dσ/dpT | stat | fit model | fit obciazenie | stat. MC | razem bez SF | SF zmierzony (oba biny) | Pythia A3 |", "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for k, x in R["pt_shape"].items():
        sp, run, bb = k.split("|")
        e = x["rel"]
        L.append(f"| {sp} | {run} | {bb} | {x['ratio']:.3f} | {x['ratio_dsdpt']:.4f} | {pc(e['stat'])} | {pc(e['fit_model'])} | {pc(e['fit_bias'])} | {pc(e['mc_stat'])} | {pc(x['tot_noSF'])} | "
                 f"{'± ' + pc(x['sf_measured_rel']) + ' %' if x['sf_measured_rel'] else '—'} | {x['pythia']:.3f} ± {x['pythia_err']:.3f} |")
    L += ["", "## T4. Test spojnosci runow: N(427929)/N(428770) (ta sama wydajnosc MC; blad = stat ⊕ model tla)",
          "", "| bin | D⁰ | D⁰ ze zmierzonymi SF obu runow | D⁺ | podwojny stosunek (D⁺/D⁰) | uzysk D* 427929/428770 (nota §2) |", "|---|---|---|---|---|---|"]
    ds = {"3-5": "niemierzalny", "5-8": "0,30 ± 0,14", "8-30": "0,14 ± 0,02"}
    for bn, x in R["run_ratio"].items():
        z, p, d = x["D0"], x["Dplus"], x["double_ratio_Dp_over_D0"]
        zs = "—" if z["ratio_with_measured_SF"] is None else f"{z['ratio_with_measured_SF']:.3f} (SF ± {pc(z['sf_rel'])} %)"
        L.append(f"| {bn} | {z['ratio']:.3f} ± {z['ratio'] * z['rel']:.3f} | {zs} | {p['ratio']:.3f} ± {p['ratio'] * p['rel']:.3f} | {d['ratio']:.2f} ± {d['ratio'] * d['rel']:.2f} | {ds[bn]} |")
    L += ["", "## T5. Przewidywanie odniesienia: Pythia 8 A3 (prawda mb_nd × σ_ND = 55,587 mb z samodzielnej Pythii 8.317; σ_gen z AMI NIEDOSTEPNE)",
          "", "| gatunek | bin | σ_ND [µb] (D + anty-D, prompt + nieprompt, \\|η\\| < 2,2) | σ wszystkie procesy [µb] (÷ f_ND) | dσ/dpT [µb/GeV] | D na zdarzenie ND | dane 427929: (σ × L)/N_zdarzen(GRL) | dane / Pythia ND (na zdarzenie) |",
          "|---|---|---|---|---|---|---|---|"]
    for sp, m in R["mc"].items():
        for bn, x in m["bins"].items():
            pe = R["per_event_427929"][f"{sp}|{bn}"]
            w = float(bn.split("-")[1]) - float(bn.split("-")[0])
            L.append(f"| {sp} | {bn} | {x['sigma_pythia_ND_ub']:.1f} ± {x['sigma_pythia_ND_ub_err']:.1f} | {x['sigma_pythia_inel_ub']:.1f} | {x['sigma_pythia_inel_ub'] / w:.1f} | {x['per_event_ND']:.3e} | "
                     f"{pe['data']:.3e} ± {pc(pe['rel_tot'])} % | {pe['data_over_pythia_ND']:.2f} |")
    L.append("")
    open(path, "w").write("\n".join(L))


if __name__ == "__main__":
    sys.exit(main())
