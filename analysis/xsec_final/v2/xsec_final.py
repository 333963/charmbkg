#!/usr/bin/env python3
"""XSEC_final v2 (20.11.2026, projekt B; kopia ../xsec_final.py ze zmianami — lista: ZMIANY.md; pliki v1 NIE sa zmieniane).
Zmiany v2: (1) D+: uzysk, blad stat., model fitu i domkniecie z FITBIAS_Dplus_2026-11-20/out/dplus_yields.json (model A2); (2) SF D+ 3-5 GeV: blad 0,08 (DECYZJA);
(3) SF D0 (tag D*) 5-8, 8-30: (+) polowa rozpietosci 4 estymatorow Dm (AUDYT_noty pkt 6 b; DECYZJA); (4) L_rel = 0,0103/0,4648 (lumi.json blad_wzgl; v1: (+) kompletnosc plikow 0,3 % z xsec.py -> 2,24 %).

XSEC_final (20.11.2026, projekt B) — koncowy rachunek sigma(D0), sigma(D+) w runie 427929. Definicje: DEFINICJA.md (zapisane przed fitami).

  sigma(bin) = N / (C_ekskl * SF * BR_PDG * L),  C_ekskl = N_reco^prawda(mb_nd, WP70, bin pT reco) / N_gen^ekskl(mb_nd, K-pi+ / K-pi+pi+ (+FSR), bin pT prawdy, |eta| < 2.2)

Wlasny rachunek zastepujacy XSEC_2026-11-19/xsec.py (kopia oryginalu: skrypty/xsec_orig_kopia_2026-11-19.py; oryginal i jego wyniki NIE sa zmieniane).
Wejscia (tylko odczyt): out/kappa_A.json (fity D0 tego zadania, skrypty/kappa_fit.py A), XSEC_2026-11-19/out/xsec_results.json (D+: uzysk, model tla, domkniecie;
wiersze nieprompt / odbicie / wybor zdarzen; Pythia; SF D0 z tagu D*), XSEC_atlas_2026-11-20/out/{truth_ngen,dstar,gen_model,fonll}.json, SFABS_wynik_2026-11-20/out/final.json,
GEN_mcstat/runs/lxplus_v20/truthD_mb.parquet (ulamek |eta| < 2.1), LUMI (przez xsec_results).
  python -B xsec_final.py   ->  out/final.json, tabele/final.md"""
import json
import math
import os
import sys

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
V1 = os.path.dirname(HERE)                                                    # XSEC_final_2026-11-20 (wejscia kappa_*.json, tylko odczyt)
CW = os.path.dirname(V1)
J = lambda *p: json.load(open(os.path.join(CW, *p)))  # noqa: E731
BINS = ["3-5", "5-8", "8-30"]
EDG = {"3-5": (3, 5), "5-8": (5, 8), "8-30": (8, 30)}
SPS = ["D0", "Dplus"]
NM = {"D0": "D⁰", "Dplus": "D⁺"}
q = lambda *x: math.sqrt(sum(float(v) ** 2 for v in x))  # noqa: E731

# ---------------------------------------------------------------- stale (zrodla w DEFINICJA.md)
BR = {"D0": (0.03936, 0.00030), "Dplus": (0.0931, 0.0014)}                     # PDG 2026; D0 bez DCS
TRK_PER_TRACK, NTRK = 0.026, {"D0": 2, "Dplus": 3}                             # 1512.02913 (XSEC_atlas C)
IP_MAT = 0.015                                                                # 2302.00336 (XSEC_atlas C)
VTX3 = {"D0": 0.0, "Dplus": 0.03}
MULT_LIT = 0.020                                                              # krotnosc (lit.), XSEC_atlas C
KAPPA_NORM = 0.058
HESSE_SCALE = 1.3
# korelacje D+/D0: czesc nieskorelowana wierszy wspolnych
RATIO_UNC = {"trk": 0.026, "ip": 0.010, "vtx": 0.030, "gen": 0.015}
FONLL136 = {"D0": {"3-5": (508.4, 218.5, 915.5), "5-8": (181.2, 113.4, 282.2), "8-30": (64.61, 49.23, 87.66)},      # ub, czastka, f = 1; XSEC2 (xsec2/fonll)
            "Dplus": {"3-5": (513.6, 216.6, 930.9), "5-8": (186.0, 114.9, 292.0), "8-30": (67.28, 50.94, 91.92)}}
K_E = {"3-5": 1.0285, "5-8": 1.036, "8-30": 1.0455}                            # FONLL sigma(13,6)/sigma(13), prompt D, |eta| < 2,2 (XSEC2)
FF = {"D0": 0.609, "Dplus": 0.240, "Dst": 0.243}                              # uzyte przez CMS (2107.01476 §6)
# CMS 2107.01476, tabela dsigma/dpT [ub/GeV], |eta| < 2,1, prompt, czastka + antyczastka (odczyt z arXiv HTML v2, 20.11.2026): (lo, hi, val, stat, syst)
CMS_T = {"Dst": [(4, 5, 166, 7, 24), (5, 6, 96, 6, 11), (6, 7, 47, 3, 5), (7, 8, 25.6, 1.7, 2.9), (8, 12, 8.8, 0.2, 1.0), (12, 16, 1.70, 0.09, 0.20), (16, 24, 0.352, 0.015, 0.042), (24, 40, 0.041, 0.003, 0.006)],
         "D0": [(4, 5, 430, 58, 38), (5, 6, 230, 17, 21), (6, 7, 136, 11, 12), (7, 8, 66, 8, 6), (8, 12, 21.0, 0.7, 1.9), (12, 16, 3.93, 0.19, 0.44), (16, 24, 0.81, 0.05, 0.09), (24, 40, 0.097, 0.009, 0.013)],
         "Dplus": [(4, 5, 250, 61, 34), (5, 6, 79, 11, 11), (6, 7, 50, 2, 7), (7, 8, 29.5, 1.4, 4.3), (8, 12, 9.4, 0.3, 1.2), (12, 16, 2.05, 0.07, 0.24), (16, 24, 0.406, 0.017, 0.051), (24, 40, 0.060, 0.006, 0.010)]}
F_30_40 = 0.32      # udzial 30-40 GeV w binie CMS 24-40 przy dsigma/dpT ~ pT^-5 (wlasne): (30^-4 - 40^-4)/(24^-4 - 40^-4)
# ATLAS 7 TeV 1512.02913 tab. 4, |eta| < 2,1, c + b, suma ladunkow: D+- dsigma/dpT [ub/GeV] (notatka w bazie wiedzy), D*+- calki z XSEC_atlas
ATLAS7_DP = [(3.5, 5, 127, 13, 12), (5, 6.5, 51.9, 4.3, 4.2), (6.5, 8, 20.0, 2.3, 1.6), (8, 12, 6.29, 0.56, 0.51), (12, 20, 0.583, 0.088, 0.050), (20, 30, 0.0736, 0.0, 0.0)]
ATLAS7_LUMI_BR = q(0.035, 0.021)
ATLAS7_DST = {"5-8": (96.3, 10.5), "8-30": (32.2, 3.4)}
# stosunki odniesienia
ALICE_R = {"3-5": (0.458, 0.029), "5-8": (0.470, 0.029), "8-30": (0.430, 0.027)}   # 2308.04877, |y| < 0,5, prompt; stat (+) ~6 % syst (XSEC WYNIKI §2)
LHCB_R = (0.405, 0.017)
EE_R = (0.395, 0.015)
# ALICE 13 TeV 2308.04877, prompt, srednia czastka/antyczastka, |y| < 0,5: dsigma/dy [ub] w 3-5 / 5-8 / 8-24 GeV (XSEC WYNIKI §5: sumy z HEPData rys. 4, NIE weryfikowane tutaj); blad przyjety 8 %
ALICE_XS = {"D0": {"3-5": 119.3, "5-8": 36.1, "8-30": 11.74}, "Dplus": {"3-5": 54.6, "5-8": 17.0, "8-30": 5.04}}
ALICE_XS_REL = 0.08
DST_TRUTH = {"5-8": 0.416, "8-30": 0.437}                                      # prawda mb_nd, XSEC_atlas RAPORT §A
DST_EE = 0.399
L_NB_ROUND, L_ERR_NB = 0.4648, 0.0103                                           # LUMI_2026-11-19/lumi.json (blad_wzgl 0,0222 = kalibracja 2,2 % (+) live 0,3 %), DEFINICJA §1
SF_DP35_ERR = 0.08                                                              # DECYZJA 20.11 (AUDYT_noty pkt 5)
CUT_SF = {"Dplus": (0.97, 0.05, 0.05), "D0": (0.95, 0.06, 0.08)}             # sprawdzian cieciowy XSEC_atlas (1/R, 5-30 GeV)


def cms_sum(sp, lo, hi):
    v = s2 = sy = 0.0
    for a, b, x, st, ss in CMS_T[sp]:
        if a >= lo and b <= hi:
            w = b - a
            v += x * w; s2 += (st * w) ** 2; sy += ss * w
    return v, math.sqrt(s2), sy


def eta_frac():
    """Ulamek prawdziwych D (wszystkie rozpady) z |eta| < 2,1 wsrod |eta| < 2,2, per gatunek i bin (prawda mb_nd)."""
    import pandas as pd
    import numpy as np
    T = pd.read_parquet(os.path.join(CW, "GEN_mcstat", "runs", "lxplus_v20", "truthD_mb.parquet"), columns=["pdg", "pt", "eta", "fromB"])
    o = {}
    for sp, pdg, m in (("D0", 421, 1.86484), ("Dplus", 411, 1.86966)):
        s = T[T.pdg.abs() == pdg]
        s = s.assign(y=np.arcsinh(s.pt * np.sinh(s.eta) / np.sqrt(s.pt ** 2 + m ** 2)))
        for bn, (lo, hi) in EDG.items():
            x = s[(s.pt >= lo) & (s.pt < hi)]
            a, b = x[x.eta.abs() < 2.2], x[x.eta.abs() < 2.1]
            hiA = 24 if bn == "8-30" else hi                                   # ALICE: ostatni bin 8-24 GeV
            ya = s[(s.pt >= lo) & (s.pt < hiA) & (s.y.abs() < 0.5) & (s.fromB == 0)]
            o[f"{sp}|{bn}"] = {"f_eta": len(b) / len(a), "n22": len(a), "f_np_22": float(a.fromB.mean()), "f_np_21": float(b.fromB.mean()),
                               "f_y05_prompt_over_eta22_prompt": len(ya) / int((a.fromB == 0).sum()), "f_y05_prompt_over_eta21_prompt": len(ya) / int((b.fromB == 0).sum())}
    return o


def main():
    KA = json.load(open(os.path.join(V1, "out", "kappa_A.json")))
    for bn_, v_ in json.load(open(os.path.join(V1, "out", "kappa_altmin.json"))).items():
        KA[bn_]["altmin"] = v_["altmin"]
    XS = J("XSEC_2026-11-19", "out", "xsec_results.json")
    NG = J("XSEC_atlas_2026-11-20", "out", "truth_ngen.json")["bins"]
    GM = J("XSEC_atlas_2026-11-20", "out", "gen_model.json")
    SFA = J("SFABS_wynik_2026-11-20", "out", "final.json")
    DS = J("XSEC_atlas_2026-11-20", "out", "dstar.json")
    FO = J("XSEC_atlas_2026-11-20", "out", "fonll.json")
    Lum, Lrel_v1 = XS["lumi"]["427929"]
    Lrel = L_ERR_NB / L_NB_ROUND                                                   # v2 (4): 2,216 %; v1: 2,2 % kalibracja (+) 0,3 % live (+) 0,3 % kompletnosc plikow = 2,24 %
    DPY = {c_["bin"]: c_ for c_ in J("FITBIAS_Dplus_2026-11-20", "out", "dplus_yields.json")["cells"] if c_["run"] == "427929"}
    AUD = J("AUDYT_noty_2026-11-20", "out", "przeliczenia.json")["6_dstar_spread"]
    EF = eta_frac()
    R = {"version": "v2", "L_nb": Lum, "L_rel": Lrel, "L_rel_v1": Lrel_v1, "cells": {}, "eta_frac": EF}
    CVER_EXCL = {"D0|3-5": 15223 + 1143, "D0|5-8": 5942 + 452, "D0|8-30": 1965 + 178, "Dplus|3-5": 18962 + 1179, "Dplus|5-8": 7486 + 437, "Dplus|8-30": 2562 + 172}  # CVER RAPORT §2

    for sp in SPS:
        for bn in BINS:
            x = XS["cells"][f"{sp}|427929|{bn}"]
            mc = XS["mc"][sp]["bins"][bn]
            ng = NG[f"{sp}|{bn}"]["N_excl"]
            assert ng == CVER_EXCL[f"{sp}|{bn}"], "mianowniki CVER i XSEC_atlas rozne"
            nreco = x["N_reco_truth_mb"]
            C = nreco / ng
            rel = {}
            c = {"N_gen_excl": ng, "N_gen_decayType": x["N_gen_sig_mb"], "N_reco_truth": nreco, "C": C, "C_XSEC": x["C"], "BR": BR[sp][0]}
            if sp == "D0":
                k = KA[bn]
                nom = k["M8|dir|1.0"]
                N = nom["N"]
                c.update({"N": N, "hesse": nom["err"], "boot": k["boot"]["rms"], "p_chi2": nom["p"], "sc": nom["sc"], "dmu": nom["dmu"], "valid": nom["valid"], "at_limit": nom["at_limit"],
                          "N_FITBIAS": k["M8|ref|1.364"]["N"], "N_ref_r1": k["M8|ref|1.0"]["N"], "N_dir_rFB": k["M8|dir|1.364"]["N"], "N_M4": k["M4|dir|1.0"]["N"], "N_M9": k["M9|dir|1.0"]["N"],
                          "N_M0": k["M0|ref|1.0"]["N"], "closure_b": k["closure|M8|dir"]["b"], "closure_b_err": k["closure|M8|dir"]["b_err"], "N_altmin": (k.get("altmin") or {}).get("N")})
                rel["stat"] = max(k["boot"]["rms"], HESSE_SCALE * nom["err"]) / N
                rel["fit_model"] = max(abs(k["M4|dir|1.0"]["N"] / N - 1), abs(k["M9|dir|1.0"]["N"] / N - 1))
                rel["fit_bias"] = q(k["closure|M8|dir"]["b"], k["closure|M8|dir"]["b_err"])
                rel["peak_tmpl"] = abs(k["M8|ref|1.0"]["N"] / N - 1)
                rel["peak_norm"] = 0.5 * abs(k["M8|dir|k+5.8"]["N"] - k["M8|dir|k-5.8"]["N"]) / N
                c["peak_norm20"] = 0.5 * abs(k["M8|dir|k+20"]["N"] - k["M8|dir|k-20"]["N"]) / N
                rel["fit_altmin"] = abs(k["altmin"]["N"] / N - 1) if k.get("altmin") else 0.0     # POST HOC (odstepstwo, RAPORT §1)
            else:
                dy = DPY[bn]                                                           # v2 (1): model A2 (FITBIAS_Dplus); zastepuje M0 z XSEC (N, stat, model, domkniecie) — bez dublowania
                assert dy["model"] == "A2"
                N = dy["N"]
                c.update({"N": N, "hesse": dy["sigma_hesse"], "boot": dy["sigma_boot"], "closure_b": dy["closure_b"], "closure_b_err": dy["closure_b_err"], "p_chi2": dy["p"], "sc": dy["sc"],
                          "N_M0_v1": x["N"], "N_A0": dy["N_A0"], "N_A4": dy["N_A4"], "fit_src": "FITBIAS_Dplus_2026-11-20 A2"})
                assert abs(dy["sigma_stat"] - max(HESSE_SCALE * dy["sigma_hesse"], dy["sigma_boot"])) < 0.2
                rel["stat"] = dy["sigma_stat"] / N                                     # max(1,3 x HESSE, bootstrap) — jak D0
                rel["fit_model"] = dy["sigma_model"] / N                               # max(|N - N_A0|, |N - N_A4|)
                rel["fit_bias"] = q(dy["closure_b"], dy["closure_b_err"])              # |b| (+) sigma_b, bez poprawki — ta sama definicja co D0 (v1 i v2)
                assert abs(rel["fit_bias"] - dy["fit_bias_rel"]) < 6e-4
            # SF
            if sp == "D0" and bn != "3-5":
                SF, sfe0 = x["SF"], x["SF"] * x["rel"]["SF"]
                spread = AUD[f"clf|{bn}"]["half_range"]                               # v2 (3): polowa rozpietosci 4 estymatorow Dm (AUDYT_noty pkt 6 b)
                sfe, src = q(sfe0, spread), "tag D* 427929 (dnet4) ⊕ rozrzut estymatorów Δm"
                c.update({"SF_err_v1": sfe0, "SF_spread": spread})
            else:
                d = SFA["D"][f"{sp}|{bn}"]
                SF, sfe, src = d["SF"], d["err"], "SFABS łączny"
                if sp == "Dplus" and bn == "3-5":                                      # v2 (2): DECYZJA ± 0,08 (AUDYT_noty pkt 5)
                    c["SF_err_v1"] = sfe
                    sfe, src = SF_DP35_ERR, "SFABS łączny, błąd ± 0,08 (decyzja)"
            c.update({"SF": SF, "SF_err": sfe, "SF_src": src})
            rel["SF"] = sfe / SF
            rel["mc_stat"] = math.sqrt((1 - C) / nreco)
            rel["BR"] = BR[sp][1] / BR[sp][0]
            rel["refl"] = x["rel"]["refl"]
            rel["nonprompt"] = x["rel"]["nonprompt"]
            rel["evsel"] = x["rel"]["evsel"]
            rel["tracking"] = TRK_PER_TRACK * NTRK[sp]
            rel["ip_material"] = IP_MAT
            rel["vertex3"] = VTX3[sp]
            rel["generator"] = q(GM[f"{sp}|WP70|{bn}"]["tot"], MULT_LIT)
            syst = q(*[v for k_, v in rel.items() if k_ != "stat"])
            sig = N / (C * SF * BR[sp][0] * Lum) / 1e3          # ub
            c.update({"rel": rel, "syst_rel": syst, "sigma_ub": sig, "stat_ub": sig * rel["stat"], "syst_ub": sig * syst, "lumi_ub": sig * Lrel, "tot_rel": q(rel["stat"], syst, Lrel),
                      "dsdpt_ub_GeV": sig / (EDG[bn][1] - EDG[bn][0]), "f_np": mc["f_np_all"], "pythia_ub": mc["sigma_pythia_inel_ub"],
                      "sigma_XSEC_ub": x["sigma_nb"] / 1e3, "syst_trk_run2": q(*[v for k_, v in rel.items() if k_ not in ("stat", "tracking")], 0.0135 * NTRK[sp])})
            if sp == "D0":
                c["sigma_closure_corr_info"] = sig / (1 + c["closure_b"])
            if sp == "Dplus":
                base = [v for k_, v in rel.items() if k_ not in ("stat", "SF")]
                c["variants"] = {"SF1_pm10": {"sigma_ub": sig * SF, "syst_rel": q(*base, 0.10)}, "SFABS_pm008": {"sigma_ub": sig, "syst_rel": q(*base, 0.08 / SF)},
                                 "SFABS_err": {"sigma_ub": sig, "syst_rel": q(*base, SFA["D"][f"{sp}|{bn}"]["err"] / SF)}}
            R["cells"][f"{sp}|{bn}"] = c

    # ---------------------------------------------------------------- D+/D0
    cells = SFA["SF_ref"]["cells"]
    wsum = sum(1 / q(u["stat"], u["model"]) ** 2 for u in cells)
    R["ratio"] = {}
    for bn in BINS:
        a, b = R["cells"][f"Dplus|{bn}"], R["cells"][f"D0|{bn}"]
        ra, rb = a["rel"], b["rel"]
        r = a["sigma_ub"] / b["sigma_ub"]
        rel = {"stat": q(ra["stat"], rb["stat"]), "fit_model": q(ra["fit_model"], rb["fit_model"]), "fit_bias": q(ra["fit_bias"], rb["fit_bias"]),
               "peak": q(rb["peak_tmpl"], rb["peak_norm"], rb["fit_altmin"]), "mc_stat": q(ra["mc_stat"], rb["mc_stat"]), "BR": q(ra["BR"], rb["BR"]), "refl": rb["refl"],
               "nonprompt": q(ra["nonprompt"], rb["nonprompt"]), "tracking": RATIO_UNC["trk"], "ip_material": RATIO_UNC["ip"], "vertex3": RATIO_UNC["vtx"], "generator": RATIO_UNC["gen"]}
        if bn == "3-5":
            extra = math.sqrt(max(ra["SF"] ** 2 - (a["SF_err_v1"] / a["SF"]) ** 2, 0))   # v2 (2): nadmiar bledu SF+ (0,08 wobec 0,048) dodany kwadratowo (jak AUDYT_noty pkt 5)
            rel["SF"] = q(SFA["ratio_3-5"]["sf_ratio_rel"], extra); sfnote = "SF⁺/SF⁰ = 0,981 ± 0,042 (SFABS; wspólny SF_ref skrócony) ⊕ nadmiar błędu SF⁺ (± 0,08)"
        else:
            cell = [u for u in cells if u["k"] == f"427929|{bn}"][0]
            w = (1 / q(cell["stat"], cell["model"]) ** 2) / wsum                       # waga komorki 427929 w SF_ref
            wC = SFA["D"][f"Dplus|{bn}"]["wC"]                                          # waga metody C (przeniesienie) w SF+ lacznym
            cov = wC * w * (b["SF_err_v1"]) ** 2 if bn == "5-8" else 0.0                   # 5-8: ta sama liczba dnet4; 8-30: w SF_ref fit 2D (inny estymator) -> korelacje pomijam
            rel["SF"] = math.sqrt(max(ra["SF"] ** 2 + rb["SF"] ** 2 - 2 * cov / (a["SF"] * b["SF"]), 0))
            sfnote = f"SF⁰ tag D*, SF⁺ SFABS; waga komórki 427929 w SF_ref {w:.3f}"
        syst = q(*[v for k_, v in rel.items() if k_ != "stat"])
        py = XS["ratio_DpD0"][f"427929|{bn}"]["pythia_A3_mbnd"]
        R["ratio"][bn] = {"R": r, "rel": rel, "stat": r * rel["stat"], "syst": r * syst, "tot": r * q(rel["stat"], syst), "syst_noSF": r * q(*[v for k_, v in rel.items() if k_ not in ("stat", "SF")]),
                          "pythia": py, "R_XSEC": XS["ratio_DpD0"][f"427929|{bn}"]["R_sf_equal"], "sf_note": sfnote}
        R["ratio"][bn]["syst_rel"] = syst

    # ---------------------------------------------------------------- D*/D0
    R["dstar"] = {}
    DST_SYST = {"5-8": 0.172, "8-30": 0.113}                                           # XSEC_atlas tabele/wyniki.md A (razem syst)
    for bn in ("5-8", "8-30"):
        d = DS[f"nom|{bn}"]
        b = R["cells"][f"D0|{bn}"]
        sd = d["sigma_ub"]; std = sd * d["N_err"] / d["N"]
        syd_full = sd * DST_SYST[bn]
        syd = q(math.sqrt(max(DST_SYST[bn] ** 2 - 0.077 ** 2 - 0.015 ** 2, 0)), 0.026)  # bez sledzenia 3 torow i IP; zostaje 1 tor
        sy0 = q(*[v for k_, v in b["rel"].items() if k_ not in ("stat", "tracking", "ip_material", "evsel")])
        r = sd / b["sigma_ub"]
        R["dstar"][bn] = {"sigma_ub": sd, "stat_ub": std, "syst_ub": syd_full, "lumi_ub": sd * Lrel, "R": r, "stat": r * q(std / sd, b["rel"]["stat"]), "syst": r * q(syd, sy0),
                          "tot": r * q(std / sd, b["rel"]["stat"], syd, sy0), "truth_mb": DST_TRUTH[bn], "ee": DST_EE}

    # ---------------------------------------------------------------- porownania
    def sig_dev(x, ex, y, ey):
        return (x - y) / q(ex, ey)
    CMP = {}
    for sp in SPS:
        for bn in BINS:
            c = R["cells"][f"{sp}|{bn}"]
            s, et = c["sigma_ub"], c["sigma_ub"] * c["tot_rel"]
            o = {}
            # Pythia A3
            o["pythia"] = {"ref": c["pythia_ub"], "ratio": s / c["pythia_ub"], "ratio_err": et / c["pythia_ub"], "Z": (s - c["pythia_ub"]) / et}
            # FONLL prompt
            fo = [2 * FF[sp] * v for v in FONLL136[sp][bn]]
            pr = s * (1 - c["f_np"]); epr = pr * q(c["tot_rel"], 0.5 * c["f_np"] / (1 - c["f_np"]))
            zb = 0.0 if fo[1] <= pr <= fo[2] else ((pr - fo[2]) / epr if pr > fo[2] else (pr - fo[1]) / epr)
            o["fonll"] = {"prompt_data": pr, "err": epr, "centr": fo[0], "min": fo[1], "max": fo[2], "ratio_centr": pr / fo[0], "ratio_err": epr / fo[0], "Z_band": zb, "ratio_to_max": pr / fo[2]}
            # CMS
            if bn != "3-5":
                lo, hi = (5, 8) if bn == "5-8" else (8, 40)
                v, st, sy = cms_sum(sp, lo, hi)
                if bn == "8-30":
                    v30 = v - F_30_40 * CMS_T[sp][-1][2] * 16
                else:
                    v30 = v
                ef = EF[f"{sp}|{bn}"]
                fac = {"f_eta": ef["f_eta"], "f_eta_flat": 2.1 / 2.2, "one_minus_fnp": 1 - c["f_np"], "inv_KE": 1 / K_E[bn]}
                conv = s * ef["f_eta"] * (1 - c["f_np"]) / K_E[bn]
                ec = conv * q(c["tot_rel"], 0.5 * c["f_np"] / (1 - c["f_np"]))
                ecms = q(st, sy)
                o["cms"] = {"cms": v, "cms_stat": st, "cms_syst": sy, "cms_8_30": v30, "factors": fac, "data_conv": conv, "data_conv_err": ec, "ratio": conv / v30, "ratio_err": (conv / v30) * q(ec / conv, ecms / v),
                            "Z": sig_dev(conv, ec, v30, ecms), "ratio_flat_eta": s * (2.1 / 2.2) * (1 - c["f_np"]) / K_E[bn] / v30}
            # ALICE (POST HOC, sprawdzian czynnika 2 i poziomu CMS; ksztalt y z Pythii A3): dane prompt -> |y| < 0,5, srednia czastka/antyczastka, 13 TeV
            ef = EF[f"{sp}|{bn}"]
            al = ALICE_XS[sp][bn]
            cv = pr * ef["f_y05_prompt_over_eta22_prompt"] / 2 / K_E[bn]
            o["alice"] = {"alice": al, "alice_err": al * ALICE_XS_REL, "f_y": ef["f_y05_prompt_over_eta22_prompt"], "data_conv": cv, "data_conv_err": cv * epr / pr, "ratio": cv / al,
                          "ratio_err": cv / al * q(epr / pr, ALICE_XS_REL), "Z": sig_dev(cv, cv * epr / pr, al, al * ALICE_XS_REL)}
            if "cms" in o:
                cc = o["cms"]["cms_8_30"] * ef["f_y05_prompt_over_eta21_prompt"] / 2
                o["alice"]["cms_conv"] = cc; o["alice"]["cms_over_alice"] = cc / al
                o["alice"]["cms_over_alice_err"] = cc / al * q(q(o["cms"]["cms_stat"], o["cms"]["cms_syst"]) / o["cms"]["cms"], ALICE_XS_REL)
            CMP[f"{sp}|{bn}"] = o
    # ATLAS 7 TeV D+ (skalowanie FONLL, sam prompt) i D*
    I = FO["integrals_ub"]
    for bn in ("5-8", "8-30"):
        lo, hi = EDG[bn]
        v = s2 = sy = 0.0
        for a, b, x, st, ss in ATLAS7_DP:
            if a >= lo and b <= hi:
                v += x * (b - a); s2 += (st * (b - a)) ** 2; sy += ss * (b - a)
        e7 = q(math.sqrt(s2), sy, v * ATLAS7_LUMI_BR)
        c19, c7 = I[f"prompt_Dp_c19_eta2.2_pt{bn}"], I[f"prompt_Dp_c7_eta2.1_pt{bn}"]
        K, Kmin = c19[0] / c7[0], c19[1] / c7[1]
        c = R["cells"][f"Dplus|{bn}"]
        s, et = c["sigma_ub"], c["sigma_ub"] * c["tot_rel"]
        ref, eref = v * K, v * K * q(e7 / v, abs(K - Kmin) / K)
        CMP[f"Dplus|{bn}"]["atlas7"] = {"atlas7": v, "atlas7_err": e7, "K": K, "K_min": Kmin, "scaled": ref, "scaled_err": eref, "ratio": s / ref, "ratio_err": (s / ref) * q(et / s, eref / ref), "Z": sig_dev(s, et, ref, eref)}
        kd = FO["K"][bn]
        d = R["dstar"][bn]
        a7, ea7 = ATLAS7_DST[bn]
        ref, eref = a7 * kd["K"], a7 * kd["K"] * q(ea7 / a7, abs(kd["K"] - kd["K_min"]) / kd["K"])
        ed = q(d["stat_ub"], d["syst_ub"], d["lumi_ub"])
        d["atlas7"] = {"scaled": ref, "scaled_err": eref, "ratio": d["sigma_ub"] / ref, "ratio_err": (d["sigma_ub"] / ref) * q(ed / d["sigma_ub"], eref / ref), "Z": sig_dev(d["sigma_ub"], ed, ref, eref)}
        # D* vs CMS (ta sama konwersja co D0)
        lo2, hi2 = (5, 8) if bn == "5-8" else (8, 40)
        v, st, sy = cms_sum("Dst", lo2, hi2)
        v30 = v - (F_30_40 * CMS_T["Dst"][-1][2] * 16 if bn == "8-30" else 0)
        c0 = R["cells"][f"D0|{bn}"]
        conv = d["sigma_ub"] * EF[f"D0|{bn}"]["f_eta"] * (1 - c0["f_np"]) / K_E[bn]
        ec = conv * q(ed / d["sigma_ub"], 0.5 * c0["f_np"] / (1 - c0["f_np"]))
        d["cms"] = {"cms": v30, "cms_err": q(st, sy), "data_conv": conv, "data_conv_err": ec, "ratio": conv / v30, "ratio_err": conv / v30 * q(ec / conv, q(st, sy) / v), "Z": sig_dev(conv, ec, v30, q(st, sy))}
        # Pythia D*
        py = {"5-8": 200.0, "8-30": 69.8}[bn]                                           # XSEC_atlas tabele/wyniki.md A
        d["pythia"] = {"ref": py, "ratio": d["sigma_ub"] / py, "ratio_err": ed / py, "Z": (d["sigma_ub"] - py) / ed}
        # CMS stosunki (wlasna arytmetyka z tabeli; syst: suma kwadratowa calkowitych syst minus wspolna swietlnosc 2,5 %)
        v0, st0, sy0 = cms_sum("D0", lo2, hi2); vp, stp, syp = cms_sum("Dplus", lo2, hi2)
        un = lambda sy_, v_: math.sqrt(max((sy_ / v_) ** 2 - 0.025 ** 2, 0))  # noqa: E731
        R["ratio"][bn]["cms"] = (vp / v0, vp / v0 * q(stp / vp, st0 / v0, un(syp, vp), un(sy0, v0)))
        d["cms_ratio"] = (v / v0, v / v0 * q(st / v, st0 / v0, un(sy, v), un(sy0, v0)))
    for bn in BINS:
        r = R["ratio"][bn]
        refs = {"ALICE": ALICE_R[bn], "LHCb": LHCB_R, "e+e-": EE_R, "Pythia A3": (r["pythia"], 0.0)}
        if "cms" in r:
            refs["CMS"] = r["cms"]
        r["cmp"] = {k: {"ref": v[0], "ref_err": v[1], "Z": sig_dev(r["R"], r["tot"], v[0], v[1])} for k, v in refs.items()}
    for bn in ("5-8", "8-30"):
        d = R["dstar"][bn]
        refs = {"prawda mb_nd": (d["truth_mb"], 0.0), "e+e-": (DST_EE, 0.0), "CMS": d["cms_ratio"]}
        d["cmp"] = {k: {"ref": v[0], "ref_err": v[1], "Z": sig_dev(d["R"], d["tot"], v[0], v[1])} for k, v in refs.items()}
    R["cmp"] = CMP
    R["cut_check_SF"] = CUT_SF
    os.makedirs(os.path.join(HERE, "out"), exist_ok=True)
    json.dump(R, open(os.path.join(HERE, "out", "final.json"), "w"), indent=1, default=float, ensure_ascii=False)
    write_md(R)


def pl(x, n=1):
    return f"{x:.{n}f}".replace(".", ",")


def sgn(x, n=1):
    return f"{x:+.{n}f}".replace(".", ",").replace("-", "−")


def nd(x):
    return 0 if x >= 10 else 1


def write_md(R):
    T = ["# XSEC_final v2 — tabele końcowe, run 427929 (generowane przez `v2/xsec_final.py`; ZASTĘPUJĄ `../tabele/final.md` (v1) — zmiany: `ZMIANY.md`; jak v1 ZASTĘPUJĄ `XSEC_2026-11-19/out/tabele.md` T1/T1b/T2, CVER §5, SFABS „Przekrój”, XSEC_atlas §C)", "",
         f"L = {pl(R['L_nb'], 4)} nb⁻¹ ± {pl(100 * R['L_rel'], 2)} % (0,0103/0,4648; v1: {pl(100 * R['L_rel_v1'], 2)} % = to samo ⊕ kompletność plików 0,3 %). D⁺: model fitu A2 (FITBIAS_Dplus). SF D⁺ 3–5 GeV ± 0,08; SF D⁰ ≥ 5 GeV ⊕ rozrzut estymatorów Δm. D + anty-D, prompt + nieprompt, |η(D)| < 2,2. σ = N / (C_ekskl × SF × BR_PDG × L). POST HOC; PRELIMINARY.", "",
         "## F1. Przekroje", "", "| gatunek | bin [GeV] | N | C_ekskl [%] | SF | **σ [µb] ± stat ± syst ± L** | dσ/dpT [µb/GeV] | razem [%] | było: XSEC 19.11 [µb] | σ_final/σ_XSEC |", "|---|---|---|---|---|---|---|---|---|---|"]
    for k, c in R["cells"].items():
        sp, bn = k.split("|")
        s = c["sigma_ub"]; n = nd(c["stat_ub"])
        T.append(f"| {NM[sp]} | {bn.replace('-', '–')} | {c['N']:.0f} ± {c['N'] * c['rel']['stat']:.0f} | {pl(100 * c['C'], 2)} | {pl(c['SF'], 3)} ± {pl(c['SF_err'], 3)} | **{pl(s, n)} ± {pl(c['stat_ub'], n)} ± {pl(c['syst_ub'], n)} ± {pl(c['lumi_ub'], n)}** | "
                 f"{pl(c['dsdpt_ub_GeV'], 1 if c['dsdpt_ub_GeV'] > 10 else 2)} | {pl(100 * c['tot_rel'])} | {pl(c['sigma_XSEC_ub'], nd(c['sigma_XSEC_ub'] / 10))} | {pl(s / c['sigma_XSEC_ub'], 3)} |")
    T += ["", "D⁺, warianty SF (nominał v2 = SFABS łączny; 3–5 GeV z błędem ± 0,08):", "", "| bin | nominał v2 | SF = 1 ± 0,10 (dawny scenariusz) | SFABS ± 0,08 | SFABS z własnym błędem (± 0,05; nominał v1) |", "|---|---|---|---|---|"]
    for bn in BINS:
        c = R["cells"][f"Dplus|{bn}"]; v = c["variants"]; n = nd(c["stat_ub"])
        f = lambda s, sy: f"{pl(s, n)} ± {pl(s * c['rel']['stat'], n)} ± {pl(s * sy, n)} ± {pl(s * R['L_rel'], n)}"  # noqa: E731
        T.append(f"| {bn.replace('-', '–')} | {f(c['sigma_ub'], c['syst_rel'])} | {f(v['SF1_pm10']['sigma_ub'], v['SF1_pm10']['syst_rel'])} | {f(v['SFABS_pm008']['sigma_ub'], v['SFABS_pm008']['syst_rel'])} | {f(v['SFABS_err']['sigma_ub'], v['SFABS_err']['syst_rel'])} |")
    cs = R["cut_check_SF"]
    T += ["", f"Niezależny sprawdzian SF (selekcja cięciowa XSEC_atlas, 1/R, 5–30 GeV; nie wchodzi do nominału): D⁺ {pl(cs['Dplus'][0], 2)} ± {pl(cs['Dplus'][1], 2)} ± {pl(cs['Dplus'][2], 2)}; D⁰ {pl(cs['D0'][0], 2)} ± {pl(cs['D0'][1], 2)} ± {pl(cs['D0'][2], 2)}.",
          "", "## F2. Budżet [%]", ""]
    rows = [("stat", "statystyczny (D⁰ i D⁺: max(bootstrap, 1,3 × HESSE))"), ("fit_model", "fit: model (D⁰: M4, M9; D⁺: A0, A4 wobec A2)"), ("fit_bias", "fit: obciążenie (domknięcie mb_nd, |b| ⊕ σ_b)"),
            ("peak_tmpl", "szablon pikujący: wprost vs interpolowany"), ("peak_norm", "normalizacja κ ± 5,8 %"), ("fit_altmin", "drugie minimum skali szerokości (POST HOC)"), ("SF", "SF klasyfikatora (D⁰ ≥ 5 GeV: ⊕ rozrzut estymatorów Δm; D⁺ 3–5: ± 0,08)"),
            ("mc_stat", "stat. MC (C)"), ("BR", "BR"), ("refl", "odbicie K↔π"), ("nonprompt", "nieprompt"), ("evsel", "wybór zdarzeń + tryger"), ("tracking", "śledzenie (2,6 %/tor) [nowy]"),
            ("ip_material", "d₀/z₀ + materiał [nowy]"), ("vertex3", "wierzchołek 3 torów [nowy]"), ("generator", "model generatora [nowy]")]
    ks = list(R["cells"])
    T += ["| źródło | " + " | ".join(f"{NM[k.split('|')[0]]} {k.split('|')[1].replace('-', '–')}" for k in ks) + " |", "|---|" + "---|" * len(ks)]
    for key, name in rows:
        T.append(f"| {name} | " + " | ".join((pl(100 * R['cells'][k]['rel'][key]) if key in R['cells'][k]['rel'] else "—") for k in ks) + " |")
    T.append("| **razem syst** | " + " | ".join(f"**{pl(100 * R['cells'][k]['syst_rel'])}**" for k in ks) + " |")
    T.append("| świetlność | " + " | ".join(pl(100 * R["L_rel"], 2) for k in ks) + " |")
    T.append("| **razem** | " + " | ".join(f"**{pl(100 * R['cells'][k]['tot_rel'])}**" for k in ks) + " |")
    T.append("| wariant: śledzenie Run 2 (1,35 %/tor) — razem syst | " + " | ".join(pl(100 * R["cells"][k]["syst_trk_run2"]) for k in ks) + " |")
    T += ["", "Pominięty świadomie: wiersz dE/dx XSEC_atlas (0,1–0,5 %) — zawarty w SF mierzonym na danych.", "", "## F3. Stosunki", "",
          "| bin | **D⁺/D⁰ ± stat ± syst** | razem | było (XSEC, SF⁺ = SF⁰) | Pythia A3 | ALICE 13 TeV \\|y\\| < 0,5 | LHCb 13 TeV | e⁺e⁻ | CMS 13 TeV | odchylenia [σ]: Pythia / ALICE / LHCb / e⁺e⁻ / CMS |", "|---|---|---|---|---|---|---|---|---|---|"]
    for bn in BINS:
        r = R["ratio"][bn]; cm = r["cmp"]
        z = " / ".join(sgn(cm[k]["Z"]) if k in cm else "—" for k in ("Pythia A3", "ALICE", "LHCb", "e+e-", "CMS"))
        T.append(f"| {bn.replace('-', '–')} | **{pl(r['R'], 3)} ± {pl(r['stat'], 3)} ± {pl(r['syst'], 3)}** | {pl(100 * r['tot'] / r['R'])} % | {pl(r['R_XSEC'], 3)} | {pl(r['pythia'], 3)} | {pl(ALICE_R[bn][0], 3)} ± {pl(ALICE_R[bn][1], 3)} | "
                 f"{pl(LHCB_R[0], 3)} ± {pl(LHCB_R[1], 3)} | {pl(EE_R[0], 3)} ± {pl(EE_R[1], 3)} | " + (f"{pl(r['cms'][0], 3)} ± {pl(r['cms'][1], 3)}" if "cms" in r else "—") + f" | {z} |")
    T += ["", "Budżet D⁺/D⁰ [%]:", "", "| bin | stat | fit model | fit obciążenie | pikująca + minimum | SF⁺/SF⁰ | stat MC | BR | nieprompt | śledzenie (1 tor) | d₀/materiał | wierzchołek | generator |", "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for bn in BINS:
        e = R["ratio"][bn]["rel"]
        T.append(f"| {bn.replace('-', '–')} | " + " | ".join(pl(100 * e[k]) for k in ("stat", "fit_model", "fit_bias", "peak", "SF", "mc_stat", "BR", "nonprompt", "tracking", "ip_material", "vertex3", "generator")) + " |")
    T += ["", "| bin | σ(D*±) [µb] (XSEC_atlas) | **D*/D⁰ ± stat ± syst** | prawda mb_nd | e⁺e⁻ | CMS 13 TeV | odchylenia [σ]: mb_nd / e⁺e⁻ / CMS |", "|---|---|---|---|---|---|---|"]
    for bn in ("5-8", "8-30"):
        d = R["dstar"][bn]; cm = d["cmp"]
        T.append(f"| {bn.replace('-', '–')} | {pl(d['sigma_ub'], 1)} ± {pl(d['stat_ub'], 1)} ± {pl(d['syst_ub'], 1)} ± {pl(d['lumi_ub'], 1)} | **{pl(d['R'], 3)} ± {pl(d['stat'], 3)} ± {pl(d['syst'], 3)}** | {pl(d['truth_mb'], 3)} | {pl(d['ee'], 3)} | "
                 f"{pl(d['cms_ratio'][0], 3)} ± {pl(d['cms_ratio'][1], 3)} | {sgn(cm['prawda mb_nd']['Z'])} / {sgn(cm['e+e-']['Z'])} / {sgn(cm['CMS']['Z'])} |")
    T += ["", "## F4. Porównania przekrojów (błąd danych = stat ⊕ syst ⊕ L; Z w σ)", "",
          "| | bin | σ [µb] | Pythia A3 [µb] | dane/Pythia (Z) | FONLL prompt 2fσ: centr. (pasmo) [µb] | dane prompt [µb] | dane/FONLL centr. | Z do krawędzi pasma | CMS 13 TeV [µb] | dane w def. CMS [µb] | dane/CMS (Z) | ATLAS 7 TeV × K_FONLL [µb] | dane/ATLAS skal. (Z) |", "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for k, c in R["cells"].items():
        sp, bn = k.split("|"); o = R["cmp"][k]; n = nd(c["stat_ub"])
        p, f = o["pythia"], o["fonll"]
        cm = f"{pl(o['cms']['cms_8_30'], 1)} ± {pl(q(o['cms']['cms_stat'], o['cms']['cms_syst']), 1)} | {pl(o['cms']['data_conv'], n)} ± {pl(o['cms']['data_conv_err'], n)} | **{pl(o['cms']['ratio'], 2)} ± {pl(o['cms']['ratio_err'], 2)}** ({sgn(o['cms']['Z'])})" if "cms" in o else "— | — | —"
        a7 = f"{pl(o['atlas7']['scaled'], 1)} ± {pl(o['atlas7']['scaled_err'], 1)} (K = {pl(o['atlas7']['K'], 2)}) | {pl(o['atlas7']['ratio'], 2)} ± {pl(o['atlas7']['ratio_err'], 2)} ({sgn(o['atlas7']['Z'])})" if "atlas7" in o else "— | —"
        T.append(f"| {NM[sp]} | {bn.replace('-', '–')} | {pl(c['sigma_ub'], n)} ± {pl(c['sigma_ub'] * c['tot_rel'], n)} | {pl(p['ref'], 0)} | {pl(p['ratio'], 2)} ± {pl(p['ratio_err'], 2)} ({sgn(p['Z'])}) | {pl(f['centr'], 0)} ({pl(f['min'], 0)}–{pl(f['max'], 0)}) | "
                 f"{pl(f['prompt_data'], n)} ± {pl(f['err'], n)} | {pl(f['ratio_centr'], 2)} ± {pl(f['ratio_err'], 2)} | {sgn(f['Z_band'])} | {cm} | {a7} |")
    for bn in ("5-8", "8-30"):
        d = R["dstar"][bn]
        T.append(f"| D*± | {bn.replace('-', '–')} | {pl(d['sigma_ub'], 1)} ± {pl(q(d['stat_ub'], d['syst_ub'], d['lumi_ub']), 1)} | {pl(d['pythia']['ref'], 0)} | {pl(d['pythia']['ratio'], 2)} ± {pl(d['pythia']['ratio_err'], 2)} ({sgn(d['pythia']['Z'])}) | (XSEC_atlas) | | | | "
                 f"{pl(d['cms']['cms'], 1)} ± {pl(d['cms']['cms_err'], 1)} | {pl(d['cms']['data_conv'], 1)} ± {pl(d['cms']['data_conv_err'], 1)} | {pl(d['cms']['ratio'], 2)} ± {pl(d['cms']['ratio_err'], 2)} ({sgn(d['cms']['Z'])}) | "
                 f"{pl(d['atlas7']['scaled'], 1)} ± {pl(d['atlas7']['scaled_err'], 1)} | {pl(d['atlas7']['ratio'], 2)} ± {pl(d['atlas7']['ratio_err'], 2)} ({sgn(d['atlas7']['Z'])}) |")
    T += ["", "Czynniki przeliczenia do definicji CMS (prompt, |η| < 2,1, 13 TeV):", "", "| | bin | f_η = N(\\|η\\| < 2,1)/N(\\|η\\| < 2,2) (prawda mb_nd) | 2,1/2,2 | 1 − f_np (mb_nd) | 1/K_E | łącznie | dane/CMS przy płaskim η |", "|---|---|---|---|---|---|---|---|"]
    for k, o in R["cmp"].items():
        if "cms" in o:
            sp, bn = k.split("|"); f = o["cms"]["factors"]
            T.append(f"| {NM[sp]} | {bn.replace('-', '–')} | {pl(f['f_eta'], 4)} | {pl(f['f_eta_flat'], 4)} | {pl(f['one_minus_fnp'], 3)} | {pl(f['inv_KE'], 3)} | {pl(f['f_eta'] * f['one_minus_fnp'] * f['inv_KE'], 3)} | {pl(o['cms']['ratio_flat_eta'], 2)} |")
    T += ["", "Sprawdzian wobec ALICE 13 TeV (POST HOC; prompt, |y| < 0,5, średnia cząstka/antycząstka; przeliczenie kształtem rapidity Pythii A3: f_y = N_prompt(|y| < 0,5)/N_prompt(|η| < 2,2), ÷ 2, ÷ K_E; ALICE: sumy binów z XSEC WYNIKI §5, błąd przyjęty 8 %):", "",
          "| | bin | f_y (płasko: 0,227) | dane w def. ALICE dσ/dy [µb] | ALICE [µb] | dane/ALICE (Z) | CMS w def. ALICE [µb] | CMS/ALICE |", "|---|---|---|---|---|---|---|---|"]
    for k, o in R["cmp"].items():
        sp, bn = k.split("|"); a = o["alice"]
        T.append(f"| {NM[sp]} | {bn.replace('-', '–')}" + (" (ALICE 8–24)" if bn == "8-30" else "") + f" | {pl(a['f_y'], 3)} | {pl(a['data_conv'], 1)} ± {pl(a['data_conv_err'], 1)} | {pl(a['alice'], 1 if a['alice'] > 10 else 2)} | {pl(a['ratio'], 2)} ± {pl(a['ratio_err'], 2)} ({sgn(a['Z'])}) | "
                 + (f"{pl(a['cms_conv'], 1)} | {pl(a['cms_over_alice'], 2)} ± {pl(a['cms_over_alice_err'], 2)}" if "cms_conv" in a else "— | —") + " |")
    T += ["", "## F5. Uzyski D⁰ (WP70, M8) — warianty (bez zmian względem v1)", "", "| bin | FITBIAS (ref, κ × 1,364) | ref, r = 1 | dir, κ × 1,364 | **dir, r = 1 (nominał)** | HESSE | bootstrap | p(χ²) | M4 | M9 | κ ± 5,8 % | drugie minimum | domknięcie mb_nd b [%] |", "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for bn in BINS:
        c = R["cells"][f"D0|{bn}"]
        T.append(f"| {bn.replace('-', '–')} | {c['N_FITBIAS']:.0f} | {c['N_ref_r1']:.0f} | {c['N_dir_rFB']:.0f} | **{c['N']:.0f} ± {c['N'] * c['rel']['stat']:.0f}** | {c['hesse']:.0f} | {c['boot']:.0f} | {pl(c['p_chi2'], 3)} | {c['N_M4']:.0f} | {c['N_M9']:.0f} | "
                 f"± {pl(100 * c['rel']['peak_norm'], 2)} % | " + (f"{c['N_altmin']:.0f}" if c["N_altmin"] else "brak") + f" | {sgn(100 * c['closure_b'])} ± {pl(100 * c['closure_b_err'])} |")
    T += ["", "## F6. Uzyski D⁺ (WP70) — model A2 (FITBIAS_Dplus_2026-11-20)", "", "| bin | M0 = A0 (v1) | **A2 (nominał v2) ± stat** | 1,3 × HESSE | bootstrap | p(χ²) | sc | A4 | model [%] | domknięcie mb_nd b [%] | \\|b\\| ⊕ σ_b [%] |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for bn in BINS:
        c = R["cells"][f"Dplus|{bn}"]
        T.append(f"| {bn.replace('-', '–')} | {c['N_M0_v1']:.0f} | **{c['N']:.0f} ± {c['N'] * c['rel']['stat']:.0f}** | {HESSE_SCALE * c['hesse']:.0f} | {c['boot']:.0f} | {pl(c['p_chi2'], 3)} | {pl(c['sc'], 2)} | {c['N_A4']:.0f} | "
                 f"{pl(100 * c['rel']['fit_model'])} | {sgn(100 * c['closure_b'])} ± {pl(100 * c['closure_b_err'])} | {pl(100 * c['rel']['fit_bias'])} |")
    T += ["", "Żadne z 9 ustawień fitu D⁺ nie spełnia kryterium domknięcia w 5–8 GeV; A2 wybrane wg rankingu Q; POST HOC. Uzysk nie jest poprawiany o b."]
    os.makedirs(os.path.join(HERE, "tabele"), exist_ok=True)
    open(os.path.join(HERE, "tabele", "final.md"), "w").write("\n".join(T) + "\n")
    print("\n".join(T))


if __name__ == "__main__":
    main()
