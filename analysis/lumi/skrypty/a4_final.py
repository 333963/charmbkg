"""a4: liczby końcowe -> lumi.json, tabele/lumi_final_perLB_<run>.csv.
L_del: OflLumi-Run3-003 (oficjalny ilumicalc cvmfs 20231005; własny iLumiCalc daje tę samą sumę, tag 004 identyczny).
live nominalny: TAV/TAP L1_RD0_FILLED per LB (COOL LVL1COUNTERS); wariant: własna pozycja L1 łańcucha.
preskale: drzewo (= L1 x HLT COOL). kompletność plików: (N_drzewo/TAV pozycji odniesienia) / plateau, <= 1."""
import pandas as pd, numpy as np, json
T = "/home/k/Documents/charm/charm_work/TAGREG_2026-11-18/runs/lxplus"
TAG = "OflLumi-Run3-003"
SEED = {"HLT_mb_mbts_L1MBTS_1": "L1_MBTS_1", "HLT_noalg_mb_L1MBTS_1": "L1_MBTS_1", "HLT_noalg_L1AFP_A_OR_C": "L1_AFP_A_OR_C", "HLT_j20_L1AFP_A_OR_C": "L1_AFP_A_OR_C", "HLT_noalg_L1J12": "L1_J12"}
CFG = {427929: dict(ref="L1_MBTS_1", chains=["HLT_mb_mbts_L1MBTS_1", "HLT_noalg_mb_L1MBTS_1"]), 428770: dict(ref="L1_J12", chains=["HLT_noalg_L1AFP_A_OR_C", "HLT_j20_L1AFP_A_OR_C", "HLT_noalg_L1J12"])}
REL_CAL = 0.022      # ZAŁOŻENIE: niepewność kalibracji 2022 (wstępna, ATL-DAPR-PUB-2023-001) przeniesiona na runy niskiego mu
out = {"_meta": {"tag": TAG, "tag_uwaga": "003 = tag oficjalnych plików GRL AFP low-mu (cvmfs 20231005/20231024); 004 daje identyczne L_del; 001/002 różnią się o -3/+7 % (427929) i +5/+7 % (428770); 005-007 zwracają 0",
                 "jednostka": "nb^-1", "blad": "względny; kalibracja 2,2 % = ZAŁOŻENIE (wartość dla danych fizycznych 2022), live i kompletność osobno; konserwatywnie 7 % (rozrzut tagów 001-003)"}}
for run, c in CFG.items():
    d = pd.read_csv(f"tabele/lumi_perLB_{run}.csv", index_col=0)
    d["live_RD0"] = d.TAV_L1_RD0_FILLED / d.TAP_L1_RD0_FILLED if "TAV_L1_RD0_FILLED" in d else np.nan
    if d.live_RD0.isna().all():
        cn = pd.read_csv(f"surowe/lxplus/cool/l1counters_{run}.csv"); g = cn[cn.item == "L1_RD0_FILLED"].set_index("LB").reindex(d.index); d["live_RD0"] = g.TAV / g.TAP
    d["Lrec"] = d.Ldel * d.live_RD0
    plateau = float(d["compl"][d["compl"] >= d["compl"].quantile(0.25)].median())
    d["fcompl"] = np.minimum(d["compl"] / plateau, 1.0)
    Ldel, Lrec = d.Ldel.sum() / 1e3, d.Lrec.sum() / 1e3
    fc = float((d.Lrec * d.fcompl).sum() / d.Lrec.sum())
    print(f"\n== {run}: n_LB {len(d)}, L_del {Ldel:.4f}, live_RD0 {Lrec/Ldel:.5f}, L_rec {Lrec:.4f} nb-1; plateau N_drzewo/TAV({c['ref']}) = {plateau:.4f}, kompletność plików {fc:.4f}, surowa {d.n_ref_tree.sum()/d['TAV_'+c['ref']].sum():.4f}")
    out[str(run)] = {}
    ps = pd.read_csv(f"{T}/ps_perLB_{run}.csv"); ps = ps[ps.kind == "HLT"]
    for name in c["chains"]:
        p = ps[ps.chain == name].set_index("lumiBlock").reindex(d.index).ps_mode
        own = d[f"TAV_{SEED[name]}"] / d[f"TAP_{SEED[name]}"]
        L = float((d.Lrec / p).sum() / 1e3); Lown = float((d.Ldel * own / p).sum() / 1e3); Lcopy = float((d.Lrec * d.fcompl / p).sum() / 1e3)
        d[f"L_{name}"] = d.Lrec / p
        live_sys = abs(Lown / L - 1)
        err_live = max(live_sys, 0.003)
        rel = float(np.hypot(REL_CAL, err_live))
        out[str(run)][name] = {"L": round(L, 5), "jednostka": "nb^-1", "błąd": round(L * rel, 5), "błąd_wzgl": round(rel, 4), "błąd_wzgl_kalibracja_ZAŁOŻENIE": REL_CAL, "błąd_wzgl_live": round(err_live, 4),
                               "tag": TAG, "LB_min": int(d.index.min()), "LB_max": int(d.index.max()), "n_LB": int(len(d)), "preskala": sorted(set(float(x) for x in np.round(p, 2))),
                               "L_wariant_live_wlasna_pozycja_L1": round(Lown, 5), "L_nasza_kopia_plikow": round(Lcopy, 5), "kompletnosc_plikow_ZAŁOŻENIE": round(fc, 4),
                               "źródło": "L_del: iLumiCalc/cvmfs GoodRunsLists data22_13p6TeV/20231005 (Inst_m_Lumi*DeltaT per LB); live: COOL /TRIGGER/LUMI/LVL1COUNTERS L1_RD0_FILLED TAV/TAP; preskala: drzewo = COOL L1xHLT (TAGREG v2); LB: GRL_MinBiasAFP_lowmu2022.xml ∩ nasze pliki"}
        print(f"  {name}: ps {out[str(run)][name]['preskala']}, L = {L:.5f} ± {L*rel:.5f} nb-1 (wariant live własny {Lown:.5f}; kopia plików {Lcopy:.5f})")
    if run == 428770:
        a, j = "HLT_noalg_L1AFP_A_OR_C", "HLT_j20_L1AFP_A_OR_C"
        pa = ps[ps.chain == a].set_index("lumiBlock").reindex(d.index).ps_mode; pj = ps[ps.chain == j].set_index("lumiBlock").reindex(d.index).ps_mode
        hlt_a = pa / pj      # czynnik HLT noalg (60)
        d["L_R1"] = d.Lrec / pa; d["L_R2"] = d.Lrec / pj * (1 - 1 / hlt_a); d["L_R3"] = d.Lrec
        early = d.index <= 16
        reg = {"R1": ("HLT_noalg_L1AFP_A_OR_C", d.L_R1), "R2": ("HLT_j20_L1AFP_A_OR_C ∧ ¬R1 (czynnik 1−1/60 z weta)", d.L_R2), "R3": ("HLT_noalg_L1J12 ∧ ¬R1 ∧ ¬R2 (weto ¬R2 = cięcie fizyczne; L bez poprawki)", d.L_R3)}
        out[str(run)]["_regiony"] = {}
        for k, (desc, s) in reg.items():
            out[str(run)]["_regiony"][k] = {"definicja": desc, "L": round(float(s.sum() / 1e3), 5), "jednostka": "nb^-1", "udział_LB_4-16": round(float(s[early].sum() / s.sum()), 4), "tag": TAG, "n_LB": int(len(d))}
            print(f"  {k}: {desc}: L = {s.sum()/1e3:.5f} nb-1; udział LB 4–16: {s[early].sum()/s.sum():.4f}")
        print(f"  L_rec LB 4–16: {d.Lrec[early].sum()/1e3:.4f}; LB 17–25: {d.Lrec[(d.index>16)&(d.index<100)].sum()/1e3:.4f}; LB 210–323: {d.Lrec[d.index>=210].sum()/1e3:.4f}")
    if run == 427929:
        tbp = d.TBP_L1_MBTS_1.sum(); sv = tbp / (d.Ldel.sum() * 1e3)
        exp_lo, exp_hi = 69.0, 72.0
        out[str(run)]["_sprawdzian_MBTS"] = {"N_TBP_L1_MBTS_1": int(tbp), "sigma_vis_zmierzone_mb": round(float(sv), 2), "sigma_vis_oczekiwane_mb": [exp_lo, exp_hi],
            "stosunek": [round(float(sv / exp_hi), 3), round(float(sv / exp_lo), 3)], "L_z_MBTS_nb": [round(float(tbp / exp_hi / 1e6), 4), round(float(tbp / exp_lo / 1e6), 4)],
            "frakcja_przeciec_z_MBTS_1": round(float(tbp / d.TBP_L1_RD0_FILLED.sum()), 6) if "TBP_L1_RD0_FILLED" in d else None, "MBTS_2/MBTS_1": round(float(d.TBP_L1_MBTS_2.sum() / tbp), 4)}
        print("  sprawdzian MBTS:", out[str(run)]["_sprawdzian_MBTS"])
    out[str(run)]["_run"] = {"L_del": round(Ldel, 5), "L_rec": round(Lrec, 5), "live_RD0_FILLED": round(Lrec / Ldel, 5), "n_LB": int(len(d)), "czas_s": round(float(d.DeltaT.sum()), 0), "mu_srednie_oficjalne": round(float(np.average(d.AvergeInteractionPerXing, weights=d.Ldel)), 5),
                            "plateau_N_drzewo/TAV": round(plateau, 4), "kompletnosc_plikow": round(fc, 4), "LB_z_ubytkiem_plikow(>3%)": [int(i) for i in d.index[d.fcompl < 0.97]]}
    keep = ["Ldel", "DeltaT", "Inst_m_Lumi", "AvergeInteractionPerXing", "live_RD0", "live", "Lrec", "n_grl", "n_ref_tree", "compl", "fcompl"] + [k for k in d if k.startswith("L_")] + [k for k in d if k.startswith(("TBP_", "TAP_", "TAV_"))]
    d[keep].to_csv(f"tabele/lumi_final_perLB_{run}.csv")
json.dump(out, open("lumi.json", "w"), indent=1, ensure_ascii=False)
