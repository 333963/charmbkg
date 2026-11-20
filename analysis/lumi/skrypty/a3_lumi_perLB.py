"""a3: świetlność per LB i per łańcuch dla LB obecnych w naszych plikach.
L_del: oficjalny ilumicalc cvmfs (OflLumi-Run3-003): Inst_m_Lumi * DeltaT [ub^-1].
live: TAV/TAP pozycji L1 o dużej częstości z COOL /TRIGGER/LUMI/LVL1COUNTERS (427929: L1_MBTS_1; 428770: L1_J12).
preskale: drzewo (= L1 x HLT z COOL, TAGREG v2). kompletność kopii: zdarzenia w drzewie z bitem / TAV."""
import pandas as pd, numpy as np, json
T = "/home/k/Documents/charm/charm_work/TAGREG_2026-11-18/runs/lxplus"
off = pd.read_csv("tabele/cvmfs_20231005__ilumicalc_histograms_None_427929-428770_OflLumi-Run3-003__LumiMetaData.csv")
CFG = {427929: dict(live="L1_MBTS_1", ref="HLT_noalg_mb_L1MBTS_1", chains=["HLT_mb_mbts_L1MBTS_1", "HLT_noalg_mb_L1MBTS_1"]),
       428770: dict(live="L1_J12", ref="HLT_noalg_L1J12", chains=["HLT_noalg_L1AFP_A_OR_C", "HLT_j20_L1AFP_A_OR_C", "HLT_noalg_L1J12", "HLT_mb_sptrk_pt4_L1AFP_A_OR_C", "HLT_noalg_mb_L1MBTS_1", "HLT_mb_mbts_L1MBTS_1"])}
res = {}
for run, c in CFG.items():
    o = off[off.RunNbr == run].set_index("LBStart"); o["Ldel"] = o.Inst_m_Lumi * o.DeltaT
    cn = pd.read_csv(f"surowe/lxplus/cool/l1counters_{run}.csv"); lbs = pd.read_csv(f"{T}/lb_summary_{run}.csv").set_index("lumiBlock")
    ps = pd.read_csv(f"{T}/ps_perLB_{run}.csv"); ps = ps[ps.kind == "HLT"]
    d = lbs[lbs.n_grl > 0][["n_grl", "mean_mu_grl"]].join(o[["Ldel", "DeltaT", "Inst_m_Lumi", "Livefrac", "IntLumi", "AvergeInteractionPerXing"]])
    for it, g in cn.groupby("item"):
        g = g.set_index("LB"); d[f"TBP_{it}"] = g.TBP; d[f"TAP_{it}"] = g.TAP; d[f"TAV_{it}"] = g.TAV
    L = c["live"]; d["live"] = d[f"TAV_{L}"] / d[f"TAP_{L}"]
    # kontrola live inną pozycją
    for alt in ("L1_MBTS_1", "L1_J12", "L1_AFP_A_OR_C", "L1_MBTS_2", "L1_EM22VHI", "L1_RD0_FILLED"):
        if alt != L and f"TAP_{alt}" in d: d[f"live_{alt}"] = d[f"TAV_{alt}"] / d[f"TAP_{alt}"].where(d[f"TAP_{alt}"] > 0)
    ref = ps[ps.chain == c["ref"]].set_index("lumiBlock"); d["n_ref_tree"] = ref.n_pass; d["compl"] = d.n_ref_tree / d[f"TAV_{L}"]
    d["Lrec"] = d.Ldel * d.live
    tot = dict(n_LB=len(d), LB_min=int(d.index.min()), LB_max=int(d.index.max()), L_del_nb=d.Ldel.sum() / 1e3, L_rec_nb=d.Lrec.sum() / 1e3, live=d.Lrec.sum() / d.Ldel.sum(),
               live_min=d.live.min(), L_rec_EM22VHI_official_nb=d.IntLumi.sum() / 1e3, compl=d.n_ref_tree.sum() / d[f"TAV_{L}"].sum(), compl_min=d["compl"].min(), compl_max=d["compl"].max(),
               mu_official=float(np.average(d.AvergeInteractionPerXing, weights=d.Ldel)), mu_tree=float(np.average(d.mean_mu_grl, weights=d.n_grl)), time_s=d.DeltaT.sum(),
               TBP=int(d[f"TBP_{L}"].sum()), TAP=int(d[f"TAP_{L}"].sum()), TAV=int(d[f"TAV_{L}"].sum()), n_ref_tree=int(d.n_ref_tree.sum()))
    tot["sigma_vis_TBP_mb"] = tot["TBP"] / (d.Ldel.sum() * 1e3) if run == 427929 else None   # ub^-1 -> mb^-1: 1 ub^-1 = 1e3 mb^-1
    for alt in [k for k in d if k.startswith("live_")]:
        tot[alt] = float((d.Ldel * d[alt]).sum() / d.Ldel[d[alt].notna()].sum())
    print(f"\n===== {run}: " + ", ".join(f"{k}={v:.5g}" if isinstance(v, float) else f"{k}={v}" for k, v in tot.items()))
    ch = {}
    for name in c["chains"]:
        p = ps[ps.chain == name].set_index("lumiBlock").reindex(d.index)
        d[f"ps_{name}"] = p.ps_mode; d[f"n_{name}"] = p.n_pass
        inv = np.where(p.ps_mode > 0, 1 / p.ps_mode, 0)
        d[f"Leff_{name}"] = d.Lrec * inv; d[f"Leff_copy_{name}"] = d.Lrec * inv * d["compl"]
        ch[name] = dict(ps=sorted(set(np.round(p.ps_mode, 2))), L_rec_ps_nb=float(d[f"Leff_{name}"].sum() / 1e3), L_copy_nb=float(d[f"Leff_copy_{name}"].sum() / 1e3), F=float(d[f"Leff_{name}"].sum() / d.Lrec.sum()), n_tree=int(p.n_pass.sum()))
        print(f"  {name}: ps {ch[name]['ps']}, F={ch[name]['F']:.5f}, L(rec/ps) = {ch[name]['L_rec_ps_nb']:.5f} nb-1, L(nasza kopia) = {ch[name]['L_copy_nb']:.5f} nb-1, N_tree = {ch[name]['n_tree']}")
    print("  kompletność per LB: kwantyle", d["compl"].quantile([0, .05, .25, .5, .75, .95, 1]).round(4).tolist(), "; LB z compl<0.5:", d.index[d["compl"] < 0.5].tolist()[:40])
    print("  live per LB: kwantyle", d.live.quantile([0, .05, .5, .95, 1]).round(4).tolist())
    d.to_csv(f"tabele/lumi_perLB_{run}.csv"); res[run] = dict(tot=tot, chains=ch)
json.dump(res, open("tabele/a3.json", "w"), indent=1, default=float)
