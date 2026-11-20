"""FORCED f3 — tabela nominał vs warianty z result_v2.json (tylko frakcje tagu przypadkowego, liczniki, bramki)."""
import json, os, math, pandas as pd
H = os.path.dirname(os.path.abspath(__file__))
ref = json.load(open(f"{H}/../runs/lxplus/out_v2/result_v2.json"))["results"]; ref = next(r for r in ref if r["sample"] == "data_428770")
R = {"raport_v2": ref}
for v in ("nominal", "W1a", "W1", "W2"):
    R[v] = json.load(open(f"{H}/{v}/result_v2.json"))["results"][0]
n = R["nominal"]; f0 = n["estimate"]["f_acc"]; s0 = f0 * n["syst_rel"]["total_rel"]; st0 = n["estimate"]["f_acc_err_stat"]
rows = []
for k, r in R.items():
    e, sy, bi, bii = r["estimate"], r["syst_rel"], r["method_Bi"], r["method_Bii"]; nt = e["n_tag"]
    g = {x["gate"]: x for x in r["gates"]}
    row = dict(wariant=k, dawcy=r["inputs"]["n_donors"], biorcy=r["inputs"]["n_targets"], N_tag=nt, N_acc_A=e["n_acc"], N_acc_A_stat=e["n_acc_err_stat"],
               f_A=e["f_acc"], f_A_stat=e["f_acc_err_stat"], f_A_tot=e["f_acc"] * sy["total_rel"], rel_stat=sy["stat"], rel_pula=sy["pool"], rel_ekstr=sy["extrapolation"],
               rel_metod=sy["half_diff_methods"], rel_razem=sy["total_rel"],
               d_f=e["f_acc"] - f0, d_sigma_tot=(e["f_acc"] - f0) / s0, d_sigma_stat=(e["f_acc"] - f0) / st0,
               N_Bi=bi["n_acc"], N_Bi_err=bi["n_acc_err"], f_Bi=bi["n_acc"] / nt, f_Bi_err=bi["n_acc_err"] / nt, Bi_tagdaw=bi["n_tagged_donors"], Bi_pula=bi["n_pool"],
               N_Bii=bii["n_acc"], N_Bii_err=bii["n_acc_err"], f_Bii=bii["n_acc"] / nt, f_Bii_err=bii["n_acc_err"] / nt, CR_obs=bii["n_obs_cr"], CR_pred=bii["n_pred_cr"],
               f_reg_A=r["estimates_by_region"]["A"]["f_acc"], f_reg_C=r["estimates_by_region"]["C"]["f_acc"], eps=e["eps_xical_acc"], tagdaw=e["n_tagged_donors"],
               werdykt=r["verdict"], G_status="; ".join(f"{a.split('_')[0]}:{x['status']}" for a, x in g.items()),
               G1_p=g["G1_pool_homogeneity"]["p"], G2_ratio=g["G2_disjoint_closure"]["ratio"], G2_z=g["G2_disjoint_closure"]["z"], G5_p=g["G5_cr_shape"]["p"],
               G7_ratio=g["G7_methods_A_vs_Bii"]["ratio"], G7_z=g["G7_methods_A_vs_Bii"]["z"], G8_ratio=g["G8_methods_A_vs_Bi"]["ratio"], G8_z=g["G8_methods_A_vs_Bi"]["z"],
               G9_obs=g["G9_opposite_side"]["n_obs"], G9_pred=g["G9_opposite_side"]["n_pred"], G0_cov=g["G0_inputs"]["coverage"], I2_p=g["I2_hasD_vs_noD"]["p"])
    row["d_Bi_sigma"] = (row["f_Bi"] - R["nominal"]["method_Bi"]["n_acc"] / n["estimate"]["n_tag"]) / (R["nominal"]["method_Bi"]["n_acc_err"] / n["estimate"]["n_tag"])
    row["d_Bii_sigma"] = (row["f_Bii"] - R["nominal"]["method_Bii"]["n_acc"] / n["estimate"]["n_tag"]) / (R["nominal"]["method_Bii"]["n_acc_err"] / n["estimate"]["n_tag"])
    row["tot_z_dodatkiem"] = math.hypot(s0, row["d_f"]) / s0 - 1
    rows.append(row)
t = pd.DataFrame(rows).set_index("wariant"); t.to_csv(f"{H}/out/t3_facc_warianty.csv")
pd.set_option("display.width", 200, "display.float_format", lambda x: f"{x:.5g}"); print(t.T.to_string())
# gatunki × pT (pik): f_acc, bez fitów
sp = []
for name in ("Dplus", "D0", "Dstar"):
    for i, r0 in enumerate(n["species"][name]["rows"]):
        d = dict(gatunek=name, pT=f"{r0['pt_lo']:g}-{r0['pt_hi']:g}")
        for k in ("nominal", "W1a", "W1", "W2"):
            p = R[k]["species"][name]["rows"][i]["peak"]; d[f"Ntag_{k}"] = p["n_tag"]; d[f"f_{k}"] = p["f_acc"]
        d["stat_nom"] = r0["peak"]["f_acc_err_stat"]; d["maxd_w_stat"] = max(abs(d[f"f_{k}"] - d["f_nominal"]) for k in ("W1a", "W1", "W2")) / d["stat_nom"]
        sp.append(d)
s = pd.DataFrame(sp); s.to_csv(f"{H}/out/t4_gatunki_pik.csv", index=False); print(s.to_string())
