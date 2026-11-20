#!/usr/bin/env python3
"""Krok 1: opis próbki strumienia AFP → describe.json (+ stdout)."""
import json, numpy as np, lib
OUT = {}
MB = ("HLT_noalg_mb_L1MBTS_1", "HLT_noalg_mb_L1MBTS_2")
for run in lib.RUNS:
    d = lib.load(run); O = OUT[str(run)] = {"info": json.loads(str(d["info"]))}
    afp, non, has = d["afp"], d["non"], d["hasD"]
    ncand = d["cnt"].sum(1)
    ch = {}
    streams = {"noalg_L1AFP": d["T_HLT_noalg_L1AFP_A_OR_C"], "sptrk_pt4_L1AFP": d["T_HLT_mb_sptrk_pt4_L1AFP_A_OR_C"], "AFP_any": afp, "AFP_only": afp & ~non,
               "nonAFP": non, "MB_noalg_L1MBTS": d["T_HLT_noalg_mb_L1MBTS_1"] | d["T_HLT_noalg_mb_L1MBTS_2"], "sptrk_pt4_L1MBTS_2": d["T_HLT_mb_sptrk_pt4_L1MBTS_2"],
               "sptrk_pt2_L1MBTS_2": d["T_HLT_mb_sptrk_pt2_L1MBTS_2"], "L1J12": d["T_HLT_noalg_L1J12"], "j20_L1MBTS_2": d["T_HLT_j20_L1MBTS_2"]}
    for nm, m in streams.items():
        r = {"n_ev": int(m.sum()), "n_ev_hasD": int((m & has).sum()), "n_cand_windows": int(ncand[m].sum()), "mean_nsel": float(d["nsel"][m].mean()) if m.any() else None,
             "mean_nsel_hasD": float(d["nsel"][m & has].mean()) if (m & has).any() else None}
        for s in lib.SIDES:
            t = d[f"tag_{s}"]; dl = d[f"xical_{s}"] - d[f"xi1_{s}"]
            for lab, mm in (("all", m), ("hasD", m & has)):
                r[f"{lab}_tag_{s}"] = int((mm & t).sum())
                r[f"{lab}_SR_{s}"] = int((mm & t & (dl <= lib.DELTA)).sum())
                r[f"{lab}_gap_{s}"] = int((mm & t & (dl > lib.DELTA) & (dl <= lib.DELTA + lib.M_NOM)).sum())
                r[f"{lab}_CR_{s}"] = int((mm & t & (dl > lib.DELTA + lib.M_NOM)).sum())
                r[f"{lab}_CR_xical>0.14_{s}"] = int((mm & t & (d[f"xical_{s}"] > 0.14)).sum())
            r[f"anyproton_{s}"] = int((m & (d[f"np_{s}"] > 0)).sum())
        r["tag_AandC"] = int((m & d["tag_A"] & d["tag_C"]).sum())
        r["hasD_tag_AandC"] = int((m & has & d["tag_A"] & d["tag_C"]).sum())
        r["no_proton"] = int((m & (d["np_A"] == 0) & (d["np_C"] == 0)).sum())
        ch[nm] = r
    O["streams"] = ch
    # efektywna preskala łańcuchów AFP: wśród zdarzeń zapisanych przez nie-AFP z L1_AFP TBP (i dla pt4: z bliźniaczym HLT pt4_L1MBTS_2) ułamek z odpalonym łańcuchem AFP
    tbp = d["trigTBP_L1_AFP_A_OR_C"]
    ps = {}
    for nm, chn, extra in (("noalg_L1AFP", "T_HLT_noalg_L1AFP_A_OR_C", np.ones(len(tbp), bool)), ("sptrk_pt4_L1AFP", "T_HLT_mb_sptrk_pt4_L1AFP_A_OR_C", d["T_HLT_mb_sptrk_pt4_L1MBTS_2"])):
        den = non & tbp & extra; num = den & d[chn]
        ps[nm] = {"n_den": int(den.sum()), "n_num": int(num.sum()), "prescale_eff": float(den.sum() / max(num.sum(), 1))}
        # per blok LB (8)
        g = lib.lb_groups(d["lb"], den, 8)
        ps[nm]["per_lb_block"] = [float((den & (g == i)).sum() / max((num & (g == i)).sum(), 1)) for i in range(8)]
    # L1 AFP TBP bez protonu / z protonem
    anyp = (d["np_A"] > 0) | (d["np_C"] > 0)
    ps["P(L1AFP TBP | proton, nonAFP)"] = float((non & anyp & tbp).sum() / max((non & anyp).sum(), 1))
    ps["P(proton | L1AFP TBP, nonAFP)"] = float((non & anyp & tbp).sum() / max((non & tbp).sum(), 1))
    for s in lib.SIDES:
        one = non & (d[f"np_{s}"] == 1) & (d["np_C" if s == "A" else "np_A"] == 0)
        ed = np.array([0.0, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.12, 0.3])
        b = np.digitize(d[f"xi1_{s}"], ed) - 1
        ps[f"L1AFP_TBP_eff_vs_xi_{s}"] = {"edges": ed.tolist(), "n": [int((one & (b == i)).sum()) for i in range(len(ed) - 1)],
                                          "eff": [float((one & tbp & (b == i)).sum() / max((one & (b == i)).sum(), 1)) for i in range(len(ed) - 1)]}
    O["prescale"] = ps
    # histogramy
    H = {}
    ex = np.linspace(0, 0.14, 29); ec = np.concatenate([np.linspace(0, 0.3, 31), [0.4, 0.6, 1.0, 3.0]]); edl = np.concatenate([[-0.2], np.linspace(-0.1, 0.3, 41), [0.4, 0.6, 1.0, 3.0]])
    for nm in ("noalg_L1AFP", "sptrk_pt4_L1AFP", "AFP_any", "nonAFP"):
        m = streams[nm]
        for lab, mm in (("all", m), ("hasD", m & has)):
            for s in lib.SIDES:
                t = mm & d[f"tag_{s}"]
                H[f"{nm}|{lab}|{s}"] = {"xi": np.histogram(d[f"xi1_{s}"][t], ex)[0].tolist(), "xical": np.histogram(d[f"xical_{s}"][t], ec)[0].tolist(),
                                        "delta": np.histogram((d[f"xical_{s}"] - d[f"xi1_{s}"])[t], edl)[0].tolist(),
                                        "xical_notag": np.histogram(d[f"xical_{s}"][mm & (d[f"np_{s}"] == 0)], ec)[0].tolist()}
    O["hist"] = {"edges_xi": ex.tolist(), "edges_xical": ec.tolist(), "edges_delta": edl.tolist(), "h": H}
    # kandydaci per kategoria w strumieniu AFP_any: tag, SR, CR (A + C)
    cat = {}
    m = afp
    for i, c in enumerate(lib.CATS):
        w = d["cnt"][:, i].astype(float); r = {}
        for s in lib.SIDES:
            t = m & d[f"tag_{s}"]; dl = d[f"xical_{s}"] - d[f"xi1_{s}"]
            r[s] = {"tag": float(w[t].sum()), "SR": float(w[t & (dl <= lib.DELTA)].sum()), "CR": float(w[t & (dl > lib.DELTA + lib.M_NOM)].sum())}
        cat[c] = r
    O["cands_AFP_any"] = cat
    print("=====", run)
    for nm, r in ch.items():
        print(nm, {k: v for k, v in r.items()})
    print(json.dumps(ps, indent=0)[:3000])
    for c, r in cat.items(): print(c, r)
json.dump(OUT, open("describe.json", "w"), indent=1)
