#!/usr/bin/env python3
"""AFPTRIG/make_report.py — składa wynik: facc_afptrig.json, rysunki PNG, RAPORT.md (z results_raw_*.json, describe.json,
extra.json, leak.json, diag2d.json, strat_var.json, toy_results.json). Bez uzysków sygnału; tylko frakcje tła i liczności."""
import json, math, time
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import lib

RUNS = ("427929", "428770")
RAW = {r: json.load(open(f"results_raw_{r}.json"))[r] for r in RUNS}
DESC = json.load(open("describe.json")); EXTRA = json.load(open("extra.json")); LEAK = json.load(open("leak.json"))
D2 = json.load(open("diag2d.json")); SV = json.load(open("strat_var.json")); TOY = json.load(open("toy_results.json"))
V2 = {"427929": dict(f=0.0087, e=0.0044, n_tag=18871, n_acc=164.7, A=(9706, 91.7, 16.9), C=(9229, 73.9, 12.6)),
      "428770": dict(f=0.0477, e=0.0135, n_tag=16065, n_acc=766.1, A=(8271, 376.5, 18.0), C=(7863, 393.4, 19.6))}
STREAM_PL = {"afp_any": "AFP (dowolny łańcuch L1AFP)", "pt4_L1AFP": "HLT_mb_sptrk_pt4_L1AFP_A_OR_C", "noalg_L1AFP": "HLT_noalg_L1AFP_A_OR_C", "nonAFP": "nie-AFP (zakres v2)"}
POOL_PL = {"twin": "bliźniak HLT_mb_sptrk_pt4_L1MBTS_2", "mb": "MB: noalg_mb_L1MBTS_1/2 ∪ sptrk(_pt2)_L1MBTS", "mbts": "MBTS = bliźniak ∪ MB", "nonAFP": "wszystkie nie-AFP", "j12": "HLT_noalg_L1J12", "self": "sam strumień (tasowanie)"}
ALT = {("427929", "afp_any"): "twin", ("427929", "noalg_L1AFP"): "twin", ("427929", "pt4_L1AFP"): "nonAFP", ("427929", "nonAFP"): "twin",
       ("428770", "afp_any"): "nonAFP", ("428770", "pt4_L1AFP"): "nonAFP", ("428770", "noalg_L1AFP"): "mbts", ("428770", "nonAFP"): "mbts"}


def comb(es):
    """suma stron: N liniowo, błąd stat. w kwadracie."""
    return {"n_acc": sum(e["n_acc"] for e in es), "err": math.sqrt(sum(e["err"] ** 2 for e in es)), "n_sr": sum(e["n_sr"] for e in es), "n_cr": sum(e["n_cr"] for e in es)}


def AC(o, path):
    es = []
    for s in "AC":
        x = o[s]
        for p in path:
            x = x[p]
        es.append(x)
    return comb(es)


def sd_nested(a, b):
    """Błąd różnicy estymatorów z zagnieżdżonych obszarów normalizacji (CR_b ⊂ CR_a lub odwrotnie); przybliżenie:
    Δ = N_b(R_b − R_a) − n·R_a, n = N_CRa − N_CRb niezależne; czynnik κ = err/(√N_CR·R_eff) uwzględnia rozrzut R między komórkami."""
    if a["n_cr"] < b["n_cr"]:
        a, b = b, a
    Ra, Rb = a["n_acc"] / max(a["n_cr"], 1), b["n_acc"] / max(b["n_cr"], 1); n = a["n_cr"] - b["n_cr"]
    kap = max(a["err"] / (math.sqrt(max(a["n_cr"], 1)) * Ra), b["err"] / (math.sqrt(max(b["n_cr"], 1)) * Rb), 1.0)
    return kap * math.sqrt(b["n_cr"] * (Rb - Ra) ** 2 + n * Ra ** 2)


def sig_part(d, sd):
    return abs(d) if abs(d) > 2 * sd else max(abs(d) - sd, 0.0)


def budget(run, stream):
    o = RAW[run][stream]; nom = AC(o, ["nominal"]); N = nom["n_acc"]; b = {}
    b["stat"] = nom["err"] / N
    alt = AC(o, ["pool", ALT[(run, stream)]]); b["pool_half_diff"] = 0.5 * abs(alt["n_acc"] - N) / N
    b["pool_alt"] = ALT[(run, stream)]; b["pool_alt_ratio"] = alt["n_acc"] / N
    pv = AC(o, ["protons", "nonafp_veto"]); b["protons"] = abs(pv["n_acc"] - N) / N
    b["lb_mu"] = max(abs(AC(o, ["strat", k])["n_acc"] - N) / N for k in ("nlb1", "nlb4", "nlb16", "mu8"))
    sv = SV[run][stream]; b["nsel_gran"] = max(abs(sv[k]["n_acc"] - sv["nominal"]["n_acc"]) / sv["nominal"]["n_acc"] for k in ("coarse", "fine", "vfine"))
    b["no_nsel_ratio_info"] = AC(o, ["strat", "no_nsel"])["n_acc"] / N
    far = AC(o, ["margin", "0.08"]); d = (far["n_acc"] - N) / N; sd = sd_nested(nom, far) / N
    b["far_norm_ratio"] = far["n_acc"] / N; b["far_norm_ratio_err"] = sd
    b["leak_down"] = sig_part(d, sd) if d < 0 else 0.0
    b["margin_up"] = sig_part(d, sd) if d > 0 else 0.0
    c = AC(o, ["closure", "norm>0.20_pred_0.10-0.20", o["nominal_pool"]]); r = c["n_sr"] / c["n_acc"]; er = r * math.hypot(1 / math.sqrt(c["n_sr"]), c["err"] / c["n_acc"])
    b["closure_far_ratio"] = r; b["closure_far_err"] = er
    b["closure_up"] = sig_part(r - 1, er) if r > 1 else 0.0; b["closure_down"] = sig_part(r - 1, er) if r < 1 else 0.0
    if stream in EXTRA[run]:
        ob = sum(EXTRA[run][stream][s]["mirror_low"]["obs"] for s in "AC"); pr = sum(EXTRA[run][stream][s]["mirror_low"]["pred"] for s in "AC")
        ke = math.sqrt(sum((EXTRA[run][stream][s]["mirror_low"]["ratio_err"] * EXTRA[run][stream][s]["mirror_low"]["pred"]) ** 2 for s in "AC")) / pr
        b["mirror_kappa"] = ob / pr; b["mirror_kappa_err"] = ke; b["mirror_down"] = math.sqrt(max(1 - ob / pr, 0) ** 2 + ke ** 2)
    else:
        b["mirror_kappa"] = None; b["mirror_down"] = None
    sym = [b["pool_half_diff"], b["protons"], b["lb_mu"], b["nsel_gran"]]
    b["syst_up"] = math.sqrt(sum(x * x for x in sym + [b["closure_up"], b["margin_up"]]))
    b["syst_down_no_mirror"] = math.sqrt(sum(x * x for x in sym + [b["leak_down"], b["closure_down"]]))
    b["syst_down"] = math.sqrt(b["syst_down_no_mirror"] ** 2 + (b["mirror_down"] or 0) ** 2)
    b["total_up"] = math.hypot(b["stat"], b["syst_up"]); b["total_down"] = math.hypot(b["stat"], b["syst_down"]); b["total_down_no_mirror"] = math.hypot(b["stat"], b["syst_down_no_mirror"])
    return nom, b


def kcheck(run, stream):
    """Kryteria K1, K4, K5 dla strumienia (A + C)."""
    o = RAW[run][stream]; np_ = o["nominal_pool"]; K = {}
    for nm in ("norm>0.10_pred_0.045-0.10", "norm>0.10_pred_0.045-0.07", "norm>0.10_pred_0.07-0.10", "norm>0.20_pred_0.10-0.20", "norm>0.045_pred_0.02-0.03", "norm>0.045_pred_0.03-0.045"):
        c = AC(o, ["closure", nm, np_]); z = (c["n_sr"] - c["n_acc"]) / math.sqrt(c["n_sr"] + c["err"] ** 2); r = c["n_sr"] / c["n_acc"]
        K[nm] = {"obs": c["n_sr"], "pred": c["n_acc"], "pred_err": c["err"], "ratio": r, "ratio_err": r * math.hypot(1 / math.sqrt(max(c["n_sr"], 1)), c["err"] / c["n_acc"]), "z": z,
                 "pass": bool(abs(z) < 2.5 or abs(r - 1) < 0.05), "per_side": {s: o[s]["closure"][nm][np_]["n_sr"] / o[s]["closure"][nm][np_]["n_acc"] for s in "AC"}}
    K["K1"] = "PASS" if K["norm>0.10_pred_0.045-0.10"]["pass"] else "FAIL"
    # K4: kształt w CR
    ed = np.array(o["A"]["hist_delta"]["edges"]); dat = sum(np.array(o[s]["hist_delta"]["data"]) for s in "AC"); pr = sum(np.array(o[s]["hist_delta"]["pred"][np_]) for s in "AC")
    m = (ed[:-1] >= 0.045 - 1e-9) & (ed[1:] <= 0.30 + 1e-9)
    d, p = dat[m], pr[m]
    # łączenie do ≥ 10 oczekiwanych
    gd, gp, cd, cp = [], [], 0.0, 0.0
    for a, b_ in zip(d, p):
        cd += a; cp += b_
        if cp >= 10:
            gd.append(cd); gp.append(cp); cd = cp = 0.0
    if cp > 0 and gp:
        gd[-1] += cd; gp[-1] += cp
    gd, gp = np.array(gd), np.array(gp)
    scale = gd.sum() / gp.sum()                      # test kształtu: normalizacja w badanym zakresie
    chi2 = float(((gd - gp * scale) ** 2 / (gp * scale)).sum()); ndf = len(gd) - 1
    pval = float(stats.chi2.sf(chi2, ndf)) if ndf > 0 else None
    # odchylenie zgrupowane: 3 równe części zakresu
    th = np.array_split(np.arange(len(gd)), 3); dev = max(abs(gd[i].sum() / (gp[i].sum() * scale) - 1) for i in th if len(i))
    K["K4"] = {"chi2": chi2, "ndf": ndf, "p": pval, "max_grouped_dev": float(dev), "status": "PASS" if (pval is not None and pval >= 0.01) or dev < 0.10 else "FAIL"}
    nom = AC(o, ["nominal"]); rows = {}; ok = True
    for mm in ("0.01", "0.04", "0.06"):
        x = AC(o, ["margin", mm]); dlt = x["n_acc"] - nom["n_acc"]; sd = sd_nested(nom, x)
        rows[mm] = {"ratio": x["n_acc"] / nom["n_acc"], "ratio_err": sd / nom["n_acc"], "z": dlt / sd if sd > 0 else 0.0}
        ok &= abs(dlt) / nom["n_acc"] < 0.05 or abs(dlt) < 2.5 * sd
    rows["0.08"] = {"ratio": AC(o, ["margin", "0.08"])["n_acc"] / nom["n_acc"]}; rows["0.12"] = {"ratio": AC(o, ["margin", "0.12"])["n_acc"] / nom["n_acc"]}; rows["0"] = {"ratio": AC(o, ["margin", "0"])["n_acc"] / nom["n_acc"]}
    K["K5"] = {"rows": rows, "status": "PASS" if ok else "FAIL"}
    return K


# ------------------------------------------------------------------------------------------------ zbierz wynik
RES = {"meta": {"generated": time.strftime("%Y-%m-%d %H:%M"), "status": "WSTĘPNE", "method": "kształt Δ = ξ_cal − ξ_AFP: N_acc(SR) = Σ_komórki N_CR(dane) · T(SR)/T(CR); SR: Δ ≤ 0,02; CR: Δ > 0,045",
                 "selection": "GRL ∧ nPileupVertices == 0 ∧ ≥ 1 kandydat D (slim, bez score); tag = dokładnie 1 proton po stronie, 0,03 ≤ ξ ≤ 0,10",
                 "note": "tylko frakcje tła przypadkowego i liczności; f_acc(A+C) = suma stron (zdarzenia z tagiem po obu stronach liczone dwukrotnie)",
                 "f_acc_far_norm": "wariant po obejrzeniu walidacji: normalizacja w Δ > 0,10 (mniej przecieku własnych protonów)"}, "runs": {}}
for run in RUNS:
    RR = RES["runs"][run] = {"streams": {}}
    for stream in ("afp_any", "pt4_L1AFP", "noalg_L1AFP", "nonAFP"):
        o = RAW[run][stream]; nom, b = budget(run, stream); K = kcheck(run, stream)
        f = nom["n_acc"] / nom["n_sr"]; far = AC(o, ["margin", "0.08"])
        S = {"central_pool": o["nominal_pool"], "n_tag_SR": nom["n_sr"], "n_CR": nom["n_cr"], "n_acc": nom["n_acc"], "n_acc_err_stat": nom["err"], "f_acc": f, "f_acc_err_stat": f * b["stat"],
             "f_acc_err_syst_up": f * b["syst_up"], "f_acc_err_syst_down": f * b["syst_down"], "f_acc_err_syst_down_no_mirror": f * b["syst_down_no_mirror"],
             "f_acc_err_total_up": f * b["total_up"], "f_acc_err_total_down": f * b["total_down"], "f_acc_err_total_down_no_mirror": f * b["total_down_no_mirror"],
             "f_acc_far_norm": far["n_acc"] / nom["n_sr"], "f_acc_far_norm_err_stat": far["err"] / nom["n_sr"],
             "budget_rel": b, "validation": K,
             "sides": {s: {k: o[s]["nominal"][k] for k in ("n_sr", "n_cr", "n_acc", "err", "f_acc", "f_acc_err_stat", "n_pool", "n_protons")} for s in "AC"}}
        if "cats" in o["A"]:
            C = {}
            for c in lib.CATS:
                e = AC(o, ["cats", c, "nominal"]); ff = e["n_acc"] / e["n_sr"] if e["n_sr"] > 0 else None
                C[c] = {"n_tag_SR": e["n_sr"], "n_CR": e["n_cr"], "n_acc": e["n_acc"], "n_acc_err_stat": e["err"], "f_acc": ff, "f_acc_err_stat": e["err"] / e["n_sr"] if e["n_sr"] > 0 else None,
                        "f_acc_err_syst_up": ff * b["syst_up"] if ff is not None else None, "f_acc_err_syst_down": ff * b["syst_down"] if ff is not None else None,
                        "n_acc_alt_pool": AC(o, ["cats", c, "alt_nonAFP"])["n_acc"],
                        "sides": {s: {"n_tag_SR": o[s]["cats"][c]["nominal"]["n_sr"], "n_CR": o[s]["cats"][c]["nominal"]["n_cr"], "n_acc": o[s]["cats"][c]["nominal"]["n_acc"], "n_acc_err_stat": o[s]["cats"][c]["nominal"]["err"]} for s in "AC"}}
            S["categories"] = C
        RR["streams"][stream] = S
    # K2
    nn = RR["streams"]["nonAFP"]; v = V2[run]; es = math.hypot(nn["f_acc_err_stat"], nn["f_acc"] * max(nn["budget_rel"]["syst_up"], nn["budget_rel"]["syst_down_no_mirror"]))
    z = (nn["f_acc"] - v["f"]) / math.hypot(es, v["e"])
    RR["K2"] = {"f_shape": nn["f_acc"], "f_shape_err": es, "f_v2": v["f"], "f_v2_err": v["e"], "z": z, "status": "PASS" if abs(z) < 2.5 else "FAIL",
                "n_acc_shape": nn["n_acc"], "n_acc_v2": v["n_acc"], "per_side": {s: {"shape": nn["sides"][s]["n_acc"], "shape_err": nn["sides"][s]["err"], "v2": v[s][1], "v2_err": v[s][2]} for s in "AC"}}
# K3 toy
TT = {}
for sc, runs in TOY["scenarios"].items():
    row = {}
    for nm in ("nominal", "no_nsel", "pool_J12_strat", "pool_J12_nostrat", "pool_self", "prot_alltags", "m0.01", "m0.06"):
        # stosunek sum po stronach, per seed
        r = np.array([sum(o[s][nm][0] for s in "AC") / sum(o[s]["true"] for s in "AC") for o in runs]); row[nm] = [float(r.mean()), float(r.std(ddof=1) / math.sqrt(len(r)))]
    c = np.array([sum(o[s]["closure"][0] for s in "AC") / sum(o[s]["closure"][1] for s in "AC") for o in runs]); row["closure"] = [float(c.mean()), float(c.std(ddof=1) / math.sqrt(len(c)))]
    row["f_acc_true"] = float(np.mean([sum(o[s]["true"] for s in "AC") / sum(o[s]["n_sr"] for s in "AC") for o in runs]))
    TT[sc] = row
m, e = TT["nominal"]["nominal"]
RES["K3_toy"] = {"table": TT, "nominal_bias": m, "nominal_bias_err": e, "status": "PASS" if abs(m - 1) < 0.05 or abs(m - 1) < 2 * e else "FAIL", "n_seeds": len(TOY["seeds"])}
for run in RUNS:
    for st, S in RES["runs"][run]["streams"].items():
        K = S["validation"]
        S["K6_usable"] = bool(K["K1"] == "PASS" and RES["runs"][run]["K2"]["status"] == "PASS" and RES["K3_toy"]["status"] == "PASS" and max(S["budget_rel"]["total_up"], S["budget_rel"]["total_down"]) < 0.5)
json.dump(RES, open("facc_afptrig.json", "w"), indent=1, ensure_ascii=False)

# ------------------------------------------------------------------------------------------------ rysunki
C1, C2, C3, C4 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"; INK = "#222222"; MUTE = "#6b6b68"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTE, "axes.labelcolor": INK, "xtick.color": MUTE, "ytick.color": MUTE, "axes.grid": True, "grid.color": "#e6e5e0", "grid.linewidth": 0.6,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": "white", "axes.facecolor": "white", "legend.frameon": False})


def atlas(ax, txt, x=0.04, y=0.97):
    from matplotlib.offsetbox import AnchoredOffsetbox, HPacker, VPacker, TextArea
    top = HPacker(children=[TextArea("ATLAS", textprops=dict(fontweight="bold", fontstyle="italic", fontsize=12, color=INK)), TextArea("", textprops=dict(fontsize=12, color=INK))], pad=0, sep=5, align="baseline")
    box = VPacker(children=[top, TextArea(txt, textprops=dict(fontsize=8.5, color=INK))], pad=0, sep=3, align="left")
    ax.add_artist(AnchoredOffsetbox(loc="upper left", child=box, frameon=False, bbox_to_anchor=(x, y), bbox_transform=ax.transAxes, borderpad=0, pad=0))


# Rys. 1: rozkład Δ, dane vs szablon przypadkowych
fig, axs = plt.subplots(2, 2, figsize=(11, 6.6), gridspec_kw={"height_ratios": [3, 1.25], "hspace": 0.07, "wspace": 0.22}, sharex="col")
for j, run in enumerate(("428770", "427929")):
    o = RAW[run]["afp_any"]; np_ = o["nominal_pool"]; ed = np.array(o["A"]["hist_delta"]["edges"])
    dat = sum(np.array(o[s]["hist_delta"]["data"]) for s in "AC").astype(float); pr = sum(np.array(o[s]["hist_delta"]["pred"][np_]) for s in "AC")
    m = ed[1:] <= 0.3001; w = np.diff(ed)[m] / 0.01; x = 0.5 * (ed[:-1] + ed[1:])[m]; d, p = dat[m] / w, pr[m] / w; de = np.sqrt(dat[m]) / w
    ax = axs[0, j]
    ax.stairs(p, ed[: m.sum() + 1], color=C2, lw=2, label="szablon przypadkowych (mieszanie,\nnormalizacja w Δ > 0,045)")
    ax.errorbar(x, d, de, fmt="o", ms=4, color=C1, lw=1, label="dane: zdarzenia z D i tagiem (A + C)")
    ax.set_yscale("log"); ax.set_ylabel("zdarzenia / 0,01"); ax.set_ylim(max(0.5, p[p > 0].min() * 0.3), d.max() * 60)
    ax.axvline(0.02, color=MUTE, lw=1, ls="--"); ax.axvline(0.045, color=MUTE, lw=1, ls=":")
    ax.text(0.012, 0.03, "obszar sygnału Δ ≤ 0,02", rotation=90, fontsize=8, color=MUTE, ha="right", va="bottom", transform=ax.get_xaxis_transform())
    ax.text(0.052, 0.03, "obszar kontrolny Δ > 0,045", rotation=90, fontsize=8, color=MUTE, ha="left", va="bottom", transform=ax.get_xaxis_transform())
    atlas(ax, f"dane 2022, run {run}\nstrumień trygera AFP, kandydaci D bez score\nWSTĘPNE")
    ax.legend(loc="upper right", fontsize=8, bbox_to_anchor=(1.0, 0.80))
    ax2 = axs[1, j]; ok = (p > 0) & (x > 0.02)
    ax2.axhline(1, color=MUTE, lw=1); ax2.errorbar(x[ok], d[ok] / p[ok], de[ok] / p[ok], fmt="o", ms=4, color=C1, lw=1)
    ax2.axvline(0.02, color=MUTE, lw=1, ls="--"); ax2.axvline(0.045, color=MUTE, lw=1, ls=":")
    ax2.set_ylim(0.4, 2.1 if run == "428770" else 2.6); ax2.set_ylabel("dane / szablon"); ax2.set_xlabel("Δ = ξ_cal − ξ_AFP"); ax2.set_xlim(-0.10, 0.30)
fig.savefig("fig1_delta.png", dpi=150, bbox_inches="tight"); plt.close(fig)

# Rys. 2: diagnostyka (ξ_p, ξ_cal) + skan marginesu
fig, axs = plt.subplots(1, 3, figsize=(15, 4.6), gridspec_kw={"wspace": 0.25})
for ax, st, ttl in ((axs[0], "pt4_L1AFP", "run 428770, łańcuch sptrk_pt4_L1AFP"), (axs[1], "nonAFP", "run 428770, strumień nie-AFP")):
    g = D2["428770"][st]; xc = np.array(g["xc_edges"]); xp = np.array(g["xp_edges"]); r = np.array(g["ratio"]); e = np.array(g["ratio_err"])
    xm = np.array([0.5 * (xc[i] + xc[i + 1]) for i in range(len(xc) - 2)] + [0.36])
    for i, (col, mk) in enumerate(zip((C1, C2, C3, C4), ("o", "s", "^", "D"))):
        ok = np.array([xc[k] - xp[i + 1] >= 0.02 - 1e-9 for k in range(len(xc) - 1)])
        ax.errorbar(xm[ok] + (i - 1.5) * 0.002, r[i][ok], e[i][ok], fmt=mk + "-", ms=5, lw=1.5, color=col, label=f"ξ_AFP {xp[i]:.3f}–{min(xp[i + 1], 0.1):.3f}".replace(".", ","))
    ax.axhline(1, color=MUTE, lw=1); ax.set_ylim(0.6, 2.1); ax.set_xlabel("ξ_cal po stronie protonu (ostatni punkt: > 0,30)"); ax.set_ylabel("dane / szablon (norm. w ξ_cal > 0,20)")
    atlas(ax, ttl + "\nzdarzenia z D, A + C; WSTĘPNE"); ax.legend(loc="upper right", fontsize=8, bbox_to_anchor=(1.0, 0.80))
ax = axs[2]; ms = ["0", "0.01", "0.025", "0.04", "0.06", "0.08", "0.12"]
for (run, st, lab), col, mk in zip((("428770", "afp_any", "428770 strumień AFP"), ("428770", "nonAFP", "428770 nie-AFP"), ("427929", "afp_any", "427929 strumień AFP")), (C1, C2, C3), ("o", "s", "^")):
    o = RAW[run][st]; n0 = AC(o, ["nominal"]); y = []; ye = []
    for mm in ms:
        x = AC(o, ["margin", mm]); y.append(x["n_acc"] / n0["n_acc"]); ye.append(sd_nested(n0, x) / n0["n_acc"] if mm != "0.025" else 0.0)
    ax.errorbar(np.array([float(v) for v in ms]) + 0.02 + (0.0015 if st == "nonAFP" else 0) - (0.0015 if run == "427929" else 0), y, ye, fmt=mk + "-", ms=5, lw=1.5, color=col, label=lab)
ax.axhline(1, color=MUTE, lw=1); ax.set_xlabel("dolna granica obszaru normalizacji Δ_min = 0,02 + margines"); ax.set_ylabel("N_acc / N_acc(nominał: Δ_min = 0,045)"); ax.set_ylim(0.55, 1.45)
atlas(ax, "skan marginesu; błąd = błąd różnicy\nwzględem nominału; WSTĘPNE"); ax.legend(loc="lower left", fontsize=8)
fig.savefig("fig2_walidacja.png", dpi=150, bbox_inches="tight"); plt.close(fig)

# Rys. 3: toy
fig, ax = plt.subplots(figsize=(8.6, 4.4)); names = list(TT.keys()); y = np.arange(len(names))[::-1]
ax.axvline(1, color=MUTE, lw=1); ax.axvspan(0.95, 1.05, color="#e6e5e0", lw=0)
ax.errorbar([TT[n]["nominal"][0] for n in names], y + 0.12, xerr=[TT[n]["nominal"][1] for n in names], fmt="o", ms=6, color=C1, lw=1.5, label="estymator nominalny (pula MB, warstwy LB × nSelBL)")
ax.errorbar([TT[n]["pool_J12_nostrat"][0] for n in names], y - 0.12, xerr=[TT[n]["pool_J12_nostrat"][1] for n in names], fmt="s", ms=5, color=C2, lw=1.5, label="pula obciążona aktywnością, bez warstw nSelBL")
ax.set_yticks(y); ax.set_yticklabels([n.split(" (")[0] for n in names], color=INK); ax.set_xlabel("N_acc(estymator) / N_acc(prawda), strony A + C"); ax.set_xlim(0.2, 2.3)
ax.grid(axis="y", visible=False); ax.legend(loc="upper right", fontsize=8, bbox_to_anchor=(1.0, 0.72))
atlas(ax, f"toy (symulacja uproszczona), {len(TOY['seeds'])} seedów\ntryger wymagający protonu", x=0.66, y=0.97)
fig.savefig("fig3_toy.png", dpi=150, bbox_inches="tight"); plt.close(fig)
print("OK: facc_afptrig.json, fig1_delta.png, fig2_walidacja.png, fig3_toy.png")
for run in RUNS:
    for st, S in RES["runs"][run]["streams"].items():
        b = S["budget_rel"]; K = S["validation"]
        print(run, st, f"f={S['f_acc']:.4f} ±{S['f_acc_err_stat']:.4f} +{S['f_acc_err_syst_up']:.4f} −{S['f_acc_err_syst_down']:.4f} (bez lustra −{S['f_acc_err_syst_down_no_mirror']:.4f}) far {S['f_acc_far_norm']:.4f}",
              "K1", K["K1"], f"{K['norm>0.10_pred_0.045-0.10']['ratio']:.3f}±{K['norm>0.10_pred_0.045-0.10']['ratio_err']:.3f} z={K['norm>0.10_pred_0.045-0.10']['z']:.1f}", "K4", K["K4"]["status"], f"p={K['K4']['p']:.2g} dev={K['K4']['max_grouped_dev']:.3f}", "K5", K["K5"]["status"], "K6", S["K6_usable"])
        print("    ", {k: (round(v, 4) if isinstance(v, float) else v) for k, v in b.items()})
    print(run, "K2", RES["runs"][run]["K2"])
print("K3", RES["K3_toy"]["status"], RES["K3_toy"]["nominal_bias"], RES["K3_toy"]["nominal_bias_err"])
