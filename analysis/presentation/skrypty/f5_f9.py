#!/usr/bin/env python3
"""F5 izolacja (B_trig btrig_iso.json), F6 tryger (REG428 summary.json), F7 budzet (XSEC_final v2 final.json), F8 sklad tla (BKGCAT comp.json), F9 siec (DNET lxplus6; od 20.11 takze zysk z fitu piku G_fit z dnet6_R.json)."""
import math

import numpy as np

import pstyle as P
from pstyle import plt

q = lambda *a: math.sqrt(sum(x * x for x in a))  # noqa: E731
XB = {"3-5": 0, "5-8": 1, "8-30": 2}
SPL = {"D0": r"$D^0$", "Dplus": r"$D^+$"}


def binax(ax, bins=P.BINS, lab=r"$p_T(D)$ bin [GeV]"):
    ax.set_xticks(range(len(bins))); ax.set_xticklabels([P.BINLAB.get(b, b.replace("-", "–")) for b in bins]); ax.set_xlabel(lab)
    ax.tick_params(axis="x", which="both", top=False, bottom=False)
    for xs in np.arange(0.5, len(bins) - 1):
        ax.axvline(xs, color=P.GRID, lw=0.8)
    ax.set_xlim(-0.55, len(bins) - 0.45)


# ---------------------------------------------------------------------------------------------- F5
def f5():
    rows = P.J("B_trig", "runs", "lxplus", "out", "btrig_iso.json")["rows"]
    get = lambda run, cls, b: next(r for r in rows if r["run"] == run and r["cls"] == cls and r["bin"] == b and r["cut"] == "ptcone40<2.0")  # noqa: E731
    bins = ["3-30", "3-5", "5-8", "8-30"]
    ser = [(427929, "MB427929", "Run 427929, minimum bias (MBTS)", P.C["blue"], "o", -0.2),
           (428770, "all", "Run 428770, all triggers (jet-trigger dominated)", P.C["orange"], "s", 0.0),
           (428770, "J12only", "Run 428770, L1_J12-only class", P.C["violet"], "D", 0.2)]
    fig, ax = plt.subplots(figsize=(9.6, 6.3), gridspec_kw=dict(left=0.09, right=0.98, top=0.80, bottom=0.11))
    ax.axhline(1, color=P.INK2, lw=0.9)
    for run, cls, lab, col, mk, dx in ser:
        for i, b in enumerate(bins):
            r = get(run, cls, b)
            if not r["determined"]:
                continue
            ax.errorbar(i + dx, r["SF"], r["SF_err"], fmt=mk, ms=7.5, color=col, ecolor=col, elinewidth=1.8, label=lab if i == 0 else None, zorder=4)
            ax.annotate(f"{r['SF']:.2f}\n±{r['SF_err']:.2f}", (i + dx, r["SF"] - r["SF_err"] - 0.035), ha="center", va="top", fontsize=8.8, linespacing=1.05)
    a, b_ = get(427929, "MB427929", "3-30"), get(428770, "all", "3-30")
    ax.axvspan(-0.55, 0.5, color=P.GRID, alpha=0.35, lw=0, zorder=0)
    binax(ax, bins)
    ax.set_ylim(0.0, 1.45); ax.set_ylabel("Isolation scale factor  ε(data) / ε(MC)")
    ax.legend(loc="upper right", fontsize=9.5)
    ax.set_title(r"$p_T^{cone40}$ < 2 GeV, true $D^0$ from the $D^*$ tag (pass/fail $\Delta m$ fit)", loc="left")
    ax.text(0.015, 0.97, f"3–30 GeV:  ε(data) = {a['eps_data']:.2f} vs ε(MC) = {a['eps_mc']:.2f} in 427929;\n"
            f"ε(data) = {b_['eps_data']:.2f} vs ε(MC) = {b_['eps_mc']:.2f} in 428770.\n"
            f"Effect: with the minimum-bias MC efficiency the $D^0$ rate after\nthe isolation cut in 428770 would be low by ×{b_['SF']:.2f}.",
            transform=ax.transAxes, ha="left", va="top", fontsize=9.5, color=P.INK)
    P.header(fig, "runs 427929 and 428770, statistical errors only")
    P.save(fig, "F5_isolation_sf")


# ---------------------------------------------------------------------------------------------- F6
def f6():
    S = P.J("REG428_2026-11-20", "out", "summary.json")
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(13.4, 6.0), gridspec_kw=dict(left=0.065, right=0.985, top=0.80, bottom=0.11, wspace=0.2, width_ratios=[1, 1.25]))
    for sp, col, mk, dx in (("D0", P.C["blue"], "o", -0.12), ("Dplus", P.C["orange"], "s", 0.12)):
        for bn in P.BINS:
            e = S["eps"][f"{sp}|{bn}"]
            a0.errorbar(XB[bn] + dx, 100 * e["eps"], 100 * e["err"], fmt=mk, ms=7.5, color=col, ecolor=col, elinewidth=1.8, label=SPL[sp] if bn == "3-5" else None)
            a0.annotate(f"{100 * e['eps']:.1f} ± {100 * e['err']:.1f}", (XB[bn] + dx, 100 * (e["eps"] + e["err"]) + 0.5), ha="center", va="bottom", fontsize=9, rotation=90)
    binax(a0); a0.set_ylim(0, 24); a0.set_ylabel(r"$\varepsilon_{J12}$ for events with a D at WP70 [%]")
    a0.set_title("Efficiency of HLT_noalg_L1J12, measured in run 427929", loc="left", fontsize=11.5)
    a0.legend(loc="upper left", title="N(MB ∧ J12) / N(MB), peak fits; HESSE errors", title_fontsize=9)
    regs = [("R3x", "L1_J12 only (R3x)", P.C["blue"]), ("rest", "none of the three chains", P.C["yellow"]), ("R2", "j20_L1AFP_A_OR_C (R2)", P.C["aqua"]), ("R1", "noalg_L1AFP_A_OR_C (R1)", P.C["magenta"])]
    xs, labs = [], []
    for isp, sp in enumerate(("D0", "Dplus")):
        for ib, bn in enumerate(P.BINS):
            x = isp * 3.6 + ib
            xs.append(x); labs.append(P.BINLAB[bn])
            tot = S["yields"][f"{sp}|{bn}|allGRL"]["N"]
            b0 = 0.0
            for rk, rl, col in regs:
                f = 100 * max(S["yields"][f"{sp}|{bn}|{rk}"]["N"], 0) / tot
                a1.bar(x, f, 0.72, bottom=b0, color=col, edgecolor="white", lw=1.2, label=rl if (isp == 0 and ib == 0) else None)
                if f > 6:
                    a1.text(x, b0 + f / 2, f"{f:.0f} %", ha="center", va="center", fontsize=9.5, color="white" if rk == "R3x" else P.INK)
                b0 += f
            a1.text(x, max(b0, 100) + 2.5, f"Σ {b0:.0f} %", ha="center", va="bottom", fontsize=8.8, color=P.INK2)
        a1.text(isp * 3.6 + 1, -13.5, SPL[sp] + r",  $p_T$ bin [GeV]", ha="center", va="top", fontsize=11)
    a1.axhline(100, color=P.INK2, lw=0.8, ls=":")
    a1.set_xticks(xs); a1.set_xticklabels(labs); a1.tick_params(axis="x", which="both", top=False, bottom=False)
    a1.set_ylim(0, 138); a1.set_yticks([0, 25, 50, 75, 100]); a1.set_ylabel("Share of the WP70 yield of run 428770 [%]")
    a1.set_title("Disjoint trigger regions of run 428770 (fit per region / fit of all events)", loc="left", fontsize=11.5)
    a1.legend(loc="upper center", ncol=2, fontsize=9.3, columnspacing=1.3)
    P.header(fig, "runs 427929 (left; minimum bias) and 428770 (right)", extra="tree H at WP70; statistical errors only")
    P.save(fig, "F6_trigger")


# ---------------------------------------------------------------------------------------------- F7
def f7():
    R = P.J("XSEC_final_2026-11-20", "v2", "out", "final.json")
    grp = [("Statistical", ["stat"], P.INK2), ("Fit bias (MC closure, not corrected)", ["fit_bias"], P.C["orange"]), ("Fit model", ["fit_model"], P.C["yellow"]),
           ("Peaking template, second minimum", ["peak_tmpl", "peak_norm", "fit_altmin"], P.C["magenta"]), ("Classifier scale factor", ["SF"], P.C["blue"]),
           ("Tracking, impact parameter, vertex (literature)", ["tracking", "ip_material", "vertex3"], P.C["aqua"]),
           ("Other (MC stat, BR, non-prompt, generator, selection, reflection)", ["mc_stat", "BR", "refl", "nonprompt", "evsel", "generator"], P.C["violet"]),
           ("Luminosity", ["lumi"], P.C["green"])]
    fig, ax = plt.subplots(figsize=(12.4, 6.9), gridspec_kw=dict(left=0.075, right=0.985, top=0.745, bottom=0.15))
    xs, labs = [], []
    for isp, sp in enumerate(("D0", "Dplus")):
        for ib, bn in enumerate(P.BINS):
            c = R["cells"][f"{sp}|{bn}"]
            rel = dict(c["rel"]); rel["lumi"] = R["L_rel"]
            used = set()
            x = isp * 3.6 + ib
            xs.append(x); labs.append(P.BINLAB[bn])
            tot2 = sum(v * v for v in rel.values())
            assert abs(math.sqrt(tot2) - c["tot_rel"]) < 2e-4, (sp, bn, math.sqrt(tot2), c["tot_rel"])
            b0 = 0.0
            for gl, keys, col in grp:
                v = q(*[rel.get(k, 0.0) for k in keys]); used |= set(keys)
                h = 100 * v * v / tot2
                ax.bar(x, h, 0.74, bottom=b0, color=col, edgecolor="white", lw=1.2, label=gl if (isp == 0 and ib == 0) else None)
                if h > 4.5:
                    ax.text(x, b0 + h / 2, f"{100 * v:.1f} %", ha="center", va="center", fontsize=9.3, color="white" if col in (P.INK2, P.C["blue"], P.C["violet"], P.C["green"]) else P.INK)
                b0 += h
            assert set(rel) <= used, set(rel) - used
            ax.text(x, 102, f"total {100 * c['tot_rel']:.1f} %\nstat {100 * rel['stat']:.1f} %\nsyst {100 * c['syst_rel']:.1f} %", ha="center", va="bottom", fontsize=9.3, linespacing=1.2)
        ax.text(isp * 3.6 + 1, -9.5, {"D0": r"$\sigma(D^0)$", "Dplus": r"$\sigma(D^+)$"}[sp] + r",  $p_T$ bin [GeV]", ha="center", va="top", fontsize=11.5)
    ax.set_xticks(xs); ax.set_xticklabels(labs); ax.tick_params(axis="x", which="both", top=False, bottom=False)
    ax.set_ylim(0, 122); ax.set_yticks([0, 20, 40, 60, 80, 100]); ax.set_ylabel("Share of the total variance [%]")
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    ax.tick_params(top=False, right=False)
    fig.legend(*ax.get_legend_handles_labels(), loc="upper left", bbox_to_anchor=(0.07, 0.915), ncol=2, fontsize=9.5, columnspacing=1.6)
    P.header(fig, P.run427_text(), extra="uncertainty budget of the inclusive cross-sections")
    fig.text(0.985, 0.015, "Bar segments: quadratic share of each group in the total uncertainty; number inside a segment: its relative uncertainty (shown if the share exceeds 4.5 %).",
             ha="right", va="bottom", fontsize=9, color=P.INK2)
    P.save(fig, "F7_uncertainty_budget")


# ---------------------------------------------------------------------------------------------- F8
def f8():
    CO = P.J("BKGCAT_2026-11-20", "out", "comp.json")
    cats = [("K↔π reflection (same $D^0$)", ["Refl"], P.C["yellow"]), ("Other decay mode of the same D species", ["DMM"], P.C["aqua"]), ("Other charm hadron", ["CM"], P.C["magenta"]),
            ("Partially combinatorial (some tracks from one charm hadron)", ["PC"], P.C["violet"]), ("Combinatorial with b-hadron tracks", ["MB", "OtherB"], P.C["orange"]),
            ("Combinatorial, no charm (incl. light resonances)", ["NoCharm", "ML"], P.C["blue"])]
    fig, ax = plt.subplots(figsize=(12.4, 6.7), gridspec_kw=dict(left=0.075, right=0.985, top=0.745, bottom=0.15))
    xs, labs = [], []
    for isp, sp in enumerate(("D0", "Dplus")):
        for ib, bn in enumerate(P.BINS):
            c = next(x for x in CO if x["sp"] == sp and x["bin"] == bn and x["cat"] == "WP70")["comp"]
            nw = lambda k: c[k]["fit"]["n_win"] if k in c else 0.0  # noqa: E731
            sig, allw = nw("Sig"), nw("all")
            bkg = allw - sig
            parts = [sum(nw(k) for k in keys) for _, keys, _ in cats]
            assert abs(sum(parts) - bkg) < 0.5, (sp, bn, sum(parts), bkg)
            x = isp * 3.6 + ib
            xs.append(x); labs.append(P.BINLAB[bn])
            b0 = 0.0
            for (cl, keys, col), n in zip(cats, parts):
                h = 100 * n / bkg
                ax.bar(x, h, 0.74, bottom=b0, color=col, edgecolor="white", lw=1.2, label=cl if (isp == 0 and ib == 0) else None)
                if h > 4.5:
                    ax.text(x, b0 + h / 2, f"{h:.0f} %", ha="center", va="center", fontsize=9.5, color="white" if col in (P.C["blue"], P.C["violet"]) else P.INK)
                b0 += h
            ax.text(x, 102, f"B/S = {bkg / sig:.2f}\nS = {sig:.0f}", ha="center", va="bottom", fontsize=9.3, linespacing=1.2)
        ax.text(isp * 3.6 + 1, -9.5, SPL[sp] + r",  $p_T$ bin [GeV]", ha="center", va="top", fontsize=11.5)
    ax.set_xticks(xs); ax.set_xticklabels(labs); ax.tick_params(axis="x", which="both", top=False, bottom=False, right=False)
    ax.set_ylim(0, 118); ax.set_yticks([0, 20, 40, 60, 80, 100]); ax.set_ylabel("Share of the background in the mass window [%]")
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    fig.legend(*ax.get_legend_handles_labels(), loc="upper left", bbox_to_anchor=(0.07, 0.915), ncol=2, fontsize=9.5, columnspacing=1.6)
    fig.text(0.01, 0.985, P.LABEL, ha="left", va="top", fontsize=10.5, fontweight="bold")
    fig.text(0.01, 0.985 - 0.03 * 6.0 / 6.7, r"$\sqrt{s}$ = 13.6 TeV, minimum-bias MC (Pythia 8 A3, non-diffractive; unweighted MC candidates, no luminosity), tree H at WP70",
             ha="left", va="top", fontsize=10, color=P.INK2)
    fig.text(0.985, 0.015, r"Mass window: |m − m$_0$| < 55 MeV ($D^0$), 45 MeV ($D^+$). S = true signal candidates in the window; B = everything else.", ha="right", va="bottom", fontsize=9, color=P.INK2)
    P.save(fig, "F8_background_composition")


# ---------------------------------------------------------------------------------------------- F9
def f9():
    P1 = P.J("DNET", "runs", "lxplus6", "dnet6_p1.json")
    RS = P.J("DNET", "runs", "lxplus6", "dnet6_results.json")
    C = RS["species"]["D0"]["cfg"]
    RJ = P.J("DNET", "runs", "lxplus6", "dnet6_R.json")                      # zysk realny z fitu piku (krok anar; jedno ziarno, G_fit bez bledu)
    GTHR, GNEED = 1.15, 2                                                    # prog etapu 4 (KRYTERIA_etap6 §5): G_fit >= 1.15 w >= 2 z 3 binow
    RD = RJ["data"]["D0"]["D0.tr.H+e04np.mc.d5.s1"]
    RKEY = {"sb427929": "427929|MB427929", "sb428770": "428770|all"}
    assert RJ["blocked"] == [] and all(RD[k][bn]["0.7"]["determined"] for k in RKEY.values() for bn in P.BINS)
    gate = [RJ["closure"]["D0"]["D0.tr.H+e04np.mc.d5.s1"][bn]["0.7"]["pass"] for bn in P.BINS]
    met = {nm: sum(RD[k][bn]["0.7"]["G_fit"] >= GTHR for bn in P.BINS) >= GNEED for nm, k in RKEY.items()}
    assert met == {"sb427929": False, "sb428770": True} and not any(gate), (met, gate)   # tekst na rysunku zaklada ten odczyt (WERDYKT_etap6 §5)
    r428 = [RD["428770|all"][bn]["0.7"]["R"] for bn in P.BINS]
    full = [f"D0.tr.H+e04.mc.d5.s{s}" for s in (1, 2, 3)]
    mask = [f"D0.tr.H+e04np.mc.d5.s{s}" for s in (1, 2, 3)]
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(13.4, 6.1), gridspec_kw=dict(left=0.065, right=0.985, top=0.80, bottom=0.11, wspace=0.2, width_ratios=[1, 1.35]))
    a0.axhspan(0.97, 1.03, color=P.GRID, alpha=0.7, lw=0, label="Criterion |DR − 1| ≤ 0.03")
    a0.axhline(1, color=P.INK2, lw=0.9)
    for tag, lab, col, mk, dx in ((full[0], "Network with surrounding tracks, no mask", P.C["orange"], "s", -0.1), (mask[0], r"Same, soft-pion candidates masked ($\Delta m$ < 155 MeV)", P.C["blue"], "o", 0.1)):
        for bn in P.BINS:
            dr = P1["models"][tag]["bins"][bn]["0.7"]["double_ratio"]
            a0.errorbar(XB[bn] + dx, dr[0], dr[1], fmt=mk, ms=7.5, color=col, ecolor=col, elinewidth=1.8, label=lab if bn == "3-5" else None)
            a0.annotate(f"{dr[0]:.3f} ± {dr[1]:.3f}", (XB[bn] + dx, dr[0] + dr[1] + 0.012), ha="center", va="bottom", fontsize=9)
    binax(a0); a0.set_ylim(0.9, 1.72)
    a0.set_ylabel(r"DR = [ε($D^0$ from $D^*$) / ε(other $D^0$)]$_{net}$ / [same]$_{tree}$")
    a0.set_title("Soft-pion effect (MC signal, WP70, seed 1)", loc="left", fontsize=11.5)
    a0.legend(loc="upper left", fontsize=9.3)
    src = [("sb427929", "Run 427929"), ("sb428770", "Run 428770"), ("mc", "MC (true background)")]
    xs, labs = [], []
    for i, (nm, ttl) in enumerate(src):
        for ib, bn in enumerate(P.BINS):
            x = i * 3.5 + ib
            xs.append(x); labs.append(P.BINLAB[bn])
            g = [C[t][f"{bn}|0.7"]["G_" + nm] for t in mask]
            ge = C[mask[0]][f"{bn}|0.7"]["G_" + nm + "_err"]
            gf = float(np.mean([C[t][f"{bn}|0.7"]["G_" + nm] for t in full]))
            gm = float(np.mean(g))
            a1.plot(x - 0.3, gf, "s", ms=7.5, color=P.C["orange"], label="Proxy Z$_{net}$/Z$_{tree}$ (data: sidebands), no mask (mean of 3 seeds)" if (i == 0 and ib == 0) else None)
            a1.errorbar(x, gm, ge, fmt="o", ms=7.5, color=P.C["blue"], ecolor=P.C["blue"], elinewidth=1.8, label="Proxy, soft pion masked (mean of 3 seeds, error of seed 1)" if (i == 0 and ib == 0) else None)
            a1.plot([x] * 3, g, "_", ms=11, color=P.INK, mew=1.1, label="Individual seeds" if (i == 0 and ib == 0) else None)
            a1.annotate(f"{gm:.2f}", (x, min(g) - max(ge, 0.0) - 0.02), ha="center", va="top", fontsize=9)
            a1.annotate(f"{gf:.2f}", (x - 0.3, gf + 0.025), ha="center", va="bottom", fontsize=9)
            if nm in RKEY:
                gfit = RD[RKEY[nm]][bn]["0.7"]["G_fit"]
                a1.plot(x + 0.3, gfit, "D", ms=7.5, color=P.C["aqua"], mec=P.INK, mew=0.6, label="Peak fit G$_{fit}$ = ratio of N$_S$/σ$_N$, soft pion masked (seed 1; no error defined)" if (i == 0 and ib == 0) else None)
                a1.annotate(f"{gfit:.2f}", (x + 0.4, gfit), ha="left", va="center", fontsize=9, color=P.C["green"], fontweight="bold")
        if nm in RKEY:
            a1.plot([i * 3.5 - 0.55, i * 3.5 + 2.55], [GTHR, GTHR], ls=(0, (4, 3)), lw=1.0, color=P.C["green"], label=f"Pre-set threshold for G$_{{fit}}$: {GTHR:.2f} in ≥ {GNEED} of 3 bins" if i == 0 else None)
            a1.text(i * 3.5 + 1, 0.795, "G$_{fit}$ threshold " + ("met" if met[nm] else "not met"), ha="center", va="bottom", fontsize=9.5, color=P.C["green"])
        else:
            a1.text(i * 3.5 + 1, 0.795, "(no peak fit)", ha="center", va="bottom", fontsize=9.5, color=P.INK2)
        a1.text(i * 3.5 + 1, 0.845, ttl, ha="center", va="bottom", fontsize=10)
        if i:
            a1.axvline(i * 3.5 - 0.75, color=P.GRID, lw=0.8)
    a1.axhline(1, color=P.INK2, lw=0.9)
    a1.text(xs[-1] + 0.5, 1.008, "tree H", ha="right", va="bottom", fontsize=9, color=P.INK2)
    a1.set_xticks(xs); a1.set_xticklabels(labs); a1.tick_params(axis="x", which="both", top=False, bottom=False)
    a1.set_xlim(-0.7, xs[-1] + 0.7); a1.set_ylim(0.78, 2.2)
    a1.set_xlabel(r"$p_T(D^0)$ bin [GeV]"); a1.set_ylabel("Gain against tree H at WP70 (sideband proxy or peak fit)")
    a1.set_title("Gain against the tree: sideband proxy and inclusive peak fit", loc="left", fontsize=11.5)
    a1.legend(loc="upper left", fontsize=9.3)
    P.header(fig, r"$D^0 \to K\pi$, runs 427929 and 428770 and MC")
    fig.text(0.985, 0.985, "The sideband proxy overestimates the gain (ε$_S$ from MC). Peak fit: inclusive $D^0$ peak in data. Where its gain passes\n"
             f"(run 428770) the network yield after the MC efficiency is {100 * (1 - max(r428)):.0f}–{100 * (1 - min(r428)):.0f} % below that of the tree; the MC closure of this ratio\n"
             "is not passed in any bin, and the $D^*$-tag scale factor of the masked network is not stable between seeds.",
             ha="right", va="top", fontsize=9, color=P.INK2)
    P.save(fig, "F9_network")


if __name__ == "__main__":
    import sys
    for f in (f5, f6, f7, f8, f9):
        if len(sys.argv) < 2 or f.__name__ in sys.argv[1:]:
            f()
