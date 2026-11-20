#!/usr/bin/env python3
"""F2 — przekroje per bin pT (D0, D+, D*±) z odniesieniami; F3 — stosunki D+/D0 i D*+/D0. Jedyne zrodlo liczb: XSEC_final_2026-11-20/v2/out/final.json."""
import math

import numpy as np

import pstyle as P
from pstyle import plt

R = P.J("XSEC_final_2026-11-20", "v2", "out", "final.json")
XB = {"3-5": 0, "5-8": 1, "8-30": 2}
q = lambda *a: math.sqrt(sum(x * x for x in a))  # noqa: E731


def databar(ax, x, y, stat, tot, **kw):
    ax.errorbar(x, y, tot, fmt="none", ecolor=P.DATA, elinewidth=1.0, capsize=0, zorder=6)
    ax.errorbar(x, y, stat, fmt="o", ms=6, color=P.DATA, ecolor=P.DATA, elinewidth=2.6, capsize=0, zorder=7, **kw)


def f2():
    fig = plt.figure(figsize=(14.2, 8.4))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.55, 1.0], hspace=0.06, wspace=0.24, left=0.065, right=0.985, top=0.80, bottom=0.085)
    sp_def = [("D0", r"$D^0 + \bar{D}^0$"), ("Dplus", r"$D^+ + D^-$"), ("Dstar", r"$D^{*+} + D^{*-}$")]
    off = {"pythia": -0.30, "fonll": -0.15, "cms": 0.0, "atlas7": 0.15, "alice": 0.30}
    sty = {"pythia": (P.C["blue"], "s", "Pythia 8 A3 (same fiducial definition)"),
           "fonll": (P.C["aqua"], "D", "FONLL 13.6 TeV, prompt only (band: scales, masses, PDF)"),
           "cms": (P.C["orange"], "^", "CMS 13 TeV (data converted to the CMS definition)"),
           "atlas7": (P.C["violet"], "v", "ATLAS 7 TeV × FONLL energy factor"),
           "alice": (P.C["magenta"], "o", "ALICE 13 TeV — post hoc, unverified reference values")}
    hand = {}
    for ic, (sp, lab) in enumerate(sp_def):
        a0, a1 = fig.add_subplot(gs[0, ic]), fig.add_subplot(gs[1, ic])
        for bn in P.BINS:
            x = XB[bn]
            if sp == "Dstar":
                c = R["dstar"].get(bn)
                if c is None:
                    a0.text(x, 0.5, "not\nmeasured", transform=a0.get_xaxis_transform(), ha="center", va="center", fontsize=9.5, color=P.INK2)
                    continue
                y, st, sy, lu = c["sigma_ub"], c["stat_ub"], c["syst_ub"], c["lumi_ub"]
                cmp_ = {k: c[k] for k in ("pythia", "cms", "atlas7") if k in c}
                py = c["pythia"]["ref"]
            else:
                c = R["cells"][f"{sp}|{bn}"]
                y, st, sy, lu = c["sigma_ub"], c["stat_ub"], c["syst_ub"], c["lumi_ub"]
                cmp_ = R["cmp"][f"{sp}|{bn}"]
                py = c["pythia_ub"]
            tot = q(st, sy, lu)
            # gora
            hand["pythia"] = a0.hlines(py, x - 0.42, x + 0.42, color=P.C["blue"], lw=2.2, zorder=4)
            if "fonll" in cmp_:
                fo = cmp_["fonll"]
                a0.fill_between([x - 0.42, x + 0.42], fo["min"], fo["max"], color=P.C["aqua"], alpha=0.30, lw=0, zorder=2)
                a0.hlines(fo["centr"], x - 0.42, x + 0.42, color=P.C["aqua"], lw=1.4, ls="--", zorder=3)
            if "atlas7" in cmp_:
                a7 = cmp_["atlas7"]
                a0.errorbar(x + 0.2, a7["scaled"], a7["scaled_err"], fmt="v", ms=6.5, mfc="white", mec=P.C["violet"], ecolor=P.C["violet"], mew=1.6, elinewidth=1.4, zorder=5)
            databar(a0, x - (0.1 if "atlas7" in cmp_ else 0), y, st, tot)
            a0.annotate(f"{y:.1f} ± {st:.1f}\n± {sy:.1f} ± {lu:.1f} µb" if y < 99.5 else f"{y:.0f} ± {st:.0f}\n± {sy:.0f} ± {lu:.0f} µb", (x, 0.045), xycoords=a0.get_xaxis_transform(),
                        ha="center", va="bottom", fontsize=9.5, color=P.INK, linespacing=1.25)
            # dol: dane / odniesienie
            for k in ("pythia", "fonll", "cms", "atlas7", "alice"):
                if k not in cmp_:
                    continue
                d = cmp_[k]
                col, mk, _ = sty[k]
                xx = x + off[k]
                if k == "fonll":
                    r, e = d["ratio_centr"], d["ratio_err"]
                    a1.fill_between([xx - 0.055, xx + 0.055], d["prompt_data"] / d["max"], d["prompt_data"] / d["min"], color=col, alpha=0.28, lw=0, zorder=2)
                else:
                    r, e = d["ratio"], d["ratio_err"]
                hand[k + "_r"] = a1.errorbar(xx, r, e, fmt=mk, ms=6.5, color=col, mfc="white" if k == "alice" else col, mec=col, mew=1.5, ecolor=col, elinewidth=1.5, zorder=5)
        a0.set_yscale("log")
        a0.set_xlim(-0.55, 2.55); a1.set_xlim(-0.55, 2.55)
        lo, hi = a0.get_ylim()
        a0.set_ylim(lo / 5.0, hi * 1.6)
        from matplotlib.ticker import NullFormatter
        a0.yaxis.set_minor_formatter(NullFormatter())
        a0.set_xticks([0, 1, 2]); a0.tick_params(labelbottom=False, which="both")
        a0.tick_params(axis="x", which="both", top=False, bottom=False)
        a0.set_title(lab + r",  $|\eta(D)| < 2.2$, prompt + non-prompt", loc="left", fontsize=11.5)
        if ic == 0:
            a0.set_ylabel(r"$\sigma$ in the $p_T$ bin [µb]")
            a1.set_ylabel("Data / reference")
        a1.axhline(1, color=P.INK2, lw=0.8)
        for xs in (0.5, 1.5):
            a0.axvline(xs, color=P.GRID, lw=0.8); a1.axvline(xs, color=P.GRID, lw=0.8)
        a1.set_ylim(0.0, 4.8 if sp != "Dstar" else 2.0)
        if sp != "Dstar":
            a1.set_yticks([0, 1, 2, 3, 4])
        a1.set_xticks([0, 1, 2]); a1.set_xticklabels([P.BINLAB[b] for b in P.BINS])
        a1.tick_params(axis="x", which="both", top=False, bottom=False)
        a1.set_xlabel(r"$p_T(D)$ bin [GeV]")
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    h = [Line2D([], [], marker="o", color=P.DATA, ms=6, lw=2.6, label="Data: thick bar = stat, thin bar = stat ⊕ syst ⊕ lumi (numbers: σ ± stat ± syst ± lumi)"),
         Line2D([], [], color=P.C["blue"], lw=2.2, marker="s", ms=6, label=sty["pythia"][2]),
         (Patch(color=P.C["aqua"], alpha=0.3), Line2D([], [], color=P.C["aqua"], lw=1.4, ls="--", marker="D", ms=5.5)),
         Line2D([], [], color=P.C["orange"], lw=0, marker="^", ms=6.5, label=sty["cms"][2]),
         Line2D([], [], color=P.C["violet"], lw=0, marker="v", ms=6.5, mfc="white", mew=1.6, label=sty["atlas7"][2]),
         Line2D([], [], color=P.C["magenta"], lw=0, marker="o", ms=6.5, mfc="white", mew=1.5, label=sty["alice"][2])]
    lab = [h[0].get_label(), sty["pythia"][2], sty["fonll"][2], sty["cms"][2], sty["atlas7"][2], sty["alice"][2]]
    fig.legend(h, lab, loc="upper left", bbox_to_anchor=(0.06, 0.925), ncol=2, fontsize=9.5, columnspacing=2.0, handlelength=2.0)
    P.header(fig, P.run427_text(), extra="no proton tag")
    fig.text(0.985, 0.985, "Lower panels: data in the definition of each reference over that reference;\nFONLL: prompt part of the data over the central value, box = FONLL band.",
             ha="right", va="top", fontsize=9, color=P.INK2)
    P.save(fig, "F2_cross_sections")


def f3():
    fig, axs = plt.subplots(1, 2, figsize=(12.6, 5.9), gridspec_kw=dict(left=0.07, right=0.985, top=0.80, bottom=0.11, wspace=0.2))
    refsty = {"Pythia A3": (P.C["blue"], "s", "Pythia 8 A3"), "ALICE": (P.C["magenta"], "o", "ALICE 13 TeV, |y| < 0.5, prompt"), "LHCb": (P.C["aqua"], "D", "LHCb 13 TeV"),
              "e+e-": (P.C["yellow"], "P", r"$e^+e^-$"), "CMS": (P.C["orange"], "^", "CMS 13 TeV"), "prawda mb_nd": (P.C["blue"], "s", "Minimum-bias MC truth (Pythia 8 A3)")}
    order = ["Pythia A3", "prawda mb_nd", "ALICE", "LHCb", "e+e-", "CMS"]
    for ax, (key, ttl) in zip(axs, (("ratio", r"$D^+ / D^0$"), ("dstar", r"$D^{*+} / D^0$"))):
        seen = {}
        for bn in P.BINS:
            x = XB[bn]
            c = R[key].get(bn)
            if c is None:
                ax.text(x, 0.5, "not measured", transform=ax.get_xaxis_transform(), ha="center", va="center", fontsize=10, color=P.INK2)
                continue
            y, st, sy = c["R"], c["stat"], c["syst"]
            refs = [k for k in order if k in c["cmp"]]
            for i, k in enumerate(refs):
                d = c["cmp"][k]
                col, mk, lb = refsty[k]
                xx = x + 0.12 + 0.085 * i
                seen[lb] = ax.errorbar(xx, d["ref"], d["ref_err"] if d["ref_err"] > 0 else None, fmt=mk, ms=7, color=col, ecolor=col, elinewidth=1.5, mec=col, zorder=4)
            databar(ax, x - 0.15, y, st, q(st, sy))
            ax.annotate(f"{y:.2f} ± {st:.2f} ± {sy:.2f}", (x, 0.04), xycoords=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=9.5)
        for xs in (0.5, 1.5):
            ax.axvline(xs, color=P.GRID, lw=0.8)
        ax.set_xlim(-0.55, 2.6); ax.set_ylim(0.0, 1.0)
        ax.set_xticks([0, 1, 2]); ax.set_xticklabels([P.BINLAB[b] for b in P.BINS])
        ax.tick_params(axis="x", which="both", top=False, bottom=False)
        ax.set_xlabel(r"$p_T(D)$ bin [GeV]"); ax.set_ylabel("Ratio of cross-sections")
        ax.set_title(ttl + r",  $|\eta(D)| < 2.2$", loc="left")
        from matplotlib.lines import Line2D
        hd = Line2D([], [], marker="o", color=P.DATA, ms=6, lw=2.6)
        ax.legend([hd] + list(seen.values()), ["Data (thick: stat, thin: stat ⊕ syst)"] + list(seen), loc="upper left", ncol=2, fontsize=9.3, columnspacing=1.2)
    P.header(fig, P.run427_text(), extra="no proton tag")
    fig.text(0.07, 0.875, "Numbers: ratio ± stat ± syst. Luminosity and event selection cancel; scale factors cancel only partly. References without error bars have no quoted uncertainty.",
             fontsize=9, color=P.INK2)
    P.save(fig, "F3_ratios")


f2(); f3()
