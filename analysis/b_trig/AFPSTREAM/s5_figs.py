#!/usr/bin/env python3
"""Rysunki: fig1_afp_classes.png (SF izolacji sygnału, kontrola na tle, ε_B per klasa), fig2_afp_activity.png (rozkłady aktywności)."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import afp_common as A

s1 = json.load(open(f"{A.OUT}/s1.json")); s3 = json.load(open(f"{A.OUT}/s3.json")); s4 = json.load(open(f"{A.OUT}/s4.json"))
BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"
INK, MUTED, GRID = "#1a1a19", "#6b6a63", "#e4e3dd"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": INK,
                     "axes.spines.top": False, "axes.spines.right": False})
ROWS = [("9:MB427929", "427929 MBTS_1 (minbias)"), ("9:MB9&j20&!J12", "427929 MBTS ∧ j20 ∧ ¬J12"), ("9:MB9&AFPnoalg", "427929 MBTS ∧ noalg_L1AFP"),
        ("8:AFPnoalg_any", "428770 noalg_L1AFP (bit, dowolne inne)"), ("8:AFPnoalg", "428770 AFP noalg-only (¬sptrk4AFP)"),
        ("8:AFPsptrk4", "428770 AFP sptrk_pt4-only"), ("8:AFP&J12", "428770 AFP ∧ J12"), ("8:AFPstream", "428770 strumień AFP (suma)"),
        ("8:none", "428770 brak zapisanego bitu"), ("8:MBTS!A", "428770 MBTS-based ¬AFP"), ("8:J12only!A", "428770 J12-only ¬AFP")]
CK = "ptcone40<2"
fig, ax = plt.subplots(1, 3, figsize=(15.5, 5.6), sharey=True, gridspec_kw={"wspace": 0.1})
y = np.arange(len(ROWS))[::-1]
for i, (k, lab) in enumerate(ROWS):
    c = s1["rows"][k + "|3-30"]["cuts"][CK]; S = s1["rows"][k + "|3-30"]["S"]
    weak = S < 120
    if "SF" in c and np.isfinite(c["SF"]):
        ax[0].errorbar(c["SF"], y[i] + 0.14, xerr=c["SF_err"], fmt="o", ms=7, color=BLUE, mfc="white" if weak else BLUE, capsize=2, lw=1.6)
        if "rw" in c:
            ax[0].errorbar(c["rw"]["SF"], y[i] - 0.14, xerr=c["rw"]["err"], fmt="s", ms=6, color=ORANGE, mfc="white" if weak else ORANGE, capsize=2, lw=1.6)
    ax[0].text(2.12, y[i], f"N={S:.0f}", va="center", ha="right", fontsize=8, color=MUTED)
    if "bkg" in c:
        b = c["bkg"]
        ax[1].errorbar(b["R"], y[i] + 0.14, xerr=b["R_err"], fmt="o", ms=7, color=BLUE, capsize=2, lw=1.6)
        ax[1].errorbar(b["R_rw"], y[i] - 0.14, xerr=b["R_rw_err"], fmt="s", ms=6, color=ORANGE, capsize=2, lw=1.6)
        ax[1].text(1.58, y[i], f"f={b['f']:.3f}", va="center", ha="right", fontsize=8, color=MUTED)
    for j, (b_, col, mk, dy) in enumerate((("3-5", BLUE, "o", 0.14), ("5-8", ORANGE, "s", -0.14))):
        r = s3["rows"].get(f"{k}|{b_}", {})
        if "R_mc_t2" in r and r["R_mc_t2_stat"][0] < 0.25:
            e = np.hypot(r["R_mc_t2_stat"][0], r["R_mc_t2"][0] * r["syst_rel"][0])
            ax[2].errorbar(r["R_mc_t2"][0], y[i] + dy, xerr=e, fmt=mk, ms=6.5, color=col, capsize=2, lw=1.6)
for a, t, xl, xr in ((ax[0], "(a) sygnał D*: SF izolacji\nptcone40 < 2 GeV, $p_T$ 3–30 GeV", "ε(dane) / ε(MC)", (-0.3, 2.15)),
                     (ax[1], "(b) tło (Δm > 150 MeV, bez score):\nfrakcja izolowana, dane / MC", "f(dane) / f(MC mb_nd)", (0.2, 1.6)),
                     (ax[2], "(c) ε$_B$ pod pikiem, WP95:\ndane / MC (truth ≥ 2)", "ε$_B$(dane) / ε$_B$(MC)", (0.6, 1.15))):
    a.axvline(1, color=MUTED, lw=1, ls="--"); a.set_title(t, fontsize=10, loc="left", color=INK); a.set_xlabel(xl); a.set_xlim(*xr)
    a.grid(axis="x", color=GRID, lw=0.8); a.set_axisbelow(True)
    for yy in (y[2] - 0.5, y[7] - 0.5):
        a.axhline(yy, color=GRID, lw=1)
ax[0].set_yticks(y); ax[0].set_yticklabels([l for _, l in ROWS]); ax[0].set_ylim(-0.6, len(ROWS) + 2.3)
h = [plt.Line2D([], [], marker="o", color=BLUE, ls="", ms=7, label="MC bez wag"), plt.Line2D([], [], marker="s", color=ORANGE, ls="", ms=6, label="MC przeważone $p_T$×lead $p_T$×nSelBL"),
     plt.Line2D([], [], marker="o", color=BLUE, mfc="white", ls="", ms=7, label="puste: N_sig(SBS) < 120 (słaba / niemiarodajna)")]
ax[0].legend(handles=h, loc="upper right", fontsize=7.5, frameon=False, bbox_to_anchor=(1.0, 0.90))
ax[1].legend(handles=h[:2], loc="upper left", fontsize=7.5, frameon=False, bbox_to_anchor=(0.0, 1.0))
ax[2].legend(handles=[plt.Line2D([], [], marker="o", color=BLUE, ls="", ms=6.5, label="$p_T$ 3–5 GeV"), plt.Line2D([], [], marker="s", color=ORANGE, ls="", ms=6.5, label="$p_T$ 5–8 GeV")],
             loc="upper left", fontsize=7.5, frameon=False, bbox_to_anchor=(0.0, 1.0))
A.atlas_label(ax[0], "$\\sqrt{s}$ = 13.6 TeV, runy 428770 / 427929; klasy z bitów trygera", x=0.02, y=0.985, fs=11)
fig.savefig(f"{A.HERE}/fig1_afp_classes.png", dpi=150, bbox_inches="tight")

# ---- fig 2
ed = s4["hist"]["edges"]
SER = [("bkg|427929", "427929 MBTS_1", BLUE, "-"), ("bkg|MC", "MC mb_nd", INK, ":"), ("bkg|J12only", "428770 J12-only", ORANGE, "-"),
       ("bkg|8:AFPnoalg_any", "428770 noalg_L1AFP (bit)", AQUA, "-"), ("bkg|8:AFPsptrk4", "428770 AFP sptrk_pt4-only", YELLOW, "-"), ("bkg|8:AFP&J12", "428770 AFP ∧ J12", MAGENTA, "--")]
VV = [("e_nSelBL", "nSelBL"), ("e_lead_pt", "lead track $p_T$ [GeV]"), ("ptcone40", "ptcone40 [GeV]")]
fig, ax = plt.subplots(2, 3, figsize=(15, 8.2), gridspec_kw={"hspace": 0.28, "wspace": 0.2})
for j, (v, xl) in enumerate(VV):
    e = np.array(ed[v]); w = np.diff(e)
    a = ax[0, j]
    for key, lab, col, ls in SER:
        hcnt = np.array(s4["hist"][key][v], float); d = hcnt / hcnt.sum() / w
        a.stairs(d, e, color=col, lw=2, ls=ls, label=lab)
    a.set_xlabel(xl); a.set_ylabel("ułamek kandydatów / jedn."); a.grid(color=GRID, lw=0.8); a.set_axisbelow(True)
    a.set_title(("(a) " if j == 0 else "") + "tło: pasma boczne m(Kπ), inkluzywne D$^0$, $p_T$ 3–30 GeV" if j == 0 else "", fontsize=10, loc="left")
    if v != "e_nSelBL":
        a.set_yscale("log")
    # sygnał: przebinowanie zgrubne
    a = ax[1, j]
    grp = {"e_nSelBL": 4, "e_lead_pt": 3, "ptcone40": 4}[v]
    idx = np.arange(0, len(e) - 1, grp); e2 = np.r_[e[idx], e[-1]]; w2 = np.diff(e2); c2 = 0.5 * (e2[1:] + e2[:-1])
    def reb(x):
        return np.add.reduceat(np.asarray(x, float), idx)
    for key, lab, col, ls in (("sig|MC", "MC fd_nd (prawdziwe D*)", INK, ":"), ("sig|9:MB427929", "427929 MBTS_1", BLUE, "-"), ("sig|8:J12only!A", "428770 J12-only", ORANGE, "-")):
        hh = s4["hist"][key][v]; hcnt = reb(hh["h"] if isinstance(hh, dict) else hh)
        a.stairs(hcnt / hcnt.sum() / w2, e2, color=col, lw=2, ls=ls, label=lab)
    hh = s4["hist"]["sig|8:AFPstream"][v]; hc, h2 = reb(hh["h"]), reb(hh["h2"])
    a.errorbar(c2, hc / hc.sum() / w2, yerr=np.sqrt(h2) / hc.sum() / w2, fmt="o", ms=7, color=GREEN, capsize=2, lw=1.6, label="428770 strumień AFP (N$_{sig}$ = 86 ± 12)")
    a.axhline(0, color=MUTED, lw=0.8)
    a.set_xlabel(xl); a.set_ylabel("ułamek sygnału / jedn."); a.grid(color=GRID, lw=0.8); a.set_axisbelow(True)
    if j == 0:
        a.set_title("(b) sygnał D*→D$^0$π (odjęte pasma Δm), score$_H$ > 0,6", fontsize=10, loc="left")
ax[0, 2].legend(fontsize=8.5, frameon=False, loc="upper right")
ax[1, 2].legend(fontsize=8.5, frameon=False, loc="upper right")
A.atlas_label(ax[0, 0], "$\\sqrt{s}$ = 13.6 TeV, runy 428770 / 427929", x=0.36, y=0.95, fs=11)
fig.savefig(f"{A.HERE}/fig2_afp_activity.png", dpi=150, bbox_inches="tight")
print("ok")
