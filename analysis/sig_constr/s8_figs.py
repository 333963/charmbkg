#!/usr/bin/env python3
"""SIG_constr — rysunki PNG z out/*.json."""
import json, os, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
COL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
INK, INK2, SURF = "#0b0b0b", "#52514e", "#fcfcfb"
plt.rcParams.update({"figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF, "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True, "grid.color": "#e4e3df", "grid.linewidth": 0.6, "axes.spines.top": False,
                     "axes.spines.right": False, "font.size": 10, "legend.frameon": False, "lines.linewidth": 2, "lines.markersize": 7})
J = lambda n: json.load(open(os.path.join(OUT, n))) if os.path.exists(os.path.join(OUT, n)) else None  # noqa: E731


def atlas(ax, extra=""):
    ax.text(0.0, 1.03, "ATLAS", transform=ax.transAxes, fontweight="bold", fontstyle="italic", fontsize=11, va="bottom")
    ax.text(0.115 if ax.get_position().width > 0.5 else 0.2, 1.03, (extra if extra else ""), transform=ax.transAxes, fontsize=10, va="bottom")


def save(fig, name):
    fig.savefig(os.path.join(HERE, name), dpi=140, bbox_inches="tight"); plt.close(fig); print("→", name)


# 1. n_eff
d = J("s1_neff.json")
if d:
    cls = ["9:MB", "8:J12only", "8:j20&!J12", "8:MBTSbased", "8:all"]
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.8), sharey=True)
    for ax, reg, ttl in zip(axs, ("PRE", "H06"), ("bez cięcia score (PRE)", "score_H > 0,6 (H06)")):
        x = np.arange(len(cls)); wd = 0.26
        for i, b in enumerate(("3-5", "5-8", "8-30")):
            v = [next((r["neff_sbs"] for r in d["rows"] if r["cls"] == c and r["reg"] == reg and r["bin"] == b), 0) for c in cls]
            ax.bar(x + (i - 1) * wd, v, wd - 0.03, color=COL[i], label=f"pT {b.replace('-', '–')} GeV")
            for xx, vv in zip(x + (i - 1) * wd, v):
                ax.text(xx, max(vv, 1) * 1.08, f"{vv:.0f}", ha="center", va="bottom", fontsize=7, color=INK2)
        ax.set_yscale("log"); ax.set_ylim(3, 6000); ax.set_xticks(x); ax.set_xticklabels(cls, fontsize=8); ax.set_title(ttl, fontsize=10, loc="right")
        atlas(ax)
    axs[0].set_ylabel("n_eff sygnału D⁰←D* (SBS)"); axs[0].legend(fontsize=8, loc="upper left")
    save(fig, "fig1_neff.png")

# 2. krzywa uczenia
d = J("s2_learn.json"); n1 = J("s1_neff.json")
if d:
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for i, (t, lab) in enumerate((("H", "zadanie H: 18 cech, PRE"), ("I", "zadanie I: 18 cech + izolacja + score, H06"))):
        s = d[t]["summary"]; N = [r["N"] for r in s]
        ax.errorbar(N, [r["Zrel0.8"] for r in s], [r["Zrel0.8_sd"] for r in s], color=COL[i], marker="o", label=lab, capsize=2)
    ax.axhline(0.95, color=INK2, lw=0.8, ls="--"); ax.text(110, 0.953, "−5 %", fontsize=8, color=INK2)
    if n1:
        for c, reg, y in (("9:MB", "H06", 0.74), ("8:J12only", "H06", 0.70), ("9:MB", "PRE", 0.66), ("8:J12only", "PRE", 0.62)):
            v = next((r["neff_sbs"] for r in n1["rows"] if r["cls"] == c and r["reg"] == reg and r["bin"] == "3-30"), None)
            if v:
                ax.annotate(f"n_eff {c} {reg}: {v:.0f}", (v, y), fontsize=8, color=INK2, ha="left", va="center", xytext=(4, 0), textcoords="offset points")
                ax.plot([v, v], [0.58, y], color=INK2, lw=0.8)
    ax.set_xscale("log"); ax.set_ylim(0.58, 1.03); ax.set_xlabel("N zdarzeń sygnału w treningu (MC)"); ax.set_ylabel("Z_rel = ε_S(N)/ε_S(pełny) przy ε_B z WP80")
    ax.legend(fontsize=8, loc="lower right"); atlas(ax, "krzywa uczenia, tło: pasma Δm 427929")
    save(fig, "fig2_krzywa_uczenia.png")

# 3. warianty na danych
rs = []
for f in sorted(glob.glob(os.path.join(OUT, "s3_variants_*.json"))):
    rs += json.load(open(f))
if rs:
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.9), sharey=True)
    for ax, task, ttl in zip(axs, ("H", "I"), ("zadanie H (PRE, 18 cech, pT 5–30)", "zadanie I (H06, + izolacja, pT 3–30)")):
        mods = ["RW", "DATA", "STACK"] + (["H"] if task == "I" else [])
        labs = {"RW": "MC przeważone w(x)", "DATA": "od zera na danych (CWoLa)", "STACK": "dostrojenie (dane + margines MC)", "H": "sam score H (bez izolacji)"}
        rr = [r for r in rs if r["task"] == task]
        x = np.arange(len(rr))
        for i, m in enumerate(mods):
            y = [r["wp"]["0.8"]["all"][m]["Zrel_vs_MC"] for r in rr]; e = [r["wp"]["0.8"]["all"][m]["Zrel_vs_MC_err"] for r in rr]
            ax.errorbar(x + (i - (len(mods) - 1) / 2) * 0.16, y, e, fmt="o", color=COL[i], label=labs[m], capsize=2)
        ax.axhline(1, color=INK2, lw=0.8); ax.axhline(1.05, color=INK2, lw=0.8, ls="--")
        ax.set_xticks(x); ax.set_xticklabels([r["cls"] for r in rr]); ax.set_xlim(-0.6, len(rr) - 0.4); ax.set_title(ttl, fontsize=10, loc="right"); atlas(ax)
    axs[0].set_ylabel("Z_rel względem modelu MC (dane OOF, ε_B z WP80)"); axs[1].legend(fontsize=8, loc="lower left")
    save(fig, "fig3_warianty.png")

# 4. domknięcie
d = J("s4_closure.json")
if d and d["iso"]:
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.0))
    ax = axs[0]
    rr = [r for r in d["iso"] if r["bin"] == "3-30" and r["cut"] == "ptcone40<2"]
    cls = list(dict.fromkeys(r["cls"] for r in rr)); sets = list(dict.fromkeys(r["set"] for r in rr))
    for i, s in enumerate(sets):
        y = [next(r["ratio"] for r in rr if r["cls"] == c and r["set"] == s) for c in cls]; e = [next(r["ratio_err"] for r in rr if r["cls"] == c and r["set"] == s) for c in cls]
        ax.errorbar(np.arange(len(cls)) + (i - 1.5) * 0.15, y, e, fmt="o", color=COL[i], label="wagi: " + s, capsize=2)
    ax.axhline(1, color=INK2, lw=0.8); ax.set_xticks(range(len(cls))); ax.set_xticklabels(cls); ax.set_ylabel("ε_dane / ε_MC,przeważone  (ptcone40 < 2 GeV)")
    ax.legend(fontsize=8); atlas(ax, "izolacja, H06, pT 3–30 GeV")
    ax = axs[1]
    rr = d["H"]
    if rr:
        keys = list(dict.fromkeys((r["cls"], r["bin"], r["cut"]) for r in rr)); sets = list(dict.fromkeys(r["set"] for r in rr))
        for i, s in enumerate(sets):
            y = [next(r["ratio"] for r in rr if (r["cls"], r["bin"], r["cut"]) == k and r["set"] == s) for k in keys]
            e = [next(r["ratio_err"] for r in rr if (r["cls"], r["bin"], r["cut"]) == k and r["set"] == s) for k in keys]
            ax.errorbar(np.arange(len(keys)) + (i - 1.5) * 0.15, y, e, fmt="o", color=COL[i], label="wagi: " + s, capsize=2, markersize=5)
        ax.axhline(1, color=INK2, lw=0.8); ax.set_xticks(range(len(keys))); ax.set_xticklabels([f"{k[0]}\n{k[1]}\n{k[2]}" for k in keys], fontsize=6.5)
        ax.set_ylim(0.5, 1.5); ax.set_ylabel("ε_dane / ε_MC,przeważone  (score H)"); ax.legend(fontsize=7, ncol=2); atlas(ax, "klasyfikator H, PRE")
    save(fig, "fig4_domkniecie.png")

# 5. przenoszenie D*-tag → inkluzywne (MC)
d = J("s5_transfer.json")
if d:
    fig, ax = plt.subplots(figsize=(6.6, 4.0))
    tests = ["ptcone40<1 (H06)", "ptcone40<2 (H06)", "ptcone40<4 (H06)", "ptcone40/pT<0.5 (H06)", "score>WP0.80 (PRE)"]
    bins = ["3-5", "5-8", "8-30"]
    for i, t in enumerate(tests):
        y = [100 * next(r["rel"] for r in d["rows"] if r["test"] == t and r["bin"] == b) for b in bins]
        e = [100 * next(r["err"] / r["eps_incl"] for r in d["rows"] if r["test"] == t and r["bin"] == b) for b in bins]
        ax.errorbar(np.arange(3) + (i - 2) * 0.12, y, e, fmt="o", color=COL[i], label=t.replace(" (H06)", "").replace(" (PRE)", ""), capsize=2)
    ax.axhline(0, color=INK2, lw=0.8); ax.axhspan(-2, 2, color="#e4e3df", alpha=0.5, lw=0)
    ax.set_xticks(range(3)); ax.set_xticklabels([b.replace("-", "–") + " GeV" for b in bins]); ax.set_ylabel("[ε(D⁰←D*, pT-rw) − ε(inkluzywne D⁰)] / ε  [%]")
    ax.legend(fontsize=8, loc="lower left"); atlas(ax, "MC fd_nd, prawda")
    save(fig, "fig5_przenoszenie.png")

# 6. stożki prostopadłe / embedding
d = J("s6_perp.json"); e7 = J("s7_embed.json") or J("s7_embed_test.json")
if d:
    rr = [r for r in d["rows"] if r["bin"] == "3-30"]
    fig, ax = plt.subplots(figsize=(8.2, 4.0)); x = np.arange(len(rr))
    ax.errorbar(x - 0.24, [r["own"] for r in rr], [r["own_err"] for r in rr], fmt="o", color=COL[0], label="własny stożek D⁰ (tag D*, bez π_s)", capsize=2)
    ax.errorbar(x - 0.08, [r["perp1"] for r in rr], [r["perp1_err"] for r in rr], fmt="s", color=COL[1], label="stożek prostopadły w zdarzeniach z D⁰ (TRKJET2)", capsize=2)
    ax.errorbar(x + 0.08, [r["r_s2"] for r in rr], [r["r_s2_err"] for r in rr], fmt="D", color=COL[2], label="stożek losowy w zdarzeniach z D⁰ (TRKJET2)", capsize=2)
    if e7:
        m = {r["cls"]: r for r in e7["rows"] if r["bin"] == "3-30"}
        xs = [i for i, r in enumerate(rr) if r["cls"] in m]
        ax.errorbar(np.array(xs) + 0.24, [m[rr[i]["cls"]]["WO (bez koła)"] for i in xs], [m[rr[i]["cls"]]["WO (bez koła)_err"] for i in xs], fmt="^", color=COL[3],
                    label="embedding „jak D⁰ wzgl. dżetów poza stożkiem” (WO" + ("; TEST 1–3 pliki)" if "test" in e7["input"] else ")"), capsize=2)
    ax.set_xticks(x); ax.set_xticklabels([r["cls"] for r in rr], fontsize=8); ax.set_ylim(0, 1.25); ax.set_ylabel("ε(Σ pT w R = 0,4 < 2 GeV)")
    ax.legend(fontsize=8, loc="upper right", ncol=1); atlas(ax, "H06, pT 3–30 GeV")
    save(fig, "fig6_embedding.png")
