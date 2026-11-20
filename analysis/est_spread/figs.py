"""EST_spread — rysunki (PNG). Wejście: out/est_new.json, out/collect.json."""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
N = json.load(open(os.path.join(OUT, "est_new.json")))
CJ = json.load(open(os.path.join(OUT, "collect.json")))
INK, MUTED = "#1a1a19", "#6b6a63"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": "white", "axes.grid": True,
                     "grid.color": "#e6e5e0", "grid.linewidth": 0.6})
# paleta kategoryczna (kolejność stała; ten sam estymator = ten sam kolor na obu rysunkach)
COL = {"pol1": "#2a78d6", "sym": "#eb6834", "lr": "#1baf7a", "pol3": "#eda100", "far": "#e87ba4", "mcV2_t2": "#008300", "mcraw_t2": "#6250d6"}
MRK = {"pol1": "o", "sym": "s", "lr": "^", "pol3": "D", "far": "v", "mcV2_t2": "P", "mcraw_t2": "X"}
LAB = {"pol1": "pasma boczne: exp(pol1)", "sym": "pasma boczne: zliczenie 75–125 MeV", "lr": "pasma boczne: interp. liniowa L/P",
       "pol3": "pasma boczne: exp(pol3)", "far": "dalekie pasma (≥ 90 MeV), exp(pol2)", "mcV2_t2": "MC przeważone ze score (typ SALAD)",
       "mcraw_t2": "MC surowe (truth ≥ 2)"}
WPS = N["meta"]["wps"]
X = np.arange(len(WPS))
BINS = ["3-5", "5-8", "8-30"]
BUD = {"0.80": 5, "0.70": 10, "0.50": 15}


def atlas(ax, extra="", x=0.04, y=0.95, fs=11):
    t = ax.text(x, y, "ATLAS", transform=ax.transAxes, fontweight="bold", fontstyle="italic", fontsize=fs, va="top")
    ax.annotate("", xy=(1, 0), xycoords=t, xytext=(4, 0), textcoords="offset points", fontsize=fs, va="bottom")
    if extra:
        ax.text(x, y - 0.085, extra, transform=ax.transAxes, fontsize=fs - 3, va="top", color=INK)


def budget(ax):
    for j, wp in enumerate(WPS):
        if wp in BUD:
            for sgn in (1, -1):
                ax.plot([j - 0.42, j + 0.42], [1 + sgn * BUD[wp] / 100] * 2, color=INK, lw=1.2, ls=(0, (4, 2)), zorder=1)


# ---------------- fig. 1: dane, estymator / pol2
fig, axs = plt.subplots(2, 3, figsize=(13.5, 7.6), sharex=True)
order = ["pol1", "sym", "lr", "far", "mcV2_t2", "mcraw_t2"]
for i, run in enumerate(("d427", "d428")):
    for k, b in enumerate(BINS):
        ax = axs[i, k]
        d = N["data"][f"{run}|{b}"]
        ref = np.array(d["pol2"]["val"])
        rerr = np.array(d["pol2"]["err"]) / ref
        for j in range(len(WPS)):
            ax.fill_between([j - 0.42, j + 0.42], 1 - rerr[j], 1 + rerr[j], color="#d9d8d2", lw=0, zorder=0)
        ax.axhline(1, color=MUTED, lw=1)
        budget(ax)
        for q, e in enumerate(order):
            if e in d:
                v, er = np.array(d[e]["val"]) / ref, np.array(d[e]["err"]) / ref
            else:
                key = e if e.startswith("mcraw") else f"{e}_{run}"
                v, er = np.array(N["mc_est"][b][key]["val"]) / ref, np.array(N["mc_est"][b][key]["err"]) / ref
            ax.errorbar(X + (q - 2.5) * 0.13, v, er, fmt=MRK[e], color=COL[e], ms=6, lw=1.2, mec="white", mew=0.6, label=LAB[e], zorder=3)
        lim = (0.78, 1.5) if b != "8-30" else (0.55, 1.95)
        ax.set_ylim(*lim)
        ax.set_xticks(X); ax.set_xticklabels([f"WP{int(float(w)*100)}" for w in WPS])
        atlas(ax, f"D⁰→Kπ, run {'427929 (MBTS)' if run == 'd427' else '428770 (82 % L1_J12)'}, {b} GeV")
        if k == 0:
            ax.set_ylabel("ε_B(SR) estymatora / ε_B(SR) exp(pol2)")
h, l = axs[0, 0].get_legend_handles_labels()
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
h += [Patch(color="#d9d8d2"), Line2D([0], [0], color=INK, lw=1.2, ls=(0, (4, 2)))]
l += ["błąd stat. nominału exp(pol2)", "przyjęty budżet D⁰: ±5 / ±10 / ±15 % (WP80/70/50)"]
fig.legend(h, l, loc="lower center", ncol=4, frameon=False, fontsize=9)
fig.suptitle("Estymatory tła pod pikiem względem nominalnej interpolacji pasm bocznych (nowe, bt_*.npz, 30 % zdarzeń)", fontsize=11, color=INK)
fig.tight_layout(rect=(0, 0.075, 1, 0.97))
fig.savefig(os.path.join(HERE, "fig1_estymatory_dane_D0.png"), dpi=140)
plt.close(fig)

# ---------------- fig. 2: domknięcie MC — estymator / prawda; D0 (bt, 3 biny) + D+ (VAL)
fig, axs = plt.subplots(1, 4, figsize=(15.5, 4.6))
order2 = ["pol1", "sym", "lr", "pol3", "far"]
for k, b in enumerate(BINS):
    ax = axs[k]
    d = N["mc_closure"][f"t2|{b}"]
    ax.axhline(1, color=MUTED, lw=1)
    budget(ax)
    te = np.array(d["true"]["err"]) / np.array(d["true"]["val"])
    for j in range(len(WPS)):
        ax.fill_between([j - 0.42, j + 0.42], 1 - te[j], 1 + te[j], color="#d9d8d2", lw=0, zorder=0)
    ax.errorbar(X - 0.3, d["pol2"]["ratio"], d["pol2"]["ratio_err"], fmt="o", color=INK, ms=6, lw=1.2, mfc="white", label="pasma boczne: exp(pol2) (nominał)", zorder=3)
    for q, e in enumerate(order2):
        ax.errorbar(X - 0.3 + (q + 1) * 0.12, d[e]["ratio"], d[e]["ratio_err"], fmt=MRK[e], color=COL[e], ms=6, lw=1.2, mec="white", mew=0.6, label=LAB[e], zorder=3)
    ax.set_ylim((0.8, 1.3) if b != "8-30" else (0.3, 2.2))
    ax.set_xticks(X); ax.set_xticklabels([f"WP{int(float(w)*100)}" for w in WPS])
    atlas(ax, f"Simulation, D⁰→Kπ tło (truth ≥ 2), {b} GeV")
    if k == 0:
        ax.set_ylabel("ε_B(SR) estymatora / prawda MC")
ax = axs[3]
ax.axhline(1, color=MUTED, lw=1)
rr = [r for r in CJ["rows"] if r["group"] == "G12_Dp_MC_VAL" and r["pt"] == "all"]
wps = [80, 70, 50]
st = {"pol2": ("o", INK, "pasma boczne: pol2 (nominał)"), "sb_sym": ("s", COL["sym"], "pasma boczne: zliczenie symetryczne"),
      "sb_asym": ("v", COL["far"], "pasma boczne: zliczenie, całe (asym.)")}
for sgn in (1, -1):
    for j, bd in enumerate((5, 10, 25)):
        ax.plot([j - 0.42, j + 0.42], [1 + sgn * bd / 100] * 2, color=INK, lw=1.2, ls=(0, (4, 2)))
for q, (e, (mk, c, lb)) in enumerate(st.items()):
    v = [next(r for r in rr if r["estimator"] == e and r["wp"] == w) for w in wps]
    ax.errorbar(np.arange(3) + (q - 1) * 0.15, [r["value"] for r in v], [r["err"] for r in v], fmt=mk, color=c, ms=6, lw=1.2,
                mfc="white" if e == "pol2" else c, mec=c if e == "pol2" else "white", mew=0.8, label=lb)
ax.set_xticks(range(3)); ax.set_xticklabels(["WP80", "WP70", "WP50"]); ax.set_ylim(0.55, 1.45)
atlas(ax, "Simulation, D⁺→Kππ tło, 3–30 GeV (VAL)")
ax.legend(loc="lower left", frameon=False, fontsize=8)
h, l = axs[0].get_legend_handles_labels()
h += [Patch(color="#d9d8d2"), Line2D([0], [0], color=INK, lw=1.2, ls=(0, (4, 2)))]
l += ["błąd stat. prawdy", "przyjęty budżet (D⁰: ±5/10/15 %, D⁺: ±5/10/25 %)"]
fig.legend(h, l, loc="lower center", ncol=4, frameon=False, fontsize=9)
fig.suptitle("Domknięcie na MC: estymatory z pasm bocznych wobec prawdziwego tła pod pikiem", fontsize=11, color=INK)
fig.tight_layout(rect=(0, 0.12, 1, 0.95))
fig.savefig(os.path.join(HERE, "fig2_domkniecie_MC.png"), dpi=140)
plt.close(fig)
print("ok")
