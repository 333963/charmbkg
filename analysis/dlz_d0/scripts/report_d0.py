"""DLZ_d0: tabele (tables.md) i rysunki z wyników węzła `d0` paczki dlz_data
(domyślnie DLZ_dplus/DATA_2026-11-17/runs/lxplus/d0).  Użycie: python report_d0.py [KATALOG]"""
import json
import os
import sys

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = sys.argv[1] if len(sys.argv) > 1 else "/home/k/Documents/charm/charm_work/DLZ_dplus/DATA_2026-11-17/runs/lxplus/d0"
FIG = os.path.join(BASE, "fig")
os.makedirs(FIG, exist_ok=True)
BINS = ["3-5", "5-8", "8-30"]
out = []
P = out.append


def atlas(ax, txt=""):
    ax.text(0.04, 0.95, "ATLAS", transform=ax.transAxes, fontweight="bold", fontstyle="italic", va="top", fontsize=11)
    ax.text(0.19, 0.95, "", transform=ax.transAxes, va="top", fontsize=11)
    if txt:
        ax.text(0.04, 0.87, txt, transform=ax.transAxes, va="top", fontsize=8)


Z = pd.read_csv(os.path.join(SRC, "D0_mc", "zproxy.csv"))
order = [v for v in ["H", "H_seed2", "H+ALL_noiso", "H+d0prod", "H+ip", "H+time", "H+NOKIN", "H+NOKIN_pid", "H+pidbb", "H+kin2", "H+kinall", "H+ALLnew"] if v in set(Z.variant)]
T = pd.read_csv(os.path.join(SRC, "transfer_D0.csv")) if os.path.exists(os.path.join(SRC, "transfer_D0.csv")) else None
SC = json.load(open(os.path.join(SRC, "D0_mc", "sculpt.json"))) if os.path.exists(os.path.join(SRC, "D0_mc", "sculpt.json")) else {}
DS = json.load(open(os.path.join(SRC, "d0_sculpt.json"))) if os.path.exists(os.path.join(SRC, "d0_sculpt.json")) else {}
MO = json.load(open(os.path.join(SRC, "modelling_D0.json"))) if os.path.exists(os.path.join(SRC, "modelling_D0.json")) else {}
ZD = pd.read_csv(os.path.join(SRC, "D0_data", "zproxy.csv")) if os.path.exists(os.path.join(SRC, "D0_data", "zproxy.csv")) else None
g = lambda df, v, b, wp, c: float(df[(df.variant == v) & (df.bin == b) & (np.isclose(df.wp, wp))][c].iloc[0])  # noqa: E731

P("## 1. Zproxy na MC (test mc), ΔZ/Z_H [%] ± bootstrap sparowany (200 replik)\n")
P("| wariant | WP70 3–5 | WP70 5–8 | WP70 8–30 | WP50 3–5 | WP80 3–5 |")
P("|---|---|---|---|---|---|")
for v in order:
    if v == "H":
        P("| H (Z) | " + " | ".join(f"{g(Z, v, b, 0.7, 'Z'):.2f} ± {g(Z, v, b, 0.7, 'Z_err'):.2f}" for b in BINS) + f" | {g(Z, v, '3-5', 0.5, 'Z'):.2f} | {g(Z, v, '3-5', 0.8, 'Z'):.2f} |")
    else:
        P(f"| {v} | " + " | ".join(f"{100*(g(Z, v, b, 0.7, 'ratio_vs_H')-1):+.1f} ± {100*g(Z, v, b, 0.7, 'ratio_err'):.1f}" for b in BINS)
          + f" | {100*(g(Z, v, '3-5', 0.5, 'ratio_vs_H')-1):+.1f} | {100*(g(Z, v, '3-5', 0.8, 'ratio_vs_H')-1):+.1f} |")
if T is not None:
    P("\n## 2. Przeniesienie na dane (model MC, ε_B z dalekich SB danych 428770), WP70\n")
    P("| wariant | ε_B SB danych 3–5 / 5–8 / 8–30 | Z_transfer/Z_H 3–5 / 5–8 / 8–30 (± bootstrap ε_B) |")
    P("|---|---|---|")
    for v in order:
        eb = " / ".join(f"{g(T, v, b, 0.7, 'epsB_dataSB'):.4f}" for b in BINS)
        zr = " / ".join(f"{g(T, v, b, 0.7, 'Z_transfer_ratio_vs_H'):.3f} ± {g(T, v, b, 0.7, 'Z_transfer_err')/g(T, v, b, 0.7, 'Z_transfer'):.3f}" for b in BINS)
        P(f"| {v} | {eb} | {zr} |")
if ZD is not None:
    P("\n## 3. Test na danych (tło trening = dalekie SB danych), WP70, ΔZ/Z_H [%]\n")
    P("| wariant | 3–5 | 5–8 | 8–30 |")
    P("|---|---|---|---|")
    for v in [x for x in order if x in set(ZD.variant) and x != "H"]:
        P(f"| {v} | " + " | ".join(f"{100*(g(ZD, v, b, 0.7, 'ratio_vs_H')-1):+.1f} ± {100*g(ZD, v, b, 0.7, 'ratio_err'):.1f}" for b in BINS) + " |")
P("\n## 4. Rzeźbienie masy tła\n")
P("| wariant | R WP70 MC truth ≥ 8 (SR) | R WP70 dane (pasmo) | nachylenie ε_B [/100 MeV] | fałszywy sygnał WP70: N_fake/σ (3–30) | per bin 3–5 / 5–8 / 8–30 | N_fake/N_tła(SR) | χ² ε_B(m) pol1 / pol2 (ndf 23 / 22) | R² regresji masy z cech |")
P("|---|---|---|---|---|---|---|---|---|")
for v in order:
    s = SC.get(v, {}).get("WP70", {})
    a, d = s.get("mc_truth8_SR"), s.get("data_band")
    fk, sh, mr = DS.get("fake", {}).get(v, {}), DS.get("shape", {}).get(v, {}).get("WP70", {}), DS.get("massreg", {}).get(v, {})
    P(f"| {v} | " + (f"{a['R']:.3f} ± {a['R_err']:.3f}" if a else "—") + " | " + (f"{d['R']:.3f} ± {d['R_err']:.3f}" if d else "—") + " | "
      + (f"{a['slope_rel_per_100MeV']:+.2f}" if a else "—") + " | "
      + (f"{fk['WP70|3-30']['N']:.0f} ± {fk['WP70|3-30']['err']:.0f} ({fk['WP70|3-30']['pull']:+.1f}σ)" if fk else "—") + " | "
      + (" / ".join(f"{fk[f'WP70|{b}']['pull']:+.1f}" for b in BINS) if fk else "—") + " | "
      + (f"{100*fk['WP70|3-30']['frac']:+.2f} %" if fk else "—") + " | "
      + (f"{sh['chi2_pol1']:.0f} / {sh['chi2_pol2']:.0f}" if sh else "—") + " | " + (f"{mr['r2']:.3f}" if mr else "—") + " |")
if DS.get("massreg"):
    P("\nRegresja masy z samych nowych grup (R², tło MC truth ≥ 8): " + ", ".join(f"{k[5:]} {v['r2']:.3f}" for k, v in DS["massreg"].items() if k.startswith("only_")) + ".")
if MO:
    P("\n## 5. Modelowanie tła (C2ST SB danych vs MC mb, 3 ziarna)\n")
    P("| zestaw | C2ST (ziarna) | średnia | Δ vs H |")
    P("|---|---|---|---|")
    h0 = np.mean(MO["c2st"]["H"])
    for k, v in MO["c2st"].items():
        P(f"| {k} | {' / '.join(f'{x:.4f}' for x in v)} | {np.mean(v):.4f} | {np.mean(v)-h0:+.4f} |")
    P("\nC2ST 1D / KS nowych cech: " + "; ".join(f"{k} {v:.3f} / {MO['c2st1d'][k + '|ks']:.3f}" for k, v in MO["c2st1d"].items() if "|" not in k) + ".")
open(os.path.join(BASE, "tables.md"), "w").write("\n".join(out) + "\n")
print("\n".join(out))

# rysunki
fig, ax = plt.subplots(figsize=(10, 5))
vs = [v for v in order if v not in ("H",)]
x = np.arange(len(vs))
for i, (b, col) in enumerate(zip(BINS, ("C0", "C1", "C2"))):
    ax.errorbar(x + 0.2 * (i - 1), [100 * (g(Z, v, b, 0.7, "ratio_vs_H") - 1) for v in vs], [100 * g(Z, v, b, 0.7, "ratio_err") for v in vs], fmt="o", color=col, capsize=2, label=f"MC, {b} GeV")
    if T is not None:
        ax.plot(x + 0.2 * (i - 1), [100 * (g(T, v, b, 0.7, "Z_transfer_ratio_vs_H") - 1) for v in vs], "x", color=col, ms=7, label=f"transfer (SB danych), {b} GeV")
ax.axhline(0, color="k", lw=0.6)
ax.axhline(10, color="0.5", ls="--", lw=0.8)
ax.axhline(4.5, color="0.7", ls=":", lw=0.8)
ax.set_xticks(x)
ax.set_xticklabels(vs, rotation=35, ha="right", fontsize=8)
ax.set_ylabel("$\\Delta Z_{proxy}/Z_H$ @WP70 [%]")
ax.legend(fontsize=7, ncol=2, loc="upper left", bbox_to_anchor=(0.3, 1.0))
atlas(ax, "D$^0\\to K\\pi$; linie: 10 % (cel), 4,5 % (3× szum)")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "fig1_zproxy_d0.png"), dpi=130)
plt.close(fig)
if DS.get("shape"):
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.5))
    c = 1.75 + (np.arange(25) + 0.5) * 0.01
    for v, col in zip([v for v in ("H", "H+NOKIN", "H+kin2", "H+kinall", "H+ALLnew") if v in DS["shape"]], ("k", "C0", "C1", "C3", "C4")):
        e, er = np.array(DS["shape"][v]["WP70"]["eps"]), np.array(DS["shape"][v]["WP70"]["err"])
        axs[0].errorbar(c, e / e.mean(), er / e.mean(), fmt="o-", ms=3, lw=0.7, color=col, label=v)
    axs[0].axvspan(1.86484 - 0.055, 1.86484 + 0.055, color="0.92")
    axs[0].set_xlabel("m(K$\\pi$) [GeV]")
    axs[0].set_ylabel("$\\varepsilon_B(m)$ / średnia, tło MC truth ≥ 8, WP70")
    axs[0].legend(fontsize=8, loc="lower left")
    axs[0].set_ylim(0.3, 2.0)
    atlas(axs[0], "D$^0$: kształt tła po cięciu")
    vv = [v for v in order if v in DS["fake"]]
    axs[1].bar(np.arange(len(vv)), [DS["fake"][v]["WP70|3-30"]["pull"] for v in vv], color=["0.5" if v in ("H", "H_seed2") else "C0" for v in vv])
    axs[1].set_xticks(np.arange(len(vv)))
    axs[1].set_xticklabels(vv, rotation=35, ha="right", fontsize=8)
    axs[1].axhline(0, color="k", lw=0.6)
    axs[1].set_ylabel("fałszywy sygnał na tle MC: $N_{fake}/\\sigma$ (WP70, 3–30 GeV)")
    atlas(axs[1], "fit Gauss + exp(pol2) do tła bez sygnału")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig2_sculpt_d0.png"), dpi=130)
    plt.close(fig)
