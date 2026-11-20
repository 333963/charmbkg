#!/usr/bin/env python3
"""TAGGED — rysunki masy (png) i tabele z out/unblind.json, out/hists.npz, out/dst_*.json. Bez nowych fitów."""
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import tg_lib as L  # noqa: E402

FIG = os.path.join(L.TOP, "fig")
C1, C2, INK, MUT = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e"
NAME = {"D0": "D⁰ → Kπ", "Dplus": "D⁺ → Kππ"}


def one(ax, e, H, mu, bkg, title, txt, xl, win=None, reb=1):
    if reb > 1:
        k = (len(H) // reb) * reb
        H, mu, bkg = (v[:k].reshape(-1, reb).sum(1) for v in (H, mu, bkg)); e = e[:k + 1:reb]
    c = 0.5 * (e[1:] + e[:-1])
    lo = np.sqrt(H + 0.25) - 0.5; hi = np.sqrt(H + 0.25) + 0.5
    ax.errorbar(c, H, yerr=[lo, hi], fmt="o", ms=3.5, color=INK, lw=0.8, label="dane (zdarzenia z tagiem AFP)", zorder=3)
    ax.stairs(mu, e, color=C1, lw=2, label="fit: sygnał (kształt stały) + tło")
    ax.stairs(bkg, e, color=C2, lw=2, ls="--", label="tło")
    if win:
        ax.axvspan(win[0], win[1], color=MUT, alpha=0.07, lw=0)
    ax.set_title(title, fontsize=9, color=INK, loc="left"); ax.set_xlabel(xl, color=MUT); ax.set_ylabel(f"kandydaci / {1000 * (e[1] - e[0]):.3g} MeV", color=MUT)
    ax.text(0.98, 0.97, txt, transform=ax.transAxes, ha="right", va="top", fontsize=8, color=INK)
    ax.set_ylim(0, None); ax.grid(alpha=0.2, lw=0.5); ax.tick_params(colors=MUT)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def main():
    os.makedirs(FIG, exist_ok=True)
    U = json.load(open(os.path.join(L.OUT, "unblind.json")))["items"]
    Hs = np.load(os.path.join(L.OUT, "hists.npz"))
    groups = {"fig1_427929_MBTS": [f"{sp}|427929|MBTS|{b}|0.7|T2" for b in ("3-30", "5-30") for sp in L.SPN],
              "fig2_428770_regiony_D0": [f"D0|428770|{r}|3-30|0.7|T2" for r in ("R1", "R2", "R3", "J12")],
              "fig3_428770_regiony_Dplus": [f"Dplus|428770|{r}|3-30|0.7|T2" for r in ("R1", "R2", "R3", "J12")],
              "fig4_kontrolki_CR_MBAFP": [f"{sp}|427929|MBTS|3-30|0.7|CR" for sp in L.SPN] + [f"{sp}|428770|J12|3-30|0.7|CR" for sp in L.SPN] + [f"{sp}|427929|MBAFP|3-30|0.7|T2" for sp in L.SPN]}
    for fn, keys in groups.items():
        keys = [k for k in keys if k in U]
        nc = 2; nr = (len(keys) + 1) // 2
        fig, axs = plt.subplots(nr, nc, figsize=(11, 3.6 * nr), squeeze=False)
        for ax in axs.flat[len(keys):]:
            ax.axis("off")
        for ax, k in zip(axs.flat, keys):
            r = U[k]; f = r["fit"]; sp = r["sp"]; kq = k.replace("|", "__")
            m0, w = L.WIN[sp]
            reb = 2 if r["n_tot"] < 3000 else 1
            txt = f"N_S = {f['N']:.1f} −{r['stat_lo']:.1f} +{r['stat_hi']:.1f}\nZ = {f['Z']:.1f}\nokno: {r['n_win']:.0f} kand., tło {r['bkg_win']:.0f}"
            one(ax, L.edges5(sp), Hs[kq + "__H"], Hs[kq + "__mu"], Hs[kq + "__bkg"], f"{NAME[sp]}, {r['run']} {r['reg']}, pT {r['bins'][0].split('-')[0]}–30 GeV, WP{int(100 * r['eps'])}, {r['xsel']}",
                txt, "m [GeV]", (m0 - w, m0 + w), reb)
        axs.flat[0].legend(fontsize=7, loc="lower left", frameon=False)
        fig.tight_layout(); fig.savefig(os.path.join(FIG, fn + ".png"), dpi=130); plt.close(fig)
    dp = os.path.join(L.OUT, "dst_unblind.json")
    if os.path.exists(dp):
        D = json.load(open(dp)); Hd = np.load(os.path.join(L.OUT, "dst_hist.npz"))
        fig, axs = plt.subplots(1, 2, figsize=(11, 3.6))
        for ax, k in zip(axs, D):
            kq = k.replace("|", "__"); f = D[k]["fit"]
            one(ax, Hd["edges"], Hd[kq + "__tag"], Hd[kq + "__mu"], Hd[kq + "__bkg"], f"D*⁺ → D⁰π (bez klasyfikatora), {k.split('|')[1]} {k.split('|')[2]}, T2",
                f"N_S = {f['N']:.1f} −{f['N'] - f['lo68']:.1f} +{f['hi68'] - f['N']:.1f}\nZ = {f['Z']:.1f}", "Δm [GeV]", (0.14345, 0.14745))
        axs[0].legend(fontsize=7, loc="lower right", frameon=False)
        fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig5_Dstar.png"), dpi=130); plt.close(fig)
    print("rysunki:", sorted(os.listdir(FIG)))


if __name__ == "__main__":
    main()
