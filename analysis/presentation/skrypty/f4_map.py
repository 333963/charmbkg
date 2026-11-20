#!/usr/bin/env python3
"""F4 — mapa m(Kpi) x Delta m kandydatow D* w danych 427929 w ukladzie wzorca `_margins`: mapa 2D + rozklady brzegowe + linie mas nominalnych + panel reszt
(obserwacja - iloczyn brzegowych / N). Wejscie: lokalne tabele par INCAL (D0__pair_*), wczytane funkcja XSEC_atlas xa_lib.load; selekcja jak XSEC_atlas xa_dstar.sel
(L_xy >= 0,25 mm, pT corek > 1 GeV, |eta| < 2,2, MBTS + GRL), ale BEZ okna masy D0. Zadnego fitu."""
import sys

import numpy as np

import pstyle as P
from pstyle import plt

sys.path.insert(0, P.CW + "/XSEC_atlas_2026-11-20/skrypty")
import xa_lib as X  # noqa: E402

PT = (5.0, 30.0)
M0, DM0 = X.M0["D0"] * 1e3, 145.4258          # PDG: m(D0); m(D*+) - m(D0) [MeV]
D = X.derived(X.load("data_427929", "D0__pair_"), "D0")
s = ((D["cls"] & 3) == 3) & (D["Lxy"] >= 0.25) & (D["ptmin"] > 1.0) & (np.abs(D["eta"]) < 2.2) & (D["ds_pt"] >= PT[0]) & (D["ds_pt"] < PT[1])
m, dm = D["m"][s] * 1e3, D["ds_dm"][s] * 1e3
ex = np.arange(1790.0, 1940.0 + 1e-6, 5.0)
ey = np.arange(140.0, 156.0 + 1e-6, 0.5)
H, _, _ = np.histogram2d(m, dm, [ex, ey])
N = H.sum()
hx, hy = H.sum(1), H.sum(0)
E = np.outer(hx, hy) / N
res = H - E
print(f"pary: {int(N)} w mapie (wybrane {int(s.sum())}); max bin {H.max():.0f}; reszty: min {res.min():.1f}, max {res.max():.1f}")

fig = plt.figure(figsize=(10.4, 10.2))
gs = fig.add_gridspec(2, 3, width_ratios=[1.0, 3.0, 0.11], height_ratios=[3.0, 1.15], wspace=0.07, hspace=0.24, left=0.085, right=0.925, top=0.875, bottom=0.125)
ax = fig.add_subplot(gs[0, 1])
im = ax.pcolormesh(ex, ey, np.ma.masked_equal(H.T, 0), cmap="Blues", vmin=0, shading="flat", rasterized=True)
cb = fig.colorbar(im, cax=fig.add_subplot(gs[0, 2])); cb.set_label("Pairs / (5 MeV × 0.5 MeV)")
ax.tick_params(labelbottom=False, labelleft=False)
lk = dict(color=P.C["orange"], lw=1.3, ls="--")
ax.axvline(M0, **lk); ax.axhline(DM0, **lk)
ax.text(M0 + 2, ey[-1] - 0.35, r"$m(D^0)$ = %.2f MeV" % M0, color=P.INK, fontsize=10, va="top", ha="left", bbox=dict(fc="white", ec="none", alpha=0.8, pad=1.5))
ax.text(ex[-1] - 2, DM0 + 1.3, r"$m(D^{*+}) - m(D^0)$ = %.2f MeV" % DM0, color=P.INK, fontsize=10, va="bottom", ha="right", bbox=dict(fc="white", ec="none", alpha=0.8, pad=1.5))
ax.set_title(r"$D^{*\pm} \to D^0\pi^\pm_s$ candidates, $D^0 \to K\pi$,  %g ≤ $p_T(D^*)$ < %g GeV,  no classifier" % PT, loc="left", fontsize=11.5)
# brzeg lewy: Delta m
al = fig.add_subplot(gs[0, 0], sharey=ax)
al.stairs(hy, ey, orientation="horizontal", fill=True, color=P.C["blue"], alpha=0.35, lw=0)
al.stairs(hy, ey, orientation="horizontal", color=P.C["blue"], lw=1.3)
al.axhline(DM0, **lk)
al.invert_xaxis(); al.set_xlabel("Pairs / 0.5 MeV"); al.set_ylabel(r"$\Delta m = m(K\pi\pi_s) - m(K\pi)$ [MeV]")
al.set_ylim(ey[0], ey[-1]); al.tick_params(axis="x", labelsize=9)
# brzeg dolny: m(Kpi)
ab = fig.add_subplot(gs[1, 1], sharex=ax)
ab.stairs(hx, ex, fill=True, color=P.C["blue"], alpha=0.35, lw=0)
ab.stairs(hx, ex, color=P.C["blue"], lw=1.3)
ab.axvline(M0, **lk)
ab.set_xlabel(r"$m(K\pi)$ [MeV]"); ab.set_ylabel("Pairs / 5 MeV"); ab.yaxis.set_label_position("right"); ab.yaxis.tick_right()
ab.set_xlim(ex[0], ex[-1]); ab.set_ylim(0, hx.max() * 1.15)
# reszty
ar = fig.add_subplot(gs[1, 0])
v = np.abs(res).max()
ir = ar.pcolormesh(ex, ey, res.T, cmap="RdBu_r", vmin=-v, vmax=v, shading="flat", rasterized=True)
ar.axvline(M0, color=P.INK2, lw=0.7, ls="--"); ar.axhline(DM0, color=P.INK2, lw=0.7, ls="--")
ar.set_xlabel(r"$m(K\pi)$ [MeV]", fontsize=10); ar.set_ylabel(r"$\Delta m$ [MeV]", fontsize=10)
ar.tick_params(labelsize=8.5); ar.set_xticks([1800, 1850, 1900])
ar.set_title("Corr. excess: obs − (margin × margin)/N", fontsize=9, loc="left")
cax = ar.inset_axes([0.0, -0.42, 1.0, 0.06])
cb2 = fig.colorbar(ir, cax=cax, orientation="horizontal"); cb2.ax.tick_params(labelsize=8.5); cb2.set_label("Pairs per cell", fontsize=9)
P.header(fig, P.run427_text(), extra=f"{int(N)} pairs")
fig.text(0.925, 0.985, "Selection: $L_{xy}$ ≥ 0.25 mm, daughter $p_T$ > 1 GeV, |η($D^0$)| < 2.2.\nDashed lines: nominal masses. No fit is drawn.", ha="right", va="top", fontsize=9, color=P.INK2)
P.save(fig, "F4_map_mKpi_dm")
