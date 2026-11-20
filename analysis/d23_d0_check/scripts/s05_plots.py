"""Rysunki: widma masy D0 2-3 GeV (oba runy, oba klasyfikatory, punkty z werdyktu etapu 4) + kontrolne (tlo MC 3-5 GeV po cieciu, eps_B(m))."""
import numpy as np, json
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from vlib import *
F = load("fd_d23"); F35 = load("fd_nd_d150"); MB = load("mb_nd_d150")
D = {r: load(f"data_{r}_d23") for r in ("427929", "428770")}
HEAD = {("427929", 0): 0.1, ("427929", 1): 0.2, ("428770", 0): 0.1, ("428770", 1): 0.1}
FN = {0: "drzewo", 1: "siec"}
J = json.load(open(os.path.join(BASE, "out", "s04.json")))
for (run, j), wp in HEAD.items():
    d = D[run]; th = mc_wp(F, j, wp); q = templ(F, j, th); s = d["lg"][:, j] >= th
    fig, ax = plt.subplots(2, 2, figsize=(12, 7), sharex="col", gridspec_kw=dict(height_ratios=[3, 1], hspace=0.05))
    for k, (rng, nb, lab) in enumerate((((1.75, 2.0), 2, "fit nominalny: exp(pol2), 1.75–2.00 GeV"), ((1.70, 2.05), 3, "wariant: exp(pol3), 1.70–2.05 GeV"))):
        H, e = hist(d["m"][s], rng); f = fit2(H, e, q, nbkg=nb); c = f["ctr"]; S, Rf, B = f["comps"]
        a = ax[0, k]
        a.errorbar(c, H, np.sqrt(np.maximum(H, 1)), fmt="o", ms=3, color="k", lw=0.8, label=f"dane {run}, {int(s.sum())} kand. po cięciu")
        a.step(e, np.r_[S + Rf + B, (S + Rf + B)[-1]], where="post", color="#c0392b", lw=1.6, label=f"sygnał + odbicie + tło (χ² = {f['chi2']:.0f}/{f['ndf']})")
        a.step(e, np.r_[Rf + B, (Rf + B)[-1]], where="post", color="#2471a3", lw=1.2, ls="--", label=f"tło + odbicie K↔π (f_refl = {q['frefl']:.2f}, z MC)")
        a.step(e, np.r_[B, B[-1]], where="post", color="#7f8c8d", lw=1.2, ls=":", label="tło")
        a.axvspan(M0 - WIN, M0 + WIN, color="#f4d03f", alpha=0.12, lw=0)
        a.set_title(f"{lab}\nN = {f['N']:.0f} ± {f['err']:.0f} (HESSE), N/ε = {f['N']/q['eps']:.0f}", fontsize=10)
        a.set_ylabel("kandydaci / 5 MeV"); a.legend(fontsize=8, frameon=False); a.set_ylim(0, None)
        pl = (H - (S + Rf + B)) / np.sqrt(np.maximum(S + Rf + B, 1))
        ax[1, k].bar(c, pl, width=0.005, color="#566573"); ax[1, k].axhline(0, color="k", lw=0.5); ax[1, k].set_ylim(-4, 4)
        ax[1, k].set_xlabel("m(Kπ) [GeV]"); ax[1, k].set_ylabel("(dane − fit)/σ")
    fig.suptitle(f"D⁰→Kπ, pT 2–3 GeV, run {run}, {'drzewo D0.xgb.H.sb.n.s1' if j == 0 else 'sieć D0.tr.H+e04.mc.d5.s1'}, ε_MC = {q['eps']:.3f} (WP{int(wp*100)}); kształt piku ustalony z MC", fontsize=11)
    fig.savefig(os.path.join(BASE, "fig", f"masa_D0_2-3_{run}_{FN[j]}_WP{int(wp*100)}.png"), dpi=130, bbox_inches="tight"); plt.close(fig)

# kontrolne: tlo MC (mb_nd 3-5 GeV, prawda: nie-D0) po cieciu WP10/WP20 + eps_B(m) w danych 2-3
fig, ax = plt.subplots(2, 2, figsize=(12, 8))
for j in (0, 1):
    a = ax[0, j]
    for wp, col in ((0.2, "#2471a3"), (0.1, "#c0392b")):
        th = mc_wp(F35, j, wp); q = templ(F35, j, th)
        b = (MB["cat"] >= 2) & (MB["lg"][:, j] >= th)
        H, e = hist(MB["m"][b], (1.70, 2.05)); c = 0.5 * (e[1:] + e[:-1])
        H2 = H.reshape(-1, 2).sum(1); c2 = c.reshape(-1, 2).mean(1)
        a.errorbar(c2, H2, np.sqrt(np.maximum(H2, 1)), fmt="o", ms=3, color=col, lw=0.8, label=f"WP{int(wp*100)}: tło MC (truth ≠ D⁰→Kπ), {int(H.sum())} kand.")
        sg = (MB["cat"] == 0) & (MB["lg"][:, j] >= th)
        Hs, _ = hist(MB["m"][sg], (1.70, 2.05)); a.step(e[::2], np.r_[Hs.reshape(-1, 2).sum(1), 0], where="post", color=col, lw=0.8, alpha=0.5, label=f"WP{int(wp*100)}: sygnał MC (dla skali)")
    a.axvspan(M0 - WIN, M0 + WIN, color="#f4d03f", alpha=0.12, lw=0); a.set_ylim(0, None)
    a.set_title(f"MC minimum-bias (mb_nd), pT 3–5 GeV, {NM[j]}: tło z prawdą po cięciu", fontsize=10); a.set_xlabel("m(Kπ) [GeV]"); a.set_ylabel("kandydaci / 10 MeV"); a.legend(fontsize=7, frameon=False)
ME = D["427929"]["ME"]; ctr = 0.5 * (ME[1:] + ME[:-1])
for k, run in enumerate(("427929", "428770")):
    a = ax[1, k]
    for j, col in ((0, "#c0392b"), (1, "#2471a3")):
        wp = HEAD[(run, j)]; th = mc_wp(F, j, wp); d = D[run]
        allm = d["H"][0, j].sum(1); npass = np.histogram(d["m"][d["lg"][:, j] >= th], ME)[0].astype(float)
        r = npass / allm; a.errorbar(ctr, 1e4 * r, 1e4 * np.sqrt(np.maximum(npass, 1)) / allm, fmt="o", ms=3, color=col, lw=0.8, label=f"{NM[j]} WP{int(wp*100)}")
    a.axvspan(M0 - WIN, M0 + WIN, color="#f4d03f", alpha=0.12, lw=0); a.set_yscale("log")
    a.set_title(f"dane {run}, pT 2–3 GeV: ułamek kandydatów przechodzących cięcie vs masa", fontsize=10); a.set_xlabel("m(Kπ) [GeV]"); a.set_ylabel("N_po cięciu / N_wszystkie  [10⁻⁴]"); a.legend(fontsize=8, frameon=False)
fig.tight_layout(); fig.savefig(os.path.join(BASE, "fig", "kontrola_tloMC_3-5_i_epsB_2-3.png"), dpi=130); plt.close(fig)
for k, v in J.items():
    if k.startswith("y35"): print(k, f"dmu {1e3*v['dmu']:+.1f} MeV sc {v['sc']:.3f}")
