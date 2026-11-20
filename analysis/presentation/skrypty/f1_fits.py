#!/usr/bin/env python3
"""F1 — widma mas z fitami nominalnymi: D0 (M8, szablon pikujacy wprost z drzewa H, r = 1) i D+ (A2), run 427929, WP70, 3 biny pT.
Fit: kod FITBIAS/BKGCAT/FITBIAS_Dplus importowany bez zmian (zadnego dostrajania); uzysk, blad stat. i p(chi2) czytane z XSEC_final v2 final.json
i porownywane z fitem odtworzonym tutaj (asercja)."""
import sys

import numpy as np

import pstyle as P
from pstyle import plt

sys.path.insert(0, P.CW + "/FITBIAS_Dplus_2026-11-20/skrypty")
sys.path.insert(0, P.CW + "/BKGCAT_2026-11-20/skrypty")
import bc_m8 as M  # noqa: E402
import fd_lib as DP  # noqa: E402
V, F, L = M.V, M.F, M.L
M._init(L.IN)
FIN = P.J("XSEC_final_2026-11-20", "v2", "out", "final.json")


def decomp(f, sigT, tied, edges, m0, nbkg):
    x = np.array(f["x"])
    ctr = 0.5 * (edges[1:] + edges[:-1])
    u = 2 * (ctr - edges[0]) / (edges[-1] - edges[0]) - 1
    sig = x[2] * sigT.bins(edges, m0, x[0], x[1])
    pk = [x[2] * kap * T.bins(edges, m0, x[0], 1.0) for T, kap in tied]
    sh = F.bshape(u, x[4:4 + nbkg], "exp")
    bkg = x[3] * sh / sh.sum()
    return ctr, sig, pk, bkg


def fit_d0(bn):
    c, S = F.SP["D0"], V.SET["M8"]
    k0 = V.wcat("D0", bn, 0.7)
    H = V.wdata(k0, V.SRC["427929"])[0]
    k = M.patch(k0, "D0", bn, 0.7, None, "dir", 1.0)
    rng = c[S["rng"]]
    edges = np.arange(rng[0], rng[1] + F.BW / 2, F.BW)
    tied = [(k["refT"], k["frefl"]), (k["peakT"], k["kap"])]
    h = V.cutr(H, "D0", rng)
    f = F.fit(h[None], [k["sigT"]], [tied], edges, c["m0"], c["nbkg"] + S["dn"], shape=S["shape"], gc=V.GC["D0"], win=c["win"])
    return h, f, edges, decomp(f, k["sigT"], tied, edges, c["m0"], c["nbkg"]), ["K↔π reflection", "other charm decays (peaking)"]


def fit_dp(bn):
    S = DP.SET["A2"]
    sigT, k = DP.sig_tmpl(bn)
    T = DP.tmpls(DP.comps("full", bn))
    H = V.wdata(k, V.SRC["427929"])[0]
    rng = S["rng"]
    edges = np.arange(rng[0], rng[1] + F.BW / 2, F.BW)
    f = DP.fitA("A2", H, sigT, T, keep_mu=True)
    tied = [(T[t][0], T[t][1]) for t in S["tied"]]
    return DP.cutr(H, rng), f, edges, decomp(f, sigT, tied, edges, DP.C["m0"], S["nbkg"]), ["other charm decays (peaking)"]


def draw(sp, name, fitf, xlab, model, m0, mlab):
    fig = plt.figure(figsize=(14.2, 8.6))
    gs = fig.add_gridspec(3, 3, height_ratios=[3.0, 2.3, 1.0], hspace=0.07, wspace=0.27, left=0.065, right=0.985, top=0.845, bottom=0.075)
    pkcol = [P.C["yellow"], P.C["aqua"]] if sp == "D0" else [P.C["aqua"]]
    for ib, bn in enumerate(P.BINS):
        h, f, edges, (ctr, sig, pk, bkg), pkn = fitf(bn)
        cell = FIN["cells"][f"{sp}|{bn}"]
        assert abs(f["N"][0] - cell["N"]) < 0.5, (sp, bn, f["N"][0], cell["N"])
        p = F.pval(f["chi2"], f["ndf"])
        assert abs(p - cell["p_chi2"]) < 2e-3, (p, cell["p_chi2"])
        estat = cell["rel"]["stat"] * cell["N"]
        mu = f["mu"][0]
        x = ctr * 1e3
        w = (edges[1] - edges[0]) * 1e3
        a0, a1, a2 = fig.add_subplot(gs[0, ib]), fig.add_subplot(gs[1, ib]), fig.add_subplot(gs[2, ib])
        # gora: widmo + fit
        a0.errorbar(x, h, np.sqrt(h), fmt="o", ms=2.6, color=P.DATA, lw=0.8, label="Data", zorder=5)
        a0.stairs(mu, edges * 1e3, color=P.C["blue"], lw=1.8, label=f"Fit ({model})", zorder=4)
        a0.stairs(bkg, edges * 1e3, color=P.C["orange"], lw=1.5, ls="--", label="Combinatorial background", zorder=3)
        lo = min(h.min(), bkg.min())
        a0.set_ylim(lo * 0.90, h.max() * 1.38)
        a0.set_ylabel(f"Candidates / {w:.0f} MeV")
        a0.set_title(f"{mlab},  {P.BINLAB[bn]} GeV", loc="left", fontsize=11.5)
        a0.text(0.97, 0.955, f"N = {cell['N']:.0f} ± {estat:.0f} (stat)\n" + rf"$\chi^2$/ndf = {f['chi2']:.1f}/{f['ndf']},  p = {cell['p_chi2']:.2f}",
                transform=a0.transAxes, ha="right", va="top", fontsize=10)
        if ib == 0:
            hl0 = a0.get_legend_handles_labels()
        # srodek: dane - tlo kombinatoryczne, skladowe
        d = h - bkg
        a1.axhline(0, color=P.INK2, lw=0.7)
        base = np.zeros_like(sig)
        for kk, (pp, nm) in enumerate(zip(pk, pkn)):
            a1.stairs(base + pp, edges * 1e3, baseline=base, fill=True, color=pkcol[kk], alpha=0.75, lw=0, label=nm.capitalize() if nm[0] != "K" else nm)
            base = base + pp
        a1.stairs(base + sig, edges * 1e3, baseline=base, fill=True, color=P.C["blue"], alpha=0.45, lw=0, label="Signal")
        a1.stairs(mu - bkg, edges * 1e3, color=P.C["blue"], lw=1.6, label="_nolegend_")
        a1.errorbar(x, d, np.sqrt(h), fmt="o", ms=2.6, color=P.DATA, lw=0.8, label="Data − comb. background", zorder=5)
        top = max((mu - bkg).max(), (d + np.sqrt(h)).max())
        a1.set_ylim(min((d - np.sqrt(h)).min() * 1.1, -0.1 * top), top * 1.15)
        a1.set_ylabel("Data − comb. bkg")
        if ib == 0:
            hl1 = a1.get_legend_handles_labels()
            hh = [hl0[0][2], hl0[0][0], hl0[0][1]] + hl1[0][:-1] + [hl1[0][-1]]
            ll = [hl0[1][2], hl0[1][0], hl0[1][1]] + hl1[1][:-1] + [hl1[1][-1] + " (middle row)"]
            fig.legend(hh, ll, loc="upper center", bbox_to_anchor=(0.525, 0.925), ncol=len(hh), fontsize=9.5, handlelength=1.6, columnspacing=1.4)
        # dol: pulle
        pull = (h - mu) / np.sqrt(mu)
        a2.axhspan(-1, 1, color=P.GRID, alpha=0.6, lw=0)
        a2.axhline(0, color=P.INK2, lw=0.7)
        a2.bar(x, pull, width=w * 0.85, color=P.INK2, lw=0)
        a2.set_ylim(-4.2, 4.2); a2.set_yticks([-3, 0, 3])
        a2.set_ylabel("Pull")
        a2.set_xlabel(xlab)
        for a in (a0, a1, a2):
            a.set_xlim(edges[0] * 1e3, edges[-1] * 1e3)
            a.axvline(m0 * 1e3, color=P.INK2, lw=0.8, ls=":", zorder=1)
        a0.tick_params(labelbottom=False); a1.tick_params(labelbottom=False)
        print(f"{sp} {bn}: N {f['N'][0]:.1f} (final.json {cell['N']:.1f}), stat {estat:.1f}, p {p:.3f}")
    P.header(fig, P.run427_text(), extra="tree H at WP70, |η(D)| < 2.2")
    fig.text(0.985, 0.985, "Dotted line: nominal mass. Statistical error of N: max(event bootstrap, 1.3 × HESSE).\nPulls: (data − fit)/√fit, band ±1.",
             ha="right", va="top", fontsize=9, color=P.INK2)
    P.save(fig, name)


draw("D0", "F1a_mass_fit_D0_M8", fit_d0, r"$m(K\pi)$ [MeV]", "model M8", F.SP["D0"]["m0"], r"$D^0 \to K^-\pi^+$")
draw("Dplus", "F1b_mass_fit_Dplus_A2", fit_dp, r"$m(K\pi\pi)$ [MeV]", "model A2", F.SP["Dplus"]["m0"], r"$D^+ \to K^-\pi^+\pi^+$")
