#!/usr/bin/env python3
"""LOKALNIE (local_tagged/) — pomiar z protonem AFP, run 427929.
T1: widmo m(Kpi) D0 w zdarzeniach z tagiem: fit nominalny (TAGGED out/hists.npz + unblind.json, bez nowego fitu) i wariant z ksztaltem tla z probki nietagowanej
    (KONSTR B: kod kb_tag.py importowany bez zmian, fit odtworzony i porownany asercja z out/kb_window.json).
T2: R^fid i sigma_tag^fid wobec Pythii (zakres) i literatury (TAGXSEC out/tagxsec_excl.json; zakresy literatury parsowane z TAGXSEC RAPORT.md, tabela §5)."""
import json
import math
import re
import sys

import numpy as np

import pstyle as P
from pstyle import plt

sys.path.insert(0, P.CW + "/KONSTR_2026-11-20/skrypty")
q = lambda *a: math.sqrt(sum(x * x for x in a))  # noqa: E731
KEY = "D0|427929|MBTS|3-30|0.7|T2"
TXT = ""


def reb(v, k=2):
    n = (len(v) // k) * k
    return v[:n].reshape(-1, k).sum(1)


def t1():
    import kb_tag as K
    T = K.T
    U = json.load(open(P.CW + "/TAGGED_2026-11-19/out/unblind.json"))["items"][KEY]
    KW = P.J("KONSTR_2026-11-20", "out", "kb_window.json")
    Hs = np.load(P.CW + "/TAGGED_2026-11-19/out/hists.npz")
    kq = KEY.replace("|", "__")
    H, mu_n, bk_n = Hs[kq + "__H"], Hs[kq + "__mu"], Hs[kq + "__bkg"]
    # wariant: odtworzenie fitu KONSTR B (konstrukcja nominalna, wagi: brak)
    sp = "D0"
    R = K.load_incl_rows(); TG = K.tagged(blind=False)
    r = R["sp"] == 0; t = TG[sp]
    Rm, Rn, Ri = R["m"][r], R["nsel"][r], R["ipt"][r]
    I, _ = T.load_incl(K.RUN); fits = T.incl_fit(I, sp, "MB_npu0", K.E0)
    sb = K.sbmask_rows(t["m"], sp)
    w = K.weights(Rm, Rn, Ri, t["m"][sb], t["nsel"][sb], t["ipt"][sb], sp, "brak")
    B, _, _ = K.incl_model(Rm, Ri, w, sp, T.MODEL[sp], None)
    S = T.sig_shape(sp, fits); Hh = T.hist5(t["m"], sp)
    assert np.abs(Hh - H).max() == 0
    fv = K.CFit(S, B).full(Hh, cl95=False)
    kv = next(v for k, v in KW.items() if k.startswith("D0|konstrukcja nominalna"))
    assert abs(fv["N"] - kv["N"]) < 0.05 and abs(fv["NB"] - kv["NB"]) < 0.5, (fv["N"], kv["N"])
    bk_v = kv["NB"] * B / B.sum(); mu_v = kv["N"] * S + bk_v
    fn = U["fit"]
    assert abs((mu_n - bk_n).sum() - fn["N"] * S.sum()) < 0.5
    e = T.edges5(sp) * 1e3
    e2 = e[::2]; c2 = 0.5 * (e2[1:] + e2[:-1])
    m0, wn = T.WIN[sp]
    ws = KW["D0|okno−pasma"]
    fig = plt.figure(figsize=(13.4, 7.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[3.2, 1], hspace=0.06, wspace=0.17, left=0.065, right=0.985, top=0.79, bottom=0.085)
    from scipy.stats import chi2 as C2
    defs = [("Nominal: free background exp(Chebyshev, order from the tagged sidebands)", mu_n, bk_n,
             f"$N_{{tag}}$ = {fn['N']:.1f} $^{{+{U['stat_hi']:.1f}}}_{{-{U['stat_lo']:.1f}}}$ (stat) ± {U['syst']['razem']:.1f} (shape)\nZ = {fn['Z']:.2f} (likelihood ratio);  p(toys) = {U['p_toys']:.4f}"),
            ("Variant: background shape fixed from the untagged sample (post hoc)", mu_v, bk_v,
             f"$N_{{tag}}$ = {kv['N']:.1f} $^{{+{kv['hi68'] - kv['N']:.1f}}}_{{-{kv['N'] - kv['lo68']:.1f}}}$ (stat)\nZ = {kv['Z']:.2f};  " + rf"$\chi^2$/ndf = {kv['chi2']:.1f}/{kv['ndf']} (p = {C2.sf(kv['chi2'], kv['ndf']):.2f})")]
    for ic, (ttl, mu, bk, txt) in enumerate(defs):
        a0, a1 = fig.add_subplot(gs[0, ic]), fig.add_subplot(gs[1, ic])
        h, m_, b_ = reb(H), reb(mu), reb(bk)
        lo = np.sqrt(h + 0.25) - 0.5; hi = np.sqrt(h + 0.25) + 0.5
        a0.axvspan((m0 - wn) * 1e3, (m0 + wn) * 1e3, color=P.GRID, alpha=0.45, lw=0, label="Mass window")
        a0.stairs(m_, e2, baseline=b_, fill=True, color=P.C["blue"], alpha=0.35, lw=0, label=r"Signal + tied components ($K\leftrightarrow\pi$ reflection, other charm decays)")
        a0.stairs(m_, e2, color=P.C["blue"], lw=1.8, label="Fit")
        a0.stairs(b_, e2, color=P.C["orange"], lw=1.6, ls="--", label="Background")
        a0.errorbar(c2, h, [lo, hi], fmt="o", ms=4, color=P.DATA, lw=0.9, label="Data, events with an AFP proton tag", zorder=5)
        a0.set_ylim(0, h.max() * 1.75); a0.set_xlim(e[0], e[-1]); a0.tick_params(labelbottom=False)
        a0.set_ylabel("Candidates / 10 MeV"); a0.set_title(ttl, loc="left", fontsize=11)
        a0.text(0.975, 0.96, txt, transform=a0.transAxes, ha="right", va="top", fontsize=10.5, linespacing=1.5)
        if ic == 0:
            fig.legend(*a0.get_legend_handles_labels(), loc="upper center", bbox_to_anchor=(0.525, 0.885), ncol=5, fontsize=9.3, columnspacing=1.4, handlelength=1.6)
        else:
            a0.text(0.975, 0.79, f"Counting check (window − background extrapolated\nfrom the sidebands): {ws['okno']:.0f} − {ws['tlo']:.0f} → N = {ws['N']:.1f}",
                    transform=a0.transAxes, ha="right", va="top", fontsize=9.5, color=P.INK2)
        pull = (h - m_) / np.sqrt(m_)
        a1.axhspan(-1, 1, color=P.GRID, alpha=0.6, lw=0); a1.axhline(0, color=P.INK2, lw=0.7)
        a1.bar(c2, pull, width=8.5, color=P.INK2, lw=0)
        a1.set_ylim(-3.6, 3.6); a1.set_yticks([-2, 0, 2]); a1.set_xlim(e[0], e[-1])
        a1.set_ylabel("Pull"); a1.set_xlabel(r"$m(K\pi)$ [MeV]")
    P.header(fig, P.run427_text(), extra="exactly one AFP proton, 0.030 ≤ ξ ≤ 0.100")
    fig.text(0.985, 0.985, TXT, ha="right", va="top", fontsize=11, fontweight="bold", color=P.C["red"])
    fig.text(0.01, 0.91, r"$D^0 \to K\pi$, 3 ≤ $p_T$ < 30 GeV, tree H at WP70, no pile-up vertex, $\xi_{cal}$ cut.  Fit in 5 MeV bins (drawn in 10 MeV); signal shape fixed from the inclusive M8 fit. "
             f"{int(H.sum())} candidates, {U['n_win']:.0f} in the window.", fontsize=9.3, color=P.INK2)
    P.save(fig, "T1_tagged_D0_mass", tagged=True)


def lit_ranges():
    """Zakresy literatury z tabeli TAGXSEC RAPORT.md: {'LMS': ((Rlo, Rhi) [%], (slo, shi) [ub]), 'SS': ((Rlo, Rhi), None)}."""
    out = {}
    rng = lambda s: tuple(float(x.replace(",", ".")) for x in re.search(r"(\d+,\d+)–(\d+,\d+)", s).groups()) if re.search(r"\d+,\d+–\d+,\d+", s) else None  # noqa: E731
    for ln in open(P.CW + "/TAGXSEC_2026-11-20/RAPORT.md"):
        c = [x.strip() for x in ln.split("|")]
        if len(c) > 3 and c[1].startswith("LMS 1412.3132"):
            out["LMS"] = (rng(c[2]), rng(c[3]))
        if len(c) > 3 and c[1].startswith("Siddikov–Schmidt 2008.12446"):
            out["SS"] = (rng(c[2]), rng(c[3]))
    return out


def t2():
    X = P.J("TAGXSEC_2026-11-20", "out", "tagxsec_excl.json")
    W, PR = X["wynik"]["D0"], X["pred"]["D0"]
    LIT = lit_ranges()
    print("literatura:", LIT)
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(13.4, 5.6), gridspec_kw=dict(left=0.285, right=0.985, top=0.73, bottom=0.12, wspace=0.08))
    rows = ["Data, nominal fit", "Data, variant: window − sidebands", "Pythia 8 A3, soft diffraction\n(point: nominal; bar: envelope of settings)", "Pythia 8, hard diffraction scenario\n(bar: envelope)",
            "Łuszczak–Maciuła–Szczurek\n(arXiv:1412.3132, converted to this ξ window)", "Siddikov–Schmidt\n(arXiv:2008.12446, converted)"]
    y = np.arange(len(rows))[::-1]
    for ax in (a0, a1):
        ax.set_ylim(-0.6, len(rows) - 0.4); ax.set_yticks(y)
        ax.tick_params(axis="y", which="both", left=False, right=False)
        ax.axhspan(y[1] - 0.5, y[0] + 0.5, color=P.GRID, alpha=0.35, lw=0)
        ax.grid(axis="x", which="major")
        ax.set_axisbelow(True)
    a0.set_yticklabels(rows, fontsize=10); a1.set_yticklabels([])

    def data(ax, yy, v, slo, shi, syst, sc, unit):
        ax.errorbar(sc * v, yy, xerr=[[sc * q(slo, syst)], [sc * q(shi, syst)]], fmt="none", ecolor=P.DATA, elinewidth=1.0)
        ax.errorbar(sc * v, yy, xerr=[[sc * slo], [sc * shi]], fmt="o", ms=7, color=P.DATA, elinewidth=2.8)
        ax.annotate(f"{sc * v:.2f} $^{{+{sc * shi:.2f}}}_{{-{sc * slo:.2f}}}$ ± {sc * syst:.2f}{unit}" if sc > 1 else f"{v:.1f} $^{{+{shi:.1f}}}_{{-{slo:.1f}}}$ ± {syst:.1f}{unit}",
                    (sc * v, yy + 0.2), ha="center", va="bottom", fontsize=9.5)

    def band(ax, yy, lo, hi, col, c=None, mk="s", txt=None):
        ax.plot([lo, hi], [yy, yy], color=col, lw=7, alpha=0.45, solid_capstyle="butt")
        if c is not None:
            ax.plot(c, yy, mk, ms=8, color=col)
        ax.annotate(txt, (math.sqrt(lo * hi), yy + 0.2), ha="center", va="bottom", fontsize=9.5)
    for i, k in enumerate(("nom", "var")):
        d = W[k]
        data(a0, y[i], d["R"], d["R_stat_lo"], d["R_stat_hi"], d["R"] * d["R_syst_rel"], 100, " %")
        data(a1, y[i], d["sigma_ub"], d["sig_stat_lo"], d["sig_stat_hi"], d["sigma_ub"] * d["syst_rel"], 1, " µb")
    band(a0, y[2], 100 * PR["R_soft_lo"], 100 * PR["R_soft_hi"], P.C["blue"], 100 * PR["R_soft"], txt=f"{100 * PR['R_soft']:.2f} % ({100 * PR['R_soft_lo']:.2f}–{100 * PR['R_soft_hi']:.2f})")
    band(a0, y[3], 100 * PR["R_hard_lo"], 100 * PR["R_hard_hi"], P.C["violet"], 100 * PR["R_hard"], "D", txt=f"{100 * PR['R_hard']:.2f} % ({100 * PR['R_hard_lo']:.2f}–{100 * PR['R_hard_hi']:.2f})")
    band(a0, y[4], *LIT["LMS"][0], P.C["aqua"], txt=f"≈ {LIT['LMS'][0][0]:g}–{LIT['LMS'][0][1]:g} %")
    band(a0, y[5], *LIT["SS"][0], P.C["yellow"], txt=f"≈ {LIT['SS'][0][0]:g}–{LIT['SS'][0][1]:g} %")
    pi = PR["py_incl_ub"]
    band(a1, y[2], PR["R_soft_lo"] * pi, PR["R_soft_hi"] * pi, P.C["blue"], PR["R_soft"] * pi, txt=f"{PR['R_soft'] * pi:.1f} µb ({PR['R_soft_lo'] * pi:.1f}–{PR['R_soft_hi'] * pi:.0f}), with Pythia's own σ(incl)")
    band(a1, y[4], *LIT["LMS"][1], P.C["aqua"], txt=f"≈ {LIT['LMS'][1][0]:g}–{LIT['LMS'][1][1]:g} µb (collinear LO)")
    for yy in (y[3], y[5]):
        a1.text(0.03, yy, "not available", transform=a1.get_yaxis_transform(), va="center", fontsize=9.5, color=P.INK2)
    a0.set_xscale("log"); a1.set_xscale("log")
    a0.set_xlim(0.08, 4.0); a1.set_xlim(0.6, 60)
    from matplotlib.ticker import FuncFormatter
    for ax in (a0, a1):
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        ax.xaxis.set_minor_formatter(FuncFormatter(lambda v, _: f"{v:g}" if any(abs(v / 10 ** math.floor(math.log10(v)) - k) < 1e-6 for k in (2, 5)) else ""))
        ax.tick_params(axis="x", which="minor", labelsize=9.5)
    a0.set_xlabel(r"$R^{fid} = \sigma^{fid}_{tag} / \sigma_{incl}$  [%]"); a1.set_xlabel(r"$\sigma^{fid}_{tag}(D^0 + \bar{D}^0)$  [µb]")
    a0.set_title("Tagged fraction", loc="left"); a1.set_title("Fiducial tagged cross-section", loc="left")
    P.header(fig, P.run427_text(), extra="post hoc with respect to the read-out of the tagged peaks")
    fig.text(0.985, 0.985, TXT, ha="right", va="top", fontsize=11, fontweight="bold", color=P.C["red"])
    fig.text(0.01, 0.865, r"$D^0 + \bar{D}^0$, 3 ≤ $p_T$ < 30 GeV, |η| < 2.2, with a proton of 0.030 ≤ ξ ≤ 0.100 on side A or C (no t cut).  Data: thick bar = stat, thin bar = stat ⊕ syst; numbers: value (stat) ± syst.  "
             "\nThe difference between the nominal fit and the window − sidebands variant is a one-sided model systematic that is not included in ± syst.", fontsize=9.3, color=P.INK2, va="top")
    P.save(fig, "T2_tagged_R_sigma", tagged=True)


t1(); t2()
