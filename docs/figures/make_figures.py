#!/usr/bin/env python
"""Figures for docs/METHODS_NOTE_DRAFT.md (project B, charmbkg).

Run:  ~/Documents/charm/.venv_h5/bin/python docs/figures/make_figures.py [F1 F4 ...]
Inputs are read (never written) from ~/Documents/charm/charm_work; outputs: docs/figures/<id>_*.{pdf,png}.
Only figures whose inputs exist locally are produced; nothing is fabricated. Values that are only available as
tabulated numbers in a source document are hard-coded below with the source named next to them.
"""
import json
import sys
import traceback
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.patches  # noqa: E402,F401

CW = Path.home() / "Documents/charm/charm_work"
OUT = Path(__file__).resolve().parent
REPO = OUT.parents[1]

# categorical slots in fixed order (CVD-checked reference palette)
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
GREY = "#8a8985"
M0 = {"D0": 1.86484, "Dplus": 1.86966}  # GeV (PDG)
MK, MPI = 0.493677, 0.13957039

plt.rcParams.update({
    "figure.dpi": 110, "savefig.dpi": 160, "font.size": 9, "axes.titlesize": 9.5, "axes.labelsize": 9,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": "#e4e3df",
    "grid.linewidth": 0.6, "legend.frameon": False, "legend.fontsize": 8, "lines.linewidth": 1.6,
    "errorbar.capsize": 2,
})
TAG = ""


def save(fig, name):
    fig.text(0.995, -0.06, TAG, ha="right", va="top", fontsize=6.5, color=GREY)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"{name}.{ext}", bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


def jload(p):
    return json.loads(Path(p).read_text())


# ---------------------------------------------------------------- F1 mass spectra
def _inv_mass(pt0, eta0, phi0, m_0, pt1, eta1, phi1, m_1):
    def p4(pt, eta, phi, m):
        px, py, pz = pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)
        return np.sqrt(px**2 + py**2 + pz**2 + m**2), px, py, pz
    e0, x0, y0, z0 = p4(pt0, eta0, phi0, m_0)
    e1, x1, y1, z1 = p4(pt1, eta1, phi1, m_1)
    return np.sqrt(np.maximum((e0 + e1) ** 2 - (x0 + x1) ** 2 - (y0 + y1) ** 2 - (z0 + z1) ** 2, 0))


def _load_sub(files, dtype_code):
    keys = ["Dmeson_mass", "Dmeson_pt", "Dmeson_decayType", "Dmeson_truth_category",
            "Dmeson_daughter0_pt", "Dmeson_daughter0_eta", "Dmeson_daughter0_phi",
            "Dmeson_daughter1_pt", "Dmeson_daughter1_eta", "Dmeson_daughter1_phi"]
    cols = {k: [] for k in keys}
    for f in files:
        z = np.load(f, allow_pickle=False)
        sel = z["Dmeson_decayType"] == dtype_code
        for k in keys:
            cols[k].append(z[k][sel])
    return {k: np.concatenate(v) for k, v in cols.items()}


def F1():
    """m(D0), m(D+) — h5 v1.1 subsamples (S_bkggen/sub, S6_kand/sub); MC mb stacked by truth_category."""
    groups = [("signal (0–1)", [0, 1], C[0]), ("partial charm (2–7)", range(2, 8), C[3]),
              ("MisMatch* (8–11)", range(8, 12), C[1]), ("OtherNoB (12)", [12], "#c9c8c2"),
              ("Other*B (13–15)", range(13, 16), C[6])]
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.8))
    specs = [("D0", 10, [CW / "S_bkggen/sub/data_p0.npz"], [CW / f"S_bkggen/sub/mb_p{i}.npz" for i in range(2)],
              (1.75, 2.00), 55, r"$m(K\pi)$ [GeV]"),
             ("Dplus", 4, [CW / "S6_kand/sub/dp_data_p0.npz"], [CW / "S6_kand/sub/dp_mb_p0.npz"],
              (1.70, 2.20), 45, r"$m(K\pi\pi)$ [GeV]")]
    info = {}
    for ax, (sp, code, fdat, fmc, rng, sr, xl) in zip(axs, specs):
        d, mc = _load_sub(fdat, code), _load_sub(fmc, code)
        m0 = M0[sp]
        def sel(a):
            return (a["Dmeson_pt"] >= 3) & (a["Dmeson_pt"] < 30) & (a["Dmeson_mass"] >= rng[0]) & (a["Dmeson_mass"] < rng[1])
        veto_note = ""
        if sp == "D0":
            # swapped-hypothesis mass: recompute both K/pi assignments, keep the one not equal to stored mass
            for a in (d, mc):
                # daughter pT stored in MeV, candidate mass/pT in GeV
                args0 = (a["Dmeson_daughter0_pt"] / 1e3, a["Dmeson_daughter0_eta"], a["Dmeson_daughter0_phi"])
                args1 = (a["Dmeson_daughter1_pt"] / 1e3, a["Dmeson_daughter1_eta"], a["Dmeson_daughter1_phi"])
                mA = _inv_mass(*args0, MK, *args1, MPI)
                mB = _inv_mass(*args0, MPI, *args1, MK)
                stored = a["Dmeson_mass"]
                useB = np.abs(mA - stored) < np.abs(mB - stored)
                a["m_swap"] = np.where(useB, mB, mA)
                a["match"] = np.minimum(np.abs(mA - stored), np.abs(mB - stored)) < 2e-3
            veto_note = "reflection veto |m_swap−m₀|>25 MeV applied"
        bins = np.linspace(*rng, 101)
        sd, sm = sel(d), sel(mc)
        if sp == "D0":
            sd &= np.abs(d["m_swap"] - m0) > 0.025
            sm &= np.abs(mc["m_swap"] - m0) > 0.025
            info["D0_mass_recompute_match_frac"] = float(np.mean(d["match"]))
        hd, _ = np.histogram(d["Dmeson_mass"][sd], bins)
        tc = mc["Dmeson_truth_category"][sm]
        mm = mc["Dmeson_mass"][sm]
        far = np.abs(0.5 * (bins[1:] + bins[:-1]) - m0) > 0.070
        hs = [np.histogram(mm[np.isin(tc, list(g))], bins)[0] for _, g, _ in groups]
        tot = np.sum(hs, axis=0)
        scale = hd[far].sum() / max(tot[far].sum(), 1)
        bottom = np.zeros(len(bins) - 1)
        for (lab, _, col), h in zip(groups, hs):
            ax.stairs(bottom + h * scale, bins, baseline=bottom, fill=True, color=col, label=f"MC mb: {lab}",
                      edgecolor="white", linewidth=0.3)
            bottom = bottom + h * scale
        x = 0.5 * (bins[1:] + bins[:-1])
        ax.errorbar(x, hd, yerr=np.sqrt(hd), fmt="o", ms=2.2, color="#0b0b0b", lw=0.8, label="data (subsample)")
        ax.axvspan(m0 - sr / 1e3, m0 + sr / 1e3, color=C[0], alpha=0.07, lw=0)
        for s in (-1, 1):
            ax.axvline(m0 + s * 0.070, color=GREY, ls="--", lw=0.8)
        ax.text(m0, ax.get_ylim()[1] * 0.97, "SR", ha="center", va="top", fontsize=8)
        ax.text(m0 + 0.072, ax.get_ylim()[1] * 0.97, "far SB →", ha="left", va="top", fontsize=7.5, color=GREY)
        ax.text(m0 - 0.072, ax.get_ylim()[1] * 0.97, "← far SB", ha="right", va="top", fontsize=7.5, color=GREY)
        ax.set_xlabel(xl)
        ax.set_ylabel("candidates / bin")
        ax.set_title(f"{'D⁰→Kπ' if sp == 'D0' else 'D⁺→Kππ'}, pT 3–30 GeV, no classifier cut"
                     + (f"\n{veto_note}" if veto_note else ""))
        info[f"{sp}_MC_scale_to_data_farSB"] = float(scale)
    h, l = axs[1].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=6, fontsize=7, bbox_to_anchor=(0.5, -0.08))
    fig.suptitle("F1  Mass spectra: data vs MC minimum-bias stacked by truth category (MC normalised in far SB, |Δm|>70 MeV)",
                 fontsize=9.5, y=1.02)
    save(fig, "F1_mass_spectra")
    (OUT / "F1_info.json").write_text(json.dumps(info, indent=1))


# ---------------------------------------------------------------- F4 B1 blind quadrant
def F4():
    rows = []
    sets = [("D⁺", CW / "B1_salad_v1/out", ["smooth", "cat_global", "cat_local", "cat_step", "cuts_local"]),
            ("D⁰", CW / "B1_salad_v1/out_D0", ["smooth", "cat_global", "cat_local", "cat_local_MM", "cat_step",
                                                "cuts_local", "real_t27"])]
    label = {"smooth": "smooth in m", "cat_global": "global composition", "cat_local": "local composition",
             "cat_local_MM": "local MisMatch only", "cat_step": "step in m", "cuts_local": "local cuts",
             "real_t27": "truth 2–7 admixture"}
    for sp, d, modes in sets:
        for mode in modes:
            f = d / f"closure_{mode}_hgb.json"
            if not f.exists():
                continue
            j = jload(f)
            r = {"sp": sp, "mode": mode, "c2st": j["c2st_salad_SR"]}
            for wp in ("0.8", "0.7", "0.5"):
                r[wp] = j[f"eff_salad_SR@{wp}"]["all"] / j[f"eff_true_SR@{wp}"]["all"]
            rows.append(r)
    fig, axs = plt.subplots(1, 3, figsize=(11, 3.6), sharey=True)
    for ax, wp in zip(axs, ("0.8", "0.7", "0.5")):
        ax.add_patch(matplotlib.patches.Rectangle((0.495, 0.05), 0.025, 0.35, color=C[7], alpha=0.08, lw=0))
        ax.axvline(0.52, color=GREY, ls="--", lw=0.8)
        ax.axhline(0.05, color=GREY, ls=":", lw=0.8)
        for sp, mk, col in (("D⁺", "o", C[0]), ("D⁰", "s", C[1])):
            rr = [r for r in rows if r["sp"] == sp]
            x = [r["c2st"] for r in rr]
            y = [abs(r[wp] - 1) for r in rr]
            ax.scatter(x, y, marker=mk, s=34, color=col, edgecolor="white", linewidth=1.2, label=sp, zorder=3)
            for r, xi, yi in zip(rr, x, y):
                if r["mode"] in ("cat_local", "smooth") and wp == "0.8":
                    ax.annotate(label[r["mode"]], (xi, yi), xytext=(4, 3), textcoords="offset points", fontsize=6.5,
                                color="#52514e")
        ax.set_title(f"WP{int(float(wp) * 100)}")
        ax.set_xlabel("C2ST AUC (SALAD vs true SR background)")
        ax.set_xlim(0.495, 0.545)
        ax.set_ylim(0, 0.33)
        ax.text(0.497, 0.31, "blind quadrant:\nC2ST passes, ε_B wrong", fontsize=7, color=C[7], va="top")
    axs[0].set_ylabel("|ε_B pred/true − 1|")
    axs[0].legend(loc="upper right")
    fig.suptitle("F4  Adversarial closure (B1): C2ST does not detect background-composition distortions (SALAD v1.0, hgb)",
                 fontsize=9.5, y=1.02)
    save(fig, "F4_B1_blind_quadrant")
    (OUT / "F4_points.json").write_text(json.dumps(rows, indent=1))


# ---------------------------------------------------------------- F6 Zproxy ceiling
def F6():
    # tabulated values: docs/signal_like_background.md §3 / §3a (setup S4, D0 MC closure, WP70); no per-run file kept
    labs = ["SB (nominal)", "SALAD", "FETA (CPU)", "ideal (true bkg)"]
    vals = [4.93, 4.96, 4.83, 4.95]
    fig, ax = plt.subplots(figsize=(5.2, 2.8))
    # dot plot (not bars): the differences are ~1 %, a zero-based bar would hide them, a cut bar would mislead
    cols = [C[0], C[0], C[0], "#52514e"]
    y = np.arange(len(vals))[::-1]
    ax.scatter(vals, y, color=cols, s=40, zorder=3)
    for yi, v in zip(y, vals):
        ax.text(v + 0.006, yi + 0.18, f"{v:.2f}", fontsize=8)
    ax.set_yticks(y, labs)
    ax.axvspan(4.95 * 0.99, 4.95 * 1.01, color=GREY, alpha=0.15, lw=0, label="ideal ±1 %")
    ax.axvline(4.95, color=GREY, ls="--", lw=0.8)
    ax.set_xlim(4.75, 5.05)
    ax.legend(loc="lower left")
    ax.set_xlabel("Zproxy = ε_S/√ε_B at WP70 (D⁰, MC closure, setup S4)")
    ax.grid(axis="y", visible=False)
    ax.set_title("F6  Classifier ceiling: training background vs Zproxy (values from signal_like_background.md §3)")
    save(fig, "F6_zproxy_ceiling")


# ---------------------------------------------------------------- F8 SF vs WP
def F8():
    d = jload(CW / "NUM_audit_2026-11-12/b2b_recompute/out/b2b_sf_new.json")["bins"]
    wps = ["0.95", "0.90", "0.80", "0.70", "0.50"]
    fig, axs = plt.subplots(1, 2, figsize=(9.5, 3.4))
    for k, (pt, col, off) in enumerate((("8-30", C[0], -0.6), ("5-8", C[1], 0.6))):
        wpd = d[pt]["regions"]["nom"]["wp"]
        x = np.array([float(w) * 100 for w in wps]) + off
        for ax, key, ekey in ((axs[0], "SF_sigswap", "SF_sigswap_err"), (axs[1], "SF_rel", "SF_rel_err")):
            y, e, det = [], [], []
            for w in wps:
                pf = wpd[w]["pf"]
                y.append(pf.get(key, np.nan) if pf.get(key) is not None else np.nan)
                e.append(pf.get(ekey, np.nan) if pf.get(ekey) is not None else np.nan)
                det.append(pf.get("determined", False))
            y, e = np.array(y, float), np.array(e, float)
            ax.errorbar(x, y, yerr=e, fmt="o-", color=col, ms=4, label=f"pT {pt.replace('-', '–')} GeV")
            # sW cross-check (unchanged by the 13.11 refit) for SF_rel
    for ax in axs:
        ax.axhline(1, color=GREY, lw=0.8)
        ax.set_xlabel("working point (ε_S in MC, %)")
        ax.set_xticks([50, 70, 80, 90, 95])
        ax.invert_xaxis()
    axs[0].set_ylabel("SF = ε_S(data)/ε_S(MC)")
    axs[1].set_ylabel("SF_rel = SF(WP)/SF(WP95)")
    axs[0].set_title("Absolute SF (pass/fail Δm fit, window ±55 MeV)")
    axs[1].set_title("Relative SF")
    axs[0].legend(loc="lower left")
    axs[0].set_ylim(0.55, 1.25)
    axs[1].set_ylim(0.55, 1.45)
    fig.suptitle("F8  D*-tagged D⁰: classifier scale factors vs WP (full luminosity, refit 13.11; 3–5 GeV not determined)",
                 fontsize=9.5, y=1.03)
    save(fig, "F8_sf_vs_wp")


# ---------------------------------------------------------------- F9 closest_mass vs pdgid
def F9():
    import pandas as pd
    sys.path.insert(0, str(REPO))
    from charmbkg.samples import mark_swaps
    cols = ["run", "event", "m", "pt", "truth", "hypo_swap", "score", "d_pt0", "d_pt1"]
    L = pd.read_parquet(CW / "HYP_d0/runs/lxplus/fd_nd.parquet", columns=cols)
    L = L[(L.m >= 1.75) & (L.m < 2.00) & (L.pt >= 3) & (L.pt < 30)].reset_index(drop=True)
    cfg = {"signal": {"truth_values": [0, 1], "truth_swap_source": "closest_mass"}}
    L0 = L.copy()
    mark_swaps(L0, {"m0": M0["D0"]}, cfg)
    true = L.truth.isin([0, 1]).to_numpy()
    cm_swap = (L0.truth >= 100).to_numpy() & true
    pd_swap = (L.hypo_swap == 1).to_numpy() & true
    stats = {"n_true": int(true.sum()), "n_hypo_swap1": int(pd_swap.sum()),
             "cm_sig_true_swap": int((true & ~cm_swap & pd_swap).sum()),
             "cm_swap_true_correct": int((cm_swap & ~pd_swap & (L.hypo_swap == 0).to_numpy()).sum())}
    bins = np.linspace(-0.11, 0.135, 99)
    fig, axs = plt.subplots(1, 3, figsize=(11, 3.4), sharey=False)
    for ax, (lo, hi) in zip(axs, ((3, 5), (5, 8), (8, 30))):
        b = (L.pt >= lo).to_numpy() & (L.pt < hi).to_numpy() & true
        dm = L.m.to_numpy() - M0["D0"]
        for lab, msk, col, ls in (("signal, closest_mass", b & ~cm_swap, C[0], "-"),
                                  ("signal, PDG-ID truth", b & ~pd_swap & (L.hypo_swap == 0).to_numpy(), C[0], "--"),
                                  ("reflection, closest_mass", b & cm_swap, C[1], "-"),
                                  ("reflection, PDG-ID truth", b & pd_swap, C[1], "--")):
            h, _ = np.histogram(dm[msk], bins)
            ax.stairs(h / max(h.sum(), 1), bins * 1e3, color=col, ls=ls, lw=1.4, label=lab)
        ax.set_yscale("log")
        ax.set_ylim(1e-4, 0.2)
        ax.set_xlabel("m(Kπ) − m₀(D⁰) [MeV]")
        ax.set_title(f"fd MC true D⁰, pT {lo}–{hi} GeV, no cut")
        for s in (-55, 55):
            ax.axvline(s, color=GREY, ls=":", lw=0.8)
    axs[0].set_ylabel("fraction / bin (unit area per curve)")
    h, l = axs[2].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=4, fontsize=7.5, bbox_to_anchor=(0.5, -0.1))
    fig.suptitle("F9  K↔π swap truth: `closest_mass` heuristic narrows the MC signal shape; PDG-ID truth (v1.2) restores the tails",
                 fontsize=9.5, y=1.03)
    save(fig, "F9_refl_truth")
    (OUT / "F9_counts.json").write_text(json.dumps(stats, indent=1))


# ---------------------------------------------------------------- F10 5–8 GeV anomaly (EV)
def F10():
    ev = {"cm": jload(CW / "EV_dstar/lxplus/ev_analysis_cm.json")["bins"],
          "pdgid": jload(CW / "EV_dstar/pdgid/ev_analysis.json")["bins"]}
    pts = ["3-5", "5-8", "8-30"]
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.5))
    ax = axs[0]
    x = np.arange(3)
    for i, (k, col) in enumerate((("cm", C[6]), ("pdgid", C[0]))):
        roc = [ev[k][p]["ring_over_core"] for p in pts]
        ax.errorbar(x + (i - 0.5) * 0.18, [r["data"] for r in roc], yerr=[r["data_err"] for r in roc], fmt="o",
                    color=col, ms=4, label=f"data ({'closest_mass' if k == 'cm' else 'PDG-ID'} MC truth)")
        ax.scatter(x + (i - 0.5) * 0.18, [r["mc_true"] for r in roc], marker="_", s=180, color=col, lw=2,
                   label=f"MC true D⁰ ({'closest_mass' if k == 'cm' else 'PDG-ID'})")
    ax.set_xticks(x, [p.replace("-", "–") + " GeV" for p in pts])
    ax.set_ylabel("N(36–55 MeV) / N(<36 MeV)  [Δm signal]")
    ax.set_title("Ring/core fraction of the D*-tagged m(D⁰) peak")
    ax.legend(fontsize=6.8, loc="upper left")
    ax.set_ylim(-0.05, 1.05)
    ax = axs[1]
    regs = [("core", "core |m−m₀|<36 MeV"), ("ring", "ring 36–55 MeV"), ("nom", "nominal <55 MeV")]
    for i, (k, col) in enumerate((("cm", C[6]), ("pdgid", C[0]))):
        y = [ev[k]["5-8"]["regions"][r]["wp"]["0.95"]["SF"] for r, _ in regs]
        e = [ev[k]["5-8"]["regions"][r]["wp"]["0.95"]["SF_err"] for r, _ in regs]
        ax.errorbar(np.arange(3) + (i - 0.5) * 0.18, y, yerr=e, fmt="o", color=col, ms=4,
                    label="closest_mass" if k == "cm" else "PDG-ID")
    y8 = ev["pdgid"]["8-30"]["regions"]["nom"]["wp"]["0.95"]
    ax.errorbar([2.35], [y8["SF"]], yerr=[y8["SF_err"]], fmt="s", color=C[2], ms=4, label="8–30 GeV, nominal (PDG-ID)")
    ax.axhline(1, color=GREY, lw=0.8)
    ax.set_xticks(range(3), [b for _, b in regs])
    ax.set_ylabel("SF95 = ε95(data)/ε95(MC)")
    ax.set_title("5–8 GeV, WP95, per m(D⁰) region (EV, full luminosity)")
    ax.legend(fontsize=7, loc="lower left")
    ax.set_ylim(0.2, 1.3)
    fig.suptitle("F10  5–8 GeV anomaly: D*-tagged D⁰ (EV_dstar, closest_mass vs PDG-ID swap truth)", fontsize=9.5, y=1.03)
    save(fig, "F10_ev_anomaly")


# ---------------------------------------------------------------- F13 PEAK_pos (new)
def F13():
    import pandas as pd
    r = pd.read_csv(CW / "PEAK_pos/out/peaks_rows.csv")
    t = pd.read_csv(CW / "PEAK_pos/out/dstar_rows.csv")
    s = jload(CW / "PEAK_pos/out/summary.json")
    sel = lambda df: df[(df["var"] == "nom") & (df.wp == 0.7) & (df.q == 0)]
    r, t = sel(r), sel(t)
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.4))
    pts = ["pT3-5", "pT5-8", "pT8-30"]
    series = [("D⁺ m(Kππ), 30 % lumi", r[(r.sp == "Dplus") & (r.veto == True)], C[0], True),  # noqa: E712
              ("D⁰ m(Kπ) from D* tag, full lumi", t[t.obs == "mKpi_tag"], C[1], True),
              ("Δm(D*), full lumi (width only)", t[t.obs == "dm"], C[2], False)]
    used = []
    for n, (lab, df, col, pos) in enumerate(series):
        rr = df.set_index("bin").loc[pts]
        ok = (rr.k > 0.401).to_numpy()  # k=0.4 boundary fits are flagged in RAPORT and not used
        xs = np.arange(3) + (n - 1) * 0.14
        ratio, eratio = rr.seff / rr.seff_mc, rr.e_seff / rr.seff_mc
        if pos:
            dmu = (rr.mu - rr.mu_mc) * 1e3
            edmu = np.hypot(rr.e_mu, rr.e_mu_mc) * 1e3
            axs[0].errorbar(xs[ok], dmu.to_numpy()[ok], yerr=edmu.to_numpy()[ok], fmt="o", color=col, ms=4, label=lab)
        axs[1].errorbar(xs[ok], ratio.to_numpy()[ok], yerr=eratio.to_numpy()[ok], fmt="o", color=col, ms=4, label=lab)
        used.append(lab)
    for ax in axs:
        ax.set_xticks(range(3), ["3–5", "5–8", "8–30"])
        ax.set_xlabel("pT [GeV]")
    axs[0].axhline(0, color=GREY, lw=0.8)
    axs[0].axhspan(-0.5, 0.5, color=C[2], alpha=0.12, lw=0, label="recommended ±0.5 MeV")
    axs[0].set_ylabel("Δμ = μ_data − μ_MC [MeV]")
    axs[0].set_title("Peak position, WP70")
    axs[1].axhline(1, color=GREY, lw=0.8)
    axs[1].axhspan(0.9, 1.1, color=C[3], alpha=0.12, lw=0, label="placeholder ±10 %")
    axs[1].set_ylabel("σ_eff(data)/σ_eff(MC)")
    axs[1].set_title("Peak width, WP70")
    h, l = axs[1].get_legend_handles_labels()
    h0, l0 = axs[0].get_legend_handles_labels()
    fig.legend(h + h0[:1], l + l0[:1], loc="lower center", ncol=5, fontsize=7.5, bbox_to_anchor=(0.5, -0.1))
    fig.suptitle("F13  Mass-peak position and width, data vs MC (PEAK_pos, nominal fit model, both charges)",
                 fontsize=9.5, y=1.03)
    save(fig, "F13_peak_pos_width")
    (OUT / "F13_info.json").write_text(json.dumps({"species_rows_plotted": used, "summary_keys": list(s)}, indent=1))


# ================================================================ v2 (17.11): trigger classes, track jets, AFP mixing
# One colour per sample, the same in F14–F16 (colour follows the sample, not its position in a panel).
BT = CW / "B_trig"
CLS = [  # (run, class key in btrig_*.json, class key in trkjet.json, label, colour, marker)
    (428770, "all", "8:all", "428770, all triggers", "#52514e", "D"),
    (428770, "J12only", "8:J12only", "428770, L1_J12 only", C[1], "s"),
    (428770, "MBTSbased", "8:MBTSbased", "428770, \"MBTS-based\" (93 % HLT_j20)", C[2], "^"),
    (427929, "MB427929", "9:MB", "427929, MBTS (minimum bias)", C[0], "o"),
]
PT3 = ["3-5", "5-8", "8-30"]


def _rows(path):
    return jload(path)["rows"]


def _pick(rows, **kw):
    out = [r for r in rows if all(str(r[k]) == str(v) for k, v in kw.items())]
    if len(out) != 1:
        raise RuntimeError(f"expected one row for {kw}, got {len(out)}")
    return out[0]


def _legend_below(fig, ax, ncol=4, y=-0.07):
    h, l = ax.get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=ncol, fontsize=7.5, bbox_to_anchor=(0.5, y))


def F14():
    """Isolation SF per trigger class and run: 3 cuts x pT bins. Input: B_trig/runs/lxplus/out/btrig_iso.json."""
    rows = _rows(BT / "runs/lxplus/out/btrig_iso.json")
    cuts = [("ptcone40<2.0", "ptcone40 < 2 GeV"), ("ptcone40<4.0", "ptcone40 < 4 GeV"),
            ("rel40<0.5", "ptcone40 / pT < 0.5")]
    fig, axs = plt.subplots(1, 3, figsize=(11, 3.5), sharey=True)
    pts, dump = [], []
    for ax, (cut, title) in zip(axs, cuts):
        for n, (run, cls, _, lab, col, mk) in enumerate(CLS):
            y, e, x = [], [], []
            for i, b in enumerate(PT3):
                r = _pick(rows, run=run, cls=cls, bin=b, cut=cut)
                if not r["determined"]:
                    continue
                x.append(i + (n - 1.5) * 0.15)
                y.append(r["SF"])
                e.append(r["SF_err"])
                dump.append({"cut": cut, "run": run, "cls": cls, "bin": b, "SF": r["SF"], "SF_err": r["SF_err"]})
            ax.errorbar(x, y, yerr=e, fmt=mk, color=col, ms=4.5, lw=1.1, label=lab)
        ax.axhline(1, color=GREY, lw=0.9)
        ax.set_xticks(range(3), [b.replace("-", "–") for b in PT3])
        ax.set_xlabel("pT(D⁰) [GeV]")
        ax.set_title(title)
        ax.grid(axis="x", visible=False)
        pts.append(ax)
    axs[0].set_ylabel("SF_iso = ε(data)/ε(MC)")
    axs[0].set_ylim(0.2, 1.3)
    _legend_below(fig, axs[0])
    fig.suptitle("F14  Isolation scale factor of D*-tagged D⁰ per trigger class and run (score_H > 0.6, |m−m₀| < 40 MeV).\n"
                 "428770 class columns use trigger-emulation MC weights that do not close (model dependent); "
                 "427929 needs none", fontsize=9, y=1.07)
    save(fig, "F14_iso_sf_trigger")
    (OUT / "F14_points.json").write_text(json.dumps(dump, indent=1))


def F15():
    """Classifier-H SF (WP95/80/50) per trigger class and run, 5-8 and 8-30 GeV. Input: btrig_sf.json."""
    rows = _rows(BT / "runs/lxplus/out/btrig_sf.json")
    wps = ["0.95", "0.80", "0.50"]
    fig, axs = plt.subplots(1, 2, figsize=(9.5, 3.5), sharey=True)
    dump = []
    for ax, b in zip(axs, ("5-8", "8-30")):
        for n, (run, cls, _, lab, col, mk) in enumerate(CLS):
            x, y, e = [], [], []
            for i, w in enumerate(wps):
                r = _pick(rows, run=run, cls=cls, bin=b, wp=w)
                if not r["determined"]:
                    continue
                x.append(i + (n - 1.5) * 0.15)
                y.append(r["SF"])
                e.append(r["SF_err"])
                dump.append({"bin": b, "wp": w, "run": run, "cls": cls, "SF": r["SF"], "SF_err": r["SF_err"]})
            ax.errorbar(x, y, yerr=e, fmt=mk, color=col, ms=4.5, lw=1.1, label=lab)
        ax.axhline(1, color=GREY, lw=0.9)
        ax.set_xticks(range(3), ["WP95", "WP80", "WP50"])
        ax.set_xlabel("working point (ε_S in MC)")
        ax.set_title(f"pT {b.replace('-', '–')} GeV")
        ax.grid(axis="x", visible=False)
    axs[0].set_ylabel("SF = ε_S(data)/ε_S(MC)")
    axs[0].set_ylim(0.4, 1.6)
    _legend_below(fig, axs[0])
    fig.suptitle("F15  Classifier-H scale factor per trigger class and run (Δm pass/fail fit, v2.0 files; 3–5 GeV not determined).\n"
                 "428770 class points use emulated MC efficiencies (model dependent)", fontsize=9, y=1.07)
    save(fig, "F15_clf_sf_trigger")
    (OUT / "F15_points.json").write_text(json.dumps(dump, indent=1))


def F16():
    """Track-jet test: fraction of signal D0 in / near a hard track jet, and SF_iso for D0 with no nearby jet.
    Input: B_trig/TRKJET/out/trkjet.json (keys a, b)."""
    t = jload(BT / "TRKJET/out/trkjet.json")
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 3.5), gridspec_kw={"width_ratios": [1.5, 1]})
    ax = axs[0]
    props = [("A10", "D⁰ inside a track jet\nwith pT > 10 GeV [uses cone tracks]"),
             ("O5", "jet pT > 5 GeV built from tracks\noutside the cone, within ΔR < 0.8")]
    dump = {"fractions": [], "sf_outside": []}
    for i, (key, _) in enumerate(props):
        a = t["a"][key]
        for n, (_, _, tk, lab, col, mk) in enumerate(CLS):
            c = a["cls"][tk]
            ax.errorbar([i + (n - 1.5) * 0.15], [c["eps"]], yerr=[c["err"]], fmt=mk, color=col, ms=4.5, lw=1.1,
                        label=lab if i == 0 else None)
            dump["fractions"].append({"prop": key, "cls": tk, "f": c["eps"], "err": c["err"]})
        m = a["MC fd_nd"]
        ax.plot([i - 0.38, i + 0.38], [m["eps"], m["eps"]], color=GREY, lw=1.4, ls="--",
                label="MC fd_nd (no trigger)" if i == 0 else None)
        dump["fractions"].append({"prop": key, "cls": "MC fd_nd", "f": m["eps"], "err": m["err"]})
    ax.set_xticks(range(2), [p[1] for p in props], fontsize=7.5)
    ax.set_ylabel("fraction of D⁰←D* signal")
    ax.set_ylim(0, 0.75)
    ax.set_title("Where the D⁰ sits relative to track jets")
    ax.grid(axis="x", visible=False)
    ax = axs[1]
    for n, (_, _, tk, lab, col, mk) in enumerate(CLS):
        b = t["b"][f"O5|poza|{tk}"]
        fit = b.get("fit") or {}
        eps, err = (fit["eps"], fit["err"]) if fit.get("det") and fit.get("valid") else (b["eps"], b["err"])
        sf = eps / b["mc"]["eps"]
        esf = sf * np.hypot(err / eps, b["mc"]["err"] / b["mc"]["eps"])
        ax.errorbar([n], [sf], yerr=[esf], fmt=mk, color=col, ms=5, lw=1.1)
        ax.text(n + 0.12, sf, f"{sf:.2f}", fontsize=7.5, va="center", color="#52514e")
        dump["sf_outside"].append({"cls": tk, "SF": sf, "SF_err": esf, "method": "fit" if fit.get("det") else "SBS"})
    ax.axhline(1, color=GREY, lw=0.9)
    ax.set_xticks(range(4), ["428770\nall", "428770\nL1_J12", "428770\n\"MBTS\"", "427929\nMBTS"], fontsize=7.5)
    ax.set_xlim(-0.5, 3.6)
    ax.set_ylim(0.2, 1.3)
    ax.set_ylabel("SF_iso (ptcone40 < 2 GeV)")
    ax.set_title("D⁰ with no such nearby jet (\"outside\")")
    ax.grid(axis="x", visible=False)
    _legend_below(fig, axs[0], ncol=5)
    fig.suptitle("F16  Track-jet test (anti-kT R = 0.4 track jets, pT 3–30 GeV): the D⁰ is more often in/near a hard jet in "
                 "L1_J12 events,\nbut the isolation deficit persists for D⁰ without a nearby jet", fontsize=9, y=1.07)
    save(fig, "F16_trkjet")
    (OUT / "F16_points.json").write_text(json.dumps(dump, indent=1))


def F17():
    """AFP proton mixing v2 (PRELIMINARY): f_acc per run and side, and relative error budget of N_acc.
    Input: AFP_mix/runs/lxplus/out_v2/result_v2.json. Background fractions only; no signal yields."""
    res = jload(CW / "AFP_mix/runs/lxplus/out_v2/result_v2.json")
    if res.get("status") != "WSTEPNE":
        raise RuntimeError("result_v2.json status changed; revisit the PRELIMINARY label")
    fig, axs = plt.subplots(1, 3, figsize=(12.5, 3.5), gridspec_kw={"width_ratios": [1, 1, 1.3], "wspace": 0.55})
    regs = [("A", "side A"), ("C", "side C"), ("AorC", "A ∨ C")]
    comps = [("stat", "statistical"), ("pool", "donor pool"), ("extrapolation", "extrapolation in ξ_cal"),
             ("half_diff_methods", "½ Δ(methods)"), ("total_rel", "total")]
    runcol = {427929: C[0], 428770: C[1]}
    dump = []
    for ax, r in zip(axs[:2], res["results"]):
        run, col = r["run_number"], runcol[r["run_number"]]
        syst = r["syst_rel"]["syst_no_stat"]
        for i, (k, lab) in enumerate(regs):
            e = r["estimates_by_region"][k]
            f, es = e["f_acc"], e["f_acc_err_stat"]
            et = np.hypot(es, f * syst)
            ax.errorbar([i], [f], yerr=[et], fmt="none", ecolor=col, elinewidth=1.0, capsize=3, alpha=0.55)
            ax.errorbar([i], [f], yerr=[es], fmt="o", color=col, ms=5, elinewidth=2.2, capsize=0)
            dump.append({"run": run, "region": k, "f_acc": f, "err_stat": es, "err_total": et, "n_tag": e["n_tag"]})
        ax.set_xticks(range(3), [lab for _, lab in regs])
        ax.set_xlim(-0.5, 2.5)
        ax.set_ylim(0, None)
        ax.set_title(f"run {run}: f_acc = N_acc/N_tag")
        ax.grid(axis="x", visible=False)
        ax.text(0.03, 0.04, "thick: stat; thin: stat ⊕ common syst\n(syst from A ∨ C, relative "
                f"{100 * syst:.0f} %)", transform=ax.transAxes, fontsize=6.8, color="#52514e")
    axs[0].set_ylabel("accidental-tag fraction f_acc")
    ax = axs[2]
    for n, r in enumerate(res["results"]):
        run = r["run_number"]
        y = np.arange(len(comps))[::-1] + (0.5 - n) * 0.3
        v = [100 * r["syst_rel"][k] for k, _ in comps]
        ax.barh(y, v, height=0.26, color=runcol[run], label=f"run {run}")
        for yi, vi in zip(y, v):
            ax.text(vi + 0.8, yi, f"{vi:.1f}", va="center", fontsize=7, color="#52514e")
    ax.set_yticks(np.arange(len(comps))[::-1], [lab for _, lab in comps])
    ax.set_xlabel("relative uncertainty of N_acc (A ∨ C) [%]")
    ax.set_xlim(0, 60)
    ax.set_title("Error budget")
    ax.grid(axis="y", visible=False)
    ax.legend(loc="center right", bbox_to_anchor=(1.0, 0.72))
    fig.suptitle("F17  PRELIMINARY — accidental AFP-tag fraction from proton mixing (v2, estimator A), non-AFP triggers only\n"
                 "(in 428770 this covers 4.4 % of all tags; events from L1AFP-only chains are excluded)", fontsize=9, y=1.07)
    save(fig, "F17_afp_mix_facc")
    (OUT / "F17_points.json").write_text(json.dumps(dump, indent=1))


ALL = {"F1": F1, "F4": F4, "F6": F6, "F8": F8, "F9": F9, "F10": F10, "F13": F13,
       "F14": F14, "F15": F15, "F16": F16, "F17": F17}

if __name__ == "__main__":
    which = sys.argv[1:] or list(ALL)
    bad = []
    for k in which:
        try:
            ALL[k]()
        except Exception:
            traceback.print_exc()
            bad.append(k)
    if bad:
        print("FAILED:", bad)
        sys.exit(1)
