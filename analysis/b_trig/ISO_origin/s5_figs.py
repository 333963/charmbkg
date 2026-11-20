#!/usr/bin/env python3
"""ISO_origin: rysunki PNG z out/s1.json i out/s2.json."""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

H = os.path.dirname(os.path.abspath(__file__))
s1 = json.load(open(f"{H}/out/s1.json")); s2 = json.load(open(f"{H}/out/s2.json"))
COL = {"9:MB": "#2a78d6", "8:J12only": "#eb6834", "8:MBTSbased": "#1baf7a", "MC": "#6250d6", "8:other": "#eda100"}
LAB = {"9:MB": "427929 MB (MBTS_1)", "8:J12only": "428770 J12only", "8:MBTSbased": "428770 „MBTSbased” (= HLT_j20)", "MC": "MC fd_nd (bez wag)", "8:other": "428770 other"}
MK = {"9:MB": "o", "8:J12only": "s", "8:MBTSbased": "^", "MC": "D", "8:other": "v"}
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": "#e6e5e0", "grid.linewidth": 0.6,
                     "axes.edgecolor": "#8a897f", "legend.frameon": False, "figure.dpi": 140})
INK = "#3d3d3a"


def g(d, k):
    v = d.get(k)
    return np.nan if v is None else v


# ---- fig1: ta sama klasa triggera w obu runach (sygnał i tło)
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
pairs = [("L1_J12", "9:MB&J12", "8:J12any"), ("HLT_j20 ∧ ¬J12", "9:MB&j20&!J12", "8:j20&!J12"), ("bez j20, bez J12\n(427929) / other (428770)", "9:MB&!j20&!J12", "8:other"),
         ("cały run", "9:MB", "8:all")]
for a, (key, ttl) in zip(ax, (("sig", "Sygnał D*→D0π (odejmowanie pasm Δm), score_H > 0,6"), ("bkg", "Tło kombinatoryczne (Δm > 150 MeV), bez cięcia score"))):
    for i, (lab, k9, k8) in enumerate(pairs):
        for off, k, c, mk in ((-0.12, k9, COL["9:MB"], "o"), (0.12, k8, COL["8:J12only"], "s")):
            r = s1["classes"][k + "|3-30"]
            if key == "sig":
                y, e = g(r, "eps"), g(r, "err"); n = r["S"]
            else:
                y, e = r["bkg_noscore"]["f"], r["bkg_noscore"]["err"]; n = r["bkg_noscore"]["n"]
            weak = key == "sig" and n < 50
            a.errorbar(i + off, y, yerr=e, fmt=mk, color=c, ms=7, capsize=3, lw=1.5, mfc="white" if weak else c)
            a.annotate(f"{n:.0f}", (i + off * 1.6, (y + e + 0.025) if np.isfinite(y) else 0), ha="right" if off < 0 else "left", fontsize=7.5, color=INK)
    a.set_xticks(range(len(pairs))); a.set_xticklabels([p[0] for p in pairs], fontsize=9); a.set_title(ttl, fontsize=10, color=INK)
    a.set_ylim(0, 1.0)
ax[0].axhline(s1["mc_eps_incl"], color=COL["MC"], lw=1.2, ls="--"); ax[0].annotate("MC bez wag", (3.45, s1["mc_eps_incl"] + 0.015), color=INK, fontsize=8, ha="right")
ax[0].set_ylabel("ε(ptcone40 < 2 GeV), pT 3–30 GeV")
ax[1].errorbar([], [], fmt="o", color=COL["9:MB"], label="run 427929 (μ ≈ 0,005)"); ax[1].errorbar([], [], fmt="s", color=COL["8:J12only"], label="run 428770 (μ ≈ 0,05)")
ax[1].errorbar([], [], fmt="o", color=INK, mfc="white", label="N_sig < 50 (niemiarodajne)"); ax[1].legend(loc="upper left", fontsize=8.5)
fig.suptitle("Ta sama klasa triggera → ta sama izolacja w obu runach (liczby: N_sig / N_kand.)", fontsize=11, color=INK)
fig.tight_layout(); fig.savefig(f"{H}/fig1_trigger_vs_run.png"); plt.close(fig)

# ---- fig2: ε vs zmienna (uniwersalność) — 4 panele
names = {"e_nSelBL": "nSelBL", "e_lead_pt": "leadTrackPtBL [GeV]", "pt": "pT(D0) [GeV]", "z": "z = pT(D0)/leadTrackPtBL"}
fig, axs = plt.subplots(2, 2, figsize=(11, 7.6))
for a, var in zip(axs.ravel(), names):
    tab = s2["nominal"]["uni"][var]["table"]
    for j, nm in enumerate(["MC", "9:MB", "8:J12only", "8:MBTSbased"]):
        rows = tab[nm]; x = np.arange(len(rows)) + (j - 1.5) * 0.13
        y = np.array([np.nan if r["eps"] is None else r["eps"] for r in rows], float); e = np.array([r["err"] or 0 for r in rows], float)
        S = np.array([r["S"] for r in rows], float)
        ok = S >= 50
        a.errorbar(x[ok], y[ok], yerr=e[ok], fmt=MK[nm], color=COL[nm], ms=5.5, capsize=2, lw=1.3, label=LAB[nm])
        bad = ~ok & np.isfinite(y) & (S > 15)
        a.errorbar(x[bad], y[bad], yerr=e[bad], fmt=MK[nm], color=COL[nm], ms=5.5, capsize=2, lw=0.8, mfc="white", alpha=0.7)
    ed = [r["lo"] for r in tab["MC"]] + [tab["MC"][-1]["hi"]]
    a.set_xticks(range(len(ed) - 1)); a.set_xticklabels([f"{lo:g}–{hi:g}" if hi < 1e8 else f"≥{lo:g}" for lo, hi in zip(ed[:-1], ed[1:])], fontsize=8)
    a.set_xlabel(names[var]); a.set_ylim(-0.05, 1.05); a.set_ylabel("ε(ptcone40 < 2 GeV)")
    c = s2["nominal"]["uni"][var]["chi2"]
    a.set_title(f"χ²/ndf: 427929 vs J12only {c['9:MB vs 8:J12only']['chi2']:.0f}/{c['9:MB vs 8:J12only']['ndf']}, 427929 vs MC {c['9:MB vs MC']['chi2']:.1f}/{c['9:MB vs MC']['ndf']}", fontsize=9, color=INK)
axs[0, 0].legend(fontsize=8, loc="upper right")
fig.suptitle("ε_iso sygnału w danych nie jest uniwersalną funkcją żadnej zmiennej zdarzeniowej (puste znaczniki: N_sig < 50)", fontsize=11, color=INK)
fig.tight_layout(); fig.savefig(f"{H}/fig2_universality.png"); plt.close(fig)

# ---- fig3: SF po przeważeniu per schemat
sch = ["brak", "nSelBL", "lead_pt", "pT(D0)", "pT×nSelBL", "pT×lead_pt", "pT×lead_pt×nSelBL", "pT×pierścień(ptcone100−40)"]
fig, a = plt.subplots(figsize=(11, 4.4))
for j, nm in enumerate(["9:MB", "8:J12only", "8:MBTSbased", "8:other"]):
    y = [s2["nominal"]["sf"][f"ptcone40<2|{nm}|3-30|{s}"]["SF"] for s in sch]; e = [s2["nominal"]["sf"][f"ptcone40<2|{nm}|3-30|{s}"]["err"] for s in sch]
    a.errorbar(np.arange(len(sch)) + (j - 1.5) * 0.15, y, yerr=e, fmt=MK[nm], color=COL[nm], ms=6, capsize=2.5, lw=1.4, label=LAB[nm])
a.axhline(1, color="#8a897f", lw=1); a.set_xticks(range(len(sch))); a.set_xticklabels([s.replace("×", "\n× ").replace("pierścień(ptcone100−40)", "pierścień\n(ptcone100−40)") for s in sch], fontsize=8.5)
a.set_ylabel("SF = ε_dane / ε_MC(przeważone)"); a.set_ylim(0.3, 1.25); a.legend(ncol=2, fontsize=8.5, loc="lower right")
a.set_title("SF izolacji (ptcone40 < 2 GeV, pT 3–30) po przeważeniu MC do rozkładu sygnału w danych — zmienna przeważania na osi x", fontsize=10, color=INK)
fig.tight_layout(); fig.savefig(f"{H}/fig3_sf_reweight.png"); plt.close(fig)

# ---- fig4: pile-up — frakcja izolowanego tła vs μ (oba runy) i vs nPileupVertices
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
pu = s1["pileup_bkg"]
for run, c, mk, lab in (("427929", COL["9:MB"], "o", "427929 (wszystkie LB, także poza GRL)"), ("428770", COL["8:J12only"], "s", "428770 J12only"), ("428770j20", COL["8:MBTSbased"], "^", "428770 j20 ∧ ¬J12")):
    ks = [k for k in pu if k.startswith(run + "|e_actual_mu|")]
    x = [pu[k]["mu_mean"] for k in ks]; y = [pu[k]["f"] for k in ks]; e = [pu[k]["err"] for k in ks]
    ax[0].errorbar(x, y, yerr=e, fmt=mk, color=c, ms=6, capsize=2, lw=1.3, label=lab)
ax[0].set_xscale("log"); ax[0].set_xlabel("⟨μ⟩ (actualInteractionsPerCrossing)"); ax[0].set_ylabel("frakcja tła z ptcone40 < 2 GeV (pT 3–5)"); ax[0].legend(fontsize=8.5)
ax[0].set_title("Tło vs μ: płasko do μ ≈ 0,1; spadek dopiero przy μ ≳ 2", fontsize=10, color=INK); ax[0].set_ylim(0, 0.6)
for j, (run, c, mk, lab) in enumerate((("427929MB", COL["9:MB"], "o", "427929 MB"), ("428770", COL["8:J12only"], "s", "428770 J12only"), ("428770j20", COL["8:MBTSbased"], "^", "428770 j20 ∧ ¬J12"))):
    ks = [k for k in pu if k.startswith(run + "|e_npu|")]
    ax[1].errorbar(np.arange(len(ks)) + (j - 1) * 0.12, [pu[k]["f"] for k in ks], yerr=[pu[k]["err"] for k in ks], fmt=mk, color=c, ms=6, capsize=2, lw=1.3, label=lab)
ax[1].set_xticks([0, 1, 2, 3]); ax[1].set_xticklabels(["0", "1", "2", "≥3"]); ax[1].set_xlabel("nPileupVertices"); ax[1].set_ylim(0, 0.6); ax[1].legend(fontsize=8.5)
ax[1].set_title("Tło vs nPileupVertices (zdarzenia z npu ≥ 1: 1–6 % próbki)", fontsize=10, color=INK)
fig.tight_layout(); fig.savefig(f"{H}/fig4_pileup.png"); plt.close(fig)
print("rysunki zapisane")
