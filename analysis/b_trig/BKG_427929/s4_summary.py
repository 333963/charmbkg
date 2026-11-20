"""Krok 4: tabele zbiorcze (out/s4.md, out/s4.json) i rysunki fig1–fig3 (PNG)."""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import bkg_common as C

s1 = json.load(open(os.path.join(C.OUT, "s1.json")))
s2 = json.load(open(os.path.join(C.OUT, "s2.json")))
s3 = json.load(open(os.path.join(C.OUT, "s3.json")))
rep = np.load(os.path.join(C.OUT, "s1_rep.npz"))
T = s1["tab"]
COL = {"428": "#2a78d6", "427": "#eb6834", "427w": "#1baf7a"}
INK, MUTED = "#1a1a19", "#6b6a63"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": "white", "axes.grid": True, "grid.color": "#e6e5e0", "grid.linewidth": 0.6})
DEN = {"A_t8": "t8", "C_t2": "t2"}
RW = "V2 + score"
out = {}
L = ["# s4 — tabele zbiorcze", ""]
for den, tn in DEN.items():
    L += [f"## Mianownik MC: {'truth ≥ 8 w SR (btrig)' if tn == 't8' else 'truth ≥ 2 w SR (całe tło niepikujące)'}", "",
          "| bin | WP | R 428770 | R 427929 | R 427929 po przeważeniu (V1 bez score) | (V2 + score) | R 428770 po przeważeniu (V2) | udział triggera f |", "|---|---|---|---|---|---|---|---|"]
    for b in C.BINS:
        for j, wp in enumerate(C.WPS):
            r8, r7 = T[f"d428|{den}|{b}|{wp}"], T[f"d427|{den}|{b}|{wp}"]
            e8 = np.sqrt(r8["stat_data"] ** 2 + r8["stat_mc"] ** 2 + r8["syst"] ** 2)
            e7 = np.sqrt(r7["stat_data"] ** 2 + r7["stat_mc"] ** 2 + r7["syst"] ** 2)
            w = {v: s3[f"d427|{v}"]["bins"][b][tn] for v in ("V1 bez score", "V2 + score")}
            w8 = s3[f"d428|{RW}"]["bins"][b][tn]
            sy7 = r7["syst"] / r7["R"]
            ew = {v: np.hypot(w[v]["err_after"][j], sy7 * w[v]["R_after"][j]) for v in w}
            ew8 = np.hypot(w8["err_after"][j], r8["syst"] / r8["R"] * w8["R_after"][j])
            # udział triggera: f = (R427 − R428)/(1 − R428); błąd z replik danych (MC wspólne — skraca się w liczniku, zostaje w mianowniku)
            Dm = r7["R"] - r8["R"]
            sD = np.hypot(r7["stat_data"], r8["stat_data"])
            f = Dm / (1 - r8["R"])
            sf = abs(f) * np.hypot(sD / abs(Dm), np.hypot(r8["stat_data"], r8["stat_mc"]) / abs(1 - r8["R"]))
            out[f"{den}|{b}|{wp}"] = {"R428": r8["R"], "e428": e8, "R427": r7["R"], "e427": e7, "R427_V1": w["V1 bez score"]["R_after"][j], "e427_V1": ew["V1 bez score"],
                                      "R427_V2": w["V2 + score"]["R_after"][j], "e427_V2": ew["V2 + score"], "R428_V2": w8["R_after"][j], "e428_V2": ew8, "f": f, "ef": sf}
            o = out[f"{den}|{b}|{wp}"]
            L.append(f"| {b} | {wp} | {o['R428']:.3f} ± {o['e428']:.3f} | {o['R427']:.3f} ± {o['e427']:.3f} | {o['R427_V1']:.3f} ± {o['e427_V1']:.3f} | "
                     f"{o['R427_V2']:.3f} ± {o['e427_V2']:.3f} | {o['R428_V2']:.3f} ± {o['e428_V2']:.3f} | {f:.2f} ± {sf:.2f} |")
    L.append("")
# średnie ważone (1/σ²) po 3 binach pT per WP — orientacyjne
L += ["## Średnie po binach pT (wagi 1/σ², σ = stat⊕syst; orientacyjnie)", "", "| mianownik | WP | R 428770 | R 427929 | R 427929 po V2 | f = (R427−R428)/(1−R428) |", "|---|---|---|---|---|---|"]
for den in DEN:
    for wp in C.WPS:
        def avg(k, e):
            x = np.array([out[f"{den}|{b}|{wp}"][k] for b in C.BINS]); s = np.array([out[f"{den}|{b}|{wp}"][e] for b in C.BINS])
            return float(np.sum(x / s ** 2) / np.sum(1 / s ** 2)), float(1 / np.sqrt(np.sum(1 / s ** 2)))
        a8, a7, aw = avg("R428", "e428"), avg("R427", "e427"), avg("R427_V2", "e427_V2")
        out[f"avg|{den}|{wp}"] = {"R428": a8, "R427": a7, "R427_V2": aw, "f": (a7[0] - a8[0]) / (1 - a8[0])}
        L.append(f"| {den} | {wp} | {a8[0]:.3f} ± {a8[1]:.3f} | {a7[0]:.3f} ± {a7[1]:.3f} | {aw[0]:.3f} ± {aw[1]:.3f} | {(a7[0] - a8[0]) / (1 - a8[0]):.2f} |")
# D dane/dane estymatorem sym (bez ekstrapolacji — mały błąd stat.)
L += ["", "## D = ε_B(427929)/ε_B(428770), dane/dane: pol2 (pod pikiem) i sym (w SB, zliczenie)", "", "| bin | WP | D pol2 | z | D sym | z |", "|---|---|---|---|---|---|"]
for b in C.BINS:
    for j, wp in enumerate(C.WPS):
        d = T[f"diff|d428|{b}|{wp}"]
        e1, e2 = s1["eps"][f"d427|{b}"]["sym"][j], s1["eps"][f"d428|{b}"]["sym"][j]
        r = rep[f"d427|{b}|sym"][:, j] / rep[f"d428|{b}|sym"][:, j]
        out[f"Dsym|{b}|{wp}"] = {"D": e1 / e2, "err": float(r.std())}
        L.append(f"| {b} | {wp} | {d['D']:.3f} ± {d['stat']:.3f} ± {d['syst']:.3f} | {d['z_tot']:.1f} | {e1 / e2:.3f} ± {r.std():.3f} | {(e1 / e2 - 1) / r.std():.1f} |")
open(os.path.join(C.OUT, "s4.md"), "w").write("\n".join(L) + "\n")
json.dump(out, open(os.path.join(C.OUT, "s4.json"), "w"), indent=1, default=float)
print("\n".join(L))

# ------------------------------------------------------------------ fig1: R vs WP
fig, axs = plt.subplots(2, 3, figsize=(12, 7.2), sharey="row", sharex=True)
x = np.arange(5)
for i, (den, tn) in enumerate(DEN.items()):
    for k, b in enumerate(C.BINS):
        ax = axs[i, k]
        ax.axhline(1, color=MUTED, lw=1)
        for key, ek, off, c, mk, lab in (("R428", "e428", -0.2, COL["428"], "o", "428770 (tryger dżetowy) / MC"), ("R427", "e427", 0.0, COL["427"], "s", "427929 (MBTS) / MC"),
                                         ("R427_V2", "e427_V2", 0.2, COL["427w"], "^", "427929 / MC przeważone (V2)")):
            y = [out[f"{den}|{b}|{wp}"][key] for wp in C.WPS]; e = [out[f"{den}|{b}|{wp}"][ek] for wp in C.WPS]
            ax.errorbar(x + off, y, e, fmt=mk, color=c, ms=6, lw=1.3, capsize=2, label=lab, mec="white", mew=0.6)
        ax.set_ylim(0.35, 1.75)
        ax.set_xticks(x); ax.set_xticklabels([f"WP{int(float(w) * 100)}" for w in C.WPS])
        C.atlas_label(ax, f"D$^0$→Kπ, tło pod pikiem, {b.replace('-', '–')} GeV\nMC tła: truth ≥ {8 if tn == 't8' else 2}")
        if k == 0:
            ax.set_ylabel("ε$_B$(dane, interpolacja z SB) / ε$_B$(MC, SR)")
        if i == 0 and k == 2:
            ax.legend(loc="lower left", fontsize=8, frameon=False)
fig.suptitle("Niezgodność tła dane/MC po cięciu na score H: run dżetowy vs minbias vs minbias po przeważeniu (błędy: stat ⊕ syst interpolacji)", fontsize=10.5, color=INK)
fig.tight_layout()
fig.savefig(os.path.join(C.HERE, "fig1_epsb_ratio.png"), dpi=150)

# ------------------------------------------------------------------ fig2: stosunek rozkładu score w SB
fig, axs = plt.subplots(1, 3, figsize=(12, 4.2), sharey=True)
cen = np.linspace(0.025, 0.975, 20)
for k, b in enumerate(C.BINS):
    ax = axs[k]
    ax.axhline(1, color=MUTED, lw=1)
    for pn, c, mk, lab in (("428770 vs MC", COL["428"], "o", "428770 / MC"), ("427929 vs MC", COL["427"], "s", "427929 / MC"), ("427929 vs 428770", COL["427w"], "^", "427929 / 428770 (dane/dane)")):
        r = s2["ratio_score"][f"{pn}|{b}"]
        ax.errorbar(cen, r["ratio"], r["err"], fmt=mk + "-", color=c, ms=4.5, lw=1.2, label=lab, mec="white", mew=0.5)
    ax.set_ylim(0.45, 1.65); ax.set_xlabel("score H")
    C.atlas_label(ax, f"pasma boczne |m−m$_0$| ≥ 75 MeV, {b.replace('-', '–')} GeV\nrozkłady unormowane, MC inkluzywne mb_nd")
    if k == 0:
        ax.set_ylabel("stosunek gęstości")
        ax.legend(loc="lower left", fontsize=8, frameon=False)
fig.tight_layout()
fig.savefig(os.path.join(C.HERE, "fig2_score_ratio_sb.png"), dpi=150)

# ------------------------------------------------------------------ fig3: C2ST
sets = ["score", "kand (score+ptcone)", "zdarzenie (nSelBL, lead_pt, ntrk500)", "wszystkie bez score", "wszystkie + score"]
psets = ["18 cech H", "18 cech + ptcone + zdarzenie"]
labs = ["score H", "score + ptcone40/100", "nSelBL, lead pT, ntrk500", "6 zmiennych bez score", "6 zmiennych + score", "pary D*: 18 cech H\n(pasmo 55–75 MeV)", "pary D*: 18 cech\n+ ptcone + zdarzenie"]
fig, ax = plt.subplots(figsize=(8.6, 4.8))
yy = np.arange(len(labs))[::-1]
ax.axvline(0.5, color=MUTED, lw=1)
for pn, c, mk, off, lab in (("428770 vs MC", COL["428"], "o", 0.2, "428770 vs MC"), ("427929 vs MC", COL["427"], "s", 0.0, "427929 vs MC"), ("427929 vs 428770", COL["427w"], "^", -0.2, "427929 vs 428770 (dane–dane)")):
    v = [s2["c2st"][f"{pn}|3-30|{s}"] for s in sets] + [s2["c2st"][f"pary|{pn}|3-30|{s}"] for s in psets]
    ax.errorbar([q["auc"] for q in v], yy + off, xerr=[q["err_folds"] for q in v], fmt=mk, color=c, ms=7, label=lab, mec="white", mew=0.6)
    for q, y in zip(v, yy + off):
        ax.text(q["auc"] + 0.008, y, f"{q['auc']:.3f}", va="center", fontsize=7.5, color=INK)
ax.set_yticks(yy); ax.set_yticklabels(labs, fontsize=8.5)
ax.set_xlim(0.48, 0.99); ax.set_xlabel("C2ST AUC (HGB, 5-fold po zdarzeniach), pasma boczne, p$_T$ 3–30 GeV")
ax.grid(axis="y", visible=False)
ax.legend(loc="upper right", fontsize=8.5, frameon=False)
C.atlas_label(ax, "", x=0.30, y=0.97)
fig.tight_layout()
fig.savefig(os.path.join(C.HERE, "fig3_c2st.png"), dpi=150)
