"""FIT_scorebins — tabele (out/tables.md) i rysunki (fig1_sigma_vs_wp.png, fig2_robustness.png)."""
import os, sys, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
SPS = [("D0", r"$D^0\to K\pi$"), ("Dplus", r"$D^+\to K\pi\pi$")]
PTB = ["3–5", "5–8", "8–30"]
COL = {"cut": "#52514e", "fix3": "#2a78d6", "fix5": "#eb6834", "fix17": "#1baf7a", "c10": "#6250d6", "free": "#e34948"}
rng = np.random.default_rng(7)


def stats(t, name, ref=None, nb=400):
    names = list(t["names"]); i = names.index(name)
    N, E, NS = t["N"][:, i], t["E"][:, i], float(t["NS"])
    pull = (N - NS) / E
    o = dict(rel=N.std(ddof=1) / NS, rel_hesse=np.mean(E) / NS, bias=(N.mean() - NS) / N.std(ddof=1), pull_m=pull.mean(), pull_s=pull.std(ddof=1),
             cov=np.mean(np.abs(pull) < 1), ok=t["OK"][:, i].mean())
    o["bias_err"] = 1 / np.sqrt(len(N))
    if ref is not None:
        R = t["N"][:, names.index(ref)]
        idx = rng.integers(0, len(N), (nb, len(N)))
        g = 1 - N[idx].std(1, ddof=1) / R[idx].std(1, ddof=1)
        o["gain"], o["gain_err"] = 1 - N.std(ddof=1) / R.std(ddof=1), g.std(ddof=1)
    return o


def main():
    md = []
    fig1, ax1 = plt.subplots(2, 3, figsize=(13, 7.2), sharex=True)
    fig2, ax2 = plt.subplots(2, 3, figsize=(13, 7.2), sharey=True)
    summary = {}
    for r, (sp, lab) in enumerate(SPS):
        for ib in range(3):
            res = json.load(open(f"{OUT}/results_{sp}_{ib}.json"))
            t = np.load(f"{OUT}/toys_{sp}_{ib}.npz", allow_pickle=True)
            names = list(t["names"]); NS = res["NS"]; bw = res["best_wp"]
            ntoy = len(t["N"])
            md.append(f"\n### {sp}, pT {PTB[ib]} GeV — N_S(oczek.) = {NS:.0f}, {ntoy} toyów, najlepszy WP (Asimov) = {bw[3:]}\n")
            md.append("| konfiguracja | granice ε_S (dolne) | σ/N_S toy [%] | σ/N_S HESSE [%] | zysk vs najlepsze cięcie [%] | obciążenie [σ] | pull śr. / szer. | pokrycie 1σ | zbieżne |")
            md.append("|---|---|---|---|---|---|---|---|---|")
            show = [bw, "cut70", "cut80", "cut50", "cut100", "fix2", "fix3", "fix5", "fix8", "fix17", "fix3_t95", "free3", "free5",
                    "c05_3", "c05_5", "c10_3", "c10_5"]
            seen = set()
            for n in show:
                if n in seen:
                    continue
                seen.add(n)
                s = stats(t, n, ref=bw if n != bw else None)
                g = f"{100*s['gain']:+.1f} ± {100*s['gain_err']:.1f}" if "gain" in s else "—"
                b = ", ".join(f"{x:.2f}" for x in res["bounds"][n]) if len(res["bounds"][n]) <= 8 else "17 komórek"
                md.append(f"| {n} | {b} | {100*s['rel']:.2f} | {100*s['rel_hesse']:.2f} | {g} | {s['bias']:+.2f} ± {s['bias_err']:.2f} | "
                          f"{s['pull_m']:+.2f} / {s['pull_s']:.2f} | {100*s['cov']:.1f} % | {100*s['ok']:.1f} % |")
                summary[(sp, ib, n)] = s
            # błąd całkowity z niepewnością SF
            sb = stats(t, bw)
            md.append(f"\nUłamki ograniczone przy delta_prawd = 0 (UWAGA: odch. std. toyów = tylko część statystyczna; błąd CAŁKOWITY z losowaną wydajnością: out/toys_sf.json i RAPORT §3), {sp} {PTB[ib]}: ")
            md.append("| σ_SF | cięcie {0} ⊕ σ_SF [%] | K=3 ogr. [%] | K=5 ogr. [%] | zysk K=5 ogr. vs cięcie [%] | K=5 swobodne [%] |".format(bw[3:]))
            md.append("|---|---|---|---|---|---|")
            for sg, tag in ((0.05, "c05"), (0.10, "c10")):
                tot = np.hypot(sb["rel"], sg)
                s3, s5, sf = stats(t, f"{tag}_3"), stats(t, f"{tag}_5"), stats(t, "free5")
                md.append(f"| {100*sg:.0f} % | {100*tot:.2f} | {100*s3['rel']:.2f} | {100*s5['rel']:.2f} | {100*(1-s5['rel']/tot):+.1f} | {100*sf['rel']:.2f} |")
                summary[(sp, ib, "tot" + tag)] = dict(cut=tot, c3=s3["rel"], c5=s5["rel"])
            # odporność (Asimov): Δ obciążenia względem nominału, w jednostkach σ_stat (toy) danej metody
            md.append(f"\nOdporność (Asimov, przesunięcie N̂_S/N_S względem nominału; w nawiasie w jednostkach σ_stat metody), {sp} {PTB[ib]}:")
            meth = [bw, "cut70", "fix3", "fix5", "fix8", "c05_5", "c10_5", "c10_3", "free5"]
            meth = list(dict.fromkeys(meth))
            md.append("| scenariusz | " + " | ".join(meth) + " |")
            md.append("|---|" + "---|" * len(meth))
            rob = {}
            for sc, d in res["robust"].items():
                row = []
                for n in meth:
                    db = d[n][0] - res["asimov_fit"][n][0]
                    sg = stats(t, n)["rel"]
                    rob[(sc, n)] = db / sg
                    row.append(f"{100*db:+.1f} % ({db/sg:+.2f}σ)")
                md.append(f"| {sc} | " + " | ".join(row) + " |")
            md.append(f"\nSkan normalizacji sygnału (Asimov): " + "; ".join(
                f"S×{k}: cięcie WP{int(round(100*v['best_wp']))} {100*v['cut']:.2f} %, K=3 {100*(1-v['fix3']/v['cut']):+.1f} %, K=5 {100*(1-v['fix5']/v['cut']):+.1f} %, K=8 {100*(1-v['fix8']/v['cut']):+.1f} %"
                for k, v in res["scan_NS"].items()))
            # --- rys. 1 ---
            a = ax1[r, ib]
            wps = sorted(int(n[3:]) for n in names if n.startswith("cut") and n != "cut100")
            y = [100 * stats(t, f"cut{w}")["rel"] for w in wps]
            a.plot(wps, y, "-o", color=COL["cut"], lw=2, ms=5, label="pojedyncze cięcie (skan WP)")
            for n, c, l in (("fix3", COL["fix3"], "K = 3 kategorie"), ("fix5", COL["fix5"], "K = 5 kategorii"),
                            ("fix17", COL["fix17"], "K = 17 (≈ fit 2D masa × score)")):
                v = 100 * stats(t, n)["rel"]
                a.axhline(v, color=c, lw=2, ls={"fix3": "-", "fix5": "--", "fix17": ":"}[n], label=l)
            ymin = min(y)
            a.set_ylim(0.8 * 100 * stats(t, "fix17")["rel"], 1.9 * ymin)
            a.set_title(f"{lab},  $p_T$ {PTB[ib]} GeV", fontsize=11, loc="left")
            a.grid(alpha=0.25, lw=0.6)
            for s_ in ("top", "right"):
                a.spines[s_].set_visible(False)
            if ib == 0:
                a.set_ylabel(r"$\sigma(N_S)/N_S$ [%]")
            if r == 1:
                a.set_xlabel(r"punkt pracy WP = $\varepsilon_S$ cięcia [%]")
            g5 = stats(t, "fix5", ref=bw)
            a.text(0.03, 0.95, f"najlepszy WP{bw[3:]}: {ymin:.2f} %\nK=5: {100*g5['rel']:.2f} %  (zysk {100*g5['gain']:.1f} ± {100*g5['gain_err']:.1f} %)",
                   transform=a.transAxes, ha="left", va="top", fontsize=9, color="#0b0b0b")
            # --- rys. 2 ---
            a = ax2[r, ib]
            scs = [("SF_A_428770_5-8", "SF jak 428770 5–8 (0,96/0,87/0,84)", "#2a78d6"),
                   ("SF_B_428770_8-30", "SF jak 428770 8–30 (0,94/0,97/1,01)", "#eb6834"),
                   ("sculpt_R089", "rzeźbienie tła R = 0,89", "#1baf7a")]
            m2 = [(bw, f"cięcie\nWP{bw[3:]}"), ("fix3", "K=3\nMC"), ("fix5", "K=5\nMC"), ("c10_5", "K=5\nogr. 10 %"), ("free5", "K=5\nswob.")]
            for k, (sc, l, c) in enumerate(scs):
                v = [rob[(sc, n)] for n, _ in m2]
                a.bar(np.arange(len(m2)) + (k - 1) * 0.26, v, 0.24, color=c, label=l)
            a.axhline(0, color="#52514e", lw=0.8)
            a.axhspan(-0.2, 0.2, color="#52514e", alpha=0.08, lw=0)
            a.set_xticks(range(len(m2))); a.set_xticklabels([l for _, l in m2], fontsize=9)
            a.set_title(f"{lab},  $p_T$ {PTB[ib]} GeV", fontsize=11, loc="left")
            a.grid(axis="y", alpha=0.25, lw=0.6)
            for s_ in ("top", "right"):
                a.spines[s_].set_visible(False)
            if ib == 0:
                a.set_ylabel(r"obciążenie $N_S$ [$\sigma_{stat}$ metody]")
    for fig, ax, name, sub in ((fig1, ax1, "fig1_sigma_vs_wp.png", "pseudo-eksperymenty: sygnał MC, tło z pasm bocznych danych (run 428770)"),
                               (fig2, ax2, "fig2_robustness.png", "Asimov: obciążenie przy ε_S ≠ MC i rzeźbieniu tła; pasmo szare = ±0,2σ")):
        h, l = ax[0, 0].get_legend_handles_labels()
        fig.legend(h, l, loc="lower center", ncol=4, frameon=False, fontsize=10)
        fig.text(0.012, 0.965, "ATLAS", fontsize=14, fontweight="bold", fontstyle="italic")
        fig.text(0.066, 0.965, "", fontsize=14)
        fig.text(0.14, 0.966, sub, fontsize=10, color="#52514e")
        fig.tight_layout(rect=(0, 0.05, 1, 0.95))
        fig.savefig(os.path.join(HERE, name), dpi=140)
    open(f"{OUT}/tables.md", "w").write("\n".join(md) + "\n")
    # zwięzłe podsumowanie
    for (sp, _l) in SPS:
        for ib in range(3):
            res = json.load(open(f"{OUT}/results_{sp}_{ib}.json"))
            line = [f"{sp} {PTB[ib]} best {res['best_wp']} {100*summary[(sp, ib, res['best_wp'])]['rel']:.2f}%"]
            for n in ("fix2", "fix3", "fix5", "fix8", "fix17"):
                s = summary[(sp, ib, n)]
                line.append(f"{n} {100*s['gain']:+.1f}±{100*s['gain_err']:.1f} (b {s['bias']:+.2f}, cov {100*s['cov']:.0f})")
            for tag in ("totc05", "totc10"):
                s = summary[(sp, ib, tag)]
                line.append(f"{tag}: cut {100*s['cut']:.2f} c3 {100*s['c3']:.2f} c5 {100*s['c5']:.2f}")
            print(" | ".join(line))


if __name__ == "__main__":
    main()
