"""PEAK_lxy010: tabele z out/lxy_rows.csv, lxy_boot.csv, lxy_trend.csv, refl_rows.csv, refl_shape.csv → out/tables.md, out/summary.json.
Kryterium „mierzalne”: fit zbieżny, Z = N_S/σ(N_S) ≥ 3, skala k poza granicami (0,4; 2,5) i σ(k) < 0,5. Inaczej podana jest tylko istotność."""
import os, json
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = f"{HERE}/out"
T = pd.read_csv(f"{OUT}/lxy_rows.csv")
B = pd.read_csv(f"{OUT}/lxy_boot.csv") if os.path.exists(f"{OUT}/lxy_boot.csv") else None
U = {"dm": "keV", "mKpi": "MeV"}; SC = {"dm": 1e6, "mKpi": 1e3}
T["good"] = T.ok & (T.Z >= 3) & (T.k > 0.41) & (T.k < 2.49) & (T.e_k < 0.5)
LAB = {"428770_wp70_v12": "428770, WP70, v1.2 pełna lumi (kontrola PEAK_pos; tylko ≥ 0,25)", "427929_wp70": "427929, WP70, v2.0 226/458 plików",
       "427929_nowp": "427929, bez klasyfikatora, v2.0 226/458 plików", "428770_nowp": "428770, bez klasyfikatora, v2.0 97/97 plików (pT = pT(D*))"}
OBS = {"dm": "Δm(D*)", "mKpi": "m(Kπ) z tagu D*"}


def cell(r, obs):
    if r is None:
        return "—", "—", "—"
    z = f"{r.Ns:.0f} ± {r.e_Ns:.0f} (Z {r.Z:.1f})"
    if not r.good:
        return z, "n.m.", "n.m."
    f = "{:+.0f} ± {:.0f}" if obs == "dm" else "{:+.2f} ± {:.2f}"
    return z, f.format(r.dmu, r.e_dmu), f"{r.k:.2f} ± {r.e_k:.2f}"


out, summ = [], {}
for obs in ("dm", "mKpi"):
    out.append(f"### {OBS[obs]} — Δμ = μ_dane − μ_MC(fd, L_xy ≥ 0,25) [{U[obs]}], k = σ_dane/σ_MC; n.m. = niemierzalne (podana tylko istotność)\n")
    out.append("| próbka | bin | N_S ≥0,25 | Δμ ≥0,25 | k ≥0,25 | N_S ≥0,10 | Δμ ≥0,10 | k ≥0,10 | N_S plaster 0,10–0,25 | Δμ plaster | k plaster | δμ(≥0,10 − ≥0,25), bootstrap | k(≥0,10)/k(≥0,25), bootstrap |")
    out.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for s in LAB:
        for b in ("pT3-5", "pT5-8", "pT8-30", "pT3-30"):
            g = T[(T["sample"] == s) & (T.bin == b) & (T.obs == obs)].set_index("cut")
            c = [cell(g.loc[k] if k in g.index else None, obs) for k in ("ge025", "ge010", "slice")]
            bd, bk = "—", "—"
            if B is not None and "ge010" in g.index:
                x = B[(B["sample"] == s) & (B.bin == b) & (B.obs == obs) & B.ok]
                if len(x) >= 50 and g.loc["ge010"].good and g.loc["ge025"].good:
                    d0 = (g.loc["ge010"].mu - g.loc["ge025"].mu) * SC[obs]; r0 = g.loc["ge010"].k / g.loc["ge025"].k
                    # błąd: połowa przedziału 16–84 % (odporny na pojedyncze nieudane fity)
                    ed = 0.5 * (np.percentile(x.dmu, 84) - np.percentile(x.dmu, 16)) * SC[obs]; er = 0.5 * (np.percentile(x.rk, 84) - np.percentile(x.rk, 16))
                    bd = ("{:+.0f} ± {:.0f}" if obs == "dm" else "{:+.2f} ± {:.2f}").format(d0, ed) + f" ({len(x)})"; bk = f"{r0:.3f} ± {er:.3f}"
                    summ[f"{s}|{b}|{obs}"] = dict(delta=d0, e_delta=ed, kratio=r0, e_kratio=er, nboot=int(len(x)))
            out.append(f"| {LAB[s]} | {b} | {c[0][0]} | {c[0][1]} | {c[0][2]} | {c[1][0]} | {c[1][1]} | {c[1][2]} | {c[2][0]} | {c[2][1]} | {c[2][2]} | {bd} | {bk} |")
    out.append("")
# trend powyżej 0,25
R = pd.read_csv(f"{OUT}/lxy_trend.csv")
out.append("### Trend w plasterkach L_xy powyżej 0,25 mm — 428770 pełna lumi, WP70; MC fd dopasowane w tym samym plasterku\n")
out.append("| obs | bin | L_xy [mm] | μ_MC | σ_MC | N_S dane (Z) | Δμ dane − MC | k |"); out.append("|---|---|---|---|---|---|---|---|")
for _, r in R.iterrows():
    good = r.ok and r.Z >= 3 and 0.41 < r.k < 2.49 and r.e_k < 0.5
    f = "{:+.0f} ± {:.0f}" if r.obs == "dm" else "{:+.2f} ± {:.2f}"
    out.append(f"| {OBS[r.obs]} | {r.bin} | {r.lxy} | {r.mu_mc:.1f} ± {r.e_mu_mc:.1f} | {r.seff_mc:.1f} | {r.Ns:.0f} ({r.Z:.1f}) | "
               + (f.format(r.dmu, r.e_dmu) if good else "n.m.") + " | " + (f"{r.k:.2f} ± {r.e_k:.2f}" if good else "n.m.") + " |")
out.append("")
# odbicia
F = pd.read_csv(f"{OUT}/refl_rows.csv"); S = pd.read_csv(f"{OUT}/refl_shape.csv")
out.append("### Odbicie K↔π z tagu D*: R = N_odbicie(1,75–2,00) / N_zgodne(1,84–1,89); dane/MC (MC zawsze L_xy ≥ 0,25); χ² kształtu vs MC fd / mb (ndf 24)\n")
out.append("| score > | bin | próbka | L_xy | N_zgodne (okno) | N_odbicie | R | R dane/MC fd | R dane/MC mb | χ² kształtu vs fd | vs mb |"); out.append("|---|---|---|---|---|---|---|---|---|---|---|")
for _, r in F.iterrows():
    dat = r["sample"].startswith("data")
    meas = r.n_tag_win / r.e_n_tag_win >= 3 and r.n_refl / r.e_n_refl >= 3
    c = {}
    for ref in ("fd_nd@ge025", "mb_nd@ge025"):
        x = S[(S.score == r.score) & (S.bin == r.bin) & (S["sample"] == r["sample"]) & (S.cut == r.cut) & (S.vs == ref)]
        c[ref] = f"{x.chi2.iloc[0]:.0f}" if len(x) and x.n_refl.iloc[0] > 200 and x.n_refl_ref.iloc[0] > 200 else "—"
    cut = {"ge025": "≥ 0,25", "ge010": "≥ 0,10", "slice": "0,10–0,25"}[r.cut]
    out.append(f"| {r.score} | {r.bin} | {r['sample']} | {cut} | {r.n_tag_win:.0f} ± {r.e_n_tag_win:.0f} | {r.n_refl:.0f} ± {r.e_n_refl:.0f} | "
               + (f"{r.R:.2f} ± {r.e_R:.2f}" if meas else "n.m.") + " | "
               + (f"{r.R_over_fd_nd:.2f} ± {r.e_R_over_fd_nd:.2f}" if dat and meas else "—") + " | "
               + (f"{r.R_over_mb_nd:.2f} ± {r.e_R_over_mb_nd:.2f}" if dat and meas else "—") + f" | {c['fd_nd@ge025'] if dat else '—'} | {c['mb_nd@ge025'] if dat else '—'} |")
open(f"{OUT}/tables.md", "w").write("\n".join(out) + "\n")
json.dump(summ, open(f"{OUT}/summary.json", "w"), indent=1)
print("\n".join(out[:60]))
