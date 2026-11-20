"""PEAK_lxy010 v2 — tabele z out/{mc_fd,mc_fd_boot,rows_*,boot_*,refl_rows,refl_boot,refl_shape}.csv → out/tables_v2.md, out/summary_v2.json.
„Mierzalne” (jak w wersji wstępnej): fit zbieżny, Z = N_S/σ(N_S) ≥ 3, k poza granicami (0,4; 2,5), σ(k) < 0,5; dodatkowo kształt MC z ≥ 200 zdarzeń i zbieżny.
Błąd bootstrapu = połowa przedziału 16–84 % (repliki po zdarzeniach)."""
import os, json
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = f"{HERE}/out"
U = {"dm": "keV", "mKpi": "MeV"}; SC = {"dm": 1e6, "mKpi": 1e3}
OBS = {"dm": "Δm(D*)", "mKpi": "m(Kπ) z tagu D*"}
BINS = ("pT3-5", "pT5-8", "pT8-30", "pT3-30"); CUT = {"ge025": "≥ 0,25", "ge010": "≥ 0,10", "slice": "0,10–0,25"}
DATA = ("data_428770", "data_427929"); SELN = {"wp70": "WP70", "nowp": "bez klasyfikatora", "s857": "score > 0,857", "s387": "score > 0,387"}
hw = lambda x: 0.5 * (np.percentile(x, 84) - np.percentile(x, 16))  # noqa: E731
fm = lambda obs: "{:+.0f} ± {:.0f}" if obs == "dm" else "{:+.2f} ± {:.2f}"  # noqa: E731

MC = pd.read_csv(f"{OUT}/mc_fd.csv"); MCB = pd.read_csv(f"{OUT}/mc_fd_boot.csv")
MC["mc_good"] = MC.ok & (MC.n >= 200)
MBF = pd.read_csv(f"{OUT}/mb_truth_free.csv")
T = pd.concat([pd.read_csv(f"{OUT}/rows_{s}.csv") for s in DATA + ("mb_nd", "mb_nd_truth")], ignore_index=True)
B = pd.concat([pd.read_csv(f"{OUT}/boot_{s}.csv") for s in DATA + ("mb_nd",)], ignore_index=True)
T = T.merge(MC[["sel", "bin", "cut", "obs", "mc_good", "n"]].rename(columns={"n": "n_mc"}), on=["sel", "bin", "cut", "obs"], how="left")
T["good"] = T.ok & (T.Z >= 3) & (T.k > 0.41) & (T.k < 2.49) & (T.e_k < 0.5) & (T.mc_good | (T["shape"] == "fd_ge025"))
out, summ = [], {}


def get(s, sel, b, c, obs, shape="fd_same"):
    g = T[(T["sample"] == s) & (T.sel == sel) & (T.bin == b) & (T.cut == c) & (T.obs == obs) & (T["shape"] == shape)]
    return g.iloc[0] if len(g) else None


def cell(r, obs):
    if r is None:
        return "—", "—", "—"
    z = f"{r.Ns:.0f} ± {r.e_Ns:.0f} (Z {r.Z:.1f})"
    if not r.good:
        return z, "n.m.", "n.m."
    return z, fm(obs).format(r.dmu, r.e_dmu), f"{r.k:.2f} ± {r.e_k:.2f}"


def mc_delta(sel, b, obs):
    """δμ_MC(≥0,10 − ≥0,25) i σ_MC(≥0,10)/σ_MC(≥0,25) z błędem bootstrapu"""
    g = MC[(MC.sel == sel) & (MC.bin == b) & (MC.obs == obs)].set_index("cut")
    x = MCB[(MCB.sel == sel) & (MCB.bin == b) & (MCB.obs == obs) & MCB.ok]
    return (g.loc["ge010"].mu - g.loc["ge025"].mu) * SC[obs], hw(x.dmu) * SC[obs], g.loc["ge010"].seff / g.loc["ge025"].seff, hw(x.rs), len(x)


def data_delta(s, sel, b, obs):
    a, c = get(s, sel, b, "ge010", obs), get(s, sel, b, "ge025", obs)
    x = B[(B["sample"] == s) & (B.sel == sel) & (B.bin == b) & (B.obs == obs) & B.ok]
    if a is None or c is None or not (a.good and c.good) or len(x) < 50:
        return None
    return (a.mu - c.mu) * SC[obs], hw(x.mu10 - x.mu25) * SC[obs], a.k / c.k, hw(x.k10 / x.k25), len(x), hw(x.mu10) * SC[obs], hw(x.mu25) * SC[obs]


# ---------- 1. MC: kształt per cięcie
out.append("## 1. MC fd_nd (sygnał z prawdy): pozycja i σ_eff per cięcie L_xy; zmiana po dodaniu plastra (bootstrap po zdarzeniach)\n")
out.append("| obs | selekcja | bin | n(≥0,25) / n(≥0,10) / n(plaster) | μ ≥0,25 | μ plaster − μ ≥0,25 | **δμ_MC(≥0,10 − ≥0,25)** | σ_eff ≥0,25 | σ plaster/σ ≥0,25 | σ(≥0,10)/σ(≥0,25) | mb_nd prawda − fd: Δμ (≥0,25 / ≥0,10 / plaster) | σ_mb/σ_fd (≥0,25 / ≥0,10 / plaster; Δm: fit kształtem fd, m(Kπ): fit swobodny, n ≥ 300) |")
out.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
for obs in OBS:
    for sel in ("wp70", "nowp"):
        for b in BINS:
            g = MC[(MC.sel == sel) & (MC.bin == b) & (MC.obs == obs)].set_index("cut")
            d, ed, rs, ers, nb = mc_delta(sel, b, obs)
            sl = g.loc["slice"]; ok = bool(sl.mc_good)
            dsl = (sl.mu - g.loc["ge025"].mu) * SC[obs]; edsl = np.hypot(sl.e_mu, g.loc["ge025"].e_mu) * SC[obs]
            mb = [get("mb_nd_truth", sel, b, c, obs) for c in CUT]
            mbs = " / ".join(fm(obs).format(r.dmu, r.e_dmu) if (r is not None and r.good) else "n.m." for r in mb)
            mbk = " / ".join(f"{r.k:.2f} ± {r.e_k:.2f}" if (r is not None and r.good) else "n.m." for r in mb)
            if obs == "mKpi":   # czysty sygnał: fit swobodny (mb_truth.py), nie fit z tłem
                q = [MBF[(MBF.sel == sel) & (MBF.bin == b) & (MBF.cut == c)].iloc[0] for c in CUT]
                mbs = " / ".join(fm(obs).format(r.dmu, r.e_dmu) if r.ok else "n.m." for r in q)
                mbk = " / ".join(f"{r.rs:.2f}" if r.ok else "n.m." for r in q)
                for c, r in zip(CUT, q):
                    if r.ok:
                        summ[f"MBminusFD|{sel}|{b}|{obs}|{c}"] = dict(dmu=r.dmu, e_dmu=r.e_dmu, rs=r.rs, n=int(r.n))
            else:
                for c, r in zip(CUT, mb):
                    if r is not None and r.good:
                        summ[f"MBminusFD|{sel}|{b}|{obs}|{c}"] = dict(dmu=r.dmu, e_dmu=r.e_dmu, rs=r.k, e_rs=r.e_k)
            f1 = "{:.1f}" if obs == "dm" else "{:.2f}"
            out.append(f"| {OBS[obs]} | {SELN[sel]} | {b} | {g.loc['ge025'].n} / {g.loc['ge010'].n} / {sl.n} | " + f1.format(g.loc['ge025'].mu * SC[obs]) + f" ± " + f1.format(g.loc['ge025'].e_mu * SC[obs])
                       + " | " + (fm(obs).format(dsl, edsl) if ok else "n.m.") + " | **" + fm(obs).format(d, ed) + "** | " + f1.format(g.loc['ge025'].seff * SC[obs])
                       + " | " + (f"{sl.seff / g.loc['ge025'].seff:.3f}" if ok else "n.m.") + f" | {rs:.3f} ± {ers:.3f} | {mbs} | {mbk} |")
            summ[f"MC|{sel}|{b}|{obs}"] = dict(dmu_mc=d, e_dmu_mc=ed, rsig=rs, e_rsig=ers, dmu_slice=dsl if ok else None, e_dmu_slice=edsl if ok else None)
out.append("")
# ---------- 2. dane − MC per cięcie
for obs in OBS:
    out.append(f"## 2. {OBS[obs]} — Δμ = μ_dane − μ_MC(fd, to samo cięcie L_xy) [{U[obs]}], k = σ_dane/σ_MC(to samo cięcie); n.m. = niemierzalne\n")
    out.append("| próbka | selekcja | bin | N_S ≥0,25 | Δμ ≥0,25 | k ≥0,25 | N_S ≥0,10 | Δμ ≥0,10 | k ≥0,10 | N_S plaster | Δμ plaster | k plaster | δμ_dane(≥0,10 − ≥0,25) boot | δμ_dane − δμ_MC | k(≥0,10)/k(≥0,25) boot | wersja wstępna: Δμ ≥0,10 z MC ≥0,25 |")
    out.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for s in DATA + ("mb_nd",):
        for sel in ("wp70", "nowp"):
            for b in BINS:
                c = [cell(get(s, sel, b, k, obs), obs) for k in CUT]
                dd = data_delta(s, sel, b, obs); md = mc_delta(sel, b, obs)
                bd = bdd = bk = "—"
                if dd:
                    bd = fm(obs).format(dd[0], dd[1]) + f" ({dd[4]})"; bk = f"{dd[2]:.3f} ± {dd[3]:.3f}"
                    bdd = fm(obs).format(dd[0] - md[0], np.hypot(dd[1], md[1]))
                    summ[f"{s}|{sel}|{b}|{obs}"] = dict(delta=dd[0], e_delta=dd[1], kratio=dd[2], e_kratio=dd[3], nboot=dd[4], delta_minus_mc=dd[0] - md[0],
                                                       e_delta_minus_mc=float(np.hypot(dd[1], md[1])), e_mu10_boot=dd[5], e_mu25_boot=dd[6])
                p = get(s, sel, b, "ge010", obs, "fd_ge025")
                lab = (s if s != "mb_nd" else "MC mb_nd jako pseudo-dane (wszystkie kandydaty)")
                out.append(f"| {lab} | {SELN[sel]} | {b} | {c[0][0]} | {c[0][1]} | {c[0][2]} | {c[1][0]} | {c[1][1]} | {c[1][2]} | {c[2][0]} | {c[2][1]} | {c[2][2]} | {bd} | {bdd} | {bk} | "
                           + (fm(obs).format(p.dmu, p.e_dmu) if p is not None and p.good else "n.m.") + " |")
    out.append("")
# ---------- 3. średnie ważone runów
out.append("## 3. Średnie ważone obu runów (tylko komórki mierzalne): dane − MC fd przy tym samym cięciu\n")
out.append("| obs | selekcja | bin | Δμ ≥0,25 | Δμ ≥0,10 | Δμ plaster | k ≥0,25 | k ≥0,10 | k plaster | δμ_dane − δμ_MC (oba runy) | χ² zgodności runów (≥0,10) |")
out.append("|---|---|---|---|---|---|---|---|---|---|---|")


def wavg(v, e):
    v, e = np.asarray(v, float), np.asarray(e, float)
    if not len(v):
        return None
    w = 1 / e ** 2
    return float((v * w).sum() / w.sum()), float(1 / np.sqrt(w.sum())), float((((v - (v * w).sum() / w.sum()) / e) ** 2).sum()), len(v)


for obs in OBS:
    for sel in ("wp70", "nowp"):
        for b in BINS:
            cs, ks, chi = [], [], "—"
            for c in CUT:
                rr = [r for r in (get(s, sel, b, c, obs) for s in DATA) if r is not None and r.good]
                a = wavg([r.dmu for r in rr], [r.e_dmu for r in rr]); k = wavg([r.k for r in rr], [r.e_k for r in rr])
                cs.append((fm(obs).format(a[0], a[1]) + f" [{a[3]}]") if a else "n.m."); ks.append(f"{k[0]:.3f} ± {k[1]:.3f}" if k else "n.m.")
                if a:
                    summ[f"AVG|{sel}|{b}|{obs}|{c}"] = dict(dmu=a[0], e_dmu=a[1], k=k[0], e_k=k[1], nruns=a[3], chi2_runs=a[2])
                if c == "ge010" and a and a[3] == 2:
                    chi = f"{a[2]:.1f}"
            dl = [summ[f"{s}|{sel}|{b}|{obs}"] for s in DATA if f"{s}|{sel}|{b}|{obs}" in summ]
            a = wavg([x["delta_minus_mc"] for x in dl], [x["e_delta_minus_mc"] for x in dl])
            if a:
                summ[f"AVG|{sel}|{b}|{obs}|ddelta"] = dict(v=a[0], e=a[1], nruns=a[3])
            out.append(f"| {OBS[obs]} | {SELN[sel]} | {b} | {cs[0]} | {cs[1]} | {cs[2]} | {ks[0]} | {ks[1]} | {ks[2]} | " + ((fm(obs).format(a[0], a[1]) + f" [{a[3]}]") if a else "—") + f" | {chi} |")
out.append("")
# ---------- 4. odbicie
F = pd.read_csv(f"{OUT}/refl_rows.csv"); S = pd.read_csv(f"{OUT}/refl_shape.csv"); RB = pd.read_csv(f"{OUT}/refl_boot.csv")
F["meas"] = (F.n_tag_win / F.e_n_tag_win >= 3) & (F.n_refl / F.e_n_refl >= 3)


def rget(sel, b, s, c):
    g = F[(F.sel == sel) & (F.bin == b) & (F["sample"] == s) & (F.cut == c)]
    if not len(g):
        return None
    r = g.iloc[0].to_dict()
    x = RB[(RB.sel == sel) & (RB.bin == b) & (RB["sample"] == s)]
    v = x[f"R_{c}"].to_numpy(); v = v[np.isfinite(v)]
    r["eb"] = hw(v) if len(v) >= 50 else np.nan
    r["rel"] = r["eb"] / abs(r["R"]) if r["R"] else np.nan
    r["meas"] = bool(r["meas"]) and np.isfinite(r["eb"]) and r["rel"] < 0.5
    return r


def rratio(sel, b, s):
    """R(plaster)/R(≥0,25) w tej samej próbce (bootstrap, rozłączne)"""
    a, c = rget(sel, b, s, "slice"), rget(sel, b, s, "ge025")
    x = RB[(RB.sel == sel) & (RB.bin == b) & (RB["sample"] == s)]
    if a is None or c is None or not (a["meas"] and c["meas"]):
        return None
    v = (x.R_slice / x.R_ge025).to_numpy(); v = v[np.isfinite(v)]
    return a["R"] / c["R"], hw(v)


out.append("## 4. Odbicie K↔π z tagu D*: R = N_odbicie(1,75–2,00) / N_zgodne(1,84–1,89); MC przy TYM SAMYM cięciu L_xy; błąd R = bootstrap po zdarzeniach (w nawiasie analityczny)\n")
out.append("| selekcja | bin | L_xy | próbka | N_zgodne | N_odbicie | R | R dane/MC fd | R dane/MC mb | χ² kształtu vs fd / vs mb (ndf 24) |")
out.append("|---|---|---|---|---|---|---|---|---|---|")
for sel in SELN:
    for b in ("pT3-30", "pT3-5", "pT5-8", "pT8-30"):
        for c in CUT:
            mc = {m: rget(sel, b, m, c) for m in ("fd_nd", "mb_nd")}
            for s in DATA + ("fd_nd", "mb_nd"):
                r = rget(sel, b, s, c)
                if r is None:
                    continue
                dat = s.startswith("data"); ra = []
                for m in ("fd_nd", "mb_nd"):
                    if dat and r["meas"] and mc[m] and mc[m]["meas"]:
                        v = r["R"] / mc[m]["R"]; e = v * np.hypot(r["rel"], mc[m]["rel"]); ra.append(f"{v:.2f} ± {e:.2f}")
                        summ[f"R|{sel}|{b}|{c}|{s}|{m}"] = dict(ratio=v, e=e, R=r["R"], e_R=r["eb"], R_mc=mc[m]["R"], e_R_mc=mc[m]["eb"])
                    else:
                        ra.append("—" if not dat else "n.m.")
                ch = []
                for m in ("fd_nd", "mb_nd"):
                    x = S[(S.sel == sel) & (S.bin == b) & (S["sample"] == s) & (S.cut == c) & (S.vs == m)]
                    ch.append(f"{x.chi2.iloc[0]:.0f}" if len(x) and x.n_refl.iloc[0] > 200 and x.n_refl_ref.iloc[0] > 200 and x.n_refl.iloc[0] / x.e_n_refl.iloc[0] > 5
                              and x.n_refl_ref.iloc[0] / x.e_n_refl_ref.iloc[0] > 5 else "—")
                out.append(f"| {SELN[sel]} | {b} | {CUT[c]} | {s} | {r['n_tag_win']:.0f} ± {r['e_n_tag_win']:.0f} | {r['n_refl']:.0f} ± {r['e_n_refl']:.0f} | "
                           + (f"{r['R']:.3f} ± {r['eb']:.3f} ({r['e_R']:.3f})" if r["meas"] else "n.m.") + f" | {ra[0]} | {ra[1]} | " + (" / ".join(ch) if dat or s == "mb_nd" else "—") + " |")
out.append("")
out.append("## 5. Czy R w plasterku 0,10–0,25 jest większe niż przy ≥ 0,25 — w danych i w MC (stosunek R_plaster/R_≥0,25 w tej samej próbce, bootstrap)\n")
out.append("| selekcja | bin | data_428770 | data_427929 | MC fd_nd | MC mb_nd | podwójny stosunek dane(średnia runów)/fd | /mb |")
out.append("|---|---|---|---|---|---|---|---|")
for sel in SELN:
    for b in ("pT3-30", "pT3-5", "pT5-8", "pT8-30"):
        q = {s: rratio(sel, b, s) for s in DATA + ("fd_nd", "mb_nd")}
        dv = [q[s] for s in DATA if q[s]]
        a = wavg([x[0] for x in dv], [x[1] for x in dv]) if dv else None
        dr = []
        for m in ("fd_nd", "mb_nd"):
            if a and q[m]:
                v = a[0] / q[m][0]; e = v * np.hypot(a[1] / a[0], q[m][1] / q[m][0]); dr.append(f"{v:.2f} ± {e:.2f}")
                summ[f"RR|{sel}|{b}|{m}"] = dict(v=v, e=e)
            else:
                dr.append("n.m.")
        for s in q:
            if q[s]:
                summ[f"Rslice_over_ge025|{sel}|{b}|{s}"] = dict(v=q[s][0], e=q[s][1])
        out.append(f"| {SELN[sel]} | {b} | " + " | ".join(f"{q[s][0]:.2f} ± {q[s][1]:.2f}" if q[s] else "n.m." for s in q) + f" | {dr[0]} | {dr[1]} |")
out.append("")
out.append("## 6. Kształt odbicia: χ²/24 — MC mb vs fd, plaster vs ≥ 0,25 w tej samej próbce\n")
out.append("| selekcja | bin | próbka | cięcie | względem | χ² | N_odbicie | N_odbicie odniesienia |"); out.append("|---|---|---|---|---|---|---|---|")
for _, x in S[(S.vs.isin(["same_data@ge025", "same_mc@ge025"])) | ((S["sample"] == "mb_nd") & (S.vs == "fd_nd"))].iterrows():
    if x.n_refl > 200 and x.n_refl_ref > 200 and x.n_refl / x.e_n_refl > 5:
        out.append(f"| {SELN[x.sel]} | {x.bin} | {x['sample']} | {CUT[x.cut]} | {x.vs} | {x.chi2:.0f} | {x.n_refl:.0f} | {x.n_refl_ref:.0f} |")
open(f"{OUT}/tables_v2.md", "w").write("\n".join(out) + "\n")
json.dump(summ, open(f"{OUT}/summary_v2.json", "w"), indent=1, default=float)
print("\n".join(out))
