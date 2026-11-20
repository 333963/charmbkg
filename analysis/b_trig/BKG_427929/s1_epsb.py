"""Krok 1: ε_B pod pikiem (interpolacja z SB) dane/MC per bin pT × WP dla 427929 i 428770.
Wyjście: out/s1.json, out/s1.md. Błędy stat.: bootstrap Poissona histogramów (kategoria score × masa), 400 replik —
zachowuje korelacje między WP. Systematyka interpolacji: max odchylenie wariantów {sym, lr, pol1} od pol2."""
import json, os
import numpy as np
import bkg_common as C

NREP = 400
rng = np.random.default_rng(427929)
thr = C.thresholds()
S = {t: C.load_all(t) for t in ("d427", "d428", "mb")}
mb = S["mb"]
tm = mb["truth"]
SEL = {
    "d427": lambda D: np.ones(len(D["m"]), bool),
    "d428": lambda D: np.ones(len(D["m"]), bool),
    "d428_MBTS": lambda D: D["mbts_any"],
    "d428_J12only": lambda D: D["j12"] & ~D["mbts_any"],
    "mc_all": lambda D: np.ones(len(D["m"]), bool),                 # inkluzywne MC (jak dane) — do „like-for-like”
    "mc_t8": lambda D: D["truth"] >= 8,                             # definicja btrig/charmbkg: truth ≥ 8 (zawiera 100/101 = swap K↔π)
    "mc_t8ns": lambda D: (D["truth"] >= 8) & (D["truth"] < 100),    # bez odbić
    "mc_t2": lambda D: D["truth"] >= 2,                             # mc_truth_min: 2
}
SRC = {"d427": "d427", "d428": "d428", "d428_MBTS": "d428", "d428_J12only": "d428", "mc_all": "mb", "mc_t8": "mb", "mc_t8ns": "mb", "mc_t2": "mb"}
METH = ["pol2", "sym", "lr", "pol1", "pol3", "true"]
res = {"thr": thr, "N": {}, "eps": {}, "rep": {}}
H = {}
for k, f in SEL.items():
    D = S[SRC[k]]
    s0 = f(D)
    for b, (lo, hi) in C.BINS.items():
        s = s0 & (D["pt"] >= lo) & (D["pt"] < hi)
        h = C.cat_hist(D["m"][s], D["score"][s], thr[b])
        H[(k, b)] = h
        ha = h.sum(0)
        res["N"][f"{k}|{b}"] = {"n": int(s.sum()), "n_sb": int(ha[C.SBm].sum()), "n_sr": int(ha[C.SRm].sum()),
                                 "n_sb_pass": [int(h[j + 1:, C.SBm].sum()) for j in range(5)],
                                 "n_sr_pass": [int(h[j + 1:, C.SRm].sum()) for j in range(5)]}
        res["eps"][f"{k}|{b}"] = {m: C.all_eps(h, m).tolist() for m in METH}
# repliki
REP = {}
for (k, b), h in H.items():
    arr = {m: np.zeros((NREP, 5)) for m in METH}
    for i in range(NREP):
        hb = C.boot(h, rng)
        for m in METH:
            arr[m][i] = C.all_eps(hb, m)
    REP[(k, b)] = arr
np.savez_compressed(os.path.join(C.OUT, "s1_rep.npz"), **{f"{k}|{b}|{m}": v for (k, b), a in REP.items() for m, v in a.items()})


def chi2(d, cov):
    d = np.asarray(d)
    return float(d @ np.linalg.solve(cov, d))


from scipy import stats
rows = []
tab = {}
# definicje mianownika MC: (A) truth ≥ 8 w SR, zliczenie wprost [btrig]; (B) like-for-like: ta sama interpolacja pol2 na inkluzywnym MC
for run in ("d427", "d428", "d428_MBTS", "d428_J12only"):
    for den, (mk, mm) in {"A_t8": ("mc_t8", "true"), "B_like": ("mc_all", "pol2"), "C_t2": ("mc_t2", "true"), "D_t8ns": ("mc_t8ns", "true")}.items():
        c2s, c2t, nd = 0.0, 0.0, 0
        for b in C.BINS:
            e_d = np.array(res["eps"][f"{run}|{b}"]["pol2"])
            e_m = np.array(res["eps"][f"{mk}|{b}"][mm])
            R = e_d / e_m
            rd, rm = REP[(run, b)]["pol2"], REP[(mk, b)][mm]
            sd, sm = rd.std(0) / e_d, rm.std(0) / e_m
            # systematyka interpolacji (dane): warianty; dla like-for-like wariant liczony spójnie w danych i MC
            var = {}
            for v in ("sym", "lr", "pol1", "pol3"):
                ev = np.array(res["eps"][f"{run}|{b}"][v])
                var[v] = (ev / np.array(res["eps"][f"{mk}|{b}"][v])) if den == "B_like" else ev / e_m
            sy = np.max([np.abs(var[v] / R - 1) for v in ("sym", "lr", "pol1")], axis=0)
            Rrep = rd / rm[rng.permutation(NREP)]
            cov = np.cov(Rrep.T)
            covt = cov + np.diag((sy * R) ** 2)
            c2s += chi2(R - 1, cov); c2t += chi2(R - 1, covt); nd += 5
            for j, wp in enumerate(C.WPS):
                tab[f"{run}|{den}|{b}|{wp}"] = {"eps_data": e_d[j], "eps_mc": e_m[j], "R": R[j], "stat_data": R[j] * sd[j], "stat_mc": R[j] * sm[j],
                                               "syst": R[j] * sy[j], **{f"R_{v}": var[v][j] for v in var}}
        tab[f"{run}|{den}|chi2"] = {"chi2_stat": c2s, "chi2_tot": c2t, "ndf": nd, "p_stat": float(stats.chi2.sf(c2s, nd)), "p_tot": float(stats.chi2.sf(c2t, nd))}
# różnica runów: D = ε_dane(427929)/ε_dane(428770) (MC się skraca)
for other in ("d428", "d428_MBTS", "d428_J12only"):
    c2s, c2t, nd = 0.0, 0.0, 0
    for b in C.BINS:
        e1, e2 = np.array(res["eps"][f"d427|{b}"]["pol2"]), np.array(res["eps"][f"{other}|{b}"]["pol2"])
        Dm = e1 / e2
        rep = REP[("d427", b)]["pol2"] / REP[(other, b)]["pol2"]
        sy = np.max([np.abs(np.array(res["eps"][f"d427|{b}"][v]) / np.array(res["eps"][f"{other}|{b}"][v]) / Dm - 1) for v in ("sym", "lr", "pol1")], axis=0)
        cov = np.cov(rep.T)
        c2s += chi2(Dm - 1, cov); c2t += chi2(Dm - 1, cov + np.diag((sy * Dm) ** 2)); nd += 5
        for j, wp in enumerate(C.WPS):
            st = rep[:, j].std()
            tab[f"diff|{other}|{b}|{wp}"] = {"D": Dm[j], "stat": st, "syst": Dm[j] * sy[j], "z_stat": (Dm[j] - 1) / st,
                                            "z_tot": (Dm[j] - 1) / np.hypot(st, Dm[j] * sy[j])}
    tab[f"diff|{other}|chi2"] = {"chi2_stat": c2s, "chi2_tot": c2t, "ndf": nd, "p_stat": float(stats.chi2.sf(c2s, nd)), "p_tot": float(stats.chi2.sf(c2t, nd))}
# domknięcie estymatora na MC: interpolacja z SB (inkl. i truth≥8) vs prawda w SR
for mk in ("mc_t8", "mc_all"):
    for b in C.BINS:
        for v in ("pol2", "sym", "lr", "pol1", "pol3"):
            num = np.array(res["eps"][f"{mk}|{b}"][v]); den_ = np.array(res["eps"][f"mc_t8|{b}"]["true"])
            # SB i SR rozłączne → niezależne repliki
            err = np.hypot(REP[(mk, b)][v].std(0) / num, REP[("mc_t8", b)]["true"].std(0) / den_) * num / den_
            tab[f"mcclos|{mk}|{b}|{v}"] = {"ratio": (num / den_).tolist(), "err": err.tolist()}
res["tab"] = tab
json.dump(res, open(os.path.join(C.OUT, "s1.json"), "w"), indent=1, default=float)

f3 = lambda x: f"{x:.3f}"
L = ["# s1 — ε_B pod pikiem dane/MC", "", "## Liczności", "", "| próbka | bin | N (1,75–2,00) | N SB | N SR | N SB pass WP95/90/80/70/50 |", "|---|---|---|---|---|---|"]
for k in SEL:
    for b in C.BINS:
        n = res["N"][f"{k}|{b}"]
        L.append(f"| {k} | {b} | {n['n']} | {n['n_sb']} | {n['n_sr']} | {'/'.join(map(str, n['n_sb_pass']))} |")
for den, ttl in {"A_t8": "A: MC truth ≥ 8 w SR (zliczenie; definicja btrig)", "B_like": "B: like-for-like (pol2 z SB na inkluzywnym MC)",
                 "C_t2": "C: MC truth ≥ 2 w SR", "D_t8ns": "D: MC truth 8–13 (bez odbić 100/101)"}.items():
    L += ["", f"## R = ε_B(dane, pol2) / ε_B(MC) — mianownik {ttl}", "", "| run | bin | WP | ε dane | ε MC | R | ± stat dane | ± stat MC | ± syst interp. | R(sym) | R(lr) | R(pol1) | R(pol3) |", "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for run in ("d427", "d428", "d428_MBTS", "d428_J12only"):
        for b in C.BINS:
            for wp in C.WPS:
                r = tab[f"{run}|{den}|{b}|{wp}"]
                L.append(f"| {run} | {b} | {wp} | {r['eps_data']:.4f} | {r['eps_mc']:.4f} | {f3(r['R'])} | {f3(r['stat_data'])} | {f3(r['stat_mc'])} | {f3(r['syst'])} | "
                         f"{f3(r['R_sym'])} | {f3(r['R_lr'])} | {f3(r['R_pol1'])} | {f3(r['R_pol3'])} |")
        c = tab[f"{run}|{den}|chi2"]
        L.append(f"| {run} | χ² vs 1 | | | | stat: {c['chi2_stat']:.1f}/{c['ndf']} | p={c['p_stat']:.2g} | stat⊕syst: {c['chi2_tot']:.1f}/{c['ndf']} | p={c['p_tot']:.2g} | | | | |")
L += ["", "## Różnica runów D = ε_B(427929)/ε_B(drugi), dane/dane", "", "| drugi | bin | WP | D | ± stat | ± syst | z stat | z stat⊕syst |", "|---|---|---|---|---|---|---|---|"]
for other in ("d428", "d428_MBTS", "d428_J12only"):
    for b in C.BINS:
        for wp in C.WPS:
            r = tab[f"diff|{other}|{b}|{wp}"]
            L.append(f"| {other} | {b} | {wp} | {f3(r['D'])} | {f3(r['stat'])} | {f3(r['syst'])} | {r['z_stat']:.1f} | {r['z_tot']:.1f} |")
    c = tab[f"diff|{other}|chi2"]
    L.append(f"| {other} | χ² vs 1 | | stat: {c['chi2_stat']:.1f}/{c['ndf']} | p={c['p_stat']:.2g} | stat⊕syst: {c['chi2_tot']:.1f}/{c['ndf']} | p={c['p_tot']:.2g} | |")
L += ["", "## Domknięcie estymatorów na MC: ε(interpolacja z SB) / ε(prawda truth ≥ 8 w SR)", "", "| MC w SB | bin | metoda | WP95 | WP90 | WP80 | WP70 | WP50 |", "|---|---|---|---|---|---|---|---|"]
for mk in ("mc_t8", "mc_all"):
    for b in C.BINS:
        for v in ("pol2", "sym", "lr", "pol1", "pol3"):
            r = tab[f"mcclos|{mk}|{b}|{v}"]
            L.append(f"| {mk} | {b} | {v} | " + " | ".join(f"{x:.3f}±{e:.3f}" for x, e in zip(r["ratio"], r["err"])) + " |")
open(os.path.join(C.OUT, "s1.md"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
