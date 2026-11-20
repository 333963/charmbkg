#!/usr/bin/env python3
"""KONSTR C — pik D0 2–3 GeV: tło = rozkład masy z odwróconego cięcia × exp(c1 x [+ c2 x²]) (KRYTERIA §C).
kc_d23.py prep   → testy rzeźbienia: pasma danych (okno zamaskowane), mb_nd 3–5 GeV samo tło   (out/kc_prep.json, tabele/C_prep.md)
kc_d23.py window → fity okna: funkcyjne (odtworzenie D23) i konstrukcyjne                      (out/kc_window.json, tabele/C_wyniki.md)"""
import json, math, os, sys
import numpy as np
from scipy.stats import chi2 as C2
H_ = os.path.dirname(os.path.abspath(__file__)); T = os.path.dirname(H_); W = os.path.dirname(T)
sys.path.insert(0, os.path.join(W, "D23_D0_weryfikacja_2026-11-19", "scripts"))
import vlib as V
from iminuit import Minuit
M0, WIN = V.M0, V.WIN
FE = np.arange(1.70, 2.05 + 1e-9, 0.005)
HEAD = {("427929", 0): 0.1, ("427929", 1): 0.2, ("428770", 0): 0.1, ("428770", 1): 0.1}
RNGS = {"1,75–2,00": (1.75, 2.00), "1,70–2,05": (1.70, 2.05)}


def cut(h, rng):
    i0, i1 = int(round((rng[0] - 1.70) / 0.005)), int(round((rng[1] - 1.70) / 0.005))
    return h[i0:i1], FE[i0:i1 + 1]


def cfit(H, edges, S, Tc, order=1, fixN=None, mask=None):
    """μ = N·S + NB·Tc·exp(Σ c_k x^k)/norma. mask = biny użyte."""
    ctr = 0.5 * (edges[1:] + edges[:-1]); x = 2 * (ctr - edges[0]) / (edges[-1] - edges[0]) - 1
    use = np.ones(len(H), bool) if mask is None else mask

    def comps(p):
        sh = Tc * np.exp(np.clip(sum(p[2 + k] * x ** (k + 1) for k in range(order)), -30, 30))
        return p[0] * S, p[1] * sh / sh[use].sum()

    def cost(p):
        mu = np.maximum(sum(comps(p)), 1e-9)[use]; n = H[use]
        return 2.0 * float(np.sum(mu - n) + np.sum(np.where(n > 0, n * np.log(np.where(n > 0, n, 1.0) / mu), 0.0)))
    m = Minuit(cost, [0.0 if fixN is None else fixN, H[use].sum()] + [0.0] * order)
    m.errordef = 1.0; m.limits[1] = (0, None)
    if fixN is not None:
        m.fixed[0] = True
    m.migrad(ncall=20000)
    if not m.valid:
        m.simplex(); m.migrad(ncall=40000)
    m.hesse()
    v = np.array(m.values); c = comps(v); mu = np.maximum(sum(c), 1e-9)
    return dict(N=float(v[0]), err=float(m.errors[0]), fval=float(m.fval), valid=bool(m.valid), chi2=float(np.sum(((H - mu) ** 2 / mu)[use])), ndf=int(use.sum() - m.nfit), c=[float(t) for t in v[2:]],
                B=float(c[1][np.abs(ctr - M0) < WIN].sum()))


def sshape(q, edges):
    s = q["sigT"].bins(edges, M0)
    return s + (q["frefl"] * q["refT"].bins(edges, M0) if q["refT"] is not None else 0.0)


def ctrl_hists(d, j, thr, cat=None):
    """próbki kontrolne (hist 5 MeV pełnego zakresu): K1 ε 0,5–0,3; K2 odrzucone (score < WP50, z histogramu wszystkich kandydatów); K3 ε 0,5–0,1 (0,5–0,2 dla punktu WP20 → patrz wywołanie)."""
    lg, m = d["lg"][:, j], d["m"]
    sel = np.ones(len(m), bool) if cat is None else np.isin(d["cat"], cat)
    k = {}
    k["K1 ε 0,5–0,3"] = np.histogram(m[sel & (lg >= thr[0.5]) & (lg < thr[0.3])], FE)[0].astype(float)
    Hh = d["H"][:, j] if cat is None else d["H"][cat, j]
    nb = int(np.searchsorted(d["LE"], thr[0.5], side="right") - 1)        # biny score w całości poniżej progu WP50
    k["K2 odrzucone (ε > 0,5)"] = Hh[..., :nb].sum(-1).sum(0)
    return k, sel


def main(mode):
    F = V.load("fd_d23")
    THR = {j: {wp: V.mc_wp(F, j, wp) for wp in (0.5, 0.4, 0.3, 0.2, 0.1)} for j in (0, 1)}
    out, L = {}, []
    if mode == "prep":
        L += ["## (a) pasma danych 2–3 GeV (okno |m − m0| < 55 MeV zamaskowane), fit samego tła konstrukcyjnego", "",
              "| run | klasyfikator | próbka kontrolna | zakres | n_ctrl/n_sygn | χ²/ndf rząd 1 (p) | 2ΔlnL(2 vs 1) | status |", "|---|---|---|---|---|---|---|---|"]
        for (run, j), wp in HEAD.items():
            d = V.load(f"data_{run}_d23")
            tight = np.histogram(d["m"][d["lg"][:, j] >= THR[j][wp]], FE)[0].astype(float)
            ks, _ = ctrl_hists(d, j, THR[j])
            for kn, hk in ks.items():
                for rn, rg in RNGS.items():
                    Hc, e = cut(tight, rg); Tc, _ = cut(hk, rg)
                    ctr = 0.5 * (e[1:] + e[:-1]); sb = np.abs(ctr - M0) >= WIN
                    f1 = cfit(Hc, e, np.zeros(len(Hc)), Tc, 1, fixN=0.0, mask=sb); f2 = cfit(Hc, e, np.zeros(len(Hc)), Tc, 2, fixN=0.0, mask=sb)
                    p = float(C2.sf(f1["chi2"], f1["ndf"])); dq = f1["fval"] - f2["fval"]
                    ok = p >= 0.01 and dq < 3.84
                    out[f"a|{run}|{V.NM[j]}|{kn}|{rn}"] = dict(chi2=f1["chi2"], ndf=f1["ndf"], p=p, dq=dq, ok=bool(ok), ratio=float(Tc.sum() / Hc.sum()), c1=f1["c"])
                    L.append(f"| {run} | {V.NM[j]} WP{int(wp * 100)} | {kn} | {rn} | {Tc.sum() / Hc.sum():.0f} | {f1['chi2']:.1f}/{f1['ndf']} ({p:.2g}) | {dq:.2f} | {'OK' if ok else 'NIE'} |")
        # (b) MC mb_nd 3–5 GeV, samo tło (kat. 2, 3)
        Fd, Mb = V.load("fd_nd_d150"), V.load("mb_nd_d150")
        T35 = {j: {wp: V.mc_wp(Fd, j, wp) for wp in (0.5, 0.3, 0.2, 0.1)} for j in (0, 1)}
        L += ["", "## (b) mb_nd 3–5 GeV, samo tło z prawdy (kategorie 2, 3): sygnał pozorny przy tle konstrukcyjnym (rząd 1) i funkcyjnym exp2", "",
              "| klasyfikator | WP | próbka kontrolna | zakres | N_spur ± σ (konstr.) | Z | N_spur (exp2) | B_okno | N_spur/B_okno | status (|Z| ≤ 2) |", "|---|---|---|---|---|---|---|---|---|---|"]
        for j in (0, 1):
            ks, sel = ctrl_hists(Mb, j, T35[j], cat=[2, 3])
            for wp in (0.1, 0.2):
                q = V.templ(Fd, j, T35[j][wp])
                tight = np.histogram(Mb["m"][sel & (Mb["lg"][:, j] >= T35[j][wp])], FE)[0].astype(float)
                for kn, hk in ks.items():
                    for rn, rg in RNGS.items():
                        Hc, e = cut(tight, rg); Tc, _ = cut(hk, rg)
                        f = cfit(Hc, e, sshape(q, e), Tc, 1)
                        g = V.fit2(Hc, e, q, nbkg=2)
                        ok = abs(f["N"] / f["err"]) <= 2
                        out[f"b|{V.NM[j]}|{wp}|{kn}|{rn}"] = dict(N=f["N"], err=f["err"], N_exp2=float(g["N"]), err_exp2=float(g["err"]), B=f["B"], ok=bool(ok), chi2=f["chi2"], ndf=f["ndf"])
                        L.append(f"| {V.NM[j]} | {int(wp * 100)} | {kn} | {rn} | {f['N']:+.0f} ± {f['err']:.0f} | {f['N'] / f['err']:+.1f} | {g['N']:+.0f} ± {g['err']:.0f} | {f['B']:.0f} | {f['N'] / f['B']:+.3f} | {'OK' if ok else 'NIE'} |")
        name = "C_prep"
    else:
        L += ["| run | klasyfikator | tło | zakres | N | σ | Z_LR | χ²/ndf (p) |", "|---|---|---|---|---|---|---|---|"]
        S = ["", "## Rozstęp Z_LR", "", "| run | klasyfikator | funkcyjne exp2–4 × 2 zakresy: min–max (rozstęp) | konstrukcyjne (K1, K2 × rząd 1, 2 × 2 zakresy): min–max (rozstęp) | tylko zaliczone w prep (rząd 1): min–max | N konstr. min–max |", "|---|---|---|---|---|---|"]
        prep = json.load(open(f"{T}/out/kc_prep.json"))
        for (run, j), wp in HEAD.items():
            d = V.load(f"data_{run}_d23"); th = THR[j][wp]; q = V.templ(F, j, th)
            tight = np.histogram(d["m"][d["lg"][:, j] >= th], FE)[0].astype(float)
            ks, _ = ctrl_hists(d, j, THR[j])
            zf, zc, zok, nc = [], [], [], []
            for rn, rg in RNGS.items():
                Hc, e = cut(tight, rg)
                for nb in (2, 3, 4):
                    f = V.fit2(Hc, e, q, nbkg=nb); f0 = V.fit2(Hc, e, q, nbkg=nb, fixN=0.0)
                    z = math.sqrt(max(f0["fval"] - f["fval"], 0)); zf.append(z)
                    out[f"{run}|{V.NM[j]}|exp{nb}|{rn}"] = dict(N=float(f["N"]), err=float(f["err"]), Z=z, chi2=f["chi2"], ndf=int(f["ndf"]))
                    L.append(f"| {run} | {V.NM[j]} WP{int(wp * 100)} | exp{nb} | {rn} | {f['N']:.0f} | {f['err']:.0f} | {z:.1f} | {f['chi2']:.1f}/{f['ndf']} ({C2.sf(f['chi2'], f['ndf']):.2g}) |")
                for kn, hk in ks.items():
                    Tc0, _ = cut(hk, rg)
                    for order in (1, 2):
                        Tc, Nprev = Tc0.copy(), 0.0
                        # poprawka na sygnał w próbce kontrolnej (tylko K1): sygnał kontrolny = N/ε · Δε, kształt wycinka
                        if kn.startswith("K1"):
                            qc = V.templ(F, j, THR[j][0.5], THR[j][0.3]); sc_ = sshape(qc, e)
                            for _ in range(4):
                                f = cfit(Hc, e, sshape(q, e), Tc, order)
                                Tc = np.maximum(Tc0 - f["N"] / q["eps"] * qc["eps"] * sc_, 0.0)
                        f = cfit(Hc, e, sshape(q, e), Tc, order); f0 = cfit(Hc, e, sshape(q, e), Tc, order, fixN=0.0)
                        z = math.sqrt(max(f0["fval"] - f["fval"], 0)); zc.append(z); nc.append(f["N"])
                        if order == 1 and prep[f"a|{run}|{V.NM[j]}|{kn}|{rn}"]["ok"]:
                            zok.append(z)
                        out[f"{run}|{V.NM[j]}|{kn}|rz{order}|{rn}"] = dict(N=f["N"], err=f["err"], Z=z, chi2=f["chi2"], ndf=f["ndf"], contam=float((Tc0 - Tc).sum() / Tc0.sum()))
                        L.append(f"| {run} | {V.NM[j]} WP{int(wp * 100)} | {kn} × exp(rząd {order}) | {rn} | {f['N']:.0f} | {f['err']:.0f} | {z:.1f} | {f['chi2']:.1f}/{f['ndf']} ({C2.sf(f['chi2'], f['ndf']):.2g}) |")
            okt = f"{min(zok):.1f}–{max(zok):.1f}" if zok else "brak"
            out[f"sum|{run}|{V.NM[j]}"] = dict(zf=[min(zf), max(zf)], zc=[min(zc), max(zc)], zok=[min(zok), max(zok)] if zok else None, nc=[min(nc), max(nc)])
            S.append(f"| {run} | {V.NM[j]} WP{int(wp * 100)} | {min(zf):.1f}–{max(zf):.1f} ({max(zf) - min(zf):.1f}) | {min(zc):.1f}–{max(zc):.1f} ({max(zc) - min(zc):.1f}) | {okt} | {min(nc):.0f}–{max(nc):.0f} |")
        L += S
        name = "C_wyniki"
    json.dump(out, open(f"{T}/out/kc_{mode}.json", "w"), indent=1)
    open(f"{T}/tabele/{name}.md", "w").write(f"# C — pik D⁰ 2–3 GeV, tło z odwróconego cięcia ({mode})\n\n" + "\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main(sys.argv[1])
