#!/usr/bin/env python3
"""KONSTR B — tło fitu tagowanego z konstrukcji (KRYTERIA §B). POST HOC po odczycie TAGGED; wyniki lokalne.
kb_tag.py prep [--nrep 200]  → test kształtu pasm (wagi), zamknięcie na próbkach pseudo-tagowanych; okno tagowane NIE jest wczytywane (blind=True)
kb_tag.py window             → fity okna: warianty konstrukcji + odtworzenie TAGGED"""
import json, math, os, sys
import numpy as np
from scipy.stats import chi2 as C2
H_ = os.path.dirname(os.path.abspath(__file__)); T_ = os.path.dirname(H_); W = os.path.dirname(T_)
sys.path.insert(0, os.path.join(W, "TAGGED_2026-11-19", "skrypty"))
import tg_lib as T
F, V = T.F, T.V
RUN, REG, E0 = "427929", "MBTS", 0.7
JS = T.jwp(E0) + 1
MOD = {"D0": {"M8": "M8", "M7": "M7", "M0": "M0"}, "Dplus": {"M0": "M0", "M7": "M7"}}
NSB = [0, 9, 13, 17, 21, 26, 33, 100000]       # przedziały nSelBL ~ równoliczne w PASMACH próbki tagowanej (kwantyle 10–90 %: 8–32; nietagowana: 17–66)
NC = len(NSB) - 1


def load_incl_rows():
    z = np.load(f"{T_}/in/kx_data_{RUN}.npz")
    R = {k[2:]: z[k] for k in z.files if k.startswith("R_")}
    q = R["s"] >= JS
    return {k: v[q] for k, v in R.items()}


def cells(nsel, ipt, scheme):
    a = np.digitize(nsel, NSB[1:-1])
    if scheme == "brak":
        return np.zeros(len(nsel), int), 1
    if scheme == "nSelBL":
        return a, NC
    return a * 3 + ipt, 3 * NC


def sbmask_rows(m, sp):
    m0, w = T.WIN[sp]
    return np.abs(m - m0) >= w


def weights(Rm, Rn, Ri, tm, tn, ti, sp, scheme):
    """wagi wierszy nietagowanych: stosunek rozkładów komórek (pasma tagowane)/(pasma nietagowane)."""
    cr, nc = cells(Rn, Ri, scheme); ct, _ = cells(tn, ti, scheme)
    ht = np.bincount(ct[sbmask_rows(tm, sp)], minlength=nc).astype(float); hr = np.bincount(cr[sbmask_rows(Rm, sp)], minlength=nc).astype(float)
    w = (ht / max(ht.sum(), 1)) / np.maximum(hr / hr.sum(), 1e-12)
    return w[cr]


def incl_model(Rm, Ri, w, sp, mname, rng):
    """fit ważonego rozkładu nietagowanego per bin pT → (b̂ w binach 5 MeV zakresu rng [nieunormowane], Ŝ unormowane tak, że współczynnik = N przy szablonie sygnału, ΣN ważone)."""
    if V.Z is None:
        V._init()
    e = T.edges5(sp, rng); fe = V.full_edges(sp)
    S = V.SET[mname]; c0 = F.SP[sp]; frng = c0[S["rng"]]
    ctr = 0.5 * (e[1:] + e[:-1]); x = 2 * (ctr - frng[0]) / (frng[1] - frng[0]) - 1
    nfit = int(round((frng[1] - frng[0]) / F.BW))
    B = np.zeros(len(ctr)); Sg = np.zeros(len(ctr)); Ntot = 0.0
    for ib, bn in enumerate(T.BINS):
        q = Ri == ib
        ww = w[q] * q.sum() / w[q].sum()
        H5 = np.histogram(Rm[q], fe, weights=ww)[0].astype(float)
        cat = V.wcat(sp, bn, E0); cat["kap"] *= T.KSC[sp]
        f = V.fit_set(H5[None], [cat], sp, mname, hesse=False)
        p = f["x"]; nb = c0["nbkg"] + S["dn"]
        xf = np.linspace(-1 + 1.0 / nfit, 1 - 1.0 / nfit, nfit)
        nrm = F.bshape(xf, p[4:4 + nb], "exp").sum()
        scale = w[q].sum() / q.sum()                               # powrót do skali Σw
        B += scale * p[3] * F.bshape(x, p[4:4 + nb], "exp") / nrm
        s = cat["sigT"].bins(e, c0["m0"], p[0], p[1])
        if cat["refT"] is not None:
            s = s + cat["frefl"] * cat["refT"].bins(e, c0["m0"], p[0], 1.0)
        if S["peak"]:
            s = s + cat["kap"] * cat["peakT"].bins(e, c0["m0"], p[0], 1.0)
        Sg += scale * p[2] * s; Ntot += scale * p[2]
    return B, Sg / Ntot, Ntot


class CFit(T.TFit):
    """tło o ustalonym kształcie b̂: μ = N_S·S + N_B·b̂/Σb̂."""

    def __init__(self, S, bhat):
        super().__init__(S, len(S), 0)
        self.b = np.asarray(bhat, float) / np.sum(bhat)

    def bsh(self, a):
        return self.b


def tagged(blind):
    C = T.load_cand(RUN, blind=blind)
    sel = C["T2"] & T.region_mask(C, RUN, REG) & C["base"] & (C["s"] >= JS)
    return {sp: {k: C[k][sel & (C["sp"] == isp)] for k in ("m", "nsel", "ipt", "pt")} for isp, sp in enumerate(T.SPN)}


def sb_chi2(tm, Rm, w, sp):
    lo, hi = T.RNG[sp]; ed = np.arange(lo, hi + 1e-9, 0.025); ctr = 0.5 * (ed[1:] + ed[:-1]); m0, wn = T.WIN[sp]
    sb = (ctr + 0.0125 <= m0 - wn + 0.006) | (ctr - 0.0125 >= m0 + wn - 0.006)
    ht = np.histogram(tm, ed)[0][sb].astype(float); hr = np.histogram(Rm, ed, weights=w)[0][sb]; hr2 = np.histogram(Rm, ed, weights=w ** 2)[0][sb]
    k = ht.sum() / hr.sum(); ex = k * hr
    c2 = float(np.sum((ht - ex) ** 2 / (ex + k ** 2 * hr2)))
    # nachylenie stosunku
    xx = ctr[sb] - m0; r = ht / ex; er = np.sqrt(ex + k ** 2 * hr2) / ex
    A = np.stack([np.ones(len(xx)), xx], 1); cov = np.linalg.inv(A.T @ (A / er[:, None] ** 2)); sl = cov @ (A.T @ (r / er ** 2))
    return c2, int(sb.sum() - 1), float(sl[1]), float(math.sqrt(cov[1, 1]))


CLO = None
_PID = None


def clo_one(seed):
    Rm, Rn, Ri, pr, n, sp, proc_s, mname, S = CLO
    global _PID
    if _PID != os.getpid():
        V._init(); _PID = os.getpid()          # własne uchwyty npz w procesie potomnym
    rng = np.random.default_rng(seed)
    ix = rng.choice(len(Rm), n, replace=False, p=pr)
    pm, pn, pi = Rm[ix], Rn[ix], Ri[ix]
    sbr = sbmask_rows(pm, sp)
    w = weights(Rm, Rn, Ri, pm[sbr], pn[sbr], pi[sbr], sp, proc_s)
    B, _, _ = incl_model(Rm, Ri, w, sp, mname, None)
    Hh = T.hist5(pm, sp)
    q = CFit(S, B).full(Hh, cl95=False)
    nbk, _, _ = T.bkg_order(S, Hh, ~T.winmask(sp))
    q2 = T.TFit(S, len(Hh), nbk).full(Hh, cl95=False)
    return [q["N"], q["lo68"], q["hi68"], q["hesse"]], [q2["N"], q2["lo68"], q2["hi68"], q2["hesse"]]


def main():
    mode = sys.argv[1]
    nrep = int(sys.argv[sys.argv.index("--nrep") + 1]) if "--nrep" in sys.argv else 200
    R = load_incl_rows()
    out, L = {}, []
    TG = tagged(blind=(mode == "prep"))
    if mode == "prep":
        L += ["## Kształt pasm: próbka tagowana wobec nietagowanej (biny 25 MeV, tylko pasma)", "", "| gatunek | wagi | χ²/ndf (p) | nachylenie stosunku [1/GeV] |", "|---|---|---|---|"]
        for isp, sp in enumerate(T.SPN):
            r = R["sp"] == isp; t = TG[sp]
            for sch in ("brak", "nSelBL", "nSelBL × pT"):
                w = weights(R["m"][r], R["nsel"][r], R["ipt"][r], t["m"], t["nsel"], t["ipt"], sp, sch)
                c2, ndf, sl, esl = sb_chi2(t["m"], R["m"][r], w, sp)
                out[f"sb|{sp}|{sch}"] = dict(chi2=c2, ndf=ndf, p=float(C2.sf(c2, ndf)), slope=sl, eslope=esl, n_tag_sb=int(len(t["m"])), neff=float(w.sum() ** 2 / (w ** 2).sum()))
                L.append(f"| {sp} | {sch} | {c2:.1f}/{ndf} ({C2.sf(c2, ndf):.2g}) | {sl:+.2f} ± {esl:.2f} |")
        json.dump(out, open(f"{T_}/out/kb_prep_sb.json", "w"), indent=1)
        nom = {}
        for sp in T.SPN:
            ok = [s for s in ("brak", "nSelBL", "nSelBL × pT") if out[f"sb|{sp}|{s}"]["p"] >= 0.05]
            nom[sp] = ok[0] if ok else None
        out["nominal_scheme"] = nom
        L += ["", f"Schemat nominalny (najprostszy z p ≥ 0,05): {nom}", ""]
        print("\n".join(L), flush=True)
        # ---- zamknięcie
        I, _ = T.load_incl(RUN)
        L += ["## Zamknięcie na próbkach pseudo-tagowanych (losowanie ∝ waga z wierszy nietagowanych; pełna procedura na każdej)", "",
              "| gatunek | wagi: prawda / procedura | n | N_oczek. | konstrukcja: ⟨N̂ − N⟩ ± błąd | rms N̂ | ⟨pull⟩ | rms pull | pokrycie 68 % | TAGGED (tło wolne): ⟨N̂ − N⟩ | rms N̂ | ⟨pull⟩ | pokrycie | status |", "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
        import multiprocessing as mp; Pool = mp.get_context("fork").Pool
        for isp, sp in enumerate(T.SPN):
            r = np.where(R["sp"] == isp)[0]; t = TG[sp]
            Rm, Rn, Ri = R["m"][r], R["nsel"][r], R["ipt"][r]
            fits = T.incl_fit(I, sp, "MB_npu0", E0); S = T.sig_shape(sp, fits)
            for truth_s, proc_s in (("brak", "brak"), ("nSelBL × pT", "brak"), ("nSelBL × pT", "nSelBL × pT")):
                wt = weights(Rm, Rn, Ri, t["m"], t["nsel"], t["ipt"], sp, truth_s)
                mname = T.MODEL[sp]
                _, _, Nw = incl_model(Rm, Ri, wt, sp, mname, None)
                fsb = np.sum(wt[sbmask_rows(Rm, sp)]) / wt.sum()       # liczność w zakresie: pasma tagowane znane, okno z proporcji nietagowanej ważonej
                n = int(round(len(t["m"]) / fsb)); nexp = Nw * n / wt.sum()
                global CLO
                CLO = (Rm, Rn, Ri, wt / wt.sum(), n, sp, proc_s, mname, S)
                with Pool(7) as pool:
                    rr = pool.map(clo_one, [1000 * isp + 7 * k + 100000 * (len(truth_s) + 3 * len(proc_s)) for k in range(nrep)])
                res, rest = [x[0] for x in rr], [x[1] for x in rr]

                def summ(a):
                    a = np.array(a); N, lo_, hi_, he = a.T
                    sig = np.where(N > nexp, N - lo_, hi_ - N); sig = np.where(np.isfinite(sig) & (sig > 0), sig, he)
                    pull = (N - nexp) / sig
                    return dict(bias=float(N.mean() - nexp), ebias=float(N.std() / math.sqrt(len(N))), rms=float(N.std()), pull=float(pull.mean()), pull_rms=float(pull.std()), cov=float(np.mean((lo_ <= nexp) & (nexp <= hi_))))
                a, b = summ(res), summ(rest)
                ok = abs(a["pull"]) <= 0.2 and 0.60 <= a["cov"] <= 0.76
                out[f"clo|{sp}|{truth_s}|{proc_s}"] = dict(n=n, nexp=float(nexp), konstr=a, tagged=b, ok=bool(ok), nrep=nrep)
                L.append(f"| {sp} | {truth_s} / {proc_s} | {n} | {nexp:.1f} | {a['bias']:+.1f} ± {a['ebias']:.1f} | {a['rms']:.1f} | {a['pull']:+.2f} | {a['pull_rms']:.2f} | {a['cov']:.2f} | {b['bias']:+.1f} ± {b['ebias']:.1f} | {b['rms']:.1f} | {b['pull']:+.2f} | {b['cov']:.2f} | {'OK' if ok else 'NIE'} |")
                print(L[-1], flush=True)
        name = "B_prep"
    else:
        prep = json.load(open(f"{T_}/out/kb_prep.json"))
        I, _ = T.load_incl(RUN)
        L += ["| gatunek | wariant | N_tag | −/+ (68 %) | Z | χ²/ndf (p) | tło w oknie |", "|---|---|---|---|---|---|---|"]
        for isp, sp in enumerate(T.SPN):
            r = R["sp"] == isp; t = TG[sp]; Rm, Rn, Ri = R["m"][r], R["nsel"][r], R["ipt"][r]
            nom = prep["nominal_scheme"][sp] or "brak"
            fits = T.incl_fit(I, sp, "MB_npu0", E0)
            # odtworzenie TAGGED
            S = T.sig_shape(sp, fits); Hh = T.hist5(t["m"], sp); msb = ~T.winmask(sp)
            Hsb = Hh.copy()
            nbk, _, _ = T.bkg_order(S, Hsb, msb)
            q = T.TFit(S, len(Hh), nbk).full(Hh, cl95=False)
            wm = T.winmask(sp)
            vsb, _, _, _ = T.TFit(S, len(Hh), nbk, mask=msb).fit(Hh, fixNS=0.0)
            tsb = T.TFit(S, len(Hh), nbk, mask=msb); bsb = (tsb.bsh(vsb[2:]) * vsb[1])[wm].sum()
            nwin = Hh[wm].sum(); nsb_est = (nwin - bsb) / S[wm].sum()
            rows = [("TAGGED (tło wolne, rząd z pasm) — odtworzenie", q, None)]
            out[f"{sp}|okno−pasma"] = dict(N=float(nsb_est), okno=float(nwin), tlo=float(bsb))
            var = [("konstrukcja nominalna (M8/M0 nom., wagi: %s)" % nom, T.MODEL[sp], nom, None)]
            for sch in ("brak", "nSelBL", "nSelBL × pT"):
                if sch != nom:
                    var.append((f"wagi: {sch}", T.MODEL[sp], sch, None))
            var.append(("zakres " + "–".join(str(x) for x in T.ALT[sp]), T.MODEL[sp], nom, T.ALT[sp]))
            if sp == "D0":
                var.append(("model nietagowany M7 (bez DMM + CM; szerokie składowe w tle)", "M7", nom, None))
                var.append(("model nietagowany M0 (bez DMM + CM, 1,75–2,00)", "M0", nom, (1.75, 2.00)))
            for lab, mname, sch, rg in var:
                w = weights(Rm, Rn, Ri, t["m"][sbmask_rows(t["m"], sp)], t["nsel"][sbmask_rows(t["m"], sp)], t["ipt"][sbmask_rows(t["m"], sp)], sp, sch)
                B, Sm, Nw = incl_model(Rm, Ri, w, sp, mname, rg)
                Suse = T.sig_shape(sp, fits, rng=rg) if mname == T.MODEL[sp] else Sm
                Hh = T.hist5(t["m"], sp, rg)
                q = CFit(Suse, B).full(Hh, cl95=False)
                q["Bwin"] = float(q["NB"] * (B / B.sum())[T.winmask(sp, rg)].sum())
                rows.append((lab, q, None))
            for lab, q, _ in rows:
                out[f"{sp}|{lab}"] = {k: q[k] for k in ("N", "lo68", "hi68", "Z", "chi2", "ndf", "NB")}
                L.append(f"| {sp} | {lab} | {q['N']:.1f} | −{q['N'] - q['lo68']:.1f} +{q['hi68'] - q['N']:.1f} | {q['Z']:.2f} | {q['chi2']:.1f}/{q['ndf']} ({C2.sf(q['chi2'], q['ndf']):.2g}) | {q.get('Bwin', float('nan')):.0f} |")
            L.append(f"| {sp} | okno − pasma (arytmetyka TAGGED §7) | {nsb_est:.1f} | – | – | – | {bsb:.0f} |")
        name = "B_wyniki"
    json.dump(out, open(f"{T_}/out/kb_{mode}.json", "w"), indent=1)
    open(f"{T_}/tabele/{name}.md", "w").write(f"# B — tło fitu tagowanego z konstrukcji ({mode}); POST HOC, lokalnie\n\n" + "\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
