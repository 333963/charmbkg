#!/usr/bin/env python3
"""SD_pred/analyze.py — przewidywania Pythii 8 (POZIOM CZĄSTEK) dla D + proton w oknie AFP (zadanie #8 planu prac z 17.11).
Wejście: merged/<model>.json z Q_lxplus/sdpred (sd_merge.py). Lekkie: tylko json + numpy, bez generowania.
Użycie: ../.venv_h5/bin/python SD_pred/analyze.py [--in SD_pred/runs/lxplus/merged] [--out SD_pred/out]
Wyjście: <out>/tables.md, <out>/results.json.
Konwencje: σ w µb, D z antycząstkami, 2 ≤ pT < 30 GeV, |η| < 2,2, BR naturalne (bez BR kanału); tag = proton 0,030 ≤ ξ ≤ 0,100 po stronie A LUB C
(obie strony zsumowane; D w zdarzeniu z tagiem A i C liczone dwa razy); N na 1 nb⁻¹ = σ[µb]·1000. Błędy: tylko stat. MC (rozrzut między jobami).
"""
import argparse, glob, json, os
import numpy as np

SPL = {"D0": "D⁰", "Dp": "D⁺", "Ds": "D*⁺"}
BR = {"D0": 0.03947, "Dp": 0.0938, "Ds": 0.677 * 0.03947}     # D0→Kπ, D+→Kππ, D*+→D0π+→Kππ (PDG, zaokrąglone)
SOFT = ["ND", "SD", "DD", "CD"]
PTL = ["2–3", "3–5", "5–8", "8–30"]


class Model:
    def __init__(self, path):
        r = json.load(open(path))
        self.name = r["model"]; self.nj = r["njobs"]; self.nok = r["n_ok"]; self.sig = r["sigma_mb"] * 1e3   # µb
        self.f = self.sig / self.nok; self.h = r["hists"]; self.edges = r["edges"]; self.proc = r["proc"]
        self.comps = sorted({k.split("|")[0] for k in self.h})

    def get(self, comps, sps, name):
        """(σ, błąd) [µb] — suma po składowych i gatunkach."""
        s = None; v = None
        for c in comps:
            for sp in sps:
                h = self.h.get(f"{c}|{sp}|{name}")
                if h is None: continue
                a = np.array(h["s"], float); a2 = np.array(h["s2"], float)
                var = self.nj / (self.nj - 1.0) * np.clip(a2 - a * a / self.nj, 0, None) if self.nj >= 5 else a
                var = np.maximum(var, a)          # nie mniej niż Poisson (mało jobów / puste biny)
                s = a if s is None else s + a; v = var if v is None else v + var
        if s is None: return None, None
        return s * self.f, np.sqrt(v) * self.f


def add(*xs):
    xs = [x for x in xs if x[0] is not None]
    if not xs: return None, None
    return sum(x[0] for x in xs), np.sqrt(sum(x[1] ** 2 for x in xs))


def tot(x):
    return (None, None) if x[0] is None else (float(np.sum(x[0])), float(np.sqrt(np.sum(x[1] ** 2))))


def ratio(a, b):
    if a[0] is None or b[0] is None: return None, None
    a0, a1, b0, b1 = (np.asarray(z, float) for z in (*a, *b))
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.where(b0 > 0, a0 / b0, np.nan)
        e = np.abs(r) * np.sqrt(np.where(a0 > 0, (a1 / a0) ** 2, 0) + np.where(b0 > 0, (b1 / b0) ** 2, 0))
    return r, e


def fmt(x, e=None, p=3, pct=False):
    if x is None or (isinstance(x, float) and not np.isfinite(x)): return "–"
    k = 100.0 if pct else 1.0
    s = f"{x * k:.{p}g}" if e is None else f"{x * k:.{p}g} ± {e * k:.2g}"
    return s.replace(".", ",")


def row(cells): return "| " + " | ".join(cells) + " |"


def quant(edges, h, qs):
    """kwantyle z histogramu (środki binów; skrajne biny = under/overflow)."""
    h = np.asarray(h, float)
    if h.sum() <= 0: return [np.nan] * len(qs)
    e = np.array(edges, float); e[0] = e[1] - (e[2] - e[1]); e[-1] = e[-2] + (e[-2] - e[-3])
    c = np.cumsum(h) / h.sum()
    return [float(np.interp(q, c, e[1:])) for q in qs]


def main():
    ap = argparse.ArgumentParser()
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--in", dest="inp", default=os.path.join(here, "runs/lxplus/merged"))
    ap.add_argument("--out", default=os.path.join(here, "out"))
    a = ap.parse_args()
    M = {}
    for p in sorted(glob.glob(os.path.join(a.inp, "*.json"))):
        if os.path.basename(p) == "summary.json": continue
        m = Model(p); M[m.name] = m
    if not M: raise SystemExit(f"brak merged/*.json w {a.inp}")
    os.makedirs(a.out, exist_ok=True)
    L = []; R = {"models": {}, "scenarios": {}}
    w = L.append
    w("# SD_pred: przewidywania Pythii 8.317 dla D + proton w oknie AFP — POZIOM CZĄSTEK (przed wydajnościami detektora)\n")
    w(f"Wejście: `{a.inp}`. σ w µb; D z antycząstkami, 2 ≤ pT < 30 GeV, |η| < 2,2; tag: proton 0,030 ≤ ξ ≤ 0,100 (A lub C, strony zsumowane); błędy tylko stat. MC.\n")
    w("## 0. Próbki\n"); w(row(["model", "joby", "zdarzenia", "σ_gen [mb]", "składowe", "P(proton w oknie, strona) per składowa"])); w(row(["---"] * 6))
    for n, m in M.items():
        pr = {}
        for k, v in m.proc.items():
            c = "ND" if k.startswith("101") else "SD" if k[:3] in ("103", "104") else "DD" if k.startswith("105") else "CD" if k.startswith("106") else "HQ" + ("nd" if k.endswith("_0") else "dif")
            q = pr.setdefault(c, [0, 0]); q[0] += v["n"]; q[1] += v["tagA"] + v["tagC"]
        w(row([n, str(m.nj), f"{m.nok:.3g}", fmt(m.sig / 1e3, p=4), ", ".join(m.comps), "; ".join(f"{c} {fmt(q[1] / 2 / max(q[0], 1), p=3, pct=True)} %" for c, q in pr.items())]))
        R["models"][n] = {"njobs": m.nj, "n_ok": m.nok, "sigma_mb": m.sig / 1e3,
                          "p_tag_side": {c: q[1] / 2 / max(q[0], 1) for c, q in pr.items()}, "sigma_tag_noD_mb": {c: q[1] / max(m.nok, 1) * m.sig / 1e3 for c, q in pr.items()}}

    # ---- scenariusze: miękkie = nd_X + dif_Y; twarde = hq_tN (licznik) / hq_incl (mianownik)
    DUP_FLUX = {"dif_SaS", "dif_DL04", "dif_DL085", "dif_DL10", "dif_H1A", "dif_H1B"}
    SOFT_LABEL = {"dif_A3": "SaS-kształt (tune A3)", "dif_Monash": "SaS-kształt (Monash)", "dif_own_SaS": "strumień SaS (own)", "dif_own_DL07": "strumień DL ε=0,07 (A3)",
                  "dif_own_DL04": "strumień DL ε=0,04", "dif_own_DL085": "strumień DL ε=0,085", "dif_own_DL10": "strumień DL ε=0,10", "dif_own_H1A": "strumień H1 A", "dif_own_H1B": "strumień H1 B"}
    nd = M.get("nd_A3"); scen = []
    for n, m in M.items():
        if n.startswith("dif_"):
            # 17.11: przy SigmaDiffractive:mode = 1 PomFlux nie działa na miękką dyfrakcję (Pythia 8.3) → dif_SaS/DL*/H1* = duplikaty dif_A3; pomijane.
            if n in DUP_FLUX: continue
            ndm = M.get("nd_Monash") if n == "dif_Monash" and "nd_Monash" in M else nd
            scen.append((SOFT_LABEL.get(n, n[4:]), "soft", m, ndm))
    if nd is not None and "nd_Monash" in M and "dif_A3" in M: scen.append(("SaS-kształt+ND_Monash", "soft", M["dif_A3"], M["nd_Monash"]))
    for n, m in M.items():
        if n.startswith("hq_t") and "hq_incl" in M: scen.append((n, "hard", m, M["hq_incl"]))

    def parts(kind, m, den, sps, name):
        """słownik składowa → (σ, błąd) dla histogramu name; mianownik inkluzywny osobno."""
        if kind == "soft":
            d = {c: m.get([c], sps, name) for c in ("SD", "DD", "CD")}
            d["ND"] = den.get(["ND"], sps, name) if den is not None else (None, None)
        else:
            d = {"HQdif": m.get(["HQdif", "HQdd"], sps, name), "HQnd": den.get(["HQnd"], sps, name)}   # HQnd: proton wiodący bez dyfrakcji (z próbki inkluzywnej)
        return d

    def incl(kind, m, den, sps, name="pt_incl"):
        if kind == "soft": return add(m.get(["SD", "DD", "CD"], sps, name), den.get(["ND"], sps, name) if den is not None else (None, None))
        return den.get(["HQnd"], sps, name)

    w("\n## 1. σ fiducjalne D + proton (µb) i R = σ(D + p)/σ_incl(D), 2 ≤ pT < 30 GeV\n")
    w("Scenariusz miękki = ND (A3) + SD + DD + CD z danym strumieniem Pomeronu. Twardy: licznik = HardQCD cc̄+bb̄ z `Diffraction:doHard` (typ 3: bez tłumienia MPI, typ 4: dynamiczna przerwa), "
      "mianownik = ta sama produkcja bez doHard (`hq_incl`); kolumna „ND/HQnd” = proton wiodący bez dyfrakcji. **σ(D + p) liczy obie strony; na jedną stronę ÷ 2.**\n")
    w(row(["scenariusz", "gatunek", "σ_SD", "σ_DD", "σ_CD", "σ_ND (lub HQnd)", "σ(D + p) razem", "σ_incl(D)", "R [%]", "udział nie-SD [%]", "N/nb⁻¹ (BR nat.)", "N/nb⁻¹ × BR kanału"])); w(row(["---"] * 12))
    for sn, kind, m, den in scen:
        S = R["scenarios"].setdefault(sn, {"kind": kind})
        for sp in SPL:
            d = {c: tot(v) for c, v in parts(kind, m, den, [sp], "pt_tag").items()}
            t = add(*d.values()); i = tot(incl(kind, m, den, [sp]))
            r = ratio(t, i); main_c = "SD" if kind == "soft" else "HQdif"
            imp = ratio(add(*[v for c, v in d.items() if c != main_c]), t)
            g = lambda c: fmt(*d[c]) if c in d and d[c][0] is not None else "–"
            w(row([sn, SPL[sp], g("SD") if kind == "soft" else g("HQdif"), g("DD"), g("CD"), g("ND") if kind == "soft" else g("HQnd"), fmt(*t), fmt(*i, p=4),
                   fmt(float(r[0]), float(r[1]), pct=True) if r[0] is not None else "–", fmt(float(imp[0]), float(imp[1]), pct=True) if imp[0] is not None else "–",
                   fmt(t[0] * 1e3, t[1] * 1e3) if t[0] is not None else "–", fmt(t[0] * 1e3 * BR[sp], t[1] * 1e3 * BR[sp]) if t[0] is not None else "–"]))
            S[sp] = {"sigma_tag_ub": t[0], "err": t[1], "sigma_incl_ub": i[0], "R": None if r[0] is None else float(r[0]), "R_err": None if r[0] is None else float(r[1]),
                     "comp_ub": {c: v[0] for c, v in d.items()}, "comp_err": {c: v[1] for c, v in d.items()}, "N_per_nb": None if t[0] is None else t[0] * 1e3,
                     "N_per_nb_channel": None if t[0] is None else t[0] * 1e3 * BR[sp]}

    w("\n## 2. R [%] w binach pT(D)\n"); w(row(["scenariusz", "gatunek"] + PTL)); w(row(["---"] * 6))
    for sn, kind, m, den in scen:
        for sp in SPL:
            t = add(*parts(kind, m, den, [sp], "pt_tag").values()); i = incl(kind, m, den, [sp]); r = ratio(t, i)
            if r[0] is None: continue
            w(row([sn, SPL[sp]] + [fmt(float(x), float(e), pct=True) for x, e in zip(*r)]))
            R["scenarios"][sn][sp]["R_pt"] = [float(x) for x in r[0]]; R["scenarios"][sn][sp]["R_pt_err"] = [float(x) for x in r[1]]
            R["scenarios"][sn][sp]["sigma_tag_pt_ub"] = [float(x) for x in t[0]]
    w("\n### dσ(D + p)/dpT [µb/GeV], suma gatunków\n"); w(row(["scenariusz"] + PTL)); w(row(["---"] * 5))
    wd = np.array([1, 2, 3, 22.0])
    for sn, kind, m, den in scen:
        t = add(*parts(kind, m, den, list(SPL), "pt_tag").values())
        if t[0] is not None: w(row([sn] + [fmt(x / d_, e / d_) for x, e, d_ in zip(t[0], t[1], wd)]))

    w("\n## 3. Zależność od ξ: dσ/dξ [µb] w binach okna i czystość dyfrakcyjna (suma gatunków)\n")
    xb = ["0,030–0,040", "0,040–0,050", "0,050–0,060", "0,060–0,075", "0,075–0,100"]; xw = np.array([0.01, 0.01, 0.01, 0.015, 0.025])
    w(row(["scenariusz", "wielkość"] + xb)); w(row(["---"] * 7))
    for sn, kind, m, den in scen:
        d = parts(kind, m, den, list(SPL), "xi_tag"); t = add(*d.values())
        if t[0] is None: continue
        w(row([sn, "dσ/dξ [µb]"] + [fmt(x / d_, e / d_) for x, e, d_ in zip(t[0], t[1], xw)]))
        for c, v in d.items():
            if v[0] is None or c in ("SD", "HQdif"): continue
            f = ratio(v, t); w(row([sn, f"udział {c} [%]"] + [fmt(float(x), float(e), pct=True) for x, e in zip(*f)]))
        R["scenarios"][sn]["xi"] = {"dsigma_dxi_ub": [float(x) for x in t[0] / xw], "comp": {c: [float(x) for x in v[0]] for c, v in d.items() if v[0] is not None}}
    w("\n### Szeroki zakres ξ protonu wiodącego (0–0,20, krok 0,02): udział nie-SD wśród D z protonem [%]\n")
    w(row(["scenariusz"] + [f"{x:.2f}–{x + 0.02:.2f}".replace(".", ",") for x in np.arange(0, 0.2, 0.02)])); w(row(["---"] * 11))
    for sn, kind, m, den in scen:
        d = parts(kind, m, den, list(SPL), "xiw"); t = add(*d.values())
        if t[0] is None: continue
        mc = "SD" if kind == "soft" else "HQdif"; o = add(*[v for c, v in d.items() if c != mc])
        rb = lambda x: (np.add.reduceat(x[0], np.arange(0, 20, 2)), np.sqrt(np.add.reduceat(x[1] ** 2, np.arange(0, 20, 2))))
        f = ratio(rb(o), rb(t)); w(row([sn] + [fmt(float(x), float(e), pct=True) for x, e in zip(*f)]))

    w("\n## 4. Rapidity i asymetria względem strony protonu (ηs = η·znak strony; ηs > 0 = półkula protonu), suma gatunków\n")
    w(row(["scenariusz", "składowa", "A_η = (N₊ − N₋)/N", "⟨ys⟩", "udział |t| < 1 GeV² [%]", "nachylenie b [GeV⁻²] (0,05 < |t| < 0,5)"])); w(row(["---"] * 6))
    for sn, kind, m, den in scen:
        for c, src in ((("SD", m), ("DD", m), ("CD", m), ("ND", den)) if kind == "soft" else (("HQdif", m), ("HQnd", den))):
            if src is None: continue
            cs = [c] if c != "HQdif" else ["HQdif", "HQdd"]
            h = src.get(cs, list(SPL), "etas_tag")
            if h[0] is None or h[0].sum() <= 0: continue
            n = len(h[0]) // 2; p_, m_ = h[0][n:].sum(), h[0][:n].sum(); N = p_ + m_
            A = (p_ - m_) / N; eA = np.sqrt(max(1 - A * A, 0) / max(N / src.f, 1))
            ys = src.get(cs, list(SPL), "ys_tag"); yc = 0.5 * (np.array(src.edges["ys"][:-1]) + np.array(src.edges["ys"][1:])); ym = float((ys[0] * yc).sum() / ys[0].sum())
            tt = src.get(cs, list(SPL), "t_tag"); te = np.array(src.edges["t"]); ft = float(tt[0][:20].sum() / tt[0].sum())
            sel = (te[:-1] >= 0.05) & (te[1:] <= 0.5) & (tt[0] > 0); b = "–"
            if sel.sum() >= 4:
                x = 0.5 * (te[:-1] + te[1:])[sel]; y = np.log(tt[0][sel]); wy = (tt[0][sel] / tt[1][sel]) ** 2
                cf, cov = np.polyfit(x, y, 1, w=np.sqrt(wy), cov="unscaled"); b = fmt(-cf[0], np.sqrt(cov[0, 0]))
            w(row([sn, c, fmt(A, eA), fmt(ym, p=2), fmt(ft, p=3, pct=True), b]))
            R["scenarios"][sn].setdefault("asym", {})[c] = {"A_eta": float(A), "err": float(eA), "mean_ys": ym, "frac_t_lt1": ft}
    w("\n### dσ/dys [µb na 0,5 jednostki] — D z tagiem, ys = y·znak strony protonu, suma gatunków i składowych\n")
    w(row(["scenariusz"] + [f"{x:+.1f}…{x + 0.5:+.1f}".replace(".", ",") for x in np.arange(-2.5, 2.5, 0.5)])); w(row(["---"] * 11))
    for sn, kind, m, den in scen:
        t = add(*parts(kind, m, den, list(SPL), "ys_tag").values())
        if t[0] is None: continue
        s2 = np.add.reduceat(t[0], np.arange(0, 20, 2)); e2 = np.sqrt(np.add.reduceat(t[1] ** 2, np.arange(0, 20, 2)))
        w(row([sn] + [fmt(x, e) for x, e in zip(s2, e2)])); R["scenarios"][sn]["dsigma_dys_ub_per_0.5"] = [float(x) for x in s2]

    w("\n## 5. ξ z cząstek centralnych względem ξ protonu (D z tagiem, suma gatunków)\n")
    w("ξ_cal = Σ(E ± p_z)/√s po widocznych cząstkach |η| < 4,9; ξ_trk: naładowane |η| < 2,5, pT > 0,5 GeV (ξ_tk1: pT > 0,1). Δ = ξ_central − ξ_p. "
      "Wydajność cięcia ξ_cal głównego łańcucha: ξ_cal ≤ ξ_p + 0,02 (poziom cząstek: bez rozdzielczości kalorymetru, szumu i pile-upu).\n")
    w(row(["scenariusz", "składowa", "zmienna", "mediana Δ", "16 %", "84 %", "mediana ξ_c/ξ_p", "ε(Δ ≤ 0)", "ε(Δ ≤ 0,01)", "**ε(Δ ≤ 0,02)**", "ε(Δ ≤ 0,05)"])); w(row(["---"] * 11))
    for sn, kind, m, den in scen:
        for c, src in ((("SD", m), ("DD", m), ("CD", m), ("ND", den)) if kind == "soft" else (("HQdif", m), ("HQnd", den))):
            if src is None: continue
            cs = [c] if c != "HQdif" else ["HQdif", "HQdd"]
            for v in ("cal", "trk", "tk1"):
                h = src.get(cs, list(SPL), f"d{v}_tag")
                if h[0] is None or h[0].sum() <= 0: continue
                ed = np.array(src.edges["d"]); q = quant(ed, h[0], [0.5, 0.16, 0.84]); N = h[0].sum(); nraw = N / src.f
                eff = []
                for cut in (0.0, 0.01, 0.02, 0.05):
                    k = int(np.argmin(np.abs(ed - cut))); p_ = h[0][:k].sum() / N; eff.append(fmt(p_, np.sqrt(max(p_ * (1 - p_), 0) / max(nraw, 1)), p=3))
                hr = src.get(cs, list(SPL), f"r{v}_tag"); qr = quant(src.edges["r"], hr[0], [0.5])[0]
                w(row([sn, c, "ξ_" + v, fmt(q[0], p=2), fmt(q[1], p=2), fmt(q[2], p=2), fmt(qr, p=2)] + eff))
                if v == "cal":
                    k = int(np.argmin(np.abs(ed - 0.02))); p_ = float(h[0][:k].sum() / N)
                    R["scenarios"][sn].setdefault("xi_cal", {})[c] = {"median_delta": q[0], "q16": q[1], "q84": q[2], "median_ratio": qr, "eff_cut_0.02": p_, "n_raw": float(nraw)}
    w("\n### Wydajność cięcia ξ_cal ≤ ξ_p + 0,02 w binach ξ_p (wszystkie składowe scenariusza, suma gatunków)\n"); w(row(["scenariusz"] + xb + ["razem"])); w(row(["---"] * 7))
    for sn, kind, m, den in scen:
        t = add(*parts(kind, m, den, list(SPL), "xi_tag").values()); p_ = add(*parts(kind, m, den, list(SPL), "xi_tag_cal").values())
        if t[0] is None or p_[0] is None: continue
        f = ratio((p_[0], np.zeros_like(p_[0])), (t[0], np.zeros_like(t[0])))[0]
        w(row([sn] + [fmt(float(x), p=3) for x in f] + [fmt(float(p_[0].sum() / t[0].sum()), p=3)]))
        R["scenarios"][sn]["eff_xical_vs_xi"] = [float(x) for x in f]; R["scenarios"][sn]["eff_xical"] = float(p_[0].sum() / t[0].sum())
    w("\n### ξ_cal w zdarzeniach z D bez wymogu protonu (do oceny, ile tła przypadkowego przepuszcza cięcie): P(ξ_cal ≤ x)\n")
    w(row(["model", "składowa", "x = 0,05", "0,07", "0,10", "0,12", "0,20", "mediana ξ_cal"])); w(row(["---"] * 8))
    for n, m in M.items():
        for c in m.comps:
            h = m.get([c], list(SPL), "xc_incl")
            if h[0] is None or h[0].sum() <= 0: continue
            cum = np.cumsum(h[0]) / h[0].sum()
            w(row([n, c] + [fmt(float(cum[k - 1]), p=3) for k in (5, 7, 10, 12, 20)] + [fmt(quant([-1] + m.edges["xc"][1:], h[0], [0.5])[0], p=2)]))

    w("\n## 6. Rozrzut modeli (oczekiwana niepewność teoretyczna)\n")
    w(row(["gatunek", "R min [%]", "R max [%]", "max/min", "σ(D + p) min [µb]", "max [µb]", "max/min", "scenariusze"])); w(row(["---"] * 8))
    for grp, sel in (("miękkie (strumienie)", lambda s: R["scenarios"][s]["kind"] == "soft"), ("wszystkie (z twardą typ 4; bez typu 3)", lambda s: R["scenarios"][s]["kind"] == "soft" or "_t4" in s)):
        for sp in SPL:
            xs = [(s, R["scenarios"][s][sp]) for s in R["scenarios"] if sel(s) and sp in R["scenarios"][s] and R["scenarios"][s][sp]["R"] is not None]
            if not xs: continue
            rr = [x[1]["R"] for x in xs]; ss = [x[1]["sigma_tag_ub"] for x in xs]
            w(row([f"{SPL[sp]} — {grp}", fmt(min(rr), pct=True), fmt(max(rr), pct=True), fmt(max(rr) / max(min(rr), 1e-12), p=3), fmt(min(ss)), fmt(max(ss)), fmt(max(ss) / max(min(ss), 1e-12), p=3),
                   f"{xs[int(np.argmin(rr))][0]} … {xs[int(np.argmax(rr))][0]}"]))
            R.setdefault("spread", {}).setdefault(grp, {})[sp] = {"R_min": min(rr), "R_max": max(rr), "sigma_min_ub": min(ss), "sigma_max_ub": max(ss)}
    t3 = [s for s in R["scenarios"] if "_t3" in s]
    if t3:
        w("\n### Przeżycie przerwy w modelu dynamicznym: S² = σ(typ 4)/σ(typ 3), D + p w oknie\n"); w(row(["strumień"] + list(SPL.values()))); w(row(["---"] * 4))
        for s3 in t3:
            s4 = s3.replace("_t3", "_t4")
            if s4 not in R["scenarios"]: continue
            c = []
            for sp in SPL:
                a3, a4 = R["scenarios"][s3].get(sp), R["scenarios"][s4].get(sp)
                if not a3 or not a4 or not a3["sigma_tag_ub"]: c.append("–"); continue
                r_ = a4["sigma_tag_ub"] / a3["sigma_tag_ub"]; c.append(fmt(r_, r_ * np.hypot(a3["err"] / a3["sigma_tag_ub"], a4["err"] / max(a4["sigma_tag_ub"], 1e-30)), p=3))
            w(row([s3.replace("hq_t3", "A3-DL" if s3 == "hq_t3" else "").strip("_")] + c)); R["scenarios"][s4]["S2_vs_t3"] = c
    w("\nUwaga: wszystko powyżej to poziom cząstek — bez akceptancji AFP w (ξ, t, φ), wydajności rekonstrukcji D, CharmNet, trygera i weta wierzchołków.")
    open(os.path.join(a.out, "tables.md"), "w").write("\n".join(L) + "\n")
    json.dump(R, open(os.path.join(a.out, "results.json"), "w"), indent=1, default=float)
    print(f"OK: {len(M)} modeli, {len(scen)} scenariuszy → {a.out}/tables.md, results.json")


if __name__ == "__main__":
    main()
