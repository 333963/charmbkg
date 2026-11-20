#!/usr/bin/env python3
"""MC_templates — test szablonów tła z MC (zadanie #7 planu prac; projekt B). Lekka analiza gotowych histogramów
(`mctmpl_hist.py`); 1 rdzeń, kilka minut (~1500 małych fitów). Kryteria: MC_templates/KRYTERIA.md (zapisane przed wynikami).
Uruchamiana w węźle `merge` na condorze (wynik: $MT_OUT/ana/) oraz lokalnie przez wrapper:
  .venv_h5/bin/python charm_work/MC_templates/analyze.py [--in MC_templates/runs/lxplus] [--out MC_templates/out]
--min-n X skaluje progi minimalnej statystyki (tylko test kodu na 1–2 plikach).
Wynik: <out>/WYNIKI.md, <out>/results.json. Region sygnału głównego łańcucha (tag AFP ∧ npu = 0 ∧ ξ_cal) nie występuje w wejściu.
"""
import argparse
import json
import math
import os
import time
import warnings

import numpy as np
from iminuit import Minuit
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
WPN = ["all", "95", "80", "50"]                      # indeks osi po skumulowaniu sb
FITWP = [1, 2, 3]
DATASETS = [("427929:all", "data_427929", "all"), ("427929:MB", "data_427929", "MB"), ("428770:all", "data_428770", "all"),
            ("428770:J12only", "data_428770", "J12only"), ("428770:MBTSbased", "data_428770", "MBTSbased"),
            ("428770:AFPstream", "data_428770", "AFPstream")]
COMP = {"D0": ["Reflection", "CharmMismatch", "DecayModeMismatch", "PartComb"], "Dplus": ["CharmMismatch", "DecayModeMismatch", "PartComb"]}
PRIOR_PC = 0.25


# ------------------------------------------------------------------------------------------------------------ wejście
def load(d, s):
    z = np.load(os.path.join(d, f"mctmpl_{s}.npz"), allow_pickle=False)
    meta = json.loads(str(z["meta"]))
    H = {}
    for k in z.files:
        if k != "meta":
            h = z[k].astype(float)
            H[k] = np.flip(np.cumsum(np.flip(h, 1), 1), 1)          # oś 1: sb → WP włączne (all, 95, 80, 50)
    return H, meta


def grp(meta, name):
    return meta["groups"][name]


def reb(h, k=2):
    return h.reshape(h.shape[:-1] + (h.shape[-1] // k, k)).sum(-1)


def edges(ax):
    return np.linspace(ax[0], ax[1], int(ax[2]) + 1)


def fmt(x, e=None, p=3):
    if x is None or not np.isfinite(x):
        return "—"
    return f"{x:.{p}f}" if e is None or not np.isfinite(e) else f"{x:.{p}f} ± {e:.{p}f}"


# ------------------------------------------------------------------------------------------------------------ wagi aktywności
def act_weights(hd, hm, sbmask):
    """w(cell) = dane/MC w pasmach bocznych (WP all). hd: [cell, m] danych (klasa), hm: [cell, m] MC (suma kategorii)."""
    d, m = hd[:, sbmask].sum(1), hm[:, sbmask].sum(1)
    w = np.where(m > 0, d / np.maximum(m, 1e-9), 0.0)
    lost = d[m == 0].sum() / max(d.sum(), 1)
    pos = w[w > 0]
    if pos.size:
        w = np.minimum(w, 20 * np.median(pos))
    tot = (w * m).sum()
    w = w * (m.sum() / tot) if tot > 0 else np.ones_like(w)
    neff = (w * m).sum() ** 2 / max((w * w * m).sum(), 1e-9) / max(m.sum(), 1)
    return w, float(neff), float(lost)


# ------------------------------------------------------------------------------------------------------------ fit szablonowy
def smooth_exp2(h, x):
    """Wygładzenie rzadkiego szablonu: dopasowanie exp(a + b·x + c·x²) (Poisson) do histogramu MC."""
    if h.sum() <= 0:
        return h
    t = (x - x.mean()) / (x.max() - x.min())

    def nll(a, b, c):
        lam = np.exp(np.clip(a + b * t + c * t * t, -50, 50))
        return 2 * np.sum(lam - h * np.log(lam))
    m = Minuit(nll, a=math.log(max(h.mean(), 1e-3)), b=0.0, c=0.0)
    m.errordef = 1
    m.limits["b"], m.limits["c"] = (-20, 20), (-20, 20)
    m.migrad()
    lam = np.exp(np.clip(m.values["a"] + m.values["b"] * t + m.values["c"] * t * t, -50, 50))
    return lam * h.sum() / lam.sum()


def shift_scale(tpl, x, x0, d, k):
    """Szablon przesunięty o d i rozciągnięty ×k wokół x0 (interpolacja dystrybuanty)."""
    w = x[1] - x[0]
    ed = np.concatenate([x - w / 2, [x[-1] + w / 2]])
    cdf = np.concatenate([[0], np.cumsum(tpl)])
    src = (ed - x0 - d) / k + x0
    return np.maximum(np.diff(np.interp(src, ed, cdf, left=0, right=cdf[-1])), 0)


def tfit(data, sig, comb, comps, kappa, free=(), fixmu=None, comb_exp=False, sig_free=False, x=None, x0=None):
    """Poissonowski fit binowany: λ = N_S·(s + Σ μ_k κ_k t_k) + N_C·c. comps: {k: szablon znormalizowany}; free: uwolnione μ_k."""
    names = list(comps)
    fixmu = fixmu or {}
    s0 = sig / sig.sum()
    c0 = comb / comb.sum() if comb.sum() > 0 else np.ones_like(sig) / len(sig)
    t = (x - x.mean()) / (x.max() - x.min())
    par = ["NS", "NC"] + [f"mu_{k}" for k in names] + ["b", "c", "d", "k"]

    def lam(p):
        s = shift_scale(s0, x, x0, p[-2], p[-1]) if sig_free else s0
        if comb_exp:
            c = np.exp(np.clip(p[-4] * t + p[-3] * t * t, -50, 50))
            c = c / c.sum()
        else:
            c = c0
        peak = s.copy()
        for i, k in enumerate(names):
            peak = peak + p[2 + i] * kappa[k] * comps[k]
        return p[0] * peak + p[1] * c

    def nll(*p):
        l = np.maximum(lam(p), 1e-9)
        return 2 * np.sum(l - data * np.log(l))
    tot = data.sum()
    init = [max(0.1 * tot, 1.0), max(0.9 * tot, 1.0)] + [fixmu.get(k, 1.0) for k in names] + [0.0, 0.0, 0.0, 1.0]
    m = Minuit(nll, *init, name=par)
    m.errordef = 1
    m.limits["NS"], m.limits["NC"] = (0, 3 * tot + 10), (0, 3 * tot + 10)
    for k in names:
        m.limits[f"mu_{k}"] = (0, 20)
        m.fixed[f"mu_{k}"] = k not in free
    m.limits["b"], m.limits["c"], m.limits["d"], m.limits["k"] = (-30, 30), (-30, 30), (-0.01, 0.01), (0.7, 1.4)
    m.fixed["b"] = m.fixed["c"] = not comb_exp
    m.fixed["d"] = m.fixed["k"] = not sig_free
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m.migrad()
        if not m.valid:
            m.migrad()
        m.hesse()
    p = [m.values[n] for n in par]
    l = lam(p)
    ok = l > 5
    chi2 = float(np.sum((data[ok] - l[ok]) ** 2 / l[ok]))
    ndf = int(ok.sum() - m.nfit)
    out = {"NS": float(m.values["NS"]), "eNS": float(m.errors["NS"]), "NC": float(m.values["NC"]), "chi2": chi2, "ndf": ndf,
           "valid": bool(m.valid and m.accurate), "mu": {k: [float(m.values[f"mu_{k}"]), float(m.errors[f"mu_{k}"])] for k in names if k in free}}
    if sig_free:
        out["d"], out["k"] = [float(m.values["d"]), float(m.errors["d"])], [float(m.values["k"]), float(m.errors["k"])]
    return out


def good(f):
    return f is not None and f["valid"] and f["ndf"] > 0 and f["chi2"] / f["ndf"] < 3 and f["eNS"] > 0 and f["NS"] / f["eNS"] > 3


# ------------------------------------------------------------------------------------------------------------ pomocnicze piki
def peak(h, ed, win, sbl, sbh):
    """Liczba nad liniowym tłem z dwóch pasm: zwraca N, σ, średnia pozycja (ważona nadwyżką)."""
    c = 0.5 * (ed[1:] + ed[:-1])
    w = (c > win[0]) & (c < win[1])
    l, u = (c > sbl[0]) & (c < sbl[1]), (c > sbh[0]) & (c < sbh[1])
    if l.sum() == 0 or u.sum() == 0 or w.sum() == 0:
        return float("nan"), float("nan"), float("nan")
    yl, yu, xl, xu = h[l].mean(), h[u].mean(), c[l].mean(), c[u].mean()
    bg = yl + (yu - yl) * (c[w] - xl) / (xu - xl)
    n = h[w].sum() - bg.sum()
    var = h[w].sum() + (w.sum() ** 2) * 0.25 * (h[l].sum() / l.sum() ** 2 + h[u].sum() / u.sum() ** 2)
    ex = h[w] - bg
    pos = float((ex * c[w]).sum() / ex.sum()) if ex.sum() > 0 else float("nan")
    return float(n), float(math.sqrt(max(var, 0))), pos


def dstar_yield(h2, dm_ed):
    """h2: [dm, m] → (N_sig, σ, szablon masy po odjęciu tła Δm). Tło: (Δm−mπ)^a·exp(−b(Δm−mπ)) dopasowane poza pikiem."""
    c = 0.5 * (dm_ed[1:] + dm_ed[:-1])
    y = h2.sum(1)
    sr = (c > 0.1444) & (c < 0.1464)
    sb = (c > 0.150)
    out = ~((c > 0.1435) & (c < 0.1480))
    q = c - 0.13957
    if y[out].sum() < 20:
        return float("nan"), float("nan"), None

    def nll(A, a, b):
        lam = np.maximum(A * q[out] ** a * np.exp(-b * q[out]), 1e-9)
        return 2 * np.sum(lam - y[out] * np.log(lam))
    m = Minuit(nll, A=max(y[sb].mean(), 1) / (0.015 ** 0.5), a=0.5, b=0.0)
    m.errordef = 1
    m.limits["A"], m.limits["a"], m.limits["b"] = (0, None), (0, 3), (-200, 200)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m.migrad()
    f = m.values["A"] * q ** m.values["a"] * np.exp(-m.values["b"] * q)
    alpha = f[sr].sum() / max(f[sb].sum(), 1e-9)
    n = y[sr].sum() - alpha * y[sb].sum()
    e = math.sqrt(y[sr].sum() + alpha * alpha * y[sb].sum())
    tpl = h2[sr].sum(0) - alpha * h2[sb].sum(0)
    etpl = np.sqrt(h2[sr].sum(0) + alpha * alpha * h2[sb].sum(0))
    return float(n), float(e), (tpl, etpl)


def ratio(a, ea, b, eb):
    if not (np.isfinite(a) and np.isfinite(b)) or b <= 0 or a <= 0:
        return float("nan"), float("nan")
    r = a / b
    return r, r * math.sqrt((ea / a) ** 2 + (eb / b) ** 2)


def k2(r, e):
    """Kryterium K2 → (δ, opis)."""
    if not np.isfinite(r) or not np.isfinite(e) or e > 0.5:
        return 0.5, "niemierzalne (δ = 0,5 a priori)"
    if abs(r - 1) < 2 * e:
        return e, "zgodne"
    return math.sqrt((r - 1) ** 2 + e * e), "NIEZGODNE — szablon wymaga korekty"


# ------------------------------------------------------------------------------------------------------------ główna
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-n", type=float, default=1.0)
    a = ap.parse_args()
    MN = a.min_n
    t0 = time.time()
    os.makedirs(a.out, exist_ok=True)
    D = {s: load(a.inp, s) for s in ("data_428770", "data_427929", "mb_nd", "fd_nd")}
    MB, mmeta = D["mb_nd"]
    FD = D["fd_nd"][0]
    PT, AX, G = mmeta["pt"], mmeta["axes"], mmeta["groups"]
    R = {"input": a.inp, "counts": {s: D[s][1]["counts"] for s in D}, "nfiles": {s: len(D[s][1]["files"]) for s in D}}
    md = []
    P = md.append
    P("# MC_templates — WYNIKI: test szablonów tła z MC (zad. #7; projekt B)\n")
    P(f"Wygenerowane {time.strftime('%Y-%m-%d %H:%M')} przez `mctmpl_analyze.py` z `{a.inp}`. Pliki: "
      + ", ".join(f"{s} {R['nfiles'][s]}" for s in D) + ". Kryteria: `KRYTERIA.md` (zapisane przed wynikami). "
      "Selekcja i klasyfikator H są nasze (nie CharmNet). Dane: GRL, bez zdarzeń regionu sygnału głównego łańcucha (tag AFP ∧ npu = 0 ∧ ξ_cal).\n")
    SPI = {"D0": dict(h="d0_m", m0=mmeta["m0"]["D0"], ax=AX["d0_m"]), "Dplus": dict(h="dp_m", m0=mmeta["m0"]["Dplus"], ax=AX["dp_m"])}

    def mass_h(H, sp):
        """→ [pt, wp, K, cell, m10] (D+: po wecie D*)."""
        h = H[SPI[sp]["h"]]
        if sp == "Dplus":
            h = h[:, :, :, :, 1]
        return reb(h, 2)

    X = {sp: 0.5 * (edges(SPI[sp]["ax"])[::2][1:] + edges(SPI[sp]["ax"])[::2][:-1]) for sp in SPI}
    SBM = {sp: np.abs(X[sp] - SPI[sp]["m0"]) > 0.075 for sp in SPI}
    mbm = {sp: mass_h(MB, sp) for sp in SPI}
    fdm = {sp: mass_h(FD, sp) for sp in SPI}

    def templates(sp, ip, iw, w=None):
        """Szablony [m] per grupa i κ względem Sig. w: wagi komórek (None = bez wag). Sig i Reflection: kształt z mb+fd."""
        def proj(h, g):
            x = h[ip, iw][G[g]].sum(0)
            return (x * (w[:, None] if w is not None else 1.0)).sum(0)
        T = {g: proj(mbm[sp], g) for g in G}
        nsig = T["Sig"].sum()
        kap = {g: (T[g].sum() / nsig if nsig > 0 else 0.0) for g in G}
        ekap = {g: (kap[g] * math.sqrt(1 / max(T[g].sum(), 1e-9) + 1 / max(nsig, 1e-9)) if T[g].sum() > 0 else 0.0) for g in G}
        Ts = dict(T)
        Ts["Sig"] = T["Sig"] + proj(fdm[sp], "Sig")
        if sp == "D0":
            rf, sf = proj(fdm[sp], "Reflection"), proj(fdm[sp], "Sig")
            Ts["Reflection"] = T["Reflection"] + rf
            kap["Reflection"] = Ts["Reflection"].sum() / max(Ts["Sig"].sum(), 1e-9)
            ekap["Reflection"] = kap["Reflection"] / math.sqrt(max(Ts["Reflection"].sum(), 1))
        return Ts, kap, ekap, T

    # ---------------------------------------------------------------- §1 skład i stabilność
    P("## 1. Skład MC wg kategorii głównego łańcucha (TRExFitter) i stabilność szablonów (K1)\n")
    P("Kategorie = próbki TRExFitter z `HistogramManager.cxx:555-556` (Sig 0,1; Reflection = Sig z zamienionymi K/π, tylko D⁰; CharmMismatch 2,3; "
      "DecayModeMismatch 6,7; PartComb 8–11; Comb 4,5,12–15). Bez podziału Sig/SigNonFid (wymaga tagu prawdy AFP — poza zakresem). "
      "Udziały w pełnym oknie (D⁰ 1,70–2,05; D⁺ 1,70–2,20 GeV, po wecie D*), mb_nd bez wag; κ = N/N_Sig.\n")
    R["composition"], R["stability"] = {}, {}
    for sp in SPI:
        P(f"### {sp}\n\n| pT | WP | N_Sig | " + " | ".join(f"κ {g}" for g in G if g != "Sig" and (sp == "D0" or g != "Reflection")) + " | Sig w oknie [%] |\n|---|---|---|"
          + "---|" * (len(G) - (1 if sp == "D0" else 2)) + "---|")
        for ip, p in enumerate(PT):
            for iw in range(4):
                T = {g: mbm[sp][ip, iw][G[g]].sum() for g in G}
                tot = sum(T.values())
                R["composition"][f"{sp}|{p}|{WPN[iw]}"] = T
                P(f"| {p} | {WPN[iw]} | {T['Sig']:.0f} | " + " | ".join(fmt(T[g] / T["Sig"], math.sqrt(T[g]) / T["Sig"]) if T["Sig"] > 0 else "—"
                                                                      for g in G if g != "Sig" and (sp == "D0" or g != "Reflection"))
                  + f" | {100 * T['Sig'] / max(tot, 1):.2f} |")
        P("")
        names = mmeta["keys"]
        sub = [i for g in ("CharmMismatch", "DecayModeMismatch", "PartComb") for i in G[g]]
        P("Podskład (WP95, pT 3–30, zdarzenia w oknie): " + ", ".join(f"{names[i]} {mbm[sp][:, 1, i].sum():.0f}" for i in sub if mbm[sp][:, 1, i].sum() > 0) + ".\n")
    P("### Stabilność przy przeważeniu MC do aktywności klasy trygera (K1)\n")
    P("Wagi w(nSelBL, pT toru wiodącego) = dane/MC w pasmach bocznych (|m − m₀| > 75 MeV, bez cięcia H), per gatunek i bin pT. "
      "Δκ/κ i odległość KS kształtu (szablon przeważony vs nieprzeważony), WP80. n_eff/N < 0,2 → „nieokreślony”.\n")
    P("| próbka danych | gatunek | pT | n_eff/N | dane poza MC [%] | " + " | ".join(f"{g}: Δκ/κ [%] (KS)" for g in ("Reflection", "CharmMismatch", "DecayModeMismatch", "PartComb", "Comb")) + " | werdykt |\n|---|---|---|---|---|---|---|---|---|---|---|")
    W = {}
    for lab, s, cl in DATASETS:
        Hd, dmeta = D[s]
        ic = dmeta["keys"].index(cl)
        for sp in SPI:
            dm_ = mass_h(Hd, sp)
            for ip, p in enumerate(PT):
                w, neff, lost = act_weights(dm_[ip, 0, ic], mbm[sp][ip, 0].sum(0), SBM[sp])
                W[(lab, sp, ip)] = (w, neff, lost)
                _, k0, ek0, T0 = templates(sp, ip, 2)
                _, k1, _, T1 = templates(sp, ip, 2, w)
                cells, verdict = [], "stabilne"
                for g in ("Reflection", "CharmMismatch", "DecayModeMismatch", "PartComb", "Comb"):
                    if (sp == "Dplus" and g == "Reflection") or T0[g].sum() < 10:
                        cells.append("—")
                        continue
                    dk = k1[g] / k0[g] - 1 if k0[g] > 0 else float("nan")
                    c0, c1 = np.cumsum(T0[g]) / T0[g].sum(), np.cumsum(T1[g]) / max(T1[g].sum(), 1e-9)
                    ks = float(np.abs(c0 - c1).max())
                    lim_k = max(0.20, 2 * ek0[g] / k0[g]) if k0[g] > 0 else 0.2
                    lim_s = max(0.05, 1.36 / math.sqrt(max(neff * T0[g].sum(), 1)))
                    bad = abs(dk) > lim_k or ks > lim_s
                    if bad:
                        verdict = "NIESTABILNE: " + g if verdict == "stabilne" else verdict + ", " + g
                    cells.append(f"{100 * dk:+.1f} ({ks:.3f}){'!' if bad else ''}")
                    R["stability"][f"{lab}|{sp}|{p}|{g}"] = {"dkappa": dk, "ks": ks, "neff": neff, "bad": bool(bad)}
                if neff < 0.2:
                    verdict = "nieokreślony (n_eff/N < 0,2)" + ("; " + verdict if verdict != "stabilne" else "")
                P(f"| {lab} | {sp} | {p} | {neff:.2f} | {100 * lost:.1f} | " + " | ".join(cells) + f" | {verdict} |")
    P("")

    # ---------------------------------------------------------------- §2a refleksja z tagu D*
    P("## 2. Testy w obszarach wzbogaconych (K2)\n")
    P("### 2a. Refleksja K↔π D⁰ z tagu D* (bliźniak o zamienionej hipotezie)\n")
    P("Wiersze D* (znak pionu miękkiego ustala hipotezę) dopasowane po torach do wierszy D⁰: R = bliźniak (zamieniona hipoteza) przechodzi WP, "
      "S = poprawna hipoteza przechodzi WP; tło pod pikiem Δm odjęte (funkcja progowa dopasowana poza pikiem). MC: ta sama procedura na mb_nd + fd_nd; "
      "„prawda” = N(Reflection)/N(Sig) z kategorii. r = (R/S)_dane / (R/S)_MC.\n")
    dm_ed, m10 = edges(AX["dsdm"]), edges(AX["m10"])
    mc10 = 0.5 * (m10[1:] + m10[:-1])
    swin = np.abs(mc10 - SPI["D0"]["m0"]) < 0.045

    def rs(Hr, Ht, key, ip, iw):
        ips = range(3) if ip is None else [ip]
        hr = sum(Hr[i, iw, key] for i in ips)[:, swin]
        ht = sum(Ht[i, iw, key] for i in ips)
        ns, es, _ = dstar_yield(hr, dm_ed)
        nr, er, tp = dstar_yield(ht, dm_ed)
        return ns, es, nr, er, tp
    mcr, mct = MB["ds_right"].sum(2, keepdims=True) + FD["ds_right"].sum(2, keepdims=True), MB["ds_twin"].sum(2, keepdims=True) + FD["ds_twin"].sum(2, keepdims=True)
    R["refl"] = {}
    REFL_TPL = {}
    P("| próbka danych | pT | WP | N_S(D*) dane | (R/S) dane | (R/S) MC tag | (R/S) MC prawda | r = dane/MC | K2 | χ²/ndf kształtu (p) |\n|---|---|---|---|---|---|---|---|---|---|")
    for lab, s, cl in DATASETS:
        Hd, dmeta = D[s]
        ic = dmeta["keys_ds"].index(cl)
        for ip, p in list(enumerate(PT)) + [(None, "3-30")]:
            for iw in FITWP:
                ns, es, nr, er, tp = rs(Hd["ds_right"], Hd["ds_twin"], ic, ip, iw)
                ms, ems, mr, emr, mtp = rs(mcr, mct, 0, ip, iw)
                rd, erd = ratio(nr, er, ns, es)
                rm, erm = ratio(mr, emr, ms, ems)
                ips = range(3) if ip is None else [ip]
                tr = sum(mbm["D0"][i, iw][G["Reflection"]].sum() + fdm["D0"][i, iw][G["Reflection"]].sum() for i in ips)
                ts = sum(mbm["D0"][i, iw][G["Sig"]].sum() + fdm["D0"][i, iw][G["Sig"]].sum() for i in ips)
                r, er_ = ratio(rd, erd, rm, erm)
                dlt, txt = k2(r, er_)
                shp = "—"
                if tp is not None and nr > 30 * MN:
                    ref = sum(mbm["D0"][i, iw][G["Reflection"]].sum((0, 1)) + fdm["D0"][i, iw][G["Reflection"]].sum((0, 1)) for i in ips)
                    ok = tp[1] > 0
                    exp = ref / max(ref.sum(), 1e-9) * tp[0].sum()
                    c2 = float(np.sum((tp[0][ok] - exp[ok]) ** 2 / (tp[1][ok] ** 2 + (exp[ok] ** 2 / np.maximum(ref[ok], 1)))))
                    nd = int(ok.sum() - 1)
                    shp = f"{c2:.1f}/{nd} ({stats.chi2.sf(c2, nd):.3f})"
                    REFL_TPL[(lab, ip, iw)] = (np.maximum(tp[0], 0), nr, rd)
                R["refl"][f"{lab}|{p}|{WPN[iw]}"] = {"NS": ns, "RS_data": [rd, erd], "RS_mc_tag": [rm, erm], "RS_mc_truth": tr / max(ts, 1), "r": [r, er_],
                                                    "delta": dlt, "k2": txt, "shape": shp}
                P(f"| {lab} | {p} | {WPN[iw]} | {fmt(ns, es, 0)} | {fmt(rd, erd)} | {fmt(rm, erm)} | {tr / max(ts, 1):.3f} | {fmt(r, er_)} | {txt} | {shp} |")
    P("")

    # ---------------------------------------------------------------- §3 fity (potrzebne też do 2b, 2c)
    def build(sp, ip, iw, w=None):
        Ts, kap, ekap, T = templates(sp, ip, iw, w)
        x = X[sp]
        comb = Ts["Comb"].copy()
        comps = {}
        for g in COMP[sp]:
            t = Ts[g]
            if t.sum() < 20:          # za mało MC na osobny szablon → do Comb (zdarzeń nie gubimy)
                comb = comb + t
                continue
            comps[g] = t / t.sum()
        if comb.sum() < 300:          # tylko bardzo rzadki Comb wygładzam; PartComb NIE (ma część pikującą)
            comb = smooth_exp2(comb, x)
        return Ts["Sig"], comb, comps, {g: kap[g] for g in comps}

    def run_configs(data, sp, ip, iw, lab=None, w=None):
        sig, comb, comps, kap = build(sp, ip, iw)
        x, x0 = X[sp], SPI[sp]["m0"]
        kw = dict(x=x, x0=x0)
        out = {"C0": tfit(data, sig, comb, comps, kap, **kw)}
        pr = [g for g in ("CharmMismatch", "DecayModeMismatch") if g in comps]
        if sp == "D0":
            if "Reflection" in comps:
                out["C1"] = tfit(data, sig, comb, comps, kap, free=("Reflection",), **kw)
            if pr:   # wspólny czynnik dla rozpadów częściowo zrekonstruowanych / innych modów
                c2 = {k: v for k, v in comps.items() if k not in pr}
                kk = {k: v for k, v in kap.items() if k not in pr}
                tpr = sum(kap[g] * comps[g] for g in pr)
                c2["PR"], kk["PR"] = tpr / tpr.sum(), float(tpr.sum())
                out["C2"] = tfit(data, sig, comb, c2, kk, free=("PR",), **kw)
                for v in (0.5, 1.5):
                    out[f"S_PR_{v}"] = tfit(data, sig, comb, c2, kk, fixmu={"PR": v}, **kw)
        else:
            if "CharmMismatch" in comps:
                out["C1"] = tfit(data, sig, comb, comps, kap, free=("CharmMismatch",), **kw)
            if "DecayModeMismatch" in comps:
                out["C2"] = tfit(data, sig, comb, comps, kap, free=("DecayModeMismatch",), **kw)
        if "PartComb" in comps:
            out["C3"] = tfit(data, sig, comb, comps, kap, free=("PartComb",), **kw)
        out["C5"] = tfit(data, sig, comb, comps, kap, free=tuple(comps), **kw)
        out["C6"] = tfit(data, sig, comb, comps, kap, comb_exp=True, **kw)
        out["C7"] = tfit(data, sig, comb, comps, kap, sig_free=True, **kw)
        for g in comps:
            if sp == "D0" and g in pr:
                continue
            for v in (0.5, 1.5):
                out[f"S_{g}_{v}"] = tfit(data, sig, comb, comps, kap, fixmu={g: v}, **kw)
        if w is not None:
            s2, c2_, cp2, k2_ = build(sp, ip, iw, w)
            out["C0w"] = tfit(data, s2, c2_, cp2, k2_, **kw)
        if sp == "D0" and lab is not None and "Reflection" in comps:
            for key in ((lab, ip, iw), (lab, None, iw), (lab.split(":")[0] + ":all", None, iw)):
                if key in REFL_TPL and REFL_TPL[key][1] > 50 * MN:
                    tpl, _, rsd = REFL_TPL[key]
                    cp, kp = dict(comps), dict(kap)
                    cp["Reflection"], kp["Reflection"] = tpl / tpl.sum(), rsd
                    out["C4"] = tfit(data, sig, comb, cp, kp, **kw)
                    out["C4"]["src"] = f"{key[0]}, pT {'3-30' if key[1] is None else PT[key[1]]}"
                    break
        return out

    P("## 3. Fity masy: konfiguracja „jak w głównym łańcuchu” vs warianty (K3, K4)\n")
    P("Model binowany (10 MeV, okna jak w ficie TRExFitter głównego łańcucha): λ = N_S·(s + Σ_k μ_k κ_k t_k) + N_C·c. **C0** = „jak w głównym łańcuchu”: N_S i Comb swobodne, kształty wszystkich "
      "składowych i κ_k z MC, μ_k = 1 (w głównym łańcuchu normalizacja absolutna σ·L; tu względem sygnału, bo nie mamy L_eff per tryger — równoważne tylko przy μ_sygnału = 1). "
      "C0w = C0 z MC przeważonym do aktywności klasy; C1 = wolna μ refleksji (D⁰) / CharmMismatch (D⁺); C2 = wolna μ rozpadów częściowo zrekonstruowanych "
      "(D⁰: CharmMismatch + DecayModeMismatch wspólnie; D⁺: DecayModeMismatch); C3 = wolna μ PartComb; C4 = refleksja (kształt i R/S) z danych z tagu D*; "
      "C5 = wszystkie μ wolne; C6 = Comb jako exp(pol2) zamiast kształtu MC; C7 = C0 + swobodne przesunięcie i szerokość sygnału. "
      "Szablony = surowe histogramy MC (jak w głównym łańcuchu), bez niepewności statystycznej MC w likelihoodzie (uproszczenie; w głównym łańcuchu gammy); składowa z < 20 zdarzeń MC "
      "jest dołączana do Comb, a Comb z < 300 zdarzeń MC wygładzany exp(pol2). PartComb nie jest wygładzany (ma część pikującą).\n")
    # K4: domknięcie
    clos, nbad, ncl = [], 0, 0
    for sp in SPI:
        for ip, p in enumerate(PT):
            for iw in FITWP:
                pseudo = mbm[sp][ip, iw].sum((0, 1))
                ntrue = mbm[sp][ip, iw][G["Sig"]].sum()
                if ntrue < 100 * MN:
                    continue
                sig, comb, comps, kap = build(sp, ip, iw)
                f = tfit(pseudo, sig, comb, comps, kap, x=X[sp], x0=SPI[sp]["m0"])
                pull = (f["NS"] - ntrue) / max(f["eNS"], 1e-9)
                ncl += 1
                nbad += abs(pull) > 1
                clos.append(f"{sp} {p} WP{WPN[iw]}: {f['NS']:.0f} ± {f['eNS']:.0f} vs {ntrue:.0f} (pull {pull:+.2f})")
    R["closure"] = {"n": ncl, "n_bad": int(nbad), "rows": clos}
    P(f"**K4 (domknięcie na mb_nd):** {ncl - nbad}/{ncl} komórek w 1σ" + (" — OK." if nbad == 0 else " — **NIE DOMYKA SIĘ, wyniki poniżej wstrzymane do wyjaśnienia.**")
      + " Szczegóły w results.json (`closure`).\n")
    R["fits"], R["syst"] = {}, {}
    hdr = "| próbka | pT | WP | N_S ± σ (C0) | χ²/ndf C0 | C0w | C1 | C2 | C3 | C4 | C5 | C6 | C7 | μ̂ C1 | μ̂ C2 | μ̂ C3 | σ(C1)/σ(C0) | σ(C2)/σ(C0) |\n|" + "---|" * 18
    FITS = {}
    for sp in SPI:
        P(f"### {sp}: ΔN_S/N_S [%] względem C0\n\n" + hdr)
        for lab, s, cl in DATASETS:
            Hd, dmeta = D[s]
            ic = dmeta["keys"].index(cl)
            dm_ = mass_h(Hd, sp)
            for ip, p in enumerate(PT):
                for iw in FITWP:
                    data = dm_[ip, iw, ic].sum(0)
                    if data.sum() < 200 * MN:
                        continue
                    F = run_configs(data, sp, ip, iw, lab, W[(lab, sp, ip)][0])
                    FITS[(lab, sp, ip, iw)] = F
                    R["fits"][f"{lab}|{sp}|{p}|{WPN[iw]}"] = F
                    f0 = F["C0"]
                    if not good(f0):
                        P(f"| {lab} | {p} | {WPN[iw]} | {fmt(f0['NS'], f0['eNS'], 0)} | {f0['chi2']:.0f}/{f0['ndf']} | nieokr. |" + " |" * 13)
                        continue

                    def d(c):
                        return f"{100 * (F[c]['NS'] / f0['NS'] - 1):+.1f}" + ("" if F[c]["valid"] else "?") if c in F else "—"

                    def mu(c):
                        if c not in F or not F[c]["mu"]:
                            return "—"
                        v = next(iter(F[c]["mu"].values()))
                        return fmt(v[0], v[1], 2)

                    def sg(c):
                        return f"{F[c]['eNS'] / f0['eNS']:.2f}" if c in F else "—"
                    P(f"| {lab} | {p} | {WPN[iw]} | {fmt(f0['NS'], f0['eNS'], 0)} | {f0['chi2']:.0f}/{f0['ndf']} | " + " | ".join(d(c) for c in ("C0w", "C1", "C2", "C3", "C4", "C5", "C6", "C7"))
                      + f" | {mu('C1')} | {mu('C2')} | {mu('C3')} | {sg('C1')} | {sg('C2')} |")
        P("")

    # ---------------------------------------------------------------- §2b
    P("### 2b. Rozpady częściowo zrekonstruowane / inne mody (dolne pasmo masy) — μ̂ z fitu C2\n")
    P("μ̂ = normalizacja danych względem MC (przy κ z MC). D⁰: CharmMismatch + DecayModeMismatch (Kππ⁰, KK, ππ …); D⁺: DecayModeMismatch. "
      "To ten sam fit co w §3 — nie jest to obszar niezależny.\n\n| próbka | gatunek | pT | WP | μ̂ ± σ | K2 |\n|---|---|---|---|---|---|")
    for (lab, sp, ip, iw), F in FITS.items():
        if "C2" in F and good(F["C0"]) and F["C2"]["mu"]:
            v = next(iter(F["C2"]["mu"].values()))
            P(f"| {lab} | {sp} | {PT[ip]} | {WPN[iw]} | {fmt(v[0], v[1], 2)} | {k2(v[0], v[1])[1]} |")
    P("")

    # ---------------------------------------------------------------- §2c D+ piki
    P("### 2c. D⁺: piki refleksji w danych vs MC (D*⁺→D⁰π⁺, φ→KK z D_s⁺/D⁺→KKπ, Λ_c→pKπ)\n")
    P("Liczby nad liniowym tłem, znormalizowane do N_S(D⁺) (dane: fit C0 po wecie D*; MC: prawda). D*: pik w Δm_min = m(Kππ) − m(Kπ) (przed wetem); "
      "φ: pik w m(KK) po przypisaniu masy K pionowi; D_s i D⁺→KKπ: piki w m(KKπ) w oknie φ; Λ_c: pik w m(pKπ). pT 3–30 łącznie. r = dane/MC.\n")
    P("| próbka | WP | składowa | N/N_S dane | N/N_S MC | r = dane/MC | pozycja dane [MeV] | pozycja MC [MeV] | K2 |\n|---|---|---|---|---|---|---|---|---|")
    PK = [("D* (Δm)", "dp_dm", None, edges(AX["dpdm"]), (0.1440, 0.1470), (0.1395, 0.1425), (0.1490, 0.1530)),
          ("φ (m_KK)", "dp_mkk", None, edges(AX["mkk"]), (1.0095, 1.0295), (0.990, 1.004), (1.036, 1.056)),
          ("D_s→φπ (m_KKπ)", "dp_kkpi", 0, edges(AX["kkpi"]), (1.950, 1.985), (1.915, 1.940), (1.995, 2.030)),
          ("D⁺→φπ (m_KKπ)", "dp_kkpi", 0, edges(AX["kkpi"]), (1.850, 1.890), (1.810, 1.840), (1.900, 1.930)),
          ("Λ_c (m_pKπ)", "dp_pkpi", None, edges(AX["pkpi"]), (2.270, 2.300), (2.235, 2.260), (2.310, 2.340))]
    R["dplus_peaks"] = {}
    DELTA_CM = {}
    for lab, s, cl in DATASETS:
        Hd, dmeta = D[s]
        ic = dmeta["keys"].index(cl)
        for iw in (1, 2):
            fs = [FITS.get((lab, "Dplus", ip, iw)) for ip in range(3)]
            fs = [f["C0"] for f in fs if f and good(f["C0"])]
            if not fs:
                continue
            nsd, ensd = sum(f["NS"] for f in fs), math.sqrt(sum(f["eNS"] ** 2 for f in fs))
            ipg = [ip for ip in range(3) if FITS.get((lab, "Dplus", ip, iw)) and good(FITS[(lab, "Dplus", ip, iw)]["C0"])]
            nsm = sum(mbm["Dplus"][ip, iw][G["Sig"]].sum() for ip in ipg)
            best = []
            for nm, hn, sub, ed, win, sbl, sbh in PK:
                hd = sum(Hd[hn][ip, iw, ic] for ip in ipg)
                hm = sum(MB[hn][ip, iw].sum(0) for ip in ipg)
                if sub is not None:
                    hd, hm = hd[sub], hm[sub]
                nd_, ed_, pd_ = peak(hd, ed, win, sbl, sbh)
                nm_, em_, pm_ = peak(hm, ed, win, sbl, sbh)
                a_, ea_ = ratio(nd_, ed_, nsd, ensd)
                b_, eb_ = ratio(nm_, em_, nsm, math.sqrt(max(nsm, 1)))
                r, er_ = ratio(a_, ea_, b_, eb_)
                dlt, txt = k2(r, er_)
                best.append((dlt, txt))
                R["dplus_peaks"][f"{lab}|{WPN[iw]}|{nm}"] = {"data": [a_, ea_], "mc": [b_, eb_], "r": [r, er_], "pos": [pd_, pm_], "k2": txt}
                P(f"| {lab} | {WPN[iw]} | {nm} | {fmt(a_, ea_, 4)} | {fmt(b_, eb_, 4)} | {fmt(r, er_, 2)} | {fmt(1e3 * pd_, None, 1)} | {fmt(1e3 * pm_, None, 1)} | {txt} |")
            meas = [d_ for d_, t_ in best if not t_.startswith("niemierz")]
            DELTA_CM[(lab, iw)] = max(meas) if meas else 0.5
    P("")

    # ---------------------------------------------------------------- §2d prompt / non-prompt
    P("### 2d. Składowa z b: |parametr zderzenia D| w oknie sygnału (po odjęciu pasm bocznych)\n")
    P("f_B = udział szablonu MC SigB (MatchedB) w dopasowaniu χ² rozkładu |impact| (0–0,2 mm, 40 binów; ostatni bin = przepełnienie) danych "
      "SR − (2/3)·(SB dolne + górne) szablonami SigC(+Reflection) i SigB z mb_nd + fd_nd. f_B MC = prawda mb_nd. WP80. Uwaga: H używa cech wierzchołkowych, "
      "więc cięcie rzeźbi ten rozkład jednakowo w danych i MC tylko o tyle, o ile MC je modeluje.\n")
    P("| próbka | gatunek | pT | f_B dane | f_B MC (prawda mb_nd) | χ²/ndf | K2 (różnica w σ) |\n|---|---|---|---|---|---|---|")
    R["fB"] = {}
    for lab, s, cl in DATASETS:
        Hd, dmeta = D[s]
        ic = dmeta["keys"].index(cl)
        for sp, hn in (("D0", "d0_imp"), ("Dplus", "dp_imp")):
            for ip, p in enumerate(PT):
                h = Hd[hn][ip, 2, ic]
                y = h[0] - (2 / 3) * (h[1] + h[2])
                ey = np.sqrt(h[0] + (4 / 9) * (h[1] + h[2]))
                mc = MB[hn][ip, 2] + FD[hn][ip, 2]
                tc = mc[[0, 2, 3], 0].sum(0) if sp == "D0" else mc[0, 0]
                tb = mc[1, 0]
                if y.sum() < 200 * MN or tb.sum() < 50 * MN or tc.sum() < 50 * MN:
                    continue
                tc, tb = tc / tc.sum(), tb / tb.sum()
                ok = ey > 0

                def c2f(N, f):
                    return float(np.sum((y[ok] - N * ((1 - f) * tc[ok] + f * tb[ok])) ** 2 / ey[ok] ** 2))
                m = Minuit(c2f, N=y.sum(), f=0.1)
                m.errordef = 1
                m.limits["f"] = (-0.5, 1.0)
                m.migrad()
                m.hesse()
                nb_, nc_ = MB[hn][ip, 2][1, 0].sum(), MB[hn][ip, 2][0, 0].sum()
                fmc = nb_ / max(nb_ + nc_, 1)
                efmc = math.sqrt(fmc * (1 - fmc) / max(nb_ + nc_, 1))
                z = (m.values["f"] - fmc) / math.hypot(m.errors["f"], efmc)
                R["fB"][f"{lab}|{sp}|{p}"] = {"f_data": [m.values["f"], m.errors["f"]], "f_mc": [fmc, efmc], "chi2": m.fval, "ndf": int(ok.sum() - 2), "z": z}
                P(f"| {lab} | {sp} | {p} | {fmt(m.values['f'], m.errors['f'])} | {fmt(fmc, efmc)} | {m.fval:.0f}/{int(ok.sum() - 2)} | {'zgodne' if abs(z) < 2 else 'NIEZGODNE'} ({z:+.1f}) |")
    P("")

    # ---------------------------------------------------------------- systematyka
    P("## 4. Proponowana systematyka na składowe z MC (w % N_S; K3)\n")
    P("syst_k = |S_k|·δ_k, S_k = [N_S(μ_k = 1,5) − N_S(μ_k = 0,5)]/N_S(C0). δ: Reflection — tag D* (§2a, ta sama próbka; przy „niemierzalne” pT 3–30); "
      "CharmMismatch D⁺ — największe δ z pików §2c; PR (D⁰) i DecayModeMismatch (D⁺) — √((μ̂−1)² + σ²) z C2 (**ten sam fit**); PartComb — założenie 25 %. "
      "„kształt R” = |N_S(C4) − N_S(C0)|/N_S. Pogrubione: > σ_stat/3 („istotna”).\n")
    P("| próbka | gatunek | pT | WP | σ_stat [%] | " + " | ".join(["Reflection", "kształt R", "CharmMismatch", "DMM / PR", "PartComb", "**suma**", "C0w−C0 (tryger)", "C6−C0 (kształt Comb)"]) + " |\n|" + "---|" * 13)
    for (lab, sp, ip, iw), F in FITS.items():
        f0 = F["C0"]
        if not good(f0):
            continue
        st = 100 * f0["eNS"] / f0["NS"]

        def S(g):
            a_, b_ = F.get(f"S_{g}_1.5"), F.get(f"S_{g}_0.5")
            return (a_["NS"] - b_["NS"]) / f0["NS"] if a_ and b_ else None
        row, tot = {}, 0.0
        if sp == "D0":
            rr = R["refl"].get(f"{lab}|{PT[ip]}|{WPN[iw]}")
            if rr is None or rr["k2"].startswith("niemierz"):
                rr = R["refl"].get(f"{lab}|3-30|{WPN[iw]}", rr)
            if S("Reflection") is not None and rr:
                row["Reflection"] = 100 * abs(S("Reflection")) * rr["delta"]
            if "C4" in F:
                row["kształt R"] = 100 * abs(F["C4"]["NS"] / f0["NS"] - 1)
            if S("PR") is not None and "C2" in F and F["C2"]["mu"]:
                v = F["C2"]["mu"]["PR"]
                row["DMM / PR"] = 100 * abs(S("PR")) * k2(v[0], v[1])[0]
        else:
            if S("CharmMismatch") is not None:
                row["CharmMismatch"] = 100 * abs(S("CharmMismatch")) * DELTA_CM.get((lab, min(iw, 2)), 0.5)
            if S("DecayModeMismatch") is not None and "C2" in F and F["C2"]["mu"]:
                v = F["C2"]["mu"]["DecayModeMismatch"]
                row["DMM / PR"] = 100 * abs(S("DecayModeMismatch")) * k2(v[0], v[1])[0]
        if S("PartComb") is not None:
            row["PartComb"] = 100 * abs(S("PartComb")) * PRIOR_PC
        tot = math.sqrt(sum(v * v for v in row.values()))
        R["syst"][f"{lab}|{sp}|{PT[ip]}|{WPN[iw]}"] = {"stat": st, "comp": row, "total": tot,
                                                       "sens": {g: S(g) for g in ("Reflection", "PR", "CharmMismatch", "DecayModeMismatch", "PartComb") if S(g) is not None}}

        def c(k):
            return "—" if k not in row else (f"**{row[k]:.2f}**" if row[k] > st / 3 else f"{row[k]:.2f}")

        def dd(k):
            return f"{100 * (F[k]['NS'] / f0['NS'] - 1):+.1f}" if k in F else "—"
        P(f"| {lab} | {sp} | {PT[ip]} | {WPN[iw]} | {st:.1f} | " + " | ".join(c(k) for k in ("Reflection", "kształt R", "CharmMismatch", "DMM / PR", "PartComb"))
          + f" | **{tot:.2f}** | {dd('C0w')} | {dd('C6')} |")
    P("")
    P("## 5. Co z tego wynika dla głównego łańcucha analizy\n")
    P("- Przenosi się (zależy od fizyki i rekonstrukcji, nie od klasyfikatora): definicje kategorii, R/S refleksji i jej kształt przy danym L_xy, "
      "pozycje i względne wielkości pików D*/φ/Λ_c przed cięciem klasyfikatora, f_B, wrażliwość S_k przy podobnym S/B.\n"
      "- Nie przenosi się 1:1: wartości κ po cięciu (CharmNet ma inne wejścia, m.in. ptcone i mKπ/m_D dla D⁺, ciaśniejsze WP), L_xy ≥ 0,25 mm dla D⁰ (w głównym łańcuchu 0,10), "
      "nasze weto D* dla D⁺ (w głównym łańcuchu brak weta w preselekcji), bin 2–3 GeV (nie liczony), region z tagiem AFP (nie oglądany — skład tła może być tam inny), "
      "normalizacja absolutna σ·L (tu względna do sygnału).\n- Hipotezy, nie wyniki: przyczyny ewentualnych niezgodności (modelowanie zmiennych wierzchołkowych, tryger).\n")
    open(os.path.join(a.out, "WYNIKI.md"), "w").write("\n".join(md) + "\n")
    json.dump(R, open(os.path.join(a.out, "results.json"), "w"), indent=1, default=lambda o: float(o) if isinstance(o, (np.floating, np.integer)) else str(o))
    print(f"zapisano {a.out}/WYNIKI.md, results.json ({time.time() - t0:.0f} s); K4: {ncl - nbad}/{ncl}")


if __name__ == "__main__":
    main()
