#!/usr/bin/env python3
"""BKGCAT (1): skład tła mb_nd per kategoria prawdy dla drzewa H (WP90/70/50, kategoria luźna [70, 50)) w oknie i pasmach + „sygnał pozorny” każdej składowej:
fit (szablon sygnału fd o kształcie stałym + exp(Czebyszew n), jak M2 z FITBIAS, bez składowych związanych) do histogramu SAMEJ składowej → N_spur/N_prawda ± HESSE.
POST HOC względem FITBIAS. Wyjście: out/comp.json, tabele/sklad.md, tabele/spur.md, tabele/ksztalt.md.   bc_comp.py [--in bkgcat.npz] [--nproc 6]"""
import argparse
import json
import os
import re
import sys
import zlib
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bc_lib as L  # noqa: E402
F = L.F

CATS = {"WP90": (0.9, None), "WP70": (0.7, None), "WP50": (0.5, None), "B[70,50)": (0.7, 0.5), "B[90,70)": (0.9, 0.7), "B[95,90)": (0.95, 0.9)}
MAIN = ["reszta (wszystko − Sig − Refl)", "PC+Comb", "PC", "Comb", "NoCharm", "DMM+CM", "DMM", "CM"]
NBOOT = 60
B = Z = None


def _init(fn):
    global B, Z
    B = L.BK(fn)
    Z = np.load(L.WPF, allow_pickle=False)


def sub(key):
    """Podskładowe: klucz → etykieta."""
    g = L.group(key)
    p = key.split("|")
    if g == "PC":
        return f"PC|{p[1]}"
    if g == "CM":
        return f"CM|{p[1]}"
    if g == "DMM":
        return "DMM|" + p[1]
    if g == "ML":
        return key
    return g


def comps(sp, bn, elo, ehi):
    c = B.cat(sp, bn, elo, ehi)
    o = {}
    add = lambda k, h: o.__setitem__(k, o.get(k, 0) + h)  # noqa: E731
    for k, h in c.items():
        g = L.group(k)
        add(g, h); add("all", h)
        if g not in ("Sig", "Refl"):
            add("reszta (wszystko − Sig − Refl)", h)
        if g in ("PC", "MB", "ML", "OtherB", "NoCharm"):
            add("PC+Comb", h)
        if g in ("MB", "ML", "OtherB", "NoCharm"):
            add("Comb", h)
        if g in ("DMM", "CM"):
            add("DMM+CM", h)
        if g in ("PC", "CM", "DMM", "ML"):
            add(sub(k), h)
        if g == "PC":
            add("PC|" + k.split("|")[2], h)
    return o


def task(t):
    sp, bn, cn = t
    c = F.SP[sp]
    elo, ehi = CATS[cn]
    jlo, jhi = B.jj(elo, ehi)
    key = f"{sp}|{L.TAG[sp]}|{bn}"
    sw = Z[key + "|fd|sig_w"][jlo + 1:(None if jhi is None else jhi + 1)].sum(0)
    sigT = F.tmpl1(sw, c["mr"][0])
    C = comps(sp, bn, elo, ehi)
    ntrue = float(C["Sig"].sum())
    res = {"sp": sp, "bin": bn, "cat": cn, "N_true": ntrue, "comp": {}}
    # kontrola z WPSCAN
    wa = Z[key + "|mb|all"][jlo + 1:(None if jhi is None else jhi + 1)].sum(0)
    ws = Z[key + "|mb|sig"][jlo + 1:(None if jhi is None else jhi + 1)].sum(0)
    res["check"] = [float(np.abs(C["all"] - wa).sum()), float(np.abs(C["Sig"] - ws).sum()), float(wa.sum())]
    for name, h1 in C.items():
        hh = L.h5(h1)
        o = {"n_full": float(hh.sum())}
        for rn in ("fit", "alt"):
            rng = c[rn]
            wm, edges = L.win_mask(sp, rng)
            hr = L.cutr(hh, sp, rng)
            o[rn] = {"n_win": float(hr[wm].sum()), "n_sb": float(hr[~wm].sum())}
            if name in ("all", "Sig", "Refl") or hr.sum() < 30:
                continue
            # Sygnał pozorny: do składowej dodany sygnał Asimova (N_prawda × szablon fd), N_spur = N_fit − N_wstrzyknięte.
            # (Fit samej składowej opiera się o dolną granicę N w fb_fit: −5√N_okno(±30 MeV) − 10 — przy D0 3–5 GeV granica była osiągana.)
            sA = ntrue * sigT.bins(edges, c["m0"])
            ninj = float(sA.sum())
            f = F.fit((hr + sA)[None], [sigT], [[]], edges, c["m0"], c["nbkg"], shape="fix", win=c["win"])
            o[rn].update({"N": f["N"][0] - ninj, "err_hesse": (f.get("err") or [None])[0], "valid": f["valid"], "chi2": f["chi2"], "ndf": f["ndf"], "err": None})
            full = cn in ("WP70", "B[70,50)", "WP50", "WP90") or name in MAIN
            if full:
                rng_ = np.random.default_rng(zlib.crc32("|".join((sp, bn, cn, rn, name)).encode()))
                bN = []
                for _ in range(NBOOT if name in MAIN else NBOOT // 2):
                    fbt = F.fit((rng_.poisson(hr) + sA)[None], [sigT], [[]], edges, c["m0"], c["nbkg"], shape="fix", win=c["win"], start=f["x"], hesse=False)
                    if fbt["valid"]:
                        bN.append(fbt["N"][0])
                if len(bN) >= 0.7 * (NBOOT if name in MAIN else NBOOT // 2):
                    o[rn]["err"] = float(np.std(bN)); o[rn]["valid"] = True; o[rn]["boot"] = True; o[rn]["nboot"] = len(bN)
            if name in MAIN and rn == "fit":
                f1 = F.fit((hr + sA)[None], [sigT], [[]], edges, c["m0"], c["nbkg"] + 1, shape="fix", win=c["win"])
                o["fit_p1"] = {"N": f1["N"][0] - ninj, "err": o[rn]["err"], "valid": f1["valid"]}
            # przewidywanie okna z samych pasm (okno zaślepione)
            fb = F.fit(hr[None], [sigT], [[]], edges, c["m0"], c["nbkg"], shape="fix", win=c["win"], blind=wm, hesse=False)
            o[rn]["pred_win_sb"] = float(fb["mu"][0][wm].sum())
        o["h25"] = hh.reshape(-1, 5).sum(1).tolist()                 # kształt: biny 25 MeV, pełny zakres
        res["comp"][name] = o
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", default=L.IN); ap.add_argument("--nproc", type=int, default=6)
    a = ap.parse_args()
    tasks = [(sp, bn, cn) for sp in L.SPS for bn in L.BINS for cn in CATS]
    with Pool(a.nproc, initializer=_init, initargs=(a.inp,)) as p:
        R = p.map(task, tasks)
    os.makedirs(os.path.join(L.BASE, "out"), exist_ok=True)
    json.dump(R, open(os.path.join(L.BASE, "out", "comp.json"), "w"))
    D = {(r["sp"], r["bin"], r["cat"]): r for r in R}
    bad = [(k, r["check"]) for k, r in D.items() if r["check"][0] or r["check"][1]]
    print("kontrola sum vs WPSCAN (mb|all, mb|sig):", "ZGODNE co do zdarzenia we wszystkich komórkach" if not bad else bad)
    T = os.path.join(L.BASE, "tabele")
    # ---- skład
    out = ["# BKGCAT — skład mb_nd per kategoria prawdy, drzewo H (POST HOC względem FITBIAS)", "",
           "Zakres nominalny fitu (D⁰ 1,75–2,00; D⁺ 1,70–2,20 GeV); okno |m − m₀| < 55 / 45 MeV; pasma = zakres − okno. Liczby kandydatów mb_nd (78 plików, wagi ≡ 1).",
           "„% tła” = udział w tle (wszystko − Sig) w oknie; „/Sig” = względem sygnału w oknie.", ""]
    for sp in L.SPS:
        for cn in ("WP70", "WP90", "WP50", "B[70,50)"):
            out += [f"## {sp} — {cn}", "", "| składowa | " + " | ".join(f"{b}: okno | pasma | % tła | /Sig" for b in L.BINS) + " |", "|---|" + "---|---|---|---|" * 3]
            names = L.GORD + ["Comb", "PC+Comb"]
            for g in names:
                row = []
                for bn in L.BINS:
                    r = D[(sp, bn, cn)]["comp"]
                    x = r.get(g, {"fit": {"n_win": 0, "n_sb": 0}})["fit"]
                    bw = r["all"]["fit"]["n_win"] - r["Sig"]["fit"]["n_win"]
                    sgw = r["Sig"]["fit"]["n_win"]
                    row.append(f"{x['n_win']:.0f} | {x['n_sb']:.0f} | {('%.1f' % (100 * x['n_win'] / bw)) if g != 'Sig' and bw else '—'} | {x['n_win'] / sgw:.3f}")
                out.append(f"| {L.GNAME.get(g, g)} | " + " | ".join(row) + " |")
            out.append("")
    out += ["## Podskładowe przy WP70 (okno / pasma; ≥ 5 kandydatów w zakresie w którymś binie)", ""]
    for sp in L.SPS:
        out += [f"### {sp}", "", "| podskładowa | " + " | ".join(f"{b}: okno | pasma" for b in L.BINS) + " |", "|---|" + "---|---|" * 3]
        names = sorted({n for bn in L.BINS for n in D[(sp, bn, "WP70")]["comp"] if re.match(r"(PC|CM|DMM|ML)\|", n)},
                       key=lambda n: (n.split("|")[0], -sum(D[(sp, bn, "WP70")]["comp"].get(n, {"n_full": 0})["n_full"] for bn in L.BINS)))
        for n in names:
            xs = [D[(sp, bn, "WP70")]["comp"].get(n, {"fit": {"n_win": 0, "n_sb": 0}})["fit"] for bn in L.BINS]
            if max(x["n_win"] + x["n_sb"] for x in xs) < 5:
                continue
            out.append(f"| {n} | " + " | ".join(f"{x['n_win']:.0f} | {x['n_sb']:.0f}" for x in xs) + " |")
        out.append("")
    open(os.path.join(T, "sklad.md"), "w").write("\n".join(out) + "\n")
    # ---- sygnał pozorny
    out = ["# BKGCAT — sygnał pozorny per składowa (POST HOC)", "",
           "Fit: szablon sygnału fd (kształt stały) + exp(Czebyszew 2 / 3) do histogramu (składowa mb_nd + sygnał Asimova N_prawda × szablon fd); N_spur = N_fit − N_wstrzyknięte, w % N_prawda; Z = N_spur/σ.",
           "Błąd: rms replik Poissona histogramu samej składowej (60 dla grup i sum, 30 dla podskładowych); sygnał Asimova bez fluktuacji. Kolumna „rząd + 1”: ten sam błąd.",
           "„obs − pasma” = (liczba w oknie − przewidywanie z fitu samych pasm)/N_prawda [%], błąd = √N_okno/N_prawda (bez błędu ekstrapolacji).",
           "Zakresy: nom = nominalny, alt = alternatywny (D⁰ 1,70–2,05 = zakres M8; D⁺ 1,75–2,10). Składowe z < 30 kandydatami w zakresie pominięte.", ""]
    ntest, nz2 = 0, 0
    for sp in L.SPS:
        for cn in ("WP70", "B[70,50)", "WP50", "WP90"):
            out += [f"## {sp} — {cn}", "", "| składowa | " + " | ".join(f"{b}: nom | Z | alt | Z | obs − pasma (nom)" for b in L.BINS) + " |", "|---|" + "---|---|---|---|---|" * 3]
            allc = {n for bn in L.BINS for n in D[(sp, bn, cn)]["comp"]} - {"all", "Sig", "Refl"}
            order = ["reszta (wszystko − Sig − Refl)", "DMM+CM", "DMM", "CM", "PC+Comb", "PC", "Comb", "NoCharm", "OtherB", "MB", "ML"]
            order += sorted(allc - set(order))
            for n in order:
                row, any_ = [], False
                for bn in L.BINS:
                    r = D[(sp, bn, cn)]
                    x = r["comp"].get(n)
                    if not x or "N" not in x["fit"] or x["fit"]["err"] is None:
                        row.append("— | — | — | — | —"); continue
                    any_ = True
                    cell = []
                    for rn in ("fit", "alt"):
                        y = x[rn]
                        if "N" not in y or y["err"] is None:
                            cell += ["—", "—"]; continue
                        if (not y["valid"]) or (not y.get("boot") and y["err"] < 0.3 * np.sqrt(max(y["n_win"], 1.0))):
                            cell += [f"{100 * y['N'] / r['N_true']:+.1f} (fit niestabilny)", "—"]; continue
                        cell += [L.pm(y["N"] / r["N_true"], y["err"] / r["N_true"], 1, True), f"{y['N'] / y['err']:+.1f}"]
                        if rn == "fit" and "|" not in n and n not in ("DMM+CM", "DMM", "CM") and cn in ("WP70", "B[70,50)"):
                            ntest += 1; nz2 += abs(y["N"] / y["err"]) > 2
                    y = x["fit"]
                    cell.append(L.pm((y["n_win"] - y["pred_win_sb"]) / r["N_true"], np.sqrt(y["n_win"]) / r["N_true"], 1, True))
                    row.append(" | ".join(cell))
                if any_:
                    out.append(f"| {L.GNAME.get(n, n)} | " + " | ".join(row) + " |")
            out.append("")
    out += ["## Średnie ważone po binach pT (N_spur/N_prawda [%], zakres nominalny) i rząd tła + 1", "",
            "| gat. | kategoria | składowa | 3–5 | 5–8 | 8–30 | średnia ważona | Z | rząd + 1: 3–5 | 5–8 | 8–30 |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    avg = {}
    for sp in L.SPS:
        for cn in ("WP70", "B[70,50)", "WP50", "WP90", "B[90,70)"):
            for n in MAIN:
                v, e, cells, c1 = [], [], [], []
                for bn in L.BINS:
                    r = D[(sp, bn, cn)]
                    x = r["comp"].get(n)
                    y = x["fit"] if x else {}
                    if "N" not in y or y["err"] is None or not y["valid"]:
                        cells.append("—"); c1.append("—"); continue
                    v.append(y["N"] / r["N_true"]); e.append(y["err"] / r["N_true"])
                    cells.append(L.pm(v[-1], e[-1], 1, True))
                    p1 = x.get("fit_p1")
                    c1.append(L.pm(p1["N"] / r["N_true"], p1["err"] / r["N_true"], 1, True) if p1 and p1["valid"] and p1["err"] else (f"{100 * p1['N'] / r['N_true']:+.1f} (?)" if p1 else "—"))
                if len(v) < 2:
                    continue
                v, e = np.array(v), np.array(e)
                w = 1 / e ** 2
                mu, se = float((w * v).sum() / w.sum()), float(1 / np.sqrt(w.sum()))
                avg[f"{sp}|{cn}|{n}"] = {"mean": mu, "err": se, "per_bin": v.tolist(), "per_bin_err": e.tolist()}
                out.append(f"| {sp} | {cn} | {L.GNAME.get(n, n)} | " + " | ".join(cells) + f" | {L.pm(mu, se, 1, True)} | {mu / se:+.1f} | " + " | ".join(c1) + " |")
    out.append("")
    json.dump(avg, open(os.path.join(L.BASE, "out", "avg.json"), "w"), indent=1)
    # ---- głębokość względna w rozłącznych przedziałach score
    out += ["## Struktura czy fluktuacja: głębokość względna d = N_spur/N_składowej w oknie [%] w ROZŁĄCZNYCH przedziałach score (zakres nominalny)", "",
            "Stała struktura kształtu składowej dawałaby to samo d w każdym przedziale; χ² = zgodność przedziałów ze średnią ważoną, Z = średnia/błąd.", "",
            "| gat. | bin | składowa | [95, 90) | [90, 70) | [70, 50) | WP50 | średnia ważona | Z | χ²/ndf zgodności |", "|---|---|---|---|---|---|---|---|---|---|"]
    SL = ["B[95,90)", "B[90,70)", "B[70,50)", "WP50"]
    summ = {}
    for sp in L.SPS:
        for bn in L.BINS:
            for n in ("reszta (wszystko − Sig − Refl)", "PC+Comb", "NoCharm", "PC", "DMM+CM"):
                v, e, cells = [], [], []
                for cn in SL:
                    x = D[(sp, bn, cn)]["comp"].get(n)
                    y = x["fit"] if x else {}
                    if "N" not in y or y["err"] is None or not y["valid"] or y["n_win"] < 20:
                        cells.append("—"); continue
                    v.append(y["N"] / y["n_win"]); e.append(y["err"] / y["n_win"])
                    cells.append(L.pm(v[-1], e[-1], 1, True))
                if len(v) < 2:
                    continue
                v, e = np.array(v), np.array(e)
                w = 1 / e ** 2
                mu, se = float((w * v).sum() / w.sum()), float(1 / np.sqrt(w.sum()))
                c2 = float((((v - mu) / e) ** 2).sum())
                summ[f"{sp}|{bn}|{n}"] = {"mean": mu, "err": se, "chi2": c2, "ndf": len(v) - 1}
                out.append(f"| {sp} | {bn} | {L.GNAME.get(n, n)} | " + " | ".join(cells) + f" | {L.pm(mu, se, 2, True)} | {mu / se:+.1f} | {c2:.1f}/{len(v) - 1} |")
    out.append("")
    json.dump(summ, open(os.path.join(L.BASE, "out", "depth.json"), "w"), indent=1)
    out.append(f"Składowe niepikujące (grupy i sumy, WP70 i [70, 50), zakres nominalny): {ntest} testów, |Z| > 2 w {nz2} (oczekiwane przy samych fluktuacjach: {0.0455 * ntest:.1f}; testy skorelowane — sumy zawierają składniki).")
    open(os.path.join(T, "spur.md"), "w").write("\n".join(out) + "\n")
    # ---- kształt 25 MeV przy WP70
    out = ["# BKGCAT — kształt składowych w masie przy WP70 (biny 25 MeV, pełny zakres histogramów; liczby kandydatów mb_nd)", ""]
    for sp in L.SPS:
        lo, hi = F.SP[sp]["mr"]
        e = np.arange(lo, hi + 1e-9, 0.025)
        for bn in L.BINS:
            r = D[(sp, bn, "WP70")]["comp"]
            out += [f"## {sp} {bn} GeV", "", "| składowa | " + " | ".join(f"{x:.3f}" for x in e[:-1]) + " |", "|---|" + "---|" * (len(e) - 1)]
            for n in ["Sig", "Refl", "DMM", "CM", "PC", "NoCharm", "OtherB", "MB", "ML"] + sorted(k for k in r if re.match(r"(PC|CM|DMM)\|", k) and r[k]["n_full"] >= 40):
                if n in r and "h25" in r[n]:
                    out.append(f"| {n} | " + " | ".join(f"{v:.0f}" for v in r[n]["h25"]) + " |")
                elif n in r:
                    hh = L.h5(B_cat_cache(sp, bn)[n]).reshape(-1, 5).sum(1)
                    out.append(f"| {n} | " + " | ".join(f"{v:.0f}" for v in hh) + " |")
            out.append("")
    open(os.path.join(T, "ksztalt.md"), "w").write("\n".join(out) + "\n")
    print("zapisano tabele/sklad.md, spur.md, ksztalt.md")


_cc = {}


def B_cat_cache(sp, bn):
    if (sp, bn) not in _cc:
        if B is None:
            _init(sys.argv[sys.argv.index("--in") + 1] if "--in" in sys.argv else L.IN)
        _cc[(sp, bn)] = comps(sp, bn, 0.7, None)
    return _cc[(sp, bn)]


if __name__ == "__main__":
    main()
