#!/usr/bin/env python3
"""SALAD_427929 — lekka analiza lokalna: składa WYNIKI.md z json pobranych z EOS (etap ORCH `salad427`).

Nic nie trenuje i nie czyta kandydatów: tylko runs/lxplus/s427_{D0,Dplus}.json i s427_{D0,Dplus}_adv.json (< 5 s, 1 rdzeń).
Kryteria i struktura: RAPORT.md (zapisane 17.11.2026 przed wynikami). Werdykty liczone tu mechanicznie wg RAPORT.md §3.
Użycie: .venv_h5/bin/python charm_work/SALAD_427929/analyze.py [--in DIR] [--out WYNIKI.md]"""
import argparse
import json
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
WPS = (0.95, 0.80, 0.70, 0.50)
BINS = ["3-5", "5-8", "8-30", "all"]
CLASSES = ["MB427929", "J12only", "MBTSbased", "AFPonly"]
BUDGET = {"Dplus": {0.80: 0.05, 0.70: 0.10, 0.50: 0.25}, "D0": {0.80: 0.05, 0.70: 0.10, 0.50: 0.15}}
NAME = {"D0": "D⁰", "Dplus": "D⁺"}
GATES = {"neff": 0.5, "nsb": 20000, "c2st": 0.53}


def f(x, n=3):
    return "—" if x is None or (isinstance(x, float) and not math.isfinite(x)) else f"{x:.{n}f}".replace(".", ",")


def ratio(ev, bn, j, tag="salad", reg="SR", train=None):
    """R = eps_B(MC, tag, reg) / eps_B(dane); dla SR dane = interpolacja pol2 z SB. Zwraca (R, stat_dane, stat_mc, syst, tot)."""
    e = ev[bn]
    num, nerr = e[f"mc_{tag}"][reg]["eps"][j], e[f"mc_{tag}"][reg]["err"][j]
    if reg == "SR":
        den, derr, syst = e["data_SR"]["pol2"][j], e["data_SR"]["pol2_stat"][j], e["data_SR"]["syst_interp"][j]
    else:
        den, derr, syst = e[f"data_{reg}"]["eps"][j], e[f"data_{reg}"]["err"][j], 0.0
    if not den or num is None or not math.isfinite(den) or not math.isfinite(num) or den <= 0:
        return (float("nan"),) * 5
    r = num / den
    sd, sm, ss = r * derr / den, (r * nerr / num if num > 0 else float("nan")), r * syst / den
    tr = train or 0.0
    return r, sd, sm, ss, math.sqrt(sd ** 2 + sm ** 2 + ss ** 2 + tr ** 2)


def verdict(r, tot, budget):
    """RAPORT.md §3: ZGODNE |R-1| <= B; NIEROZSTRZYGNIĘTE gdy przekracza B o mniej niż 2 sigma albo sigma > B; PRZEKROCZONE gdy >= 2 sigma ponad B."""
    if not math.isfinite(r) or not math.isfinite(tot):
        return "brak"
    d = abs(r - 1)
    if d <= budget:
        return "ZGODNE" if tot <= budget else "zgodne, niemiarodajne (σ > B)"
    return "PRZEKROCZONE" if (d - budget) >= 2 * tot else "nierozstrzygnięte"


def seed_spread(cl, bn, j):
    k = cl.get("kfold", {})
    if "seed0" in k and "seed1" in k:
        a, b = ratio(k["seed0"]["eval"], bn, j)[0], ratio(k["seed1"]["eval"], bn, j)[0]
        if math.isfinite(a) and math.isfinite(b):
            return abs(a - b) / math.sqrt(2)
    return 0.0


def section_species(sp, R, out):
    meta = R["meta"]
    out.append(f"\n## {NAME[sp]}\n")
    out.append(f"Wersja {meta['version']}, yaml SALAD {meta['salad_cfg_version']}, foldów {meta['k_folds']}, replik {meta['nboot']}, "
               f"ziarna {meta['seeds']}, MC tła truth ≥ {int(R['mc_truth_min'])} (N = {R['n_mc_bkg']}), tryb testowy: {'TAK' if meta['fast'] else 'nie'}.")
    out.append(f"Tagger H427: AUC OOF (fd vs SB danych MB427929) {', '.join(f(x) for x in meta['tagger_auc_oof'])}. "
               f"Liczności kandydatów: " + ", ".join(f"{k} {v}" for k, v in meta["n"].items()) + ".")
    # --- bramki v1.0
    out.append("\n### Bramki SALAD v1.0 per klasa (bieg „prod”: jeden fit na całym SB)\n")
    out.append("| klasa | N_SB danych | AUC_SB | n_eff/N (SR) | C2ST pasma | bramki (≥ 0,5 / ≥ 20 tys. / ≤ 0,53) | `ok` |")
    out.append("|---|---|---|---|---|---|---|")
    for c in CLASSES:
        cl = R["classes"].get(c)
        if not cl:
            continue
        if "skipped" in cl:
            out.append(f"| {c} | {cl['n_data_SB']} | — | — | — | pominięta: {cl['skipped']} | — |")
            continue
        p = cl["prod"]
        g = [p["neff_over_N_SR"] >= GATES["neff"], cl["n_data_SB"] >= GATES["nsb"], (p.get("c2st_band") or 0.5) <= GATES["c2st"]]
        out.append(f"| {c} | {cl['n_data_SB']} | {f(p['val_auc_SB'])} | {f(p['neff_over_N_SR'], 2)} | {f(p.get('c2st_band'))} | "
                   f"{' / '.join('tak' if x else 'NIE' for x in g)} | {'tak' if p['ok'] else 'NIE'} |")
    # --- k-fold: C2ST, n_eff, rzeźbienie
    out.append("\n### K-fold po zdarzeniach (wagi OOF): C2ST przed/po, n_eff/N, rzeźbienie masy\n")
    out.append("| klasa | AUC_SB (średnia foldów) | C2ST SB przed → po (podłoga) | C2ST pasmo przed → po | n_eff/N SR 3–5 / 5–8 / 8–30 / razem | "
               "wagi obcięte | corr(w, m) w SB | średnia w SR / pasmo / SB | χ²/ndf masy w SB przed → po | C2ST vs MB427929 (dane–dane, SB) |")
    out.append("|---|---|---|---|---|---|---|---|---|---|")
    for c in CLASSES:
        cl = R["classes"].get(c)
        if not cl or "kfold" not in cl:
            continue
        k = cl["kfold"]["seed0"]
        cc, s, ne = k["c2st"], k["sculpt"], k["neff_over_N"]["SR"]
        auc = sum(k["diag"]["val_auc_SB"]) / len(k["diag"]["val_auc_SB"])
        mr, ms = s["mass_SB_raw"], s["mass_SB_salad"]
        out.append(f"| {c} | {f(auc)} | {f(cc['raw_SB'])} → {f(cc['salad_SB'])} ({f(cc['floor_SB'])}) | {f(cc['raw_band'])} → {f(cc['salad_band'])} | "
                   f"{' / '.join(f(ne[b], 2) for b in BINS)} | {f(k['diag']['frac_clipped'], 4)} | {f(s['corr_w_m_SB'])} | "
                   f"{' / '.join(f(s['mean_w'][r]) for r in ('SR', 'band', 'SB'))} | {f(mr['chi2'] / mr['ndf'], 2)} → {f(ms['chi2'] / ms['ndf'], 2)} | "
                   f"{f(R.get('c2st_data_data_SB', {}).get(c)) if c != 'MB427929' else '—'} |")
    # --- eps_B: SB OOF, pasmo, SR
    for reg, title in (("band", "Domknięcie ε_B w paśmie walidacyjnym (dane zliczone; pasmo nie wchodzi do fitu SALAD)"),
                       ("SB", "Domknięcie ε_B w SB (wagi OOF — fold MC niewidziany w treningu)"),
                       ("SR", "Przewidywanie pod pikiem: ε_B(MC ważone, SR) / ε_B(dane, interpolacja exp(pol2) z SB)")):
        out.append(f"\n### {title}\n")
        out.append("R = ε_B(MC)/ε_B(dane); „przed” = MC surowe, „po” = MC z wagami SALAD (k-fold, ziarno 0); błąd = stat danych ⊕ stat MC"
                   + (" ⊕ syst interpolacji ⊕ rozrzut ziaren" if reg == "SR" else "") + ".\n")
        out.append("| klasa | bin pT | " + " | ".join(f"WP{int(w * 100)} przed → po" for w in WPS) + " |")
        out.append("|---|---|" + "---|" * len(WPS))
        for c in CLASSES:
            cl = R["classes"].get(c)
            if not cl or "kfold" not in cl:
                continue
            ev = cl["kfold"]["seed0"]["eval"]
            for bn in BINS:
                cells = []
                for j in range(len(WPS)):
                    r0 = ratio(ev, bn, j, "raw", reg)[0]
                    r1 = ratio(ev, bn, j, "salad", reg, seed_spread(cl, bn, j) if reg == "SR" else None)
                    cells.append(f"{f(r0)} → {f(r1[0])} ± {f(r1[4])}")
                out.append(f"| {c} | {bn} | " + " | ".join(cells) + " |")
    # --- budżet
    out.append("\n### Porównanie z budżetem (pod pikiem, k-fold) — werdykty wg RAPORT.md §3\n")
    out.append("| klasa | bin pT | WP | R − 1 | σ (stat dane / stat MC / syst interp. / ziarna) | budżet | werdykt | „prod” R |")
    out.append("|---|---|---|---|---|---|---|---|")
    for c in CLASSES:
        cl = R["classes"].get(c)
        if not cl or "kfold" not in cl:
            continue
        ev = cl["kfold"]["seed0"]["eval"]
        for bn in BINS:
            for j, wp in enumerate(WPS):
                tr = seed_spread(cl, bn, j)
                r, sd, sm, ss, tot = ratio(ev, bn, j, "salad", "SR", tr)
                b = BUDGET[sp].get(wp)
                rp = ratio(cl["prod"]["eval"], bn, j)[0]
                out.append(f"| {c} | {bn} | {int(wp * 100)} | {f(r - 1)} | {f(tot)} ({f(sd)} / {f(sm)} / {f(ss)} / {f(tr)}) | "
                           f"{'±' + f(b, 2) if b else 'brak (tylko raport)'} | {verdict(r, tot, b) if b else '—'} | {f(rp)} |")
    # --- estymator na MC, kształt masy
    out.append("\n### Kontrole estymatora (klasa MB427929)\n")
    cl = R["classes"].get("MB427929")
    if cl and "kfold" in cl:
        ev = cl["kfold"]["seed0"]["eval"]
        out.append("| bin pT | obciążenie interpolacji na MC ważonym: pol2(SB)/prawda(SR) @WP95/80/70/50 | N_SR/N_SB: dane (pol2) | MC surowe | MC ważone |")
        out.append("|---|---|---|---|---|")
        for bn in BINS:
            m = ev[bn]["mc_salad"]
            bias = [a / b if b else float("nan") for a, b in zip(m["pol2_on_mc"], m["SR"]["eps"])]
            out.append(f"| {bn} | {' / '.join(f(x) for x in bias)} | {f(ev[bn]['data_nSR_pred_over_SBbins'])} | "
                       f"{f(ev[bn]['mc_raw']['n_SR_over_SBbins'])} | {f(m['n_SR_over_SBbins'])} |")
    if "variant_truth8_MB427929" in R:
        v = R["variant_truth8_MB427929"]
        out.append("\nWariant kontrolny D⁰: MC truth ≥ 8 (bez truth 2–7 i odbić), MB427929, R pod pikiem (bin „all”): "
                   + " / ".join(f"WP{int(w * 100)} {f(ratio(v['eval'], 'all', j)[0])} ± {f(ratio(v['eval'], 'all', j)[4])}" for j, w in enumerate(WPS))
                   + f"; n_eff/N(SR) {f(v['neff_over_N']['SR']['all'], 2)}. Budżet WP50 dla tego wariantu: ±0,20.")
    # --- transfer
    out.append("\n### Transfer wag: wagi uczone w klasie a (wiersz), dane klasy b (kolumna)\n")
    out.append("Komórka: R pod pikiem (bin „all”) @WP80 / WP70 / WP50 ± σ(WP80); w nawiasie C2ST SB po przeważeniu. Przekątna = k-fold w klasie.\n")
    cs = [c for c in CLASSES if "kfold" in R["classes"].get(c, {})]
    out.append("| a \\ b | " + " | ".join(cs) + " |")
    out.append("|---|" + "---|" * len(cs))
    for a in cs:
        row = []
        for b in cs:
            if a == b:
                ev, cc = R["classes"][a]["kfold"]["seed0"]["eval"], R["classes"][a]["kfold"]["seed0"]["c2st"]["salad_SB"]
            else:
                t = R["transfer"].get(a, {}).get(b)
                if not t:
                    row.append("—")
                    continue
                ev, cc = t["eval"], t["c2st_salad_SB"]
            rr = [ratio(ev, "all", j) for j in (1, 2, 3)]
            row.append(f"{' / '.join(f(x[0]) for x in rr)} ± {f(rr[0][4])} ({f(cc)})")
        out.append(f"| {a} | " + " | ".join(row) + " |")
    out.append("\nPer bin pT, WP80 (R ± σ):\n")
    out.append("| a → b | 3–5 | 5–8 | 8–30 |")
    out.append("|---|---|---|---|")
    for a in cs:
        for b in cs:
            if a != b and R["transfer"].get(a, {}).get(b):
                ev = R["transfer"][a][b]["eval"]
                out.append(f"| {a} → {b} | " + " | ".join(f"{f(ratio(ev, bn, 1)[0])} ± {f(ratio(ev, bn, 1)[4])}" for bn in BINS[:3]) + " |")


def section_adv(sp, A, out):
    out.append(f"\n### {NAME[sp]}: closure adwersarialny B1 na MC v2.0 (pseudo-dane o statystyce SB klasy MB427929: {A['n_sb_target']})\n")
    out.append("| tryb | N pseudo SB / SR | f_MM SR pseudo / MC | C2ST surowe SB / SR | C2ST SALAD SB / pasmo / SR (podłoga) | n_eff/N | "
               "ε_B pred/prawda @WP95 / 80 / 70 / 50 (SALAD) | to samo, MC surowe | σ_stat prawdy |")
    out.append("|---|---|---|---|---|---|---|---|---|")
    for mode, r in A["modes"].items():
        e = r["eps"]["all"]
        ps = [a / b if b else float("nan") for a, b in zip(e["salad"], e["true"])]
        pr = [a / b if b else float("nan") for a, b in zip(e["mcraw"], e["true"])]
        st = [a / b if b else float("nan") for a, b in zip(e["true_err"], e["true"])]
        out.append(f"| {mode} | {r['n_pseudo_SB']} / {r['n_pseudo_SR']} | {f(r['f_MM_SR_pseudo'])} / {f(r['f_MM_SR_A'])} | "
                   f"{f(r['c2st_raw_SB'])} / {f(r['c2st_raw_SR'])} | {f(r['c2st_salad_SB'])} / {f(r['c2st_salad_band'])} / {f(r['c2st_salad_SR'])} "
                   f"({f(r['c2st_floor_SR'])}) | {f(r['neff_over_N_SR'], 2)} | {' / '.join(f(x) for x in ps)} | {' / '.join(f(x) for x in pr)} | "
                   f"{' / '.join(f(x) for x in st)} |")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", default=os.path.join(HERE, "runs", "lxplus"))
    ap.add_argument("--out", default=os.path.join(HERE, "WYNIKI.md"))
    a = ap.parse_args()
    out = ["# SALAD_427929 — wyniki (generowane przez analyze.py z json etapu `salad427`; kryteria: RAPORT.md)\n",
           "Tylko pasma boczne i tło. W SR danych liczona jest wyłącznie interpolacja z pasm bocznych; żadnych fitów piku. Pary D–D poza zakresem."]
    n = 0
    for sp in ("D0", "Dplus"):
        p = os.path.join(a.inp, f"s427_{sp}.json")
        if os.path.exists(p):
            section_species(sp, json.load(open(p)), out)
            n += 1
        else:
            out.append(f"\n## {NAME[sp]}\n\nBrak {p}.")
        p = os.path.join(a.inp, f"s427_{sp}_adv.json")
        if os.path.exists(p):
            section_adv(sp, json.load(open(p)), out)
    open(a.out, "w").write("\n".join(out) + "\n")
    print(f"zapisano {a.out} ({n} gatunków, {len(out)} wierszy)")


if __name__ == "__main__":
    main()
