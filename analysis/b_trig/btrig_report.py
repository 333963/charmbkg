#!/usr/bin/env python3
"""B_trig — RAPORT_DANE.md: bramki G1–G8 z PROCEDURA.md na wynikach btrig_analyze (effmap/sf/epsb/iso) + bt_summary.json.
Użycie: btrig_report.py --in runs/lxplus --out runs/lxplus/out   (wywoływane przez run_analysis.sh)."""
import argparse
import json
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import btrig_common as C  # noqa: E402

BINS_MAIN = ["5-8", "8-30"]
BINS_ALL = ["3-5", "5-8", "8-30"]
SHOW = [(428770, "all"), (428770, "J12only"), (428770, "MBTSbased"), (428770, "other"), (427929, "all427929"), (427929, "MB427929"),
        (427929, "sptrk4_427929"), (428770, "sptrk4")]


def jl(p):
    try:
        return json.load(open(p))
    except Exception:  # noqa: BLE001
        return None


def fin(x):
    return isinstance(x, (int, float)) and x is not None and math.isfinite(x)


def ve(v, e, nd=2):
    return f"{v:.{nd}f} ± {e:.{nd}f}" if fin(v) and fin(e) else "–"


def get(rows, run, cls, b, key, val):
    for r in rows:
        if r["run"] == run and r["cls"] == cls and r["bin"] == b and r.get(key) == val:
            return r
    return None


def ok(r):
    return bool(r) and fin(r.get("SF")) and fin(r.get("SF_err")) and r["SF_err"] > 0 and r.get("determined")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    summ = jl(os.path.join(a.inp, "bt_summary.json")) or {}
    sf, eb, iso, em = (jl(os.path.join(a.out, f"btrig_{n}.json")) for n in ("sf", "epsb", "iso", "effmap"))
    G = []          # (bramka, status, opis)
    L = ["# B_trig — RAPORT NA DANYCH (h5 v2.0, L_xy ≥ 0,25 mm): SF / ε_B / izolacja per klasa triggera i run", "",
         f"Wygenerowano {time.strftime('%Y-%m-%d %H:%M')} przez `run_analysis.sh` (`btrig_report.py`). Wejście `{a.inp}`, wyniki cząstkowe "
         f"`{a.out}/btrig_{{effmap,sf,epsb,iso}}.md`. Projekt B, pojedyncze D (nie ślepe). Bramki: `PROCEDURA.md` §2.", ""]

    # ------------------------------------------------------------------ próbki / G1
    if summ:
        L += ["## Próbki (bt_summary.json)", "", "| próbka | grupy | pary D0←D* | D0 all | zdarzenia | dup (run, event) |", "|---|---|---|---|---|---|"]
        for s, v in summ.items():
            L.append(f"| {s} | {v.get('n_groups')} | {v.get('n_pairs')} | {v.get('n_all')} | {v.get('n_events')} | {v.get('ev_dup_run_event')} |")
        L.append("")
        dup = {s: v.get("ev_dup_run_event", 0) for s, v in summ.items() if s.startswith("data_")}
        nf = max([d.get("notfound", 0.0) for d in (sf or {}).get("inputs", {}).get("data", [])] or [0.0])
        src = [d.get("src") for d in (sf or {}).get("inputs", {}).get("data", [])]
        G.append(("G1 złączenie", "ZDANE" if all(v == 0 for v in dup.values()) and nf < 1e-3 else "NIE",
                  f"duplikaty (run, event) {dup}; bity z {src}; nieznalezione {nf:.4%}"))
    else:
        G.append(("G1 złączenie", "BRAK", "brak bt_summary.json"))

    rows = (sf or {}).get("rows", [])
    tests = {(t["bin"], t["wp"]): t for t in (sf or {}).get("tests", [])}
    wps = sorted({r["wp"] for r in rows}, reverse=True)

    # ------------------------------------------------------------------ G2
    if rows:
        st, desc = "ZDANE", []
        for b in BINS_MAIN:
            r = get(rows, 428770, "all", b, "wp", "0.95")
            nv, ne = C.NOMINAL_SF95[b]
            if not ok(r):
                st = "BRAK"; desc.append(f"{b}: fit nieokreślony"); continue
            d = r["SF"] - nv
            if abs(d) > 0.02:
                st = "NIE"
            desc.append(f"{b}: SF95(all) {ve(r['SF'], r['SF_err'])} vs nominał {nv:.2f} ± {ne:.2f} (Δ {d:+.3f}), N_sig {r['N_sig']:.0f} ± {r['N_sig_err']:.0f}")
        G.append(("G2 domknięcie z nominałem (±0,02)", st, "; ".join(desc) + ". N_sig ±2 % vs EV — porównać ręcznie z `EV_dstar/EV_pdgid.md` (nie w JSON)."))
    else:
        G.append(("G2 domknięcie z nominałem", "BRAK", "brak btrig_sf.json"))

    # ------------------------------------------------------------------ G3
    cl = (sf or {}).get("classes", {})
    if cl:
        n = {k.split(":")[1]: v["n_pairs"] for k, v in cl.items() if k.startswith("428770:")}
        desc, st = [], "ZDANE"
        if n:
            s3 = n.get("J12only", 0) + n.get("MBTSbased", 0) + n.get("other", 0)
            st = "ZDANE" if s3 == n.get("all") else "NIE"
            tot = max(n.get("all", 0), 1)
            desc.append("428770 pary: " + ", ".join(f"{c} {n.get(c, 0)} ({n.get(c, 0) / tot:.1%})" for c in C.CLASSES_428770) + f"; suma rozłącznych {s3}")
        n9 = {k.split(":")[1]: v["n_pairs"] for k, v in cl.items() if k.startswith("427929:")}
        if n9:
            desc.append("427929 pary: " + ", ".join(f"{c} {v}" for c, v in n9.items()))
        p9 = os.path.join(a.inp, "bt_data_427929.npz")
        if os.path.exists(p9):
            z = np.load(p9, allow_pickle=False)
            k = "ev_ps_HLT_mb_mbts_L1MBTS_1"
            if k in z.files:
                g = z["ev_grl"] == 1 if "ev_grl" in z.files else np.ones(len(z[k]), bool)
                f1 = float((z[k][g] == 1).mean()) if g.any() else float("nan")
                desc.append(f"427929: PS(mbts_L1MBTS_1) = 1 w {f1:.2%} zdarzeń GRL z kandydatem (miękki próg 95 %{'' if f1 >= 0.95 else ' — NIE: waga PS'})")
        G.append(("G3 ułamki klas / PS", st, "; ".join(desc)))

    # ------------------------------------------------------------------ G4
    toy = [f for f in os.listdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")) if f.startswith("btrig_toy")] \
        if os.path.isdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")) else []
    G.append(("G4 toy pass/fail", "INFO" if toy else "NIE LICZONE", f"`out/{toy}`" if toy else
              "`btrig_toy.py` przy (f_S, f_B) klas — osobno; klasy×biny z det = NIE traktować jako niemierzalne"))

    # ------------------------------------------------------------------ G5
    c5 = (em or {}).get("closure_428770")
    if c5 and c5.get("ndf"):
        x = c5["chi2"] / c5["ndf"]
        mr = c5.get("max_rel_dev")
        G.append(("G5 domknięcie emulacji (miękka)", "ZDANE" if x < 2 and fin(mr) and mr < 0.20 else "NIE",
                  f"χ²/ndf = {c5['chi2']:.1f}/{c5['ndf']} = {x:.2f}, max |pred/obs − 1| = {mr if not fin(mr) else round(mr, 3)}; "
                  f"n_ref = {em.get('n_ref')}, ε_k = {({k: round(v, 4) for k, v in em.get('eps_global', {}).items() if fin(v)})}"
                  + ("" if x < 2 else " → ε_MC z metody B (gęstość), |A − B| jako syst.")))
    else:
        G.append(("G5 domknięcie emulacji", "BRAK", "brak btrig_effmap.json / tabeli 428770 (wagi MC = 1)"))

    # ------------------------------------------------------------------ tabela SF + G6
    if rows:
        for wp in wps:
            L += [f"## SF klasyfikatora H (fit Δm pass/fail), WP {wp}", "", "| run | klasa | " + " | ".join(f"{b} GeV: N_sig; SF" for b in BINS_ALL) + " |",
                  "|---|---|" + "---|" * len(BINS_ALL)]
            for run, c in SHOW:
                cells = []
                for b in BINS_ALL:
                    r = get(rows, run, c, b, "wp", wp)
                    cells.append("–" if not r else f"{r['N_sig']:.0f}; {ve(r['SF'], r['SF_err'])}" + ("" if r.get("determined") else " (nieokr.)"))
                if any(x != "–" for x in cells):
                    L.append(f"| {run} | {c} | " + " | ".join(cells) + " |")
            L.append("")
            if wp != "0.95":
                continue
            L += ["Bin 3–5 GeV tylko informacyjnie (niemierzalny w nominale). χ² klas rozłącznych {J12only, MBTSbased, other}₄₂₈₇₇₀ ∪ {MB427929}:", "",
                  "| bin | χ²/ndf (p) | człony | J12only − MBTSbased [σ] | MB427929 vs nominał [σ] | MB427929 vs 1 [σ] | 427929 vs 428770 all [σ] |", "|---|---|---|---|---|---|---|"]
            for b in BINS_ALL:
                t = tests.get((b, wp))
                if not t:
                    continue
                c = t["chi2_disjoint"]
                rj, rm = get(rows, 428770, "J12only", b, "wp", wp), get(rows, 428770, "MBTSbased", b, "wp", wp)
                pjm = C.pull_indep(rj["SF"], rj["SF_err"], rm["SF"], rm["SF_err"]) if ok(rj) and ok(rm) else None
                f = lambda x: f"{x:+.2f}" if fin(x) else "–"  # noqa: E731
                cs = f"{c['chi2']:.1f}/{c['ndf']} ({c['p']:.3f})" if fin(c.get("chi2")) else "–"
                L.append(f"| {b} | {cs} | {', '.join(t.get('members', []))} | {f(pjm)} | {f(t['pulls'].get('MB427929_vs_nominal'))} | "
                         f"{f(t.get('MB427929_vs_1'))} | {f(t['pulls'].get('MB427929_vs_428770all'))} |")
            L.append("")
        st, desc = "ZDANE", []
        for b in BINS_MAIN:
            t = tests.get((b, "0.95"))
            r9 = get(rows, 427929, "MB427929", b, "wp", "0.95")
            if not t:
                st = "BRAK"; continue
            p = t["chi2_disjoint"].get("p")
            pn = t["pulls"].get("MB427929_vs_nominal")
            bad = (fin(p) and p <= 0.01) or (fin(pn) and abs(pn) >= 2)
            if not fin(p) or not fin(pn):
                st = "NIEKONKLUZYWNE" if st == "ZDANE" else st
            if bad:
                st = "NIE"
            v1 = t.get("MB427929_vs_1")
            s = f"{b}: p(χ²) = {p:.3f}" if fin(p) else f"{b}: χ² niedostępne"
            if ok(r9):
                s += f", SF95(MB427929) = {ve(r9['SF'], r9['SF_err'])} (vs nominał {pn:+.1f}σ" + (f", vs 1 {v1:+.1f}σ)" if fin(v1) else ")")
                if b == "5-8":
                    if r9["SF_err"] > 0.12:
                        s += " → §3: NIEKONKLUZYWNE (σ > 0,12): bez zmiany nominału"
                    elif fin(v1) and abs(v1) < 1 and fin(pn) and pn >= 2:
                        s += " → §3: anomalia = obciążenie J12/aktywnością 428770; nominał minbias = MBTS/427929"
                    elif fin(pn) and abs(pn) < 1 and fin(v1) and v1 <= -2:
                        s += " → §3: anomalia w modelowaniu MC (nie trigger); syst. 0,17 zostaje"
                    else:
                        s += " → §3: nie rozróżnia hipotez „0,83” i „1”"
            else:
                s += ", MB427929 nieokreślone"
            desc.append(s)
        G.append(("G6 zgodność klas i runu (WP95)", st, "; ".join(desc)))
    else:
        G.append(("G6 zgodność klas i runu", "BRAK", "brak btrig_sf.json"))

    # ------------------------------------------------------------------ ε_B / G7
    if eb and eb.get("rows"):
        er = eb["rows"]
        ewps = sorted({r["wp"] for r in er}, reverse=True)
        L += ["## ε_B pod pikiem per klasa (symetryczne SB; pol2 ML z kształtem klasy odniesienia; R = pol2 / MC tła z emulacją)", "",
              "| run | klasa | bin | " + " | ".join(f"WP{w}: sym; pol2; R" for w in ewps) + " |", "|---|---|---|" + "---|" * len(ewps)]
        f3 = lambda x: f"{x:.3f}" if fin(x) else "–"  # noqa: E731
        for b in BINS_ALL:
            for run, c in SHOW:
                cells = []
                for w in ewps:
                    r = get(er, run, c, b, "wp", w)
                    cells.append("–" if not r else f"{f3(r.get('sym'))}; {f3(r.get('pol2'))} ± {f3((r.get('pol2') or 0) * (r.get('pol2_stat_rel') or float('nan')))}; {f3(r.get('R_data_mc'))}")
                if any(x != "–" for x in cells):
                    L.append(f"| {run} | {c} | {b} | " + " | ".join(cells) + " |")
        L.append("")
        dec = [d for d in eb.get("decision", []) if d["bin"] in BINS_MAIN]
        need = [d for d in dec if d.get("per_trigger_reweighting_needed")]
        G.append(("G7 ε_B per klasa (reweighting per trigger?)", "BRAK" if not dec else ("TAK — per trigger" if need else "nie — wspólny SALAD"),
                  "; ".join(f"{d['bin']}@{d['wp']}: max|Δ| {d['max_rel_dev']:.1%} (tol {d['tol']:.0%}), p = {d['chi2']['p']:.3f}" if fin(d['chi2'].get('p'))
                            else f"{d['bin']}@{d['wp']}: max|Δ| {d['max_rel_dev']:.1%}, χ² niedostępne" for d in dec) + f" [{dec[0]['quantity'] if dec else ''}]"))
    else:
        G.append(("G7 ε_B per klasa", "BRAK", "brak btrig_epsb.json"))

    # ------------------------------------------------------------------ izolacja / G8
    if iso and iso.get("rows"):
        ir = iso["rows"]
        cuts = sorted({r["cut"] for r in ir})
        L += ["## SF izolacji per klasa (score_H > 0,6, |m − m0| < 40 MeV)", "",
              "| cięcie | pT | 428770 all | J12only | MBTSbased | MB427929 | max |pull| (pary klas) |", "|---|---|---|---|---|---|---|"]
        worst, nbad, ntest = 0.0, 0, 0
        for cut in cuts:
            for b in ["3-30"] + BINS_MAIN:
                rr = {c: get(ir, run, c, b, "cut", cut) for run, c in [(428770, "all"), (428770, "J12only"), (428770, "MBTSbased"), (427929, "MB427929")]}
                tri = [rr[c] for c in ("J12only", "MBTSbased", "MB427929") if ok(rr[c])]
                pl = [abs(C.pull_indep(x["SF"], x["SF_err"], y["SF"], y["SF_err"])) for i, x in enumerate(tri) for y in tri[i + 1:]]
                mp = max(pl) if pl else float("nan")
                if pl:
                    ntest += 1; worst = max(worst, mp); nbad += mp >= 2
                L.append(f"| {cut} | {b} | " + " | ".join(ve(rr[c]["SF"], rr[c]["SF_err"]) + ("" if rr[c].get("determined") else " (nieokr.)") if rr[c] else "–"
                                                            for c in ("all", "J12only", "MBTSbased", "MB427929")) + f" | {mp:.1f} |" if fin(mp) else
                         f"| {cut} | {b} | " + " | ".join(ve(rr[c]["SF"], rr[c]["SF_err"]) if rr[c] else "–" for c in ("all", "J12only", "MBTSbased", "MB427929")) + " | – |")
        L.append("")
        G.append(("G8 izolacja per klasa (pull < 2)", "BRAK" if not ntest else ("ZDANE" if nbad == 0 else "NIE"),
                  f"{ntest} zestawów (cięcie × pT), max |pull| = {worst:.1f}, ≥ 2σ: {nbad}" + (f" (oczekiwane z fluktuacji przy {ntest} skorelowanych testach: kilka %)" if nbad else "")))
    else:
        G.append(("G8 izolacja per klasa", "BRAK", (iso or {}).get("skipped", "brak btrig_iso.json")))

    head = ["## Bramki", "", "| bramka | wynik | szczegóły |", "|---|---|---|"] + [f"| {g} | **{s}** | {d} |" for g, s, d in G] + [""]
    thr = (sf or {}).get("thresholds")
    tail = []
    if thr:
        tail = ["## Progi WP (kwantyle score na prawdziwych D0 z fd_nd v2.0, okno ±55 MeV)", "", "| bin | " + " | ".join(C.WPS) + " |", "|---|" + "---|" * len(C.WPS)]
        tail += [f"| {b} | " + " | ".join(f"{thr[b][w]:.4f}" if fin(thr[b].get(w)) else "–" for w in C.WPS) + " |" for b in BINS_ALL if b in thr]
        tail += ["", "Odniesienie EV (v1.1): `EV_dstar/lxplus/ev_analysis_cm.json` → thresholds.", ""]
    tail += ["## Uwagi", "", "- Błędy SF: profil (albo HESSE przy `FAST=1`) ⊕ statystyka MC (n_eff po wagach emulacji).",
             "- Klasy 428770 po prescale (trigPassed); MC bez bitów → wagi emulacji z mapy 427929 (G5) albo 1.",
             "- 427929: wymagany passGRL; kompletność kontenera (458/583 plików AOD) niepotwierdzona — nie wpływa na SF, tylko na statystykę."]
    out = os.path.join(a.out, "RAPORT_DANE.md")
    open(out, "w").write("\n".join(L[:3] + head + L[3:] + tail) + "\n")
    json.dump([{"gate": g, "status": s, "detail": d} for g, s, d in G], open(os.path.join(a.out, "gates.json"), "w"), indent=1, ensure_ascii=False)
    print(f"zapisano {out}")
    for g, s, d in G:
        print(f"  {g}: {s}")


if __name__ == "__main__":
    sys.exit(main())
