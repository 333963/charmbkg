#!/usr/bin/env python3
"""DNET etap 6 — JEDNORAZOWY ODCZYT werdyktu wg DNET/KRYTERIA_etap6.md (progi wiążące, zapisane przed wynikami).

  ~/Documents/charm/.venv_h5/bin/python DNET/verdict6.py [--dir DNET/runs/lxplus6] [--out PLIK.md]

Czyta: dnet6_p1.json (dn6_ana.py), dnet6_results.json (analyze.py compute --sf-all), dnet6_mask.json (dn6_study.py), opcjonalnie dnet6_R.json
(dn4_ana.py r --force-data), meta/, fail/, DONE, DONE_R. Wypisuje i zapisuje <dir>/WERDYKT_etap6.md.
ŻADNYCH fitów: tylko odczyt i porównanie z progami. Brak pliku / pola / fit nieokreślony → NIEROZSTRZYGNIĘTE albo NIESPEŁNIONE wg zapisu (§4: brak
pomiaru SF = S niespełnione), nigdy „spełnione”."""
import argparse
import glob
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

# ------------------------------------------------------------------------------------------------ progi (KRYTERIA_etap6.md)
P1_DR, P1_NSIG = 0.03, 2.0                         # §3: |DR − 1| ≤ 0,03 w 3/3 binów, ziarno 1; przekroczenie < 2σ → adnotacja „na granicy mocy”
S_RANGE = 0.04                                     # §4: rozstęp SF(D*) WP70 między ziarnami 1–3, biny 5–8 i 8–30, oba runy (klasy główne)
S_BINS = ["5-8", "8-30"]
K2_NSIG, K2_SLOPE = 2.0, 0.10                      # RAPORT §4 K2
G_HAVE, G_NONE, G_NBIN = 1.20, 1.10, 2             # §4: etykiety opisowe proxy G na pasmach 427929
D_GFIT, D_NBIN = 1.15, 2                           # §5 (próg etapu 4): G_fit ≥ 1,15 w ≥ 2 z 3 binów w klasie głównej

BINS = ["3-5", "5-8", "8-30"]
MAIN = {"427929": "MB427929", "428770": "all"}
REF = "D0.xgb.H.sb.n.s1"
NP = [f"D0.tr.H+e04np.mc.d5.s{s}" for s in (1, 2, 3)]
FULL = [f"D0.tr.H+e04.mc.d5.s{s}" for s in (1, 2, 3)]
OK, NO, UN, DESC = "SPEŁNIONE", "NIESPEŁNIONE", "NIEROZSTRZYGNIĘTE", "opisowo"
ST = {True: OK, False: NO, None: UN}


def fin(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def pm(x, e=None, nd=3):
    if not fin(x):
        return "–"
    return f"{x:.{nd}f}" + (f" ± {e:.{nd}f}" if fin(e) else "")


def dig(d, *ks):
    for k in ks:
        if not isinstance(d, dict) or k not in d:
            return None
        d = d[k]
    return d


def both(*s):
    return False if any(x is False for x in s) else (None if any(x is None for x in s) else True)


class Doc:
    def __init__(self):
        self.L = []

    def h(self, s, lvl=2):
        self.L += ["", "#" * lvl + " " + s, ""]

    def p(self, s=""):
        self.L.append(s)

    def table(self, rows, head=("kryterium", "wielkość", "wartość ± błąd", "próg", "status", "uwagi")):
        if not rows:
            self.L.append("(brak wierszy)")
            return
        self.L += ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
        for r in rows:
            r = list(r) + [""] * (len(head) - len(r))
            self.L.append("| " + " | ".join(str(x).replace("|", "\\|") for x in r) + " |")
        self.L.append("")


def load(d, rel, miss):
    p = os.path.join(d, rel)
    try:
        j = json.load(open(p))
        if not isinstance(j, dict):
            raise ValueError
        return j
    except Exception as e:  # noqa: BLE001
        miss.append(f"{rel}: {'brak pliku' if not os.path.isfile(p) else 'nieczytelny (' + type(e).__name__ + ')'}")
        return None


# ------------------------------------------------------------------------------------------------ kryteria
def crit_P1(P, tag):
    """→ (status True/False/None, wiersze, lista DR)."""
    rows, fl, drs = [], [], []
    for b in BINS:
        x = dig(P, "models", tag, "bins", b, "0.7")
        dr = (x or {}).get("double_ratio")
        if not dr or not fin(dr[0]):
            rows.append(("P1", f"{b}: DR = [ε_a/ε_b]_sieć / [ε_a/ε_b]_drzewo", "–", f"|… − 1| ≤ {P1_DR}", UN, "brak wartości"))
            fl.append(None); drs.append(None)
            continue
        ok = abs(dr[0] - 1) <= P1_DR
        ns = abs(dr[0] - 1) / dr[1] if fin(dr[1]) and dr[1] > 0 else float("nan")
        note = f"ε_a {pm(x['a_tag'][0], x['a_tag'][1])}, ε_b {pm(x['b_rest'][0], x['b_rest'][1])}, a/b {pm(x['a_over_b'][0], x['a_over_b'][1])}; {ns:.1f}σ od 1"
        if not ok and fin(ns) and ns < P1_NSIG:
            note += " — NA GRANICY MOCY TESTU (< 2σ)"
        rows.append(("P1", f"{b}: DR = [ε_a/ε_b]_sieć / [ε_a/ε_b]_drzewo", pm(dr[0], dr[1], 4), f"|… − 1| ≤ {P1_DR}", OK if ok else NO, note))
        fl.append(ok); drs.append(dr)
    st = False if any(f is False for f in fl) else (None if any(f is None for f in fl) else True)
    return st, rows, drs


def sf_get(R, tag, run, cls, b, wp="0.7"):
    x = dig(R, "species", "D0", "sf", f"{tag}|{run}|{cls}|{b}|{wp}")
    if not x:
        return None
    good = bool(x.get("determined")) and bool(x.get("valid")) and fin(x.get("SF")) and fin(x.get("SF_err"))
    return {"SF": x.get("SF"), "err": x.get("SF_err"), "good": good, "N": x.get("N_sig")}


def crit_S(R, tags):
    rows, fl = [], []
    for run, cls in MAIN.items():
        for b in S_BINS:
            v = [sf_get(R, t, run, cls, b) for t in tags]
            h = sf_get(R, REF, run, cls, b)
            txt = " / ".join(pm(x["SF"], x["err"]) + ("" if x["good"] else " (nieokr.)") if x else "brak" for x in v)
            if any(x is None or not x["good"] for x in v):
                fl.append(False)
                rows.append(("S", f"{run} {cls} {b}: SF(D*) WP70 ziarna 1/2/3", txt, f"określony ×3 i rozstęp ≤ {S_RANGE}", NO, "SF nieokreślony (brak pomiaru / fit nieokreślony w ≥ 1 ziarnie)"
                             + (f"; drzewo H {pm(h['SF'], h['err'])}" if h else "")))
                continue
            rg = max(x["SF"] for x in v) - min(x["SF"] for x in v)
            ok = rg <= S_RANGE
            fl.append(ok)
            mean = sum(x["SF"] for x in v) / len(v)
            rows.append(("S", f"{run} {cls} {b}: SF(D*) WP70 ziarna 1/2/3", txt, f"określony ×3 i rozstęp ≤ {S_RANGE}", OK if ok else NO,
                         f"rozstęp {rg:.3f}; średnia {mean:.3f}" + (f"; drzewo H {pm(h['SF'], h['err'])}" if h else "")))
    for run, cls, b in (("427929", "MB427929", "3-5"), ("428770", "all", "3-5"), ("428770", "J12only", "5-8"), ("428770", "J12only", "8-30")):
        v = [sf_get(R, t, run, cls, b) for t in tags]
        h = sf_get(R, REF, run, cls, b)
        rows.append(("S-opis.", f"{run} {cls} {b}: SF(D*) WP70 ziarna 1/2/3", " / ".join(pm(x["SF"], x["err"]) + ("" if x["good"] else " (nieokr.)") if x else "brak" for x in v),
                     "–", DESC, f"drzewo H {pm(h['SF'], h['err'])}" if h else ""))
    if not dig(R, "species", "D0", "sf"):
        return None if R is None else False, rows
    return all(fl), rows


def crit_K2(R, tags):
    C = dig(R, "species", "D0", "cfg") or {}
    k = "sculpt|mc_SR|0.7"
    rows, fl = [], []
    for t in tags:
        s = dig(C, t, k)
        sr = dig(C, t, "ref_same_rows", k) or dig(C, REF, k)
        if not s or not sr or not fin(s.get("R")) or not fin(sr.get("R")):
            rows.append(("K2", f"{t}: |R_sculpt − 1| vs drzewo H", "–", f"nie większe niż H w {K2_NSIG:.0f}σ", UN, "brak wartości"))
            fl.append(None)
            continue
        ok = abs(s["R"] - 1) <= abs(sr["R"] - 1) + K2_NSIG * math.hypot(s.get("R_err") or 0, sr.get("R_err") or 0)
        fl.append(ok)
        note = ""
        if fin(s.get("slope")) and fin(sr.get("slope")):
            d = s["slope"] - sr["slope"]
            note = f"nachylenie ε_B(m)/100 MeV: sieć {s['slope']:+.3f}, H {sr['slope']:+.3f}, Δ {d:+.3f} ({'w' if abs(d) <= K2_SLOPE else 'POZA'} ±{K2_SLOPE})"
        for kk, lab in (("sculpt|data427_band|0.7", "pasmo 427929"), ("sculpt|data_band|0.7", "pasmo 428770")):
            x = dig(C, t, kk)
            if x and fin(x.get("R")):
                note += f"; R {lab} {pm(x['R'], x.get('R_err'))}"
        cm = dig(C, t, "corr_score_m")
        if fin(cm):
            note += f"; corr(score, m) {cm:+.4f}"
        rows.append(("K2", f"{t}: |R_sculpt − 1| (tło MC, okno, WP70) vs drzewo H", f"sieć {pm(s['R'], s.get('R_err'))}; H {pm(sr['R'], sr.get('R_err'))}",
                     f"nie większe niż H w {K2_NSIG:.0f}σ", OK if ok else NO, note))
    return (False if any(f is False for f in fl) else (None if any(f is None for f in fl) or not fl else True)), rows


def gain_rows(R, tags, full):
    C = dig(R, "species", "D0", "cfg") or {}
    rows, lab427 = [], []
    for nm, title in (("sb427929", "pasma 427929"), ("sb428770", "pasma 428770"), ("mc", "MC (tło truth ≥ 8)")):
        for b in BINS:
            g = [dig(C, t, f"{b}|0.7", "G_" + nm) for t in tags]
            ge = [dig(C, t, f"{b}|0.7", "G_" + nm + "_err") for t in tags]
            gf = [dig(C, t, f"{b}|0.7", "G_" + nm) for t in full]
            gv, gfv = [x for x in g if fin(x)], [x for x in gf if fin(x)]
            if not gv:
                rows.append(("G-proxy", f"{title} {b}: G = Z/Z_H WP70", "–", "–", DESC, "brak wartości"))
                continue
            m = sum(gv) / len(gv)
            mf = sum(gfv) / len(gfv) if gfv else None
            keep = (m - 1) / (mf - 1) if mf is not None and abs(mf - 1) > 1e-6 else None
            e1 = next((x for x in ge if fin(x)), None)
            rows.append(("G-proxy", f"{title} {b}: G = Z/Z_H WP70 (średnia ziaren)", pm(m, e1), "–", DESC,
                         f"ziarna {' / '.join(pm(x) for x in g)}; rozstęp {max(gv) - min(gv):.3f}; pełna sieć {pm(mf)}; zachowana część zysku f = {pm(keep, nd=2)}"
                         + (" — PROXY, zawyża (etap 4)" if nm != "mc" else "")))
            if nm == "sb427929":
                lab427.append(m)
    lab = None
    if len(lab427) == 3:
        if sum(x >= G_HAVE for x in lab427) >= G_NBIN:
            lab = f"proxy G(427929) ≥ {G_HAVE} w ≥ {G_NBIN} z 3 binów → „jest co mierzyć fitem” (NIE zysk realny)"
        elif sum(x < G_NONE for x in lab427) >= G_NBIN:
            lab = f"proxy G(427929) < {G_NONE} w ≥ {G_NBIN} z 3 binów → „po usunięciu pionu sieć nie daje zysku w minimum-bias nawet w proxy”"
        else:
            lab = "proxy G(427929) pośrednie (między etykietami §4)"
    return rows, lab


def gfit_rows(RR, tag):
    rows, out = [], {}
    for run, cls in MAIN.items():
        fl = []
        for b in BINS:
            x = dig(RR, "data", "D0", tag, f"{run}|{cls}", b, "0.7")
            c = dig(RR, "closure", "D0", tag, b, "0.7") or {}
            if not x or not fin(x.get("G_fit")):
                rows.append(("G_fit", f"{run} {cls} {b}", "–", f"≥ {D_GFIT}", UN, "brak")); fl.append(None); continue
            ok = x["G_fit"] >= D_GFIT
            fl.append(ok if x.get("determined") else None)
            rows.append(("G_fit", f"{run} {cls} {b}: [N_S/σ_N]_sieć / [N_S/σ_N]_drzewo", pm(x["G_fit"], nd=2), f"≥ {D_GFIT}", (OK if ok else NO) if x.get("determined") else UN,
                         f"R {pm(x.get('R'), x.get('R_err'))} (opisowo; domknięcie MC {'zaliczone' if c.get('pass') else 'NIEzaliczone'}: R_fit − R_prawda {pm(c.get('bias'))})"
                         + ("" if x.get("determined") else "; fit nieokreślony")))
        n_ok, n_un = sum(f is True for f in fl), sum(f is None for f in fl)
        out[run] = True if n_ok >= D_NBIN else (None if n_ok + n_un >= D_NBIN else False)
    return rows, out


# ------------------------------------------------------------------------------------------------ main
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=os.path.join(HERE, "runs", "lxplus6")); ap.add_argument("--out", default="")
    a = ap.parse_args(argv)
    d = os.path.abspath(a.dir)
    miss = []
    P, R, M = load(d, "dnet6_p1.json", miss), load(d, "dnet6_results.json", miss), load(d, "dnet6_mask.json", miss)
    RR = load(d, "dnet6_R.json", []) if os.path.isfile(os.path.join(d, "dnet6_R.json")) else None
    D = Doc()
    D.L.append("# WERDYKT DNET etap 6 — odczyt wg KRYTERIA_etap6.md")
    D.p(f"Katalog: `{d}`. Skrypt `DNET/verdict6.py` — sam odczyt i porównanie z progami zapisanymi przed wynikami; żadnych fitów.")
    D.h("0. Kompletność wejścia")
    metas = sorted(os.path.basename(f)[:-5] for f in glob.glob(os.path.join(d, "meta", "*.json")))
    have = [t for t in NP if t in metas]
    for f in ("dnet6_p1.json", "dnet6_results.json", "dnet6_mask.json", "dnet6_R.json", "DONE", "DONE_R"):
        D.p(f"- `{f}`: {'jest' if os.path.exists(os.path.join(d, f)) else '**BRAK**'}")
    D.p(f"- `meta/`: {len(have)} z 3 ziaren ({', '.join(have) or 'brak'})")
    dev = {}
    for t in have:
        try:
            m = json.load(open(os.path.join(d, "meta", t + ".json")))
            dev[t] = sorted({str(f.get("device")) for f in m.get("fit", [])})
            mk = m.get("mask")
        except Exception:  # noqa: BLE001
            dev[t], mk = ["?"], None
        D.p(f"  - `{t}`: urządzenie foldów {'/'.join(dev[t])}, maska {mk}")
    if len({tuple(v) for v in dev.values()}) > 1 or any(len(v) > 1 for v in dev.values()):
        D.p("- **UWAGA: ziarna/foldy liczone różnymi urządzeniami (GPU/CPU) — rozstęp ziaren zawiera różnicę urządzeń (§2).**")
    fails = sorted(os.path.basename(f) for f in glob.glob(os.path.join(d, "fail", "*")))
    D.p(f"- węzły w `fail/`: {', '.join(fails) if fails else 'brak'}")
    for m_ in miss:
        D.p(f"- **UWAGA: {m_}**")
    if M:
        mkey = f"dm<{M.get('DM_MAX')}|sign{M.get('SIGN')}"
        pr = dig(M, "pis_removed", mkey) or {}
        D.p(f"- maska (MC + pasma): Δm < {M.get('DM_MAX')} GeV, znak {M.get('SIGN')}; π_s usunięty {pm(pr.get('of_seen_e04'), nd=4)} widzianych w e04, zostaje widoczny w "
            f"{pm(pr.get('seen_e04_left_of_all_tag'), nd=4)} próbki tagu; tory usunięte: " + ", ".join(
                f"{lab} {pm(dig(M, k, mkey, 'frac_tracks_removed'), nd=4)}"
                for k, lab in (("loss_sig_incl", "sygnał inkl."), ("loss_bkg_mc", "tło MC"), ("loss_sb427929", "pasma 427929"), ("loss_sb428770", "pasma 428770"))))

    D.h("1. P1 — czy maska usuwa efekt wolnego pionu (§3; MC, WP70)")
    p1, rows, _ = crit_P1(P, NP[0]) if P else (None, [], [])
    D.table(rows)
    seeds = {NP[0]: p1}
    for t in NP[1:]:
        st_, rws, _ = crit_P1(P, t) if P else (None, [], [])
        seeds[t] = st_
        D.p(f"- opisowo `{t}`: P1 = {ST[st_]}; DR " + " / ".join(r[2] for r in rws))
    for t in FULL[:1]:
        _, rws, _ = crit_P1(P, t) if P else (None, [], [])
        D.p(f"- odniesienie bez maski `{t}` (score etapu 2): DR " + " / ".join(r[2] for r in rws))
    x = dig(P, "models", NP[0], "bins")
    if x:
        D.p("- opisowo ε_b (ziarno 1) wiersze z usuniętym torem / bez: " + "; ".join(
            f"{b}: {pm(*(dig(x, b, '0.7', 'b_rest_masked_rows') or [None, None])[:2])} / {pm(*(dig(x, b, '0.7', 'b_rest_unmasked_rows') or [None, None])[:2])}" for b in BINS))
    arch = all(v is True for v in seeds.values())
    D.p(f"- **P1 (ziarno 1): {ST[p1]}**" + ("" if p1 is not True else (" — wszystkie 3 ziarna spełniają (własność architektury)" if arch else " — nie wszystkie ziarna spełniają: „konkretne modele”")))

    D.h("2. S — SF z tagu D* (§4; WP70)")
    S, rows = crit_S(R, NP) if R else (None, [])
    D.table(rows)
    D.p(f"- **S: {ST[S]}**" + (" („SF nieokreślony / niestabilny”)" if S is False else ""))

    D.h("3. K2 — rzeźbienie (§4)")
    K2, rows = crit_K2(R, NP) if R else (None, [])
    D.table(rows)
    D.p(f"- **K2: {ST[K2]}**")

    D.h("4. Zysk — WYŁĄCZNIE proxy na pasmach (§4); nie jest zyskiem realnym")
    rows, lab = gain_rows(R, NP, FULL) if R else ([], None)
    D.table(rows)
    D.p("Zastrzeżenie wiążące: proxy G na pasmach zawyża (etap 4, D⁺: proxy 1,26–1,36 → G_fit 1,03–1,11 w 427929); ε_S pochodzi z MC.")
    if lab:
        D.p(f"- etykieta opisowa: {lab}")

    D.h("5. Zysk realny (G_fit; krok opcjonalny `app` + `anar`, §5)")
    gf = {}
    if RR:
        rows, gf = gfit_rows(RR, NP[0])
        D.table(rows)
        for run in MAIN:
            D.p(f"- {run}: G_fit ≥ {D_GFIT} w ≥ {D_NBIN} z 3 binów: **{ST[gf.get(run)]}**")
    else:
        D.p("`dnet6_R.json` brak → **zysk realny NIEZNANY** (krok `app` + `anar` nie policzony: bramka P1 zamknięta, jeszcze liczy albo nieudany). Liczby z §4 to proxy.")

    D.h("6. WERDYKT KOŃCOWY")
    closed = both(p1, S, K2)
    if closed is True:
        head = "SIEĆ DOMKNIĘTA jako skalibrowana alternatywa dla drzewa H (P1 ∧ S ∧ K2)"
    elif closed is False:
        why = [n for n, s in (("P1 — maska nie usuwa efektu", p1), ("S — SF nieokreślony / niestabilny", S), ("K2 — rzeźbi bardziej niż drzewo", K2)) if s is False]
        head = "SIEĆ NIE DOMKNIĘTA: " + "; ".join(why)
    else:
        head = "NIEROZSTRZYGNIĘTE (brak danych): " + ", ".join(n for n, s in (("P1", p1), ("S", S), ("K2", K2)) if s is None)
    D.p(f"**{head}**")
    D.p(f"- P1 {ST[p1]}, S {ST[S]}, K2 {ST[K2]}.")
    D.p("- Znaczenie „domknięta”: SF z tagu D* określony, stabilny między ziarnami i dotyczy próbki inkluzywnej. NIE oznacza realnego zysku ani zmiany nominału (drzewo H).")
    D.p("- Zysk: " + ("G_fit — " + ", ".join(f"{r} {ST[v]}" for r, v in gf.items()) if gf else "realny nieznany; proxy z §4 z zastrzeżeniem o zawyżeniu."))
    if closed is False:
        D.p("- Wynik negatywny → notatka w `charm-bkg-notes` (podejścia odrzucone).")
    txt = "\n".join(D.L) + "\n"
    out = a.out or os.path.join(d, "WERDYKT_etap6.md")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    open(out, "w").write(txt)
    print(txt)
    print(f"zapisano {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
