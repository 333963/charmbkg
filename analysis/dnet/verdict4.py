#!/usr/bin/env python3
"""DNET etap 4 — JEDNORAZOWY ODCZYT werdyktu wg DNET/KRYTERIA_etap4.md (progi wiążące, zapisane przed wynikami).

  ~/Documents/charm/.venv_h5/bin/python DNET/verdict4.py [--dir DNET/runs/lxplus4] [--out PLIK.md]

Czyta: dnet4_R.json (dn4_ana.py r), dnet4_results.json i dnet4_ref2_results.json (analyze.py compute --sf-all), dnet4_d23.json (dn4_ana.py d23),
dnet4_pis.json (dn4_pis.py), fail/, meta/, DONE, DONE23. Wypisuje i zapisuje <dir>/WERDYKT_etap4.md.
ŻADNYCH fitów: tylko odczyt i porównanie z progami (arytmetyka: sumy kwadratowe błędów, rozstęp, χ² stałej = średnia ważona).
Brak pliku / pola / fit nieokreślony → NIEROZSTRZYGNIĘTE (nigdy „spełnione”). Każdy próg niżej ma odsyłacz do paragrafu KRYTERIA_etap4.md."""
import argparse
import glob
import json
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

# ------------------------------------------------------------------------------------------------ progi (KRYTERIA_etap4.md)
A_NSIG, A_MAXBIAS, A_MAXERR = 2.5, 0.04, 0.05      # §1.A: |R_fit − R_prawda| ≤ 2,5 σ_stat i ≤ 0,04; σ_stat ≤ 0,05; fit określony
A_GLOB = 3.0                                       # §1.A: |średni pull| po wszystkich punktach WP70 ≤ 3/√n
A_TRANSFER = 0.03                                  # §1.A: |R_prawda − 1| ≤ 0,03; przekroczenie nie blokuje, wchodzi do σ_tot binu
B_SIGTOT = {"3-5": 0.08, "5-8": 0.05, "8-30": 0.05}  # §1.B: σ_tot(R)
B_SIGBKG = 0.03                                    # §1.B: σ_tło ≤ 0,03, większa → „metoda ograniczona rzeźbieniem tła”
C1_RANGE = {"3-5": 0.06, "5-8": 0.04, "8-30": 0.04}  # §1.C1: rozstęp R między ziarnami 1–3, obie klasy główne
C1_DSTAR = 0.04                                    # §1.C1: rozstęp SF(D*) WP70 8–30 GeV 428770 all trzech ziaren
C2_NSIG = 2.0                                      # §1.C2: |R(J12only) − R(notJ12only)| ≤ 2σ (σ_stat kwadratowo), każdy bin
C3_PMIN = 0.01                                     # §1.C3: χ² stałej dla trzech binów (σ_tot), p ≥ 0,01
C4_NSIG = 2.0                                      # §1.C4: R wobec SF_sieć(D*)/SF_H(D*) w 2σ, biny 5–8 i 8–30 (tylko D0; §9: tylko gdy P1 niespełnione)
MULT_NMAX, MULT_NSIG = 1, 3.0                      # §7: C2 i C4 niespełnione dopiero przy > 1 przekroczeniu na 3 biny albo jednym > 3σ
D_GFIT, D_GSB, D_NBIN = 1.15, 1.20, 2              # §1.D: G_fit ≥ 1,15 w ≥ 2 z 3 binów; kontrolnie G(pasma tego runu) ≥ 1,20 w ≥ 2 z 3
K2_NSIG, K2_SLOPE = 2.0, 0.10                      # RAPORT.md §4 K2: |R_sculpt − 1| nie większe niż dla H w 2σ; nachylenie ±0,10/100 MeV od H
SB_DG = 0.05                                       # §4: sb2mc „lepszy niż mc”: G_fit większe o ≥ 0,05 w ≥ 2 z 3 binów w obu runach przy spełnionym B
D23_Z, D23_ZALT, D23_GAIN, D23_ZREF = 5.0, 4.0, 1.15, 3.1   # §5: Z ≥ 5 nominalnie i ≥ 4 dla obu sąsiednich rzędów tła; D+: Z_sieć/Z_drzewo ≥ 1,15
P1_CA, P1_DR, P1_NSIG, P1_NBIN = 0.03, 0.03, 2.0, 2  # §9 P1: |ε_c/ε_a − 1| > 0,03 w ≥ 1 binie albo podwójny stosunek ≠ 1 o > 0,03 i > 2σ w ≥ 2 z 3
NBOOT_REQ = 100                                    # §1: bootstrap Poissona po zdarzeniach, 100 replik

BINS = ["3-5", "5-8", "8-30"]
MAIN = {"427929": "MB427929", "428770": "all"}     # klasy główne (§1.B)
REF = {"D0": "D0.xgb.H.sb.n.s1", "Dplus": "Dplus.xgb.H.sb.n.s1"}
SPN = {"D0": "D⁰", "Dplus": "D⁺"}
OK, NO, UN, DESC = "SPEŁNIONE", "NIESPEŁNIONE", "NIEROZSTRZYGNIĘTE", "opisowo"
FILES = ["dnet4_R.json", "dnet4_results.json", "dnet4_ref2_results.json", "dnet4_d23.json", "dnet4_pis.json"]


# ------------------------------------------------------------------------------------------------ narzędzia
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


def st3(flags):
    """Lista True/False/None (None = brak) → status „wszystkie muszą spełniać”."""
    if any(f is False for f in flags):
        return NO
    return UN if (not flags or any(f is None for f in flags)) else OK


def nof(flags, need):
    """„≥ need z len(flags)”: True/False/None → status."""
    n_ok, n_un = sum(f is True for f in flags), sum(f is None for f in flags)
    if n_ok >= need:
        return OK
    return UN if n_ok + n_un >= need else NO


def mult(pulls, nsig):
    """Reguła §7 dla testów „2σ” (C2, C4): niespełnione przy > 1 przekroczeniu albo jednym > 3σ. pulls: lista |pull| albo None."""
    have = [p for p in pulls if p is not None]
    nex = sum(p > nsig for p in have)
    if nex > MULT_NMAX or any(p > MULT_NSIG for p in have):
        return NO
    return UN if len(have) < len(pulls) or not pulls else OK


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

    def text(self):
        return "\n".join(self.L) + "\n"


class Ctx:
    """Wczytane pliki + dostęp odporny na braki."""

    def __init__(self, d):
        self.dir, self.miss, self.J = d, [], {}
        for f in FILES:
            p = os.path.join(d, f)
            try:
                self.J[f] = json.load(open(p))
                if not isinstance(self.J[f], dict):
                    raise ValueError("nie słownik")
            except Exception as e:  # noqa: BLE001
                self.J[f] = None
                self.miss.append(f"{f}: {'brak pliku' if not os.path.isfile(p) else 'nieczytelny (' + type(e).__name__ + ')'}")
        self.R, self.RES, self.REF2, self.D23, self.PIS = (self.J[f] for f in FILES)

    # --- metoda R
    def tags(self, sp, kind="mc"):
        """Sieci gatunku: kind 'mc' → {ziarno: tag}, 'sb2mc' → tag|None. Źródła: dnet4_R.json, dnet4_results.json; domyślnie nazwy z KRYTERIA §2–§4."""
        names = set(dig(self.R, "tags", sp, "nets") or []) | set(dig(self.R, "closure", sp) or {}) | set(dig(self.RES, "species", sp, "cfg") or {})
        if kind == "sb2mc":
            t = sorted(x for x in names if re.match(rf"^{sp}\.tr\.H\+e04\.sb2mc\.d\d+\.s1$", x))
            return t[0] if t else (f"{sp}.tr.H+e04.sb2mc.d5.s1" if sp == "D0" else None)
        o = {}
        for x in names:
            m = re.match(rf"^{sp}\.tr\.H\+e04\.mc\.d\d+\.s([123])$", x)
            if m:
                o[int(m.group(1))] = x
        for s in (1, 2, 3):
            o.setdefault(s, f"{sp}.tr.H+e04.mc.d5.s{s}")
        return o

    def clo(self, sp, tag, b):
        return dig(self.R, "closure", sp, tag, b, "0.7")

    def transfer(self, sp, tag, b):
        """|R_prawda − 1| (kontrola przenoszenia ε fd_nd → mb_nd, §1.A) albo None."""
        o = self.clo(sp, tag, b)
        return abs(o["R_truth"] - 1) if o and fin(o.get("R_truth")) else None

    def rd(self, sp, tag, run, cls, b, wp="0.7"):
        """Pomiar R na danych: dict {R, stat, tot, bkg, G, o} albo None (brak / zablokowany / fit nieokreślony)."""
        o = dig(self.R, "data", sp, tag, f"{run}|{cls}", b, wp)
        if not o or not o.get("determined") or not fin(o.get("R")) or not fin(o.get("R_err")):
            return None
        tot = o.get("R_err_tot")
        tr = self.transfer(sp, tag, b)
        extra = tr if (tr is not None and tr > A_TRANSFER) else 0.0      # §1.A: odchylenie przenoszenia wchodzi do σ_tot (kod dn4_ana tego NIE robi)
        tot = math.hypot(tot, extra) if fin(tot) else None
        return {"R": o["R"], "stat": o["R_err"], "tot": tot, "bkg": o.get("syst_bkg"), "G": o.get("G_fit"), "o": o, "extra": extra}

    def why(self, sp, tag, run, cls, b):
        """Powód braku pomiaru R."""
        if self.R is None:
            return "brak dnet4_R.json"
        if f"{sp}|{tag}|{run}|{cls}|{b}|0.7" in (self.R.get("blocked") or []):
            return "zablokowane bramką domknięcia"
        o = dig(self.R, "data", sp, tag, f"{run}|{cls}", b, "0.7")
        if o is None:
            return "brak wpisu w JSON"
        return "fit nieokreślony" + (" (dmu/sc na granicy)" if dig(o, "fit", "at_limit") else "") + (f", warianty tła OK: {o.get('n_var_ok')}/4" if "n_var_ok" in o else "")

    # --- analyze.py
    def sf(self, tag, run, cls, b, res=None):
        x = dig(self.RES if res is None else res, "species", "D0", "sf", f"{tag}|{run}|{cls}|{b}|0.7")
        if not x or not x.get("determined") or not fin(x.get("SF")) or not fin(x.get("SF_err")):
            return None
        return x["SF"], x["SF_err"]

    def cfg(self, sp, tag):
        return dig(self.RES, "species", sp, "cfg", tag)


# ------------------------------------------------------------------------------------------------ kryteria (każda funkcja → (status, wiersze))
def crit_A(C, sp, tag):
    rows, fl = [], []
    for b in BINS:
        o = C.clo(sp, tag, b)
        if o is None:
            rows.append(("A", f"{b}: R_fit − R_prawda", "–", f"≤ {A_NSIG}σ i ≤ {A_MAXBIAS}; σ ≤ {A_MAXERR}", UN, "brak wpisu domknięcia")); fl.append(None); continue
        if not (fin(o.get("R")) and fin(o.get("R_err")) and fin(o.get("R_truth"))):
            rows.append(("A", f"{b}: R_fit − R_prawda", "–", f"≤ {A_NSIG}σ i ≤ {A_MAXBIAS}; σ ≤ {A_MAXERR}", NO, o.get("why") or "fit nieokreślony / brak R")); fl.append(False); continue
        bias, e = o["R"] - o["R_truth"], o["R_err"]
        ok = bool(o.get("determined")) and abs(bias) <= A_NSIG * e and abs(bias) <= A_MAXBIAS and e <= A_MAXERR
        note = f"R_fit {pm(o['R'], e)}, R_prawda {pm(o['R_truth'], o.get('R_truth_err'))}, pull {bias / e:+.2f}" if e > 0 else ""
        if not o.get("determined"):
            note += "; fit NIEOKREŚLONY"
        if bool(o.get("pass")) != ok:
            note += f"; UWAGA: pole pass w JSON = {o.get('pass')} ≠ odczyt"
        rows.append(("A", f"{b}: R_fit − R_prawda", pm(bias, e, 4), f"≤ {A_NSIG}σ i ≤ {A_MAXBIAS}; σ ≤ {A_MAXERR}", OK if ok else NO, note)); fl.append(ok)
        tr = abs(o["R_truth"] - 1)
        rows.append(("A-przen.", f"{b}: |R_prawda − 1|", pm(tr, o.get("R_truth_err"), 4), f"≤ {A_TRANSFER} (nie blokuje)", OK if tr <= A_TRANSFER else NO,
                     "" if tr <= A_TRANSFER else "dodane kwadratowo do σ_tot tego binu (w tym skrypcie; JSON tego nie zawiera)"))
    return st3(fl), rows, fl


def crit_Aglob(C):
    pulls = []
    for sp, v in (dig(C.R, "closure") or {}).items():
        for t, w in v.items():
            for b, x in w.items():
                o = x.get("0.7") or {}
                if fin(o.get("R")) and fin(o.get("R_truth")) and fin(o.get("R_err")) and o["R_err"] > 0:
                    pulls.append((o["R"] - o["R_truth"]) / o["R_err"])
    if not pulls:
        return UN, [("A-glob", "średni pull domknięć WP70", "–", f"|śr.| ≤ {A_GLOB}/√n", UN, "brak punktów domknięcia")]
    mean, lim = sum(pulls) / len(pulls), A_GLOB / math.sqrt(len(pulls))
    s = OK if abs(mean) <= lim else NO
    return s, [("A-glob", f"średni pull domknięć WP70 (n = {len(pulls)}; oczekiwane 21)", f"{mean:+.3f}", f"|śr.| ≤ {lim:.3f}", s,
                "wszystkie gatunki i modele razem (§1.A „dodatkowo globalnie”)")]


def crit_B(C, sp, tag, run):
    cls, rows, fl = MAIN[run], [], {}
    for b in BINS:
        x = C.rd(sp, tag, run, cls, b)
        if x is None or x["tot"] is None:
            rows.append(("B", f"{run}:{cls} {b}: σ_tot(R)", "–", f"≤ {B_SIGTOT[b]}", UN, C.why(sp, tag, run, cls, b))); fl[b] = None; continue
        ok = x["tot"] <= B_SIGTOT[b]
        o = x["o"]
        note = f"R = {pm(x['R'], x['stat'])} (stat); tło {pm(o.get('syst_bkg'))}, kształt {pm(o.get('syst_shape'))}, MC {pm(o.get('syst_mcstat'))}"
        if x["extra"]:
            note += f", przenoszenie ε {x['extra']:.3f} (dodane)"
        rows.append(("B", f"{run}:{cls} {b}: σ_tot(R)", pm(x["tot"]), f"≤ {B_SIGTOT[b]}", OK if ok else NO, note)); fl[b] = ok
        if fin(x["bkg"]):
            rows.append(("B-tło", f"{run}:{cls} {b}: σ_tło", pm(x["bkg"]), f"≤ {B_SIGBKG}", OK if x["bkg"] <= B_SIGBKG else NO,
                         ("" if x["bkg"] <= B_SIGBKG else "„metoda ograniczona rzeźbieniem tła”; ") + ("σ_tło > σ_stat" if x["bkg"] > x["stat"] else "σ_tło ≤ σ_stat")))
        if sp == "D0":
            s = C.sf(tag, run, cls, b)
            rows.append(("B-moc", f"{run}:{cls} {b}: σ_tot(R) < σ(SF sieci z tagu D*)", f"{pm(x['tot'])} vs {pm(s[1]) if s else '–'}", "σ_tot(R) mniejsze", DESC,
                         ("większa moc niż tag D*" if x["tot"] < s[1] else "NIE większa moc niż tag D*") if s else "SF z tagu D* nieokreślony / brak (3–5 GeV: niemierzalny)"))
    if fl["5-8"] is False or fl["8-30"] is False:
        s = NO
    elif fl["5-8"] is None or fl["8-30"] is None:
        s = UN
    elif fl["3-5"] is False:
        s = "PT5"                                   # §1.E: B zawodzi tylko w 3–5 → „przyjęta dla pT ≥ 5 GeV”
    else:
        s = UN if fl["3-5"] is None else OK
    return s, rows


def crit_C1(C, sp, p1):
    """Rozstęp R między ziarnami (6 punktów na gatunek) + (D0) rozstęp SF z tagu D* 8–30 GeV 428770 all."""
    seeds, rows, fl = C.tags(sp), [], []
    for run, cls in MAIN.items():
        for b in BINS:
            v = {s: C.rd(sp, seeds[s], run, cls, b) for s in (1, 2, 3)}
            have = {s: x["R"] for s, x in v.items() if x}
            rng = max(have.values()) - min(have.values()) if len(have) >= 2 else None
            txt = ", ".join(f"s{s} {pm(v[s]['R'], v[s]['stat']) if v[s] else '–'}" for s in (1, 2, 3))
            if len(have) == 3:
                ok = rng <= C1_RANGE[b]
            else:
                ok = False if (rng is not None and rng > C1_RANGE[b]) else None
                txt += "; brak: " + ", ".join(f"s{s} ({C.why(sp, seeds[s], run, cls, b)})" for s in (1, 2, 3) if not v[s])
            rows.append(("C1", f"{run}:{cls} {b}: rozstęp R ziaren 1–3", pm(rng), f"≤ {C1_RANGE[b]}", {True: OK, False: NO, None: UN}[ok], txt)); fl.append(ok)
    s_r = st3(fl)
    s_d = None
    if sp == "D0":
        for res, lab in ((C.RES, "modele etapu 4"), (C.REF2, "score etapu 2")):
            v = {s: C.sf(seeds[s], "428770", "all", "8-30", res) for s in (1, 2, 3)}
            have = [x[0] for x in v.values() if x]
            rng = max(have) - min(have) if len(have) >= 2 else None
            ok = (rng <= C1_DSTAR) if len(have) == 3 else (False if rng is not None and rng > C1_DSTAR else None)
            stt = {True: OK, False: NO, None: UN}[ok]
            if lab == "modele etapu 4":
                s_d = stt
            rows.append(("C1-D*", f"428770:all 8–30: rozstęp SF(D*) WP70 ziaren 1–3, {lab}", pm(rng), f"≤ {C1_DSTAR}", stt if lab == "modele etapu 4" else f"{stt} ({DESC})",
                         ", ".join(f"s{s} {pm(*v[s]) if v[s] else '– (fit nieokreślony/brak)'}" for s in (1, 2, 3))
                         + ("; §9(iii–iv): P1 spełnione → nie wchodzi do werdyktu E" if p1 is True else "")))
    return s_r, s_d, rows


def crit_C2(C, sp, tag):
    rows, pulls = [], []
    for b in BINS:
        a, n = C.rd(sp, tag, "428770", "J12only", b), C.rd(sp, tag, "428770", "notJ12only", b)
        if not a or not n:
            rows.append(("C2", f"428770 {b}: R(J12only) − R(notJ12only)", "–", f"≤ {C2_NSIG}σ", UN,
                         "; ".join(f"{c}: {C.why(sp, tag, '428770', c, b)}" for c, x in (("J12only", a), ("notJ12only", n)) if not x))); pulls.append(None); continue
        d, e = a["R"] - n["R"], math.hypot(a["stat"], n["stat"])
        pl = abs(d) / e if e > 0 else None
        pulls.append(pl)
        rows.append(("C2", f"428770 {b}: R(J12only) − R(notJ12only)", pm(d, e), f"≤ {C2_NSIG}σ", UN if pl is None else (OK if pl <= C2_NSIG else NO),
                     f"{pl:.2f}σ; J12only {pm(a['R'], a['stat'])}, notJ12only {pm(n['R'], n['stat'])}" if pl is not None else ""))
    s = mult(pulls, C2_NSIG)
    rows.append(("C2", "łącznie (§7: niespełnione przy > 1 przekroczeniu na 3 biny albo jednym > 3σ)", f"przekroczeń 2σ: {sum(p is not None and p > C2_NSIG for p in pulls)}/3", "≤ 1 i żadne > 3σ", s,
                 "niespełnione → SF per klasa triggera („warunkowo”)"))
    return s, rows


def crit_C3(C, sp, tag, run):
    cls = MAIN[run]
    v = [C.rd(sp, tag, run, cls, b) for b in BINS]
    if any(x is None or not x["tot"] for x in v):
        return UN, [("C3", f"{run}:{cls}: χ² stałej, 3 biny pT (σ_tot)", "–", f"p ≥ {C3_PMIN}", UN, "brak pomiaru R w: " + ", ".join(b for b, x in zip(BINS, v) if x is None or not x["tot"]))]
    w = [1 / x["tot"] ** 2 for x in v]
    mean = sum(wi * x["R"] for wi, x in zip(w, v)) / sum(w)
    chi2 = sum(wi * (x["R"] - mean) ** 2 for wi, x in zip(w, v))
    p = math.exp(-chi2 / 2)                         # ndf = 2 → p = exp(−χ²/2) dokładnie
    s = OK if p >= C3_PMIN else NO
    return s, [("C3", f"{run}:{cls}: χ² stałej, 3 biny pT (σ_tot)", f"χ² = {chi2:.2f}/2, p = {p:.3f}", f"p ≥ {C3_PMIN}", s,
                f"średnia ważona R = {pm(mean, 1 / math.sqrt(sum(w)))}; " + ", ".join(f"{b} {pm(x['R'], x['tot'])}" for b, x in zip(BINS, v))
                + ("" if s == OK else " → SF per bin pT (dopuszczalne, jeśli każdy bin spełnia B)"))]


def crit_C4(C, tag, run):
    cls, rows, pulls = MAIN[run], [], []
    for b in ("5-8", "8-30"):
        x, sn, sh = C.rd("D0", tag, run, cls, b), C.sf(tag, run, cls, b), C.sf(REF["D0"], run, cls, b)
        if not x or not x["tot"] or not sn or not sh:
            why = [w for w, c in ((f"R: {C.why('D0', tag, run, cls, b)}", not x or not x["tot"]), ("SF sieci (D*) nieokreślony/brak", not sn), ("SF_H (D*) nieokreślony/brak", not sh)) if c]
            rows.append(("C4", f"{run}:{cls} {b}: R − SF_sieć(D*)/SF_H(D*)", "–", f"≤ {C4_NSIG}σ", UN, "; ".join(why))); pulls.append(None); continue
        q = sn[0] / sh[0]
        qe = q * math.hypot(sn[1] / sn[0], sh[1] / sh[0])
        e = math.hypot(x["tot"], qe)
        pl = abs(x["R"] - q) / e
        pulls.append(pl)
        rows.append(("C4", f"{run}:{cls} {b}: R − SF_sieć(D*)/SF_H(D*)", pm(x["R"] - q, e), f"≤ {C4_NSIG}σ", OK if pl <= C4_NSIG else NO,
                     f"{pl:.2f}σ (σ = σ_tot(R) ⊕ σ_ilorazu); R {pm(x['R'], x['tot'])}, iloraz {pm(q, qe)}; z samym σ_ilorazu: {abs(x['R'] - q) / qe:.2f}σ"))
    s = mult(pulls, C4_NSIG)
    rows.append(("C4", f"{run}:{cls}: łącznie (§7: niespełnione przy > 1 przekroczeniu albo jednym > 3σ)", f"przekroczeń 2σ: {sum(p is not None and p > C4_NSIG for p in pulls)}/2", "≤ 1 i żadne > 3σ", s, ""))
    return s, rows


def crit_D(C, sp, tag, run):
    cls, rows, fl, flc = MAIN[run], [], [], []
    cf = C.cfg(sp, tag)
    for b in BINS:
        x = C.rd(sp, tag, run, cls, b)
        g = x["G"] if x and fin(x["G"]) else None
        fl.append(None if g is None else g >= D_GFIT)
        note = ""
        if g is not None:
            o = x["o"]
            note = f"σ_N/N sieć {pm(o.get('relerr_net'))}, drzewo {pm(o.get('relerr_tree'))}"
            if g > 1 and x["tot"]:
                note += f"; próg opłacalności: błąd stat. wzgl. drzewa x > {x['tot'] / math.sqrt(1 - 1 / g ** 2):.3f}"
            elif x["tot"]:
                note += "; próg opłacalności: nie istnieje (G_fit ≤ 1)"
        else:
            note = C.why(sp, tag, run, cls, b)
        rows.append(("D", f"{run}:{cls} {b}: G_fit", pm(g), f"≥ {D_GFIT}", {True: OK, False: NO, None: UN}[fl[-1]], note))
        r = dig(cf, f"{b}|0.7") or {}
        gs = r.get(f"G_sb{run}")
        flc.append(gs >= D_GSB if fin(gs) else None)
        rows.append(("D-kontr.", f"{run} {b}: G na dalekich pasmach tego runu (analyze.py)", pm(gs, r.get(f"G_sb{run}_err")), f"≥ {D_GSB}", {True: OK, False: NO, None: UN}[flc[-1]],
                     f"G_mc {pm(r.get('G_mc'), r.get('G_mc_err'))}" if r else "brak w dnet4_results.json"))
    s, sc = nof(fl, D_NBIN), nof(flc, D_NBIN)
    rows.append(("D", f"{run}:{cls}: łącznie G_fit ≥ {D_GFIT}", f"{sum(f is True for f in fl)}/3 binów", f"≥ {D_NBIN} z 3", s, ""))
    rows.append(("D-kontr.", f"{run}: łącznie G(pasma) ≥ {D_GSB}", f"{sum(f is True for f in flc)}/3 binów", f"≥ {D_NBIN} z 3", sc, "kontrola (§1.D „kontrolnie”)"))
    return s, sc, rows


def crit_K2(C, sp, tag):
    cf, ref = C.cfg(sp, tag), C.cfg(sp, dig(C.RES, "species", sp, "ref") or REF[sp])
    k = "sculpt|mc_SR|0.7"
    s = dig(cf, k)
    same = dig(cf, "ref_same_rows", k)
    sr = same or dig(ref, k)
    if not s or not sr:
        return UN, [("K2", "|R_sculpt − 1| (tło MC truth ≥ 8, okno, WP70) vs drzewo H", "–", f"nie większe niż H w {K2_NSIG:.0f}σ", UN,
                     "brak dnet4_results.json" if C.RES is None else f"brak pola {k} (sieć: {'jest' if s else 'brak'}, H: {'jest' if sr else 'brak'})")], None
    slope_ok = None
    ok = abs(s["R"] - 1) <= abs(sr["R"] - 1) + K2_NSIG * math.hypot(s["R_err"], sr["R_err"])
    rows = [("K2", "|R_sculpt − 1| (tło MC truth ≥ 8, okno, WP70) vs drzewo H", f"sieć {pm(s['R'], s['R_err'])}; H {pm(sr['R'], sr['R_err'])}", f"nie większe niż H w {K2_NSIG:.0f}σ",
             OK if ok else NO, "H na tych samych wierszach" if same else "H na wszystkich wierszach (jak analyze.py report)")]
    if fin(s.get("slope")) and fin(sr.get("slope")):
        d = s["slope"] - sr["slope"]
        slope_ok = abs(d) <= K2_SLOPE
        rows.append(("K2-kontr.", "nachylenie ε_B(m) /100 MeV: sieć − H", f"{d:+.3f}", f"|Δ| ≤ {K2_SLOPE}", OK if abs(d) <= K2_SLOPE else NO, f"sieć {s['slope']:+.3f}, H {sr['slope']:+.3f} (kontrola K2, bramka G3)"))
    for kk, lab in (("sculpt|data_band|0.7", "428770"), ("sculpt|data427_band|0.7", "427929")):
        x = dig(cf, kk)
        if x:
            rows.append(("K2-kontr.", f"R_sculpt pasmo danych {lab} [okno, dalekie SB)", pm(x["R"], x["R_err"]), "–", DESC, ""))
    return (OK if ok else NO), rows, slope_ok


def crit_P1(C, tag):
    tree = REF["D0"]
    M = dig(C.PIS, "models")
    if not M or tag not in M or tree not in M:
        why = "brak dnet4_pis.json" if C.PIS is None else f"brak modelu w pliku: {', '.join(t for t in (tag, tree) if t not in (M or {}))}"
        return None, [("P1", "diagnostyka wolnego pionu D* (MC)", "–", f"|ε_c/ε_a − 1| > {P1_CA} w ≥ 1 binie albo podw. stosunek", UN, why)]
    rows, c1, c2 = [], [], []
    for b in BINS:
        x, t = dig(M, tag, b, "0.7"), dig(M, tree, b, "0.7")
        ca = x.get("c_over_a") if x else None
        c1.append(abs(ca - 1) > P1_CA if fin(ca) else None)
        rows.append(("P1", f"{b}: ε_c/ε_a (maska wolnego pionu), sieć", pm(ca, None, 4), f"|… − 1| > {P1_CA} → efekt", UN if c1[-1] is None else ("EFEKT" if c1[-1] else "bez efektu"),
                     f"ε_a {pm(*(x['a_tag'][:2]))}, ε_c {pm(*(x['c_masked'][:2]))}, ε_b {pm(*(x['b_rest'][:2])) if x.get('b_rest') else '–'}" if x and x.get("a_tag") and x.get("c_masked") else "bin pominięty w dn4_pis (za mało wierszy)"))
        ab, tb = (x or {}).get("a_over_b"), (t or {}).get("a_over_b")
        if ab and tb and ab[0] > 0 and tb[0] > 0:
            dr = ab[0] / tb[0]
            dre = dr * math.hypot(ab[1] / ab[0], tb[1] / tb[0])
            f = abs(dr - 1) > P1_DR and abs(dr - 1) > P1_NSIG * dre
            c2.append(f)
            rows.append(("P1", f"{b}: [ε_a/ε_b]_sieć / [ε_a/ε_b]_drzewo", pm(dr, dre, 4), f"|… − 1| > {P1_DR} i > {P1_NSIG:.0f}σ → efekt", "EFEKT" if f else "bez efektu",
                         f"sieć {pm(*ab)}, drzewo {pm(*tb)}; {abs(dr - 1) / dre:.1f}σ (błędy dwumianowe MC, kwadratowo)"))
        else:
            c2.append(None)
            rows.append(("P1", f"{b}: [ε_a/ε_b]_sieć / [ε_a/ε_b]_drzewo", "–", f"|… − 1| > {P1_DR} i > {P1_NSIG:.0f}σ → efekt", UN, "brak a_over_b"))
    if any(f is True for f in c1) or sum(f is True for f in c2) >= P1_NBIN:
        p1 = True
    elif any(f is None for f in c1) or sum(f is True for f in c2) + sum(f is None for f in c2) >= P1_NBIN:
        p1 = None
    else:
        p1 = False
    rows.append(("P1", "łącznie: „efekt istotny”", f"(c/a): {sum(f is True for f in c1)}/3 binów; podw. stosunek: {sum(f is True for f in c2)}/3", f"≥ 1 bin albo ≥ {P1_NBIN} z 3",
                 {True: OK, False: NO, None: UN}[p1],
                 {True: "P1 SPEŁNIONE → R = kalibracja podstawowa; C4 NIE jest wymaganiem (opisowo); B/C1 w werdykcie E tylko z R (§9)",
                  False: "P1 niespełnione → obie kalibracje równoprawne, C4 OBOWIĄZUJE (§9)", None: "nie wiadomo, czy C4 obowiązuje"}[p1]))
    return p1, rows


def crit_d23(C, sp, run):
    """→ (status „widoczne” sieci, Z sieci, Z drzewa, wiersze)."""
    cls = MAIN[run]
    O = dig(C.D23, "species", sp)
    T = dig(O, "runs", f"{run}|{cls}", "tags")
    if not T:
        why = "brak dnet4_d23.json" if C.D23 is None else (dig(O, "why") or f"brak {sp} {run}|{cls} w pliku")
        return UN, None, None, [("2–3", f"{SPN[sp]} {run}:{cls}", "–", f"Z ≥ {D23_Z:.0f}; Z(tło ±1) ≥ {D23_ZALT:.0f}", UN, why)]
    rows, res = [], {}
    for t, o in T.items():
        kind = "drzewo" if ".xgb." in t else "sieć"
        if "wp" not in o:
            rows.append(("2–3", f"{SPN[sp]} {run}:{cls} {kind} `{t}`", "–", "", UN, o.get("why") or "brak wyboru WP")); res[kind] = (UN, None); continue
        f = dig(o, "fits", str(o["wp"])) or {}
        nom = f.get("nom") or {}
        alt = {k: v for k, v in f.items() if k.startswith("exp_pol") and isinstance(v, dict)}
        z = nom.get("Z") if nom.get("valid") else None
        za = [v.get("Z") if v.get("valid") else None for v in alt.values()]
        if not fin(z) or len(za) < 2 or any(not fin(x) for x in za):
            s = NO if (fin(z) and z < D23_Z) or any(fin(x) and x < D23_ZALT for x in za) else UN
        else:
            s = OK if z >= D23_Z and min(za) >= D23_ZALT else NO
        res[kind] = (s, z)
        nb = ", ".join(f"WP{int(round(float(w) * 100))}: y = N/ε {pm(dig(v, 'y'), None, 0)}, Z {pm(dig(v, 'nom', 'Z'), None, 1)}" for w, v in sorted((o.get("fits") or {}).items(), reverse=True))
        rows.append(("2–3", f"{SPN[sp]} {run}:{cls} {kind} `{t}` WP{int(round(o['wp'] * 100))} ({o.get('rule')})",
                     f"Z = {pm(z, None, 1)}; tło rząd −1/+1: " + " / ".join(pm(x, None, 1) for x in za) + f"; N = {pm(nom.get('N'), nom.get('N_err'), 0)}",
                     f"Z ≥ {D23_Z:.0f}; Z(tło ±1) ≥ {D23_ZALT:.0f}", s if (sp == "D0" and kind == "sieć") else f"{s} ({DESC})",
                     f"S/B {pm(f.get('S_over_B'))}; fit nom. {'ważny' if nom.get('valid') else 'NIEWAŻNY'}; bootstrap {nom.get('boot_ok', 'brak')}; stabilność: {nb}; bez klasyfikatora Z ≤ {D23_ZREF} (D23 §5)"))
    return res.get("sieć", (UN, None))[0], res.get("sieć", (UN, None))[1], res.get("drzewo", (UN, None))[1], rows


# ------------------------------------------------------------------------------------------------ werdykt E
def verdict_E(sp, run, A, nA_fail, Ag, B, C1, C2, C3, C4, D, K2, c4_req):
    """Reguła §1.E + §6 (+ §9 dla C4). Zwraca (nagłówek, lista uzasadnień)."""
    why = []
    hard_no, hard_un = [], []
    for nm, s in (("A (domknięcie MC, ziarno 1)", A), ("A-glob (średni pull)", Ag), ("B (precyzja R)", B), ("D (G_fit)", D), ("K2 (rzeźbienie)", K2)) + ((("C4 (zgodność z tagiem D*)", C4),) if c4_req is not False else ()):
        if s == NO:
            hard_no.append(nm)
        elif s == UN:
            hard_un.append(nm)
    if c4_req is None and C4 == OK:
        hard_un = [x for x in hard_un if not x.startswith("C4")]
    if hard_no:
        if A == NO and nA_fail >= 2:
            why.append("§6: A niezaliczone w ≥ 2 z 3 binów → metoda R nie działa; brak nowej kalibracji, sieć zostaje kandydatem warunkowym jak po etapie 3")
        if D == NO:
            why.append("§6: G_fit < 1,15 → zysk z pasm nie przenosi się na pik; sieć nie zastępuje drzewa")
        if B == NO:
            why.append("§6: B niespełnione w 5–8 lub 8–30 GeV → metoda nie ma wymaganej precyzji")
        head = "ODRZUCONA" if (D == NO or K2 == NO) else "NIE PRZYJĘTA jako nominał"
        return f"**{head}** — niespełnione: {', '.join(hard_no)}" + (f"; nierozstrzygnięte: {', '.join(hard_un)}" if hard_un else ""), why
    if hard_un:
        if c4_req is None and C4 != OK:
            why.append("nie wiadomo, czy C4 obowiązuje (brak P1 z dnet4_pis.json)")
        return f"**NIEROZSTRZYGNIĘTE** — brak danych / fit nieokreślony: {', '.join(hard_un)}", why
    q, un = [], []
    if B == "PT5":
        q.append("tylko dla pT ≥ 5 GeV (B zawodzi w 3–5 GeV)")
    if C1 == NO:
        q.append("przyjęty konkretny model (ziarno 1) z własnym SF, nie architektura (C1 niespełnione)")
    elif C1 == UN:
        un.append("C1 (ziarna) → nie wiadomo: architektura czy konkretny model")
    if run == "428770":
        if C2 == NO:
            q.append("SF per klasa triggera (C2 niespełnione; dopuszczalne tylko tam, gdzie analiza już rozróżnia klasy)")
        elif C2 == UN:
            un.append("C2 (klasy triggera) → nie wiadomo, czy SF per klasa")
    if C3 == NO:
        q.append("SF per bin pT (C3: p < 0,01; każdy użyty bin spełnia B)")
    elif C3 == UN:
        un.append("C3 (biny pT) → nie wiadomo, czy jeden SF na run")
    if sp == "Dplus":
        why.append("D⁺: SF bezwzględny pozostaje nieznany (jak dla drzewa); R mówi „sieć skalibrowana nie gorzej niż drzewo”")
    if c4_req is False:
        why.append("C4 nie jest wymaganiem (D⁺ albo P1 spełnione, §9)")
    if un:
        return "**PRZYJĘTA co najmniej WARUNKOWO — forma NIEROZSTRZYGNIĘTA**: " + "; ".join(q + un), why
    if q:
        return "**PRZYJĘTA WARUNKOWO**: " + "; ".join(q), why
    return "**SIEĆ PRZYJĘTA JAKO NOMINAŁ** (jeden SF na run)", why


# ------------------------------------------------------------------------------------------------ główny
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=os.path.join(HERE, "runs", "lxplus4"))
    ap.add_argument("--out", default="")
    a = ap.parse_args(argv)
    C, D = Ctx(a.dir), Doc()
    D.L.append("# WERDYKT DNET etap 4 — odczyt wg KRYTERIA_etap4.md")
    D.p(f"Katalog: `{a.dir}`. Skrypt `DNET/verdict4.py` — sam odczyt i porównanie z progami zapisanymi przed wynikami; żadnych fitów. "
        "Statusy: SPEŁNIONE / NIESPEŁNIONE / NIEROZSTRZYGNIĘTE (brak danych albo fit nieokreślony).")

    # ---- kompletność
    D.h("0. Kompletność wejścia")
    for f in FILES:
        j = C.J[f]
        D.p(f"- `{f}`: " + ("**BRAK**" if j is None else f"jest (wygenerowano {j.get('generated', '?')}" + (f", nboot {j['nboot']}" if "nboot" in j else "") + ")"))
    for f in ("DONE", "DONE23"):
        D.p(f"- `{f}`: {'jest' if os.path.isfile(os.path.join(a.dir, f)) else '**BRAK**'}")
    fails = sorted(os.path.basename(x) for x in glob.glob(os.path.join(a.dir, "fail", "*")))
    D.p(f"- węzły w `fail/`: {('**' + ', '.join(fails) + '**') if fails else 'brak'}")
    metas = sorted(os.path.basename(x)[:-5] for x in glob.glob(os.path.join(a.dir, "meta", "*.json")))
    exp = [REF["D0"], REF["Dplus"]] + [f"{sp}.tr.H+e04.mc.d5.s{s}" for sp in ("D0", "Dplus") for s in (1, 2, 3)] + ["D0.tr.H+e04.sb2mc.d5.s1"]
    D.p(f"- `meta/`: {len(metas)} z 9 oczekiwanych modeli" + (f"; **brak: {', '.join(t for t in exp if t not in metas)}**" if any(t not in metas for t in exp) else ""))
    warn = []
    if C.R:
        if C.R.get("force_data"):
            warn.append("dnet4_R.json policzony z --force-data (bramka domknięcia WYŁĄCZONA) — niezgodne z §1.A")
        if (C.R.get("nboot") or 0) < NBOOT_REQ:
            warn.append(f"bootstrap R: {C.R.get('nboot')} replik < {NBOOT_REQ} (§1)")
        bl = C.R.get("blocked") or []
        b70 = [x for x in bl if x.endswith("|0.7")]
        D.p(f"- zablokowane przez bramkę domknięcia (`blocked`): {len(bl)}, w tym WP70: {len(b70)}" + (": " + "; ".join(b70) if b70 else ""))
        for sp in ("D0", "Dplus"):
            if sp not in (C.R.get("tags") or {}):
                warn.append(f"dnet4_R.json: brak gatunku {sp} (brak plików app albo modeli)")
            else:
                D.p(f"- R {sp}: sieci {', '.join(C.R['tags'][sp].get('nets', []))}; pliki app {C.R['tags'][sp].get('files')}")
                for smp in ("data_427929", "data_428770"):
                    if smp not in (C.R["tags"][sp].get("files") or {}):
                        warn.append(f"dnet4_R.json: {sp} bez próbki {smp}")
    for nm, j in (("dnet4_results.json", C.RES), ("dnet4_ref2_results.json", C.REF2)):
        if j and j.get("fast"):
            warn.append(f"{nm}: tryb --fast (błędy SF z HESSE, nie z przedziałów)")
    for x in C.miss + warn:
        D.p(f"- **UWAGA: {x}**")

    # ---- P1 (§9) — decyduje o C4
    seeds0 = C.tags("D0")
    D.h("1. Diagnostyka wolnego pionu D* (§9, P1) — rozstrzyga, czy C4 obowiązuje")
    p1, rows = crit_P1(C, seeds0[1])
    D.table(rows)
    for s in (2, 3):
        _p, _ = crit_P1(C, seeds0[s])
        D.p(f"- opisowo ziarno {s} (`{seeds0[s]}`): P1 = " + {True: "spełnione", False: "niespełnione", None: "nierozstrzygnięte"}[_p])
    sb = C.tags("D0", "sb2mc")
    _p, _ = crit_P1(C, sb)
    D.p(f"- opisowo sb2mc (`{sb}`): P1 = " + {True: "spełnione", False: "niespełnione", None: "nierozstrzygnięte"}[_p])
    if C.PIS:
        D.p(f"- wolny pion znaleziony w tokenach: {pm(C.PIS.get('pis_found'))}, widziany przez e04: {pm(C.PIS.get('pis_seen_e04'))}; próbka tagu {C.PIS.get('n_p_tag')}, sygnał a_ {C.PIS.get('n_a')}")

    Ag, rowsAg = crit_Aglob(C)
    verdicts, nsig_bkg = [], [0, 0]
    for isp, sp in enumerate(("D0", "Dplus")):
        seeds = C.tags(sp)
        tag = seeds[1]
        D.h(f"{2 + isp}. {SPN[sp]} — sieć `{tag}` (ziarno 1) wobec drzewa `{REF[sp]}`")
        D.h("A. Domknięcie na MC (bramka)", 3)
        A, rows, flA = crit_A(C, sp, tag)
        D.table(rows + rowsAg)
        for s in (2, 3):
            _s, _, _f = crit_A(C, sp, seeds[s])
            D.p(f"- ziarno {s} (`{seeds[s]}`): A {_s} ({sum(x is True for x in _f)}/3 binów)")
        c4_req = False if sp == "Dplus" else (None if p1 is None else (not p1))
        C1r, C1d, rowsC1 = crit_C1(C, sp, p1)
        C1 = C1r
        if sp == "D0" and p1 is not True:           # §9(iv): przy P1 spełnionym C1 opiera się wyłącznie na R
            C1 = NO if NO in (C1r, C1d) else (UN if UN in (C1r, C1d) or (p1 is None and C1d != OK) else OK)
        C2, rowsC2 = crit_C2(C, sp, tag)
        K2, rowsK2, slope_ok = crit_K2(C, sp, tag)
        per_run = {}
        for run in MAIN:
            B, rowsB = crit_B(C, sp, tag, run)
            C3, rowsC3 = crit_C3(C, sp, tag, run)
            Dd, Dc, rowsD = crit_D(C, sp, tag, run)
            C4, rowsC4 = crit_C4(C, tag, run) if sp == "D0" else (None, [])
            per_run[run] = (B, C3, Dd, Dc, C4)
            D.h(f"Run {run}, klasa {MAIN[run]}: B precyzja, C3 biny pT, D zysk" + (", C4 zgodność z tagiem D*" if sp == "D0" else ""), 3)
            if sp == "D0":
                rowsC4 = [(r[0] + ("" if c4_req else (" (opisowo, §9)" if c4_req is False else " (obowiązywanie nieznane)")),) + tuple(r[1:]) for r in rowsC4]
            D.table(rowsB + rowsC3 + rowsD + rowsC4)
            for b in BINS:
                x = C.rd(sp, tag, run, MAIN[run], b)
                if x and fin(x["bkg"]):
                    nsig_bkg[1] += 1
                    nsig_bkg[0] += x["bkg"] > x["stat"]
        D.h("C1 ziarna, C2 klasy triggera (428770), K2 rzeźbienie", 3)
        D.table(rowsC1 + rowsC2 + rowsK2)
        if sp == "D0":
            D.h("SF_sieć = R × SF_H (SF_H z tagu D*, błędy kwadratowo; §1) — WP70, ziarno 1", 3)
            rows = []
            for run, cl in (("427929", "MB427929"), ("428770", "all"), ("428770", "J12only")):
                for b in BINS:
                    x, sh, sn = C.rd(sp, tag, run, cl, b), C.sf(REF[sp], run, cl, b), C.sf(tag, run, cl, b)
                    v = (x["R"] * sh[0], x["R"] * sh[0] * math.hypot(x["tot"] / x["R"], sh[1] / sh[0])) if x and x["tot"] and sh else None
                    rows.append((f"{run}:{cl}", b, pm(x["R"], x["tot"]) if x else "–", pm(*sh) if sh else "– (nieokreślony)", pm(*v) if v else "–", pm(*sn) if sn else "– (nieokreślony)"))
            D.table(rows, head=("run:klasa", "bin", "R ± σ_tot", "SF_H (D*)", "SF_sieć = R × SF_H", "SF_sieć bezpośrednio z tagu D*"))
            D.p("(klasa `notJ12only`: analyze.py nie liczy dla niej SF z tagu D* → SF_sieć = R × SF_H niedostępny.)")
        D.h("WP50 / WP80 (opisowo, HESSE)", 3)
        rows = []
        for run, cl in MAIN.items():
            for b in BINS:
                v = [dig(C.R, "data", sp, tag, f"{run}|{cl}", b, wp) or {} for wp in ("0.5", "0.7", "0.8")]
                rows.append((f"{run}:{cl}", b) + tuple(pm(o.get("R"), o.get("R_err")) for o in v))
        D.table(rows, head=("run:klasa", "bin", "R WP50", "R WP70", "R WP80"))
        D.h(f"Werdykt E — {SPN[sp]}", 3)
        for run in MAIN:
            B, C3, Dd, Dc, C4 = per_run[run]
            head, why = verdict_E(sp, run, A, sum(f is False for f in flA), Ag, B, C1, C2, C3, C4, Dd, K2, c4_req)
            comp = f"A {A}, A-glob {Ag}, B {'SPEŁNIONE tylko pT ≥ 5' if B == 'PT5' else B}, C1 {C1}, " + (f"C2 {C2}, " if run == "428770" else "") + f"C3 {C3}, " \
                   + (f"C4 {C4} ({'wymagane' if c4_req else ('opisowo' if c4_req is False else 'obowiązywanie nieznane')}), " if sp == "D0" else "") + f"D {Dd}, D-kontr. {Dc}, K2 {K2}"
            if Dc == NO and Dd == OK:
                why.append("UWAGA: kontrola G na pasmach (≥ 1,20) NIESPEŁNIONA mimo G_fit ≥ 1,15 — §1.D nie mówi, czy to zmienia D; do decyzji")
            if slope_ok is False:
                why.append("UWAGA: kontrola K2 — nachylenie ε_B(m) różni się od H o > 0,10/100 MeV (RAPORT §4 K2 „kontrola”; analyze.py nie wlicza jej do K2) — do decyzji")
            line = f"- **{SPN[sp]}, run {run}**: {head}  \n  składowe: {comp}" + ("".join(f"  \n  – {w}" for w in why))
            D.p(line)
            verdicts.append(line)

    # ---- sb2mc (§4)
    D.h("4. Wariant z tłem treningowym z danych `sb2mc` (§4, tylko D⁰, 1 ziarno)")
    if sb:
        A, rows, flA = crit_A(C, "D0", sb)
        C2, rowsC2 = crit_C2(C, "D0", sb)
        K2, rowsK2, _ = crit_K2(C, "D0", sb)
        better, parts = [], []
        for run in MAIN:
            B, rowsB = crit_B(C, "D0", sb, run)
            C3, rowsC3 = crit_C3(C, "D0", sb, run)
            Dd, Dc, rowsD = crit_D(C, "D0", sb, run)
            rows += rowsB + rowsC3 + rowsD
            fl = []
            for b in BINS:
                x, y = C.rd("D0", sb, run, MAIN[run], b), C.rd("D0", seeds0[1], run, MAIN[run], b)
                d = x["G"] - y["G"] if x and y and fin(x["G"]) and fin(y["G"]) else None
                fl.append(None if d is None else d >= SB_DG)
                rows.append(("sb2mc>mc", f"{run}:{MAIN[run]} {b}: G_fit(sb2mc) − G_fit(mc)", pm(d), f"≥ {SB_DG}", {True: OK, False: NO, None: UN}[fl[-1]], ""))
            sg = nof(fl, D_NBIN)
            better.append(NO if (sg == NO or B in (NO, "PT5")) else (UN if UN in (sg, B) else OK))
            parts.append(f"run {run}: A–D: B {'tylko pT ≥ 5' if B == 'PT5' else B}, C3 {C3}, D {Dd}; zysk o ≥ {SB_DG} w ≥ 2 z 3 binów: {sg}")
        D.table(rows + rowsC2 + rowsK2)
        bt = NO if NO in better else (UN if UN in better else OK)
        line = f"- **sb2mc**: A {A}, C2 {C2}, K2 {K2}; " + "; ".join(parts) + f" → „lepszy niż wariant mc”: **{bt}**" + (" → najwyżej **KANDYDAT** (1 ziarno, wymaga ziaren)" if bt == OK else "") \
               + ". Oczekiwanie zapisane z góry: R wyraźnie < 1 — samo nie dyskwalifikuje, o ile B i C2–C3 spełnione."
        D.p(line)
        verdicts.append(line)
    # ---- bin 2–3 (§5)
    D.h("5. Bin 2–3 GeV (§5)")
    rows, vis, lines = [], {}, []
    for sp in ("D0", "Dplus"):
        for run in MAIN:
            s, zn, zt, r = crit_d23(C, sp, run)
            rows += r
            if sp == "D0":
                vis[run] = (s, zn)
            else:
                g = zn / zt if fin(zn) and fin(zt) and zt > 0 else None
                sg = UN if g is None else (OK if g >= D23_GAIN else NO)
                rows.append(("2–3", f"D⁺ {run}: Z_sieć/Z_drzewo (każdy przy własnym WP)", pm(g, None, 2), f"≥ {D23_GAIN}", sg, f"Z sieć {pm(zn, None, 1)}, drzewo {pm(zt, None, 1)}"))
                lines.append(f"- **2–3 GeV D⁺, run {run}**: zysk sieci Z_sieć/Z_drzewo ≥ {D23_GAIN}: **{sg}**")
    D.table(rows)
    for i, (run, (s, z)) in enumerate(vis.items()):
        lines.insert(i, f"- **2–3 GeV D⁰, run {run}**: „D⁰ widoczne” (sieć; Z ≥ 5 i Z ≥ 4 dla obu sąsiednich rzędów tła): **{s}** (Z = {pm(z, None, 1)})")
    if vis and all(fin(z) and z < D23_Z for _, z in vis.values()):
        lines.append("- **§6: D⁰ Z < 5 w obu runach → „sieć nie otwiera binu 2–3 dla D⁰”** (wynik negatywny → notatka w charm-bkg-notes)")
    lines.append("- zastrzeżenia: Z wariantów rzędu tła = N/σ(HESSE) bez bootstrapu (nominalny: max(HESSE, bootstrap)); ε_S(dane) w 2–3 GeV niekalibrowane — Z dotyczy widoczności piku, nie przekroju (§5)")
    for x in lines:
        D.p(x)
    verdicts += lines[:-1]

    # ---- podsumowanie
    D.h("6. WERDYKT KOŃCOWY")
    D.p("- **P1 (wolny pion D*)**: " + {True: "SPEŁNIONE — efekt istotny: SF sieci z tagu D* nie kalibruje próbki inkluzywnej; kalibracja podstawowa = R × SF_H; C4 opisowo",
                                       False: "NIESPEŁNIONE — obie kalibracje równoprawne; C4 obowiązuje", None: "NIEROZSTRZYGNIĘTE — nie wiadomo, czy C4 obowiązuje"}[p1])
    for x in verdicts:
        D.p(x)
    if nsig_bkg[1]:
        D.p(f"- §6: σ_tło > σ_stat w {nsig_bkg[0]} z {nsig_bkg[1]} punktów (ziarno 1, klasy główne)" + (" → **kalibrację ogranicza rzeźbienie tła, nie statystyka**" if nsig_bkg[0] > nsig_bkg[1] / 2 else ""))
    if C.miss or fails:
        D.p("- **braki wejścia**: " + "; ".join(C.miss + ([f"fail/: {', '.join(fails)}"] if fails else [])))
    D.p("- każdy wynik negatywny → notatka w `charm-bkg-notes` (§6); wszystko po obejrzeniu = „post hoc” (§8)")
    D.h("Luki odczytu (wielkości, których kod etapu 4 nie zapisuje — patrz też nagłówek)", 3)
    for x in ("§1.A: odchylenie przenoszenia ε (|R_prawda − 1| > 0,03) NIE jest w `R_err_tot` z dn4_ana.py — skrypt dodaje je kwadratowo sam (kolumna σ_tot powyżej już je zawiera).",
              "§1.C4: błąd różnicy przyjęty jako σ_tot(R) ⊕ σ(ilorazu SF); w uwagach podano też liczbę σ z samym σ(ilorazu).",
              "§1.D: G_fit i σ_N/N w JSON pochodzą z kowariancji HESSE (bez bootstrapu po zdarzeniach).",
              "§5: Z dla sąsiednich rzędów tła — HESSE bez bootstrapu; „sąsiedzi WP” mają tylko fit nominalny HESSE.",
              "§9 P1: ε_c/ε_a bez błędu w JSON (próg bezwzględny — błąd niewymagany); błąd podwójnego stosunku = suma kwadratowa (próbki sieci i drzewa te same → zawyżony).",
              "SF_sieć = R × SF_H dla klasy notJ12only niedostępny (analyze.py: klasy SF = all, J12only, MB427929)."):
        D.p(f"- {x}")
    txt = D.text()
    out = a.out or os.path.join(a.dir, "WERDYKT_etap4.md")
    try:
        os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
        open(out, "w").write(txt)
        txt += f"\n[zapisano {out}]\n"
    except OSError as e:
        txt += f"\n[NIE zapisano {out}: {e}]\n"
    sys.stdout.write(txt)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception:  # noqa: BLE001 — nieoczekiwana struktura wejścia: nie ma werdyktu, a nie „spełnione”
        import traceback
        traceback.print_exc()
        print("\nBŁĄD ODCZYTU (nieoczekiwana struktura pliku wejściowego) → całość NIEROZSTRZYGNIĘTA; werdykt NIE został zapisany.")
        sys.exit(2)
