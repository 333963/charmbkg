#!/usr/bin/env python3
"""DNET etap 5a — JEDNORAZOWY ODCZYT czterech werdyktów (+ reguła GPU) wg DNET/KRYTERIA_etap5.md (progi wiążące, zapisane przed wynikami).

  ~/Documents/charm/.venv_h5/bin/python DNET/verdict5.py [--dir DNET/runs/lxplus5] [--out PLIK.md]

Czyta: dnet5_cmp.json (dn5_ana.py: pary sparowane + kaskada H × E), dnet5_results.json (analyze.py compute), gpu/probe.json, test/{bench,equiv}.json,
folds/<tag>/fold<k>.json, meta/, fail/, DONE. Wypisuje i zapisuje <dir>/WERDYKT_etap5.md.
ŻADNYCH fitów: tylko odczyt i porównanie z progami. Brak pliku / pola → NIEROZSTRZYGNIĘTE (nigdy „spełnione”)."""
import argparse
import glob
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

# ------------------------------------------------------------------------------------------------ progi (KRYTERIA_etap5.md)
LC_CLOSED_100_50, LC_CLOSED_50_25 = 0.02, 0.03     # §2 „zamknięte liczbą”: r100/50 − 1 < 0,02 w ≥ 2 z 3 binów i r50/25 − 1 < 0,03 w ≥ 2 z 3
LC_OPEN_100_50 = 0.03                              # §2 „otwiera #8”: r100/50 − 1 ≥ 0,03 w ≥ 2 z 3 binów i r50/25 ≥ r100/50 w ≥ 2 binach
PAIR_MC, PAIR_NSIG, PAIR_SB = 1.05, 3.0, 1.03      # §3: Z(trp)/Z(tr) ≥ 1,05 na MC w ≥ 2 z 3 binów, każdy ≥ 3σ; pasma 427929 ≥ 1,03 w ≥ 2 z 3
PAIR_DCORR = 0.02                                  # §3 kontrola: |corr(score, m)| na tle MC nie większe niż dla tr o > 0,02
K2_NSIG, K2_SLOPE = 2.0, 0.10                      # RAPORT.md §4 K2 (kontrola §3)
E_G, E_DG = 1.15, 0.10                             # §4: G_E(427929) ≥ 1,15 w ≥ 2 z 3 binów i G_E − G_Htight ≥ 0,10 w tych samych binach
E_CORR, E_RATIO = 0.02, (0.85, 1.0)                # §4 opisowo: |corr(score_E, m)| ≤ 0,02; (G_E − 1)/(G_tr.H+e04 − 1) oczekiwane 0,85–1,0
NBIN = 2
GPU_WAIT_S, GPU_SPEEDUP = 7200, 10.0               # §1: czas w kolejce ≤ 2 h i s/epokę(GPU) ≤ 1/10 s/epokę(CPU, 4 wątki)
CPU_SEC_EPOCH_FALLBACK = {"tr": 82.0, "trp": 107.0}  # §1 (pomiar zgrubny lxplus932, 4 wątki) — używane TYLKO gdy brak test/bench.json
FOLD_MAX_S = 7200                                  # §1: najdłuższy job foldu ≤ 2 h
NBOOT_REQ = 200                                    # §2: bootstrap sparowany, 200 replik

BINS = ["3-5", "5-8", "8-30"]
BASE = [f"D0.tr.H+e04.mc.d5.s{s}" for s in (1, 2, 3)]
LC50, LC25, PAIR = "D0.tr.H+e04.mcS50.d5.s1", "D0.tr.H+e04.mcS25.d5.s1", "D0.trp.H+e04.mc.d5.s1"
ETAGS = {"D0": ["D0.trE.E04.mc.n.s1", "D0.trE.E04h70.mc.n.s1"], "Dplus": ["Dplus.trE.E04.mc.n.s1", "Dplus.trE.E04h70.mc.n.s1"]}
REF = {"D0": "D0.xgb.H.sb.n.s1", "Dplus": "Dplus.xgb.H.sb.n.s1"}
SPN = {"D0": "D⁰", "Dplus": "D⁺"}
OK, NO, UN, DESC = "SPEŁNIONE", "NIESPEŁNIONE", "NIEROZSTRZYGNIĘTE", "opisowo"
S3 = {True: OK, False: NO, None: UN}


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


def nof(flags, need=NBIN):
    n_ok, n_un = sum(f is True for f in flags), sum(f is None for f in flags)
    if n_ok >= need:
        return True
    return None if n_ok + n_un >= need else False


def both(*s):
    """Koniunkcja True/False/None."""
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
            raise ValueError("nie słownik")
        return j
    except Exception as e:  # noqa: BLE001
        miss.append(f"{rel}: {'brak pliku' if not os.path.isfile(p) else 'nieczytelny (' + type(e).__name__ + ')'}")
        return None


def pair(CMP, name, b, bkg="mc", wp="0.7"):
    x = dig(CMP, "pairs", name, f"{b}|{wp}|{bkg}")
    return (x["r"], x.get("err")) if x and fin(x.get("r")) else None


def why_pair(CMP, name):
    if CMP is None:
        return "brak dnet5_cmp.json"
    for m in CMP.get("missing") or []:
        if m.get("pair") == name:
            return "brak score: " + ", ".join(m.get("tags", []))
    return "brak wpisu (za mało wierszy w binie?)" if name in (CMP.get("pairs") or {}) else "pary nie ma w pliku"


# ------------------------------------------------------------------------------------------------ §2 krzywa uczenia
def v_lc(CMP):
    rows, c1, c2, o1, o2 = [], [], [], [], []
    for b in BINS:
        a, c = pair(CMP, "lc_100over50", b), pair(CMP, "lc_50over25", b)
        c1.append(None if not a else a[0] - 1 < LC_CLOSED_100_50)
        o1.append(None if not a else a[0] - 1 >= LC_OPEN_100_50)
        c2.append(None if not c else c[0] - 1 < LC_CLOSED_50_25)
        o2.append(None if not (a and c) else c[0] >= a[0])
        rows.append(("§2", f"{b}: r100/50 = Z(100 %; śr. ziaren 1–3)/Z(50 %), MC WP70", pm(*a) if a else "–", f"zamknięte: − 1 < {LC_CLOSED_100_50}; otwiera: − 1 ≥ {LC_OPEN_100_50}",
                     UN if not a else ("zamknięte" if c1[-1] else ("otwiera" if o1[-1] else "pomiędzy")), "" if a else why_pair(CMP, "lc_100over50")))
        rows.append(("§2", f"{b}: r50/25 = Z(50 %)/Z(25 %), MC WP70", pm(*c) if c else "–", f"zamknięte: − 1 < {LC_CLOSED_50_25}; otwiera: r50/25 ≥ r100/50",
                     UN if not c else ("zamknięte" if c2[-1] else "–") + ("; ≥ r100/50" if o2[-1] else ("; < r100/50" if o2[-1] is False else "")), "" if c else why_pair(CMP, "lc_50over25")))
    closed, opens = both(nof(c1), nof(c2)), both(nof(o1), nof(o2))
    rows.append(("§2", "„zamknięte liczbą”", f"r100/50: {sum(f is True for f in c1)}/3; r50/25: {sum(f is True for f in c2)}/3", "≥ 2 z 3 i ≥ 2 z 3", S3[closed], ""))
    rows.append(("§2", "„otwiera #8”", f"r100/50 − 1 ≥ {LC_OPEN_100_50}: {sum(f is True for f in o1)}/3; r50/25 ≥ r100/50: {sum(f is True for f in o2)}/3", "≥ 2 z 3 i ≥ 2 binach", S3[opens], ""))
    for b in BINS:                                  # opisowo: pasma danych, WP80/50, rozrzut ziaren
        d = []
        for nm in ("lc_100over50", "lc_50over25"):
            for bk in ("sb427929", "sb428770"):
                x = pair(CMP, nm, b, bk)
                d.append(f"{nm[3:]} {bk}: {pm(*x) if x else '–'}")
        for wp in ("0.8", "0.5"):
            x = pair(CMP, "lc_100over50", b, "mc", wp)
            d.append(f"100over50 MC WP{wp[2:]}0: {pm(*x) if x else '–'}")
        for nm in ("seed2_over_seed1", "seed3_over_seed1"):
            x = pair(CMP, nm, b)
            d.append(f"{nm}: {pm(*x) if x else '–'}")
        rows.append(("§2", f"{b}: pasma danych, WP80/50, rozrzut ziaren", "; ".join(d), "–", DESC, ""))
    if closed is True:
        v = "**ZAMKNIĘTE LICZBĄ** — pretrening, wspólny trening D⁰/D⁺ (#8), sieci per bin: nie warto (więcej MC sygnału nie pomaga)"
    elif opens is True:
        v = "**OTWIERA #8** — krzywa uczenia nienasycona (wspólny trening D⁰/D⁺ / więcej sygnału ma sens)"
    elif closed is False and opens is False:
        v = "**NIEROZSTRZYGNIĘTE wg kryterium §2** (ani „zamknięte”, ani „otwiera”) → 5b: ziarna 2–3 dla wariantu 50 % (10 jobów foldów)"
    else:
        v = "**NIEROZSTRZYGNIĘTE — brak danych** (" + "; ".join(sorted({why_pair(CMP, "lc_100over50"), why_pair(CMP, "lc_50over25")})) + ")"
    return v, rows


# ------------------------------------------------------------------------------------------------ §3 cechy par
def v_pair(CMP, RES):
    rows, f_mc, f_sb = [], [], []
    for b in BINS:
        a, s = pair(CMP, "pair_over_base", b), pair(CMP, "pair_over_base", b, "sb427929")
        ns = (a[0] - 1) / a[1] if a and fin(a[1]) and a[1] > 0 else None
        f_mc.append(None if not a or ns is None else (a[0] >= PAIR_MC and ns >= PAIR_NSIG))
        f_sb.append(None if not s else s[0] >= PAIR_SB)
        rows.append(("§3", f"{b}: Z(trp)/Z(tr; śr. ziaren 1–3), MC WP70", pm(*a) if a else "–", f"≥ {PAIR_MC} i ≥ {PAIR_NSIG:.0f}σ", S3[f_mc[-1]],
                     f"istotność (r − 1)/σ = {ns:+.1f}σ" if ns is not None else why_pair(CMP, "pair_over_base")))
        x8 = pair(CMP, "pair_over_base", b, "sb428770")
        rows.append(("§3", f"{b}: ten sam stosunek, dalekie pasma 427929", pm(*s) if s else "–", f"≥ {PAIR_SB}", S3[f_sb[-1]],
                     f"pasma 428770 (opisowo): {pm(*x8) if x8 else '–'}" if s else why_pair(CMP, "pair_over_base")))
    acc = both(nof(f_mc), nof(f_sb))
    rows.append(("§3", "łącznie", f"MC: {sum(f is True for f in f_mc)}/3; pasma 427929: {sum(f is True for f in f_sb)}/3", "≥ 2 z 3 i ≥ 2 z 3", S3[acc], ""))
    # kontrole (mogą unieważnić „przyjęte”)
    ctrl = []
    cf = dig(RES, "species", "D0", "cfg", PAIR)
    k = "sculpt|mc_SR|0.7"
    s, sr = dig(cf, k), dig(cf, "ref_same_rows", k) or dig(RES, "species", "D0", "cfg", dig(RES, "species", "D0", "ref") or REF["D0"], k)
    if s and sr:
        ok = abs(s["R"] - 1) <= abs(sr["R"] - 1) + K2_NSIG * math.hypot(s["R_err"], sr["R_err"])
        ctrl.append(("K2", ok))
        rows.append(("§3 kontr.", "K2: |R_sculpt − 1| (tło MC, WP70) vs drzewo H", f"trp {pm(s['R'], s['R_err'])}; H {pm(sr['R'], sr['R_err'])}", f"nie większe niż H w {K2_NSIG:.0f}σ", S3[ok], ""))
        if fin(s.get("slope")) and fin(sr.get("slope")):
            d = s["slope"] - sr["slope"]
            ctrl.append(("nachylenie", abs(d) <= K2_SLOPE))
            rows.append(("§3 kontr.", "K2: nachylenie ε_B(m)/100 MeV, trp − H", f"{d:+.3f}", f"|Δ| ≤ {K2_SLOPE}", S3[abs(d) <= K2_SLOPE], ""))
    else:
        ctrl.append(("K2", None))
        rows.append(("§3 kontr.", "K2: |R_sculpt − 1| vs drzewo H", "–", f"nie większe niż H w {K2_NSIG:.0f}σ", UN, "brak dnet5_results.json albo wpisu trp"))
    cp = dig(cf, "corr_score_m")
    cb = [dig(RES, "species", "D0", "cfg", t, "corr_score_m") for t in BASE]
    cbv = [abs(x) for x in cb if fin(x)]
    if fin(cp) and cbv:
        d = abs(cp) - sum(cbv) / len(cbv)
        ctrl.append(("corr", d <= PAIR_DCORR))
        rows.append(("§3 kontr.", "|corr(score, m)| tło MC: trp − tr (śr. dostępnych ziaren)", f"{d:+.4f}", f"≤ {PAIR_DCORR}", S3[d <= PAIR_DCORR],
                     f"trp {cp:+.4f}; tr ziarna: " + ", ".join(pm(x, None, 4) for x in cb) + f"; względem najmniejszego |corr|: {abs(cp) - min(cbv):+.4f}"))
    else:
        ctrl.append(("corr", None))
        rows.append(("§3 kontr.", "|corr(score, m)| tło MC: trp − tr", "–", f"≤ {PAIR_DCORR}", UN, "brak corr_score_m w dnet5_results.json"))
    sf = dig(RES, "species", "D0", "sf") or {}
    d = [f"{kk.split('|', 1)[1].rsplit('|', 1)[0]}: {pm(v.get('SF'), v.get('SF_err'), 2)}{'' if v.get('determined') else ' (n)'}" for kk, v in sorted(sf.items()) if kk.startswith(PAIR + "|") and kk.endswith("|0.7")]
    rows.append(("§3", "SF z tagu D* trp WP70", "; ".join(d) if d else "–", "–", DESC, "może nie dotyczyć próbki inkluzywnej (KRYTERIA_etap4 §9)"))
    if acc is True:
        bad, unk = [n for n, f in ctrl if f is False], [n for n, f in ctrl if f is None]
        v = "**PRZYJĘTE DO DALSZYCH ZIAREN**" + (f" — ale kontrola NIESPEŁNIONA: {', '.join(bad)} (§3: może unieważnić „przyjęte” — do decyzji)" if bad else "") \
            + (f"; kontrole bez danych: {', '.join(unk)}" if unk else "")
    elif acc is False:
        v = "**PORZUCIĆ PO 1 ZIARNIE** (kontrole nie zmieniają tego werdyktu) → notatka w charm-bkg-notes"
    else:
        v = "**NIEROZSTRZYGNIĘTE — brak danych** (" + why_pair(CMP, "pair_over_base") + ")"
    return v, rows


# ------------------------------------------------------------------------------------------------ §4 tagger E
def e_variant(CMP, et):
    """→ (status True/False/None, średnie G_mc albo None, wiersze, uwagi)."""
    C = dig(CMP, "cascade", et)
    rows, fl, gmc, notes = [], [], [], []
    if not C:
        why = "brak dnet5_cmp.json" if CMP is None else next(("brak score: " + ", ".join(m.get("tags", [])) for m in CMP.get("missing") or [] if m.get("cascade") == et), "brak wpisu kaskady")
        return None, None, [("§4", f"`{et}`", "–", f"G_E(427929) ≥ {E_G} i G_E − G_Htight ≥ {E_DG}", UN, why)], [why]
    for b in BINS:
        o = C.get(b)
        if not o:
            rows.append(("§4", f"`{et}` {b}", "–", f"G_E(427929) ≥ {E_G} i G_E − G_Htight ≥ {E_DG}", UN, "bin pominięty (< 50 wierszy sygnału po H)")); fl.append(None); gmc.append(None); continue
        wp = o.get("wp_nom")
        w = o.get(f"wp{wp}") or {}
        g, ge, ht = w.get("G_sb427929"), w.get("G_sb427929_err"), w.get("G_Htight_sb427929")
        f = (g >= E_G and g - ht >= E_DG) if fin(g) and fin(ht) else (False if fin(g) and g < E_G else None)
        fl.append(f)
        gmc.append(w.get("G_mc") if fin(w.get("G_mc")) else None)
        nmc = dig(o, "n_b_H", "mc")
        note = f"WP_E {wp}; G_E(MC) {pm(w.get('G_mc'), w.get('G_mc_err'))}; G_E(428770, opisowo) {pm(w.get('G_sb428770'), w.get('G_sb428770_err'))}; tło po H: MC {nmc}, 427929 {dig(o, 'n_b_H', 'sb427929')}"
        if fin(nmc) and nmc < 50:
            note += "; UWAGA: < 50 wierszy tła MC po H → WP_E = 0,80 domyślnie (NIE wybrany na MC)"
            notes.append(f"{et} {b}: WP_E domyślny 0,80 (tło MC po cięciu H: {nmc} wierszy)")
        fn = dig(o, "G_fullnet", "sb427929")
        if fin(fn) and fin(g) and fn != 1:
            q = (g - 1) / (fn - 1)
            note += f"; (G_E − 1)/(G_tr.H+e04 − 1) = {q:.2f} (oczekiwane {E_RATIO[0]}–{E_RATIO[1]}; G pełnej sieci {fn:.3f})"
        c = o.get("corr_E_m_afterH")
        sc = o.get("sculpt_afterH")
        if fin(c):
            note += f"; corr(E, m) po H {c:+.4f} ({'≤' if abs(c) <= E_CORR else '> — 5b: E z DisCo'} {E_CORR})"
        if sc:
            note += f"; R_sculpt po H {pm(sc.get('R'), sc.get('R_err'))}" + ("" if abs(sc["R"] - 1) <= 2 * sc["R_err"] else " (≠ 1 o > 2σ)")
        rows.append(("§4", f"`{et}` {b}: G_E(pasma 427929); G_E − G_Htight", f"{pm(g, ge)}; {pm(g - ht) if fin(g) and fin(ht) else '–'} (G_Htight {pm(ht)})",
                     f"≥ {E_G} i ≥ {E_DG}", S3[f], note))
    s = nof(fl)
    mean = sum(gmc) / 3 if all(x is not None for x in gmc) else None
    rows.append(("§4", f"`{et}`: łącznie", f"{sum(f is True for f in fl)}/3 binów; średnie G_E(MC) = {pm(mean)}", f"≥ {NBIN} z 3", S3[s], ""))
    return s, mean, rows, notes


def v_E(CMP, sp):
    a, b = ETAGS[sp]
    sa, ma, ra, na = e_variant(CMP, a)
    sb, mb, rb, nb = e_variant(CMP, b)
    rows, notes = ra + rb, na + nb
    if ma is not None and mb is not None:
        sel, ssel, soth = (a, sa, sb) if ma >= mb else (b, sb, sa)
        choice = f"wybrany wariant (większe średnie G_E na MC po 3 binach): `{sel}` ({pm(max(ma, mb))} vs {pm(min(ma, mb))})"
    else:
        sel, ssel, soth = None, (sa if sa == sb else None), None
        choice = "wybór wariantu NIEMOŻLIWY (brak G_E(MC) we wszystkich 3 binach obu wariantów)" + ("; oba warianty dają ten sam status" if sa == sb and sa is not None else "")
    rows.append(("§4", f"{SPN[sp]}: wybór wariantu", choice, "max średniego G_E(MC)", DESC, ""))
    if ssel is True:
        v = f"**SPEŁNIONE** — kaskada H × E ma zysk na pasmach 427929 ({'`' + sel + '`' if sel else 'oba warianty'}) → 5b: pomiar SF_E dla {SPN[sp]}"
    elif ssel is False:
        v = f"**NIESPEŁNIONE** dla {'wybranego wariantu `' + sel + '`' if sel else 'obu wariantów'}"
        if sa is False and sb is False:
            v += " — oba warianty: kaskada nie ma z czego dawać zysku; pomiar SF_E w 5b dla tego gatunku NIE jest uruchamiany; notatka w charm-bkg-notes"
        elif soth is True:
            v += " (UWAGA: drugi wariant spełnia kryterium — §4 „wynik negatywny” wymaga niespełnienia dla OBU; do decyzji)"
        elif soth is None:
            v += " (drugi wariant: nierozstrzygnięty)"
    else:
        v = f"**NIEROZSTRZYGNIĘTE — brak danych**: `{a}` {S3[sa]}, `{b}` {S3[sb]}; " + ("wybór wariantu niemożliwy (brak G_E(MC) w 3 binach obu wariantów)" if sel is None else f"wybrany `{sel}`")
    wn = [n for n in notes if "WP_E" in n]
    if wn:
        v += "  \n  – UWAGA (punkt pracy E nie z MC): " + "; ".join(wn)
    return v, rows


# ------------------------------------------------------------------------------------------------ §1 GPU, czasy
def v_gpu(d, miss):
    P, Bn = load(d, "gpu/probe.json", miss), load(d, "test/bench.json", miss)
    rows = []
    if not P:
        return "**NIEROZSTRZYGNIĘTE — brak gpu/probe.json** (sonda wycofana po 6 h albo nie wystartowała → wg §1 zostaje CPU po foldach)", [("§1", "sonda GPU", "–", "", UN, "brak gpu/probe.json")]
    qw = P.get("queue_wait_s")
    f_q = None if not fin(qw) else qw <= GPU_WAIT_S
    rows.append(("§1", "czas sondy w kolejce [s]", pm(qw, None, 0), f"≤ {GPU_WAIT_S}", S3[f_q], f"karta: {P.get('gpu')}, host {P.get('host')}"))
    fl = []
    for t, x in (P.get("folds") or {}).items():
        arch = "trp" if ".trp." in t else "tr"
        cpu, src = dig(Bn, "sec_epoch_full", f"{arch}|4"), "test/bench.json"
        if not fin(cpu):
            cpu, src = CPU_SEC_EPOCH_FALLBACK[arch], "KRYTERIA §1 (pomiar zgrubny; brak test/bench.json)"
        se = x.get("sec_epoch")
        on_gpu = x.get("device") == "cuda"
        f = None if not fin(se) else (on_gpu and se <= cpu / GPU_SPEEDUP)
        fl.append(f)
        rows.append(("§1", f"`{t}` fold 0: s/epokę GPU vs CPU (4 wątki)", f"{pm(se, None, 2)} vs {pm(cpu, None, 1)}" + (f" (×{cpu / se:.0f})" if fin(se) and se > 0 else ""), f"GPU ≤ CPU/{GPU_SPEEDUP:.0f}", S3[f],
                     f"urządzenie {x.get('device')} ({x.get('device_name')}); CPU z: {src}; fold {pm(x.get('sec'), None, 0)} s, {x.get('epochs')} epok" + ("" if on_gpu else "; UWAGA: sonda NIE liczyła na GPU")))
    f_s = None if not fl else both(*fl)
    r = both(f_q, f_s)
    if r is True:
        v = "**SPEŁNIONE** → trening kolejnych etapów na GPU (5 foldów w jednym jobie) z awaryjnym CPU po foldach"
    elif r is False:
        v = "**NIESPEŁNIONE** → zostaje CPU po foldach"
    else:
        v = "**NIEROZSTRZYGNIĘTE** (brak czasu w kolejce albo s/epokę w probe.json)"
    return v, rows


def infra(d, miss, D):
    rows = []
    tags = [LC25, LC50, PAIR] + ETAGS["D0"] + ETAGS["Dplus"]
    for t in tags:
        fs = sorted(glob.glob(os.path.join(d, "folds", t, "fold*.json")))
        js = []
        for f in fs:
            try:
                js.append(json.load(open(f)))
            except Exception:  # noqa: BLE001
                pass
        meta = None
        mp = os.path.join(d, "meta", t + ".json")
        if os.path.isfile(mp):
            try:
                meta = json.load(open(mp))
            except Exception:  # noqa: BLE001
                meta = None
        mx = max((j.get("sec_total") or 0 for j in js), default=None)
        mx = dig(meta, "sec_wall_fold_max") or mx
        dev = sorted({str(j.get("device")) for j in js})
        st = UN if not fin(mx) else (OK if mx <= FOLD_MAX_S else NO)
        rows.append(("§1", f"`{t}`: najdłuższy fold [s]", pm(mx, None, 0), f"≤ {FOLD_MAX_S}", st,
                     f"foldów {len(js)}/5, urządzenie {', '.join(dev) or '–'}, meta {'jest' if meta else 'BRAK (wariant niezłączony)'}" + (f", łączenie {meta.get('sec_merge')} s" if meta else "")))
        if len(js) < 5 or not meta:
            miss.append(f"{t}: foldów {len(js)}/5, meta {'jest' if meta else 'brak'}")
    Eq = load(d, "test/equiv.json", [])
    if Eq:
        rows.append(("§1", "test równoważności stara/nowa ścieżka (test/equiv.json)", f"equivalent = {Eq.get('equivalent')}, mechanics_ok = {Eq.get('mechanics_ok')}", "True", OK if Eq.get("equivalent") else NO, "dotyczy ścieżki CPU (§9)"))
    else:
        rows.append(("§1", "test równoważności (test/equiv.json)", "–", "True", UN, "brak pliku (od 07:20 trening na GPU; test dotyczył ścieżki CPU)"))
    D.table(rows)


# ------------------------------------------------------------------------------------------------ główny
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=os.path.join(HERE, "runs", "lxplus5"))
    ap.add_argument("--out", default="")
    a = ap.parse_args(argv)
    miss, D = [], Doc()
    CMP, RES = load(a.dir, "dnet5_cmp.json", miss), load(a.dir, "dnet5_results.json", miss)
    D.L.append("# WERDYKT DNET etap 5a — odczyt wg KRYTERIA_etap5.md")
    D.p(f"Katalog: `{a.dir}`. Skrypt `DNET/verdict5.py` — sam odczyt i porównanie z progami zapisanymi przed wynikami; żadnych fitów. "
        "Statusy: SPEŁNIONE / NIESPEŁNIONE / NIEROZSTRZYGNIĘTE (brak danych).")
    D.h("0. Kompletność wejścia")
    for nm, j in (("dnet5_cmp.json", CMP), ("dnet5_results.json", RES)):
        D.p(f"- `{nm}`: " + ("**BRAK**" if j is None else f"jest (wygenerowano {j.get('generated', '?')}, nboot {j.get('nboot')})"))
    D.p(f"- `DONE`: {'jest' if os.path.isfile(os.path.join(a.dir, 'DONE')) else '**BRAK**'}")
    fails = sorted(os.path.basename(x) for x in glob.glob(os.path.join(a.dir, "fail", "*")))
    D.p(f"- węzły w `fail/`: {('**' + ', '.join(fails) + '**') if fails else 'brak'}")
    if CMP:
        ms = CMP.get("missing") or []
        D.p("- pole `missing` w dnet5_cmp.json: " + ("puste" if not ms else "**" + "; ".join(f"{m.get('pair') or m.get('cascade')}: brak {', '.join(m.get('tags', []))}" for m in ms) + "**"))
        src = CMP.get("sources") or {}
        D.p("- źródła score odniesienia 100 % (powinny być z etapu 2, §2): " + ("; ".join(f"`{t}` ← {src.get(t, 'BRAK')}" for t in BASE)))
        if (CMP.get("nboot") or 0) < NBOOT_REQ:
            miss.append(f"dnet5_cmp.json: bootstrap {CMP.get('nboot')} replik < {NBOOT_REQ} (§2)")
        for nm, R in (CMP.get("pairs") or {}).items():
            if nm.startswith("etap4_over_etap2"):
                D.p(f"- opisowo `{nm}` (rozrzut powtórzenia, MC WP70): " + ", ".join(f"{b} {pm(*(pair(CMP, nm, b) or (None,)))}" for b in BINS))
    if RES and RES.get("fast"):
        miss.append("dnet5_results.json: tryb --fast (błędy SF z HESSE)")

    D.h("1. Krzywa uczenia (§2, #5)")
    v1, rows = v_lc(CMP)
    D.table(rows)
    D.p(f"- werdykt: {v1}")
    D.h("2. Cechy par w uwadze (§3, #4)")
    v2, rows = v_pair(CMP, RES)
    D.table(rows)
    D.p(f"- werdykt: {v2}")
    vE = {}
    for i, sp in enumerate(("D0", "Dplus")):
        D.h(f"{3 + i}. Tagger otoczenia E, kaskada H × E — {SPN[sp]} (§4, #3)")
        vE[sp], rows = v_E(CMP, sp)
        D.table(rows)
        D.p(f"- werdykt: {vE[sp]}")
    D.h("5. Sonda GPU i infrastruktura (§1; bez progu fizycznego)")
    vg, rows = v_gpu(a.dir, miss)
    D.table(rows)
    infra(a.dir, miss, D)
    D.p(f"- reguła GPU: {vg}")

    D.h("6. CZTERY WERDYKTY (+ reguła GPU)")
    D.p(f"1. **Krzywa uczenia (§2)**: {v1}")
    D.p(f"2. **Cechy par `trp` (§3)**: {v2}")
    D.p(f"3. **Tagger E · D⁰ (§4)**: {vE['D0']}")
    D.p(f"4. **Tagger E · D⁺ (§4)**: {vE['Dplus']}")
    D.p(f"5. **Reguła GPU (§1)**: {vg}")
    if miss or fails:
        D.p("- **braki / uwagi wejścia**: " + "; ".join(miss + ([f"fail/: {', '.join(fails)}"] if fails else [])))
    D.p("- progi bezwzględne, bez poprawki na wielokrotność (§5); wszystko po obejrzeniu = „post hoc” (§7); poza 5a: SF_E, G_fit(H×E/H), domknięcie MC (5b, §6)")
    D.h("Luki odczytu", 3)
    for x in ("§4: gdy po cięciu H zostaje < 50 wierszy tła MC w binie, dn5_ana.py bierze WP_E = 0,80 zamiast wyboru na MC (flaga w uwagach wiersza); przy < 30 wierszach brak G_E(MC) → wybór wariantu niemożliwy.",
              "§4: G_Htight bez błędu w JSON (próg różnicy G_E − G_Htight ≥ 0,10 bezwzględny — porównanie bez niepewności).",
              "§3: |corr(score, m)| dla `tr` — KRYTERIA nie mówią, które ziarno; użyto średniej |corr| dostępnych ziaren 1–3 (w uwagach też różnica względem najmniejszego).",
              "§1: s/epokę CPU (4 wątki) z test/bench.json; gdy pliku brak — liczby zgrubne z KRYTERIA §1 (82 / 107 s)."):
        D.p(f"- {x}")
    txt = "\n".join(D.L) + "\n"
    out = a.out or os.path.join(a.dir, "WERDYKT_etap5.md")
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
