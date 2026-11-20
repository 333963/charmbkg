#!/usr/bin/env python3
"""Tabele do raportu z out/s0..s4.json → out/tables.md; werdykty K1–K6 (pulle), średnia po składzie strumienia AFP."""
import json
import numpy as np
import afp_common as A

s0, s1, s2, s3, s4 = (json.load(open(f"{A.OUT}/s{i}.json")) for i in range(5))
L = []
P = L.append
f2 = lambda v, e: f"{v:.2f} ± {e:.2f}" if v is not None and np.isfinite(v) else "–"   # noqa
f3 = lambda v, e: f"{v:.3f} ± {e:.3f}" if v is not None and np.isfinite(v) else "–"   # noqa
pull = lambda a, ea, b, eb: (a - b) / np.hypot(ea, eb)                                    # noqa


def grade(S, Serr):
    Sf = S / 0.82
    if Sf < 50 or S < 3 * Serr:
        return "niemiarodajna"
    return "słaba" if Sf < 150 else "miarodajna"


CLS = ["AFPnoalg", "AFPsptrk4", "AFPboth", "AFP&J12", "AFP&MBTS", "J12only!A", "MBTS!A", "none"]
EXTRA = ["AFPonly", "AFPnoalg_any", "AFPstream", "J12only", "MBTSbased", "other", "all"]
# ---- T1 liczności
P("## T1. Klasy: liczności (428770)\n")
P("| klasa | zdarzenia z kandydatem | udział | pary D⁰←D* | kandydaci all_* (30 %) | N_sig SBS (3–30 GeV, score > 0,6) | N_sig fit | miarodajność |")
P("|---|---|---|---|---|---|---|---|")
ev, pr, al = s0["428770"]["zdarzenia"], s0["428770"]["pary"], s0["428770"]["all"]
for k in CLS + EXTRA:
    r = s1["rows"][f"8:{k}|3-30"]; c = r["cuts"]["ptcone40<2"]
    fit = f"{c['fit']['N']:.0f} ± {c['fit']['N_err']:.0f}" if "fit" in c else "–"
    nev = ev["cls"].get(k, ev["bit"][A.AN] if k == "AFPnoalg_any" else None); npr = pr["cls"].get(k, pr["bit"][A.AN]); nal = al["cls"].get(k, al["bit"][A.AN])
    P(f"| {k} | {nev} | {100 * nev / ev['n']:.2f} % | {npr} | {nal} | {r['S']:.0f} ± {r['S_err']:.0f} | {fit} | {grade(r['S'], r['S_err'])} |")
P("")
# ---- T2 izolacja
ref9 = s1["rows"]["9:MB427929|3-30"]["cuts"]; refJ = s1["rows"]["8:J12only!A|3-30"]["cuts"]
ROWS = ["9:MB427929", "9:MB9&!j20&!J12", "9:MB9&j20&!J12", "9:MB9&AFPnoalg", "8:AFPnoalg_any", "8:AFPnoalg", "8:AFPsptrk4", "8:AFPboth", "8:AFP&J12", "8:AFP&MBTS",
        "8:AFPonly", "8:AFPstream", "8:none", "8:other", "8:MBTS!A", "8:J12only!A", "8:all"]
P("## T2. SF izolacji sygnału (SBS; w nawiasie fit gdzie N_sig ≥ 150), pT 3–30 GeV\n")
cuts = s1["cuts"]
P("| próbka | N_sig SBS | " + " | ".join(f"{c} bez wag / przeważone" for c in cuts) + " |")
P("|---|---|" + "---|" * len(cuts))
for k in ROWS:
    r = s1["rows"][k + "|3-30"]
    cells = []
    for c in cuts:
        x = r["cuts"][c]
        if "SF" not in x or not np.isfinite(x["SF"]):
            cells.append("–"); continue
        s = f2(x["SF"], x["SF_err"])
        if "SF_fit" in x:
            s += f" ({f2(x['SF_fit'], x['SF_fit_err'])})"
        s += " / " + (f2(x["rw"]["SF"], x["rw"]["err"]) if "rw" in x else "–")
        cells.append(s)
    P(f"| {k} | {r['S']:.0f} ± {r['S_err']:.0f} | " + " | ".join(cells) + " |")
P("")
P("## T2b. SF izolacji ptcone40 < 2 GeV w binach pT (SBS; bez wag / przeważone)\n")
P("| próbka | 3–5 (N) | 3–5 SF | 5–8 (N) | 5–8 SF | 8–30 (N) | 8–30 SF |"); P("|---|---|---|---|---|---|---|")
for k in ROWS:
    cells = []
    for pb in ("3-5", "5-8", "8-30"):
        r = s1["rows"][f"{k}|{pb}"]; x = r["cuts"]["ptcone40<2"]
        cells.append(f"{r['S']:.0f} ± {r['S_err']:.0f}")
        cells.append((f2(x["SF"], x["SF_err"]) + " / " + (f2(x["rw"]["SF"], x["rw"]["err"]) if "rw" in x else "–")) if "SF" in x and np.isfinite(x["SF"]) and r["S"] > 3 * r["S_err"] else "–")
    P(f"| {k} | " + " | ".join(cells) + " |")
P("")
# ---- T3 tło
P("## T3. Kontrola na tle (pary, Δm 150–170 MeV, |m−m₀| < 40 MeV, bez score), pT 3–30 GeV: frakcja izolowana\n")
P("| próbka | N | f(ptcone40 < 2) | f / f_MC | f / f_MC przeważone | f(ptcone40 < 1) R / R_rw | f(ptcone40 < 4) R / R_rw | f(rel40 < 0,5) R / R_rw | f(ptcone100 < 10) R / R_rw |")
P("|---|---|---|---|---|---|---|---|---|")
for k in ROWS:
    r = s1["rows"][k + "|3-30"]; b = r["cuts"]["ptcone40<2"].get("bkg")
    if not b:
        continue
    oth = []
    for c in ("ptcone40<1", "ptcone40<4", "rel40<0.5", "ptcone100<10"):
        bb = r["cuts"][c]["bkg"]; oth.append(f"{bb['R']:.2f} / {bb['R_rw']:.2f} ± {bb['R_rw_err']:.2f}")
    P(f"| {k} | {b['n']} | {b['f']:.4f} ± {b['err']:.4f} | {b['R']:.3f} ± {b['R_err']:.3f} | {b['R_rw']:.3f} ± {b['R_rw_err']:.3f} | " + " | ".join(oth) + " |")
P(f"\nMC mb_nd (pary, ta sama selekcja, N = {s1['mc']['3-30']['n_bkg']}): f(ptcone40 < 2) = {s1['mc']['3-30']['f_bkg']['ptcone40<2']:.4f}.\n")
# ---- T4 SF H
P("## T4. SF klasyfikatora H (SBS; w nawiasie fit profilowy gdzie N_sig(SBS) ≥ 120)\n")
P("| próbka | bin | N_sig SBS | miarodajność | WP95 | WP80 | WP50 |"); P("|---|---|---|---|---|---|---|")
HR = ["9:MB427929", "8:AFPstream", "8:AFPonly", "8:AFPsptrk4", "8:AFPnoalg_any", "8:AFP&J12", "8:none", "8:other", "8:MBTS!A", "8:J12only!A", "8:all"]
for b in ("5-8", "8-30"):
    for k in HR:
        cells = []; S = Se = None
        for wp in ("0.95", "0.80", "0.50"):
            r = s2["rows"][f"{k}|{b}|{wp}"]
            if "S" not in r or not np.isfinite(r.get("SF", np.nan)):
                cells.append("–"); continue
            S, Se = r["S"], r["S_err"]
            s = f2(r["SF"], r["SF_err"])
            if "SF_fit" in r and r["fit"]["det"]:
                s += f" ({f2(r['SF_fit'], r['SF_fit_err'])})"
            cells.append(s)
        P(f"| {k} | {b} | " + (f"{S:.0f} ± {Se:.0f} | {grade(S, Se)}" if S is not None else "– | niemiarodajna") + " | " + " | ".join(cells) + " |")
P("")
# ---- T5 ε_B
P("## T5. ε_B pod pikiem: R = ε_B(dane, exp(pol2) z SB) / ε_B(MC truth ≥ 2, zliczenie w SR); ± stat (bootstrap) [syst. interpolacji wzgl.]; oraz R_SB = frakcja pass w SB dane / MC inkluzywne\n")
P("| próbka | bin | N w SB | R WP95 | R WP80 | R WP50 | R_SB WP95 | R_SB WP80 | R_SB WP50 |"); P("|---|---|---|---|---|---|---|---|---|")
ER = ["9:MB427929", "9:MB9&!j20&!J12", "9:MB9&j20&!J12", "8:AFPnoalg_any", "8:AFPnoalg", "8:AFPsptrk4", "8:AFP&J12", "8:AFPstream", "8:none", "8:MBTS!A", "8:J12only!A", "8:all"]
for b in ("3-5", "5-8", "8-30"):
    for k in ER:
        r = s3["rows"][f"{k}|{b}"]
        if "pol2" not in r:
            P(f"| {k} | {b} | {r['n_sb']} | – | – | – | – | – | – |"); continue
        c1 = [f"{r['R_mc_t2'][j]:.3f} ± {r['R_mc_t2_stat'][j]:.3f} [{r['syst_rel'][j]:.2f}]" if r["R_mc_t2_stat"][j] < 0.5 else "–" for j in (0, 2, 4)]
        c2 = [f"{r['Rsb_mc'][j]:.3f} ± {r['Rsb_mc_err'][j]:.3f}" for j in (0, 2, 4)]
        P(f"| {k} | {b} | {r['n_sb']} | " + " | ".join(c1 + c2) + " |")
P("")
# ---- T6 aktywność
P("## T6. Aktywność: zdarzenia z kandydatem (ev_*)\n")
P("| próbka | N zdarzeń | ⟨nSelBL⟩ | mediana nSelBL | ⟨lead_pt⟩ [GeV] | mediana lead_pt | ułamek lead_pt < 4 GeV | ułamek nPileupVertices = 0 |"); P("|---|---|---|---|---|---|---|---|")
for k, r in s4["ev"].items():
    if not r:
        continue
    P(f"| {k} | {r['nSelBL']['n']} | {r['nSelBL']['mean']:.1f} | {r['nSelBL']['med']:.0f} | {r['lead_pt']['mean']:.2f} | {r['lead_pt']['med']:.2f} | " +
      (f"{r['frac_lead_lt4']:.3f} | {r['frac_npu0']:.3f}" if "frac_lead_lt4" in r else "– | –") + " |")
P("")
P("## T7. Kandydaci tła (pasma boczne masy, all_*, pT 3–30): średnie\n")
P("| próbka | N | ⟨pT⟩ | ⟨ptcone40⟩ | ⟨ptcone100⟩ | ⟨nSelBL⟩ | ⟨lead_pt⟩ | ⟨score⟩ |"); P("|---|---|---|---|---|---|---|---|")
for k, r in s4["bkg"].items():
    P(f"| {k} | {r['pt']['n']} | {r['pt']['mean']:.2f} | {r['ptcone40']['mean']:.2f} | {r['ptcone100']['mean']:.2f} | {r['e_nSelBL']['mean']:.1f} | {r['e_lead_pt']['mean']:.2f} | {r['score']['mean']:.4f} |")
P("")
P("## T8. Kandydaci sygnału D* (wagi SBS, score > 0,6, pT 3–30): średnie ± błąd\n")
P("| próbka | N_sig | ⟨pT⟩ | ⟨ptcone40⟩ | ⟨ptcone100⟩ | ⟨nSelBL⟩ | ⟨lead_pt⟩ | ⟨ptcone40/pT⟩ |"); P("|---|---|---|---|---|---|---|---|")
for k, r in s4["sig"].items():
    P(f"| {k} | {r['pt']['n']:.0f} | " + " | ".join(f2(r[v]["mean"], r[v]["mean_err"]) for v in ("pt", "ptcone40", "ptcone100", "e_nSelBL", "e_lead_pt", "rel40")) + " |")
P("")
P("## T9. C2ST (AUC, HGB 5-fold GroupKFold po zdarzeniu; pasma boczne, pT 3–30) i KS\n")
P("| para | n / klasę | AUC zdarzenie(3) | AUC kandydat(3) | AUC wszystkie(6) | AUC +score(7) | KS nSelBL | KS lead_pt | KS ptcone40 | KS ptcone100 | KS pT | KS score |"); P("|---|---|---|---|---|---|---|---|---|---|---|---|")
pairs = sorted({k.rsplit("|", 1)[0] for k in s4["c2st"] if not k.startswith("NULL")}, key=lambda s: (s.split("|")[1], s.split("|")[0]))
for p in pairs:
    c = [s4["c2st"].get(f"{p}|{fn}") for fn in ("zdarzenie(3)", "kandydat(3)", "wszystkie(6)", "wszystkie+score(7)")]
    ks = s4["ks"].get(p, {})
    P(f"| {p.replace('|', ' vs ')} | {c[0]['n']} | " + " | ".join(f"{x['auc']:.3f} ± {x['err']:.3f}" for x in c) + " | " +
      " | ".join(f"{ks.get(v, float('nan')):.3f}" for v in ("e_nSelBL", "e_lead_pt", "ptcone40", "ptcone100", "pt", "score")) + " |")
n = s4["c2st"]["NULL 427929 polowy|wszystkie(6)"]
P(f"\nTest zerowy (dwie losowe, rozłączne części 427929, n = {n['n']}): AUC = {n['auc']:.3f} ± {n['err']:.3f}.\n")
# ---- werdykty
P("## T10. Werdykty wg kryteriów §0\n")
m9 = ref9["ptcone40<2"]; mJ = refJ["ptcone40<2"]
P(f"Odniesienia (SBS, ten sam estymator): 427929 SF = {f2(m9['SF'], m9['SF_err'])}; J12only¬A SF = {f2(mJ['SF'], mJ['SF_err'])}; tło f: {m9['bkg']['f']:.3f} / {mJ['bkg']['f']:.3f}.\n")
P("| klasa | miarodajność | SF bez wag | pull vs 427929 | pull vs J12only | pull vs 1 | K1 | f tła | K1b | SF przeważone (pull vs 1) | tło R_rw: klasa − 427929 [σ] |"); P("|---|---|---|---|---|---|---|---|---|---|---|")
for k in ["AFPnoalg_any", "AFPnoalg", "AFPsptrk4", "AFPboth", "AFP&J12", "AFP&MBTS", "AFPonly", "AFPstream", "none", "MBTS!A"]:
    r = s1["rows"][f"8:{k}|3-30"]; x = r["cuts"]["ptcone40<2"]; g = grade(r["S"], r["S_err"])
    p9, pJ, p1 = pull(x["SF"], x["SF_err"], m9["SF"], m9["SF_err"]), pull(x["SF"], x["SF_err"], mJ["SF"], mJ["SF_err"]), (x["SF"] - 1) / x["SF_err"]
    if g == "niemiarodajna":
        v = "nierozstrzygnięte (niemiarodajna)"
    elif abs(p9) < 2 and abs(pJ) > 3:
        v = "minbias-podobna"
    elif abs(pJ) < 2 and abs(p1) > 3:
        v = "dżetowa"
    elif abs(p9) > 2 and abs(pJ) > 2:
        v = "pośrednia"
    else:
        v = "nierozstrzygnięte"
    if g == "słaba":
        v += " (słaba: tylko wskazanie)"
    b = x["bkg"]; vb = "minbias-podobne" if abs(b["f"] - 0.480) < 0.05 else ("dżetowe" if abs(b["f"] - 0.158) < 0.05 else "pośrednie")
    rw = x.get("rw"); b9 = m9["bkg"]
    P(f"| {k} | {g} | {f2(x['SF'], x['SF_err'])} | {p9:+.1f} | {pJ:+.1f} | {p1:+.1f} | {v} | {b['f']:.3f} ± {b['err']:.3f} | {vb} | " +
      (f"{f2(rw['SF'], rw['err'])} ({(rw['SF'] - 1) / rw['err']:+.1f})" if rw else "–") + f" | {pull(b['R_rw'], b['R_rw_err'], b9['R_rw'], b9['R_rw_err']):+.1f} |")
P("")
# ---- K6: średnia po składzie strumienia AFP
P("## T11. Strumień AFP: skład i średnie ważone (ptcone40 < 2 GeV, pT 3–30)\n")
P("| klasa | udział zdarzeń | udział par | udział N_sig (SBS) | SF sygnału bez wag | tło f/f_MC | tło f/f_MC przeważone |"); P("|---|---|---|---|---|---|---|")
tot_ev = ev["cls"]["AFPstream"]; tot_p = pr["cls"]["AFPstream"]
Ssum = sum(max(s1["rows"][f"8:{k}|3-30"]["S"], 0) for k in A.AFP_STREAM)
acc = []
for k in A.AFP_STREAM:
    r = s1["rows"][f"8:{k}|3-30"]; x = r["cuts"]["ptcone40<2"]; b = x["bkg"]
    P(f"| {k} | {100 * ev['cls'][k] / tot_ev:.1f} % | {100 * pr['cls'][k] / tot_p:.1f} % | {100 * max(r['S'], 0) / Ssum:.0f} % | {f2(x['SF'], x['SF_err'])} | {b['R']:.3f} ± {b['R_err']:.3f} | {b['R_rw']:.3f} ± {b['R_rw_err']:.3f} |")
    acc.append((ev["cls"][k] / tot_ev, pr["cls"][k] / tot_p, max(r["S"], 0) / Ssum, x["SF"], x["SF_err"], b["R"], b["R_err"], b["R_rw"], b["R_rw_err"], b["n"]))
a = np.array(acc)
x = s1["rows"]["8:AFPstream|3-30"]["cuts"]["ptcone40<2"]
wS = a[:, 2]; sfw = (wS * a[:, 3]).sum(); sfe = np.sqrt((wS ** 2 * a[:, 4] ** 2).sum())
chi = (((a[:, 3] - x["SF"]) / a[:, 4]) ** 2).sum(); ndf = len(a) - 1
P(f"\n- SF sygnału, próbka połączona (= wagi N_sig z konstrukcji): **{f2(x['SF'], x['SF_err'])}**; przeważone: {f2(x['rw']['SF'], x['rw']['err'])}.")
P(f"- SF sygnału, średnia klas z wagami N_sig: {f2(sfw, sfe)}; rozrzut klas wokół wartości połączonej χ²/ndf = {chi:.1f}/{ndf} → błąd × {max(1, np.sqrt(chi / ndf)):.2f} = {sfe * max(1, np.sqrt(chi / ndf)):.2f}.")
for nm, w in (("udział par", a[:, 1]), ("udział zdarzeń", a[:, 0])):
    P(f"- SF sygnału, średnia klas z wagami = {nm}: {f2((w * a[:, 3]).sum(), np.sqrt((w ** 2 * a[:, 4] ** 2).sum()))}.")
wb = a[:, 9] / a[:, 9].sum()
for j, nm in ((5, "bez wag"), (7, "przeważone")):
    mu = (wb * a[:, j]).sum(); e = np.sqrt((wb ** 2 * a[:, j + 1] ** 2).sum()); sp = np.sqrt((wb * (a[:, j] - mu) ** 2).sum())
    P(f"- Tło f/f_MC {nm}, wagi = liczba kandydatów tła: {mu:.3f} ± {e:.3f} (stat), rozrzut klas (RMS ważony) {sp:.3f}.")
P("\nSF izolacji strumienia AFP dla pozostałych cięć i binów — T2, T2b (wiersz 8:AFPstream).")
open(f"{A.OUT}/tables.md", "w").write("\n".join(L) + "\n")
print("ok", len(L))
