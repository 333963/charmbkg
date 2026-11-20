#!/usr/bin/env python3
"""SD_var/analyze_var.py — niepewność modelowa R = σ(D + p w oknie AFP)/σ_incl(D) z wariantów modelu dyfrakcji Pythii 8 (POZIOM CZĄSTEK).
Kryterium i zbiory ustalone PRZED wynikami: PLAN.md §4 + Q_lxplus/sd_var/scenarios.txt (ten skrypt je tylko wykonuje).
Wejście: merged/<model>.json z sd_var (+ xsec.json) oraz z sdpred (punkty odniesienia). Lekkie: json + numpy, bez generowania.
Użycie: .venv_h5/bin/python charm_work/SD_var_2026-11-19/analyze_var.py [--var DIR] [--ref DIR] [--scen PLIK] [--out DIR]
Wyjście: <out>/tables.md, <out>/results.json.
"""
import argparse, json, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); W = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(W, "SD_pred"))
from analyze import Model, add, tot, ratio, SPL          # ten sam kod odczytu histogramów co w SD_pred  # noqa: E402

PTL = ["2–3", "3–5", "5–8", "8–30"]
# stałe kryterium (PLAN.md §4) — NIE zmieniać po obejrzeniu wyników
STAT_MAX = 0.05                      # wariant z błędem stat. R > 5 % nie wchodzi do obwiedni (oznaczany STAT)
SD_SCALE = (0.64, 1.0)               # skala σ_SD: dane ALFA 8 TeV / A3 (SD_pred §7.4) … bez zmiany; σ_ND kompensuje (σ_inel stałe)
DDCD_SCALE = (0.5, 2.0)              # skala σ_DD i σ_CD (brak pomiaru w oknie; czynnik 2)
XSEC_REF = {"dif_A3": "ref_A3", "dif_Monash": "ref_A3", "dif_sigPp25": "ref_A3", "dif_MBR": "ref_MBR"}   # dif_own_* → ref_own
# przeliczenie na uzysk 427929 (TAGREG_2026-11-18 §3, scenariusz M): N_incl D⁰ WP80 w binach 3–5/5–8/8–30 × A·ε_AFP × ε_ξcal
N_INCL_D0 = [6773, 2010, 813]; AEPS = 0.55; EXICAL = 0.99


def load(d):
    M = {}
    if d and os.path.isdir(d):
        for f in sorted(os.listdir(d)):
            if f.endswith(".json") and f not in ("summary.json", "xsec.json"):
                m = Model(os.path.join(d, f)); M[m.name] = m
    return M


def scen_table(path):
    S = []
    for ln in open(path, encoding="utf-8"):
        if ln.lstrip().startswith("#") or not ln.strip(): continue
        c = [x.strip() for x in ln.split("|")]
        S.append(dict(name=c[0], kind=c[1], num=c[2], den=c[3], set=c[4], grp=c[5], ref=c[6], desc=c[7]))
    return S


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--var", default=os.path.join(HERE, "runs/lxplus/merged")); ap.add_argument("--ref", default=os.path.join(W, "SD_pred/runs/lxplus/merged"))
    ap.add_argument("--scen", default=os.path.join(W, "Q_lxplus/sd_var/scenarios.txt")); ap.add_argument("--out", default=os.path.join(HERE, "out"))
    ap.add_argument("--over", default=os.path.join(HERE, "runs/lxplus/merged_close"))   # domknięcie 20.11: modele z dosłanymi częściami (dv_CRqcd, nv_CRqcd) nadpisują te z --var; brak katalogu → bez zmian
    a = ap.parse_args()
    M = load(a.ref); M.update(load(a.var)); over = load(a.over); M.update(over)
    xs = json.load(open(os.path.join(a.var, "xsec.json"))) if os.path.exists(os.path.join(a.var, "xsec.json")) else {}
    S = scen_table(a.scen); os.makedirs(a.out, exist_ok=True)

    def kND(s):      # σ_ND konfiguracji licznika / σ_ND próbki ND (σ_inel − Σσ_dyfr. zależy od modelu dyfrakcji)
        key = s["num"] if s["num"] in xs else XSEC_REF.get(s["num"], "ref_own" if s["num"].startswith("dif_own_") else None)
        if key not in xs or "ND" not in xs[key]: return 1.0, False
        return xs[key]["ND"] * 1e3 / M[s["den"]].sig, True

    def comps(s, sp, name, sc=None):
        """składowe (σ, błąd) [µb] licznika `name` i mianownik inkluzywny; sc = skale składowych {SD, DD, CD} (ND kompensuje σ_inel)."""
        m, den = M[s["num"]], M[s["den"]]; sc = sc or {}
        if s["kind"] == "hard":
            d = {"HQdif": tot(m.get(["HQdif", "HQdd"], [sp], name)), "HQnd": tot(den.get(["HQnd"], [sp], name))}
            return d, tot(den.get(["HQnd"], [sp], "pt_incl"))
        k, _ = kND(s)
        if sc:       # zmiana σ_dyfr. o Δ przesuwa σ_ND o −Δ
            sig = {c: sum(v["n"] for q, v in m.proc.items() if q[:3] in cc) / m.nok * m.sig for c, cc in (("SD", ("103", "104")), ("DD", ("105",)), ("CD", ("106",)))}
            snd = den.sig * k; k *= max(snd - sum(sig[c] * (sc.get(c, 1) - 1) for c in sig), 0) / snd
        d = {}; inc = []
        for c in ("SD", "DD", "CD"):
            f = sc.get(c, 1.0); x = tot(m.get([c], [sp], name)); y = tot(m.get([c], [sp], "pt_incl"))
            d[c] = (None, None) if x[0] is None else (x[0] * f, x[1] * f); inc.append((None, None) if y[0] is None else (y[0] * f, y[1] * f))
        x = tot(den.get(["ND"], [sp], name)); y = tot(den.get(["ND"], [sp], "pt_incl"))
        d["ND"] = (None, None) if x[0] is None else (x[0] * k, x[1] * k); inc.append((None, None) if y[0] is None else (y[0] * k, y[1] * k))
        return d, add(*inc)

    def R_of(s, sp, sc=None, name="pt_tag"):
        d, i = comps(s, sp, name, sc); t = add(*d.values())
        if t[0] is None or i[0] is None or not i[0] > 0: return None
        r = t[0] / i[0]; e = r * math.hypot(t[1] / t[0] if t[0] > 0 else 0, i[1] / i[0])
        nsd = 1 - (d.get("SD", d.get("HQdif"))[0] or 0) / t[0] if t[0] > 0 else float("nan")
        return dict(R=r, err=e, sigma_tag_ub=t[0], sigma_incl_ub=i[0], non_main=nsd, comp={c: v[0] for c, v in d.items()})

    def R_pt(s, sp):
        m, den = M[s["num"]], M[s["den"]]; k = 1.0 if s["kind"] == "hard" else kND(s)[0]
        if s["kind"] == "hard": t = add(m.get(["HQdif", "HQdd"], [sp], "pt_tag"), den.get(["HQnd"], [sp], "pt_tag")); i = den.get(["HQnd"], [sp], "pt_incl")
        else:
            nt, ni = den.get(["ND"], [sp], "pt_tag"), den.get(["ND"], [sp], "pt_incl")
            t = add(m.get(["SD", "DD", "CD"], [sp], "pt_tag"), (nt[0] * k, nt[1] * k) if nt[0] is not None else nt); i = add(m.get(["SD", "DD", "CD"], [sp], "pt_incl"), (ni[0] * k, ni[1] * k))
        r = ratio(t, i); return (None, None) if r[0] is None else ([float(x) for x in r[0]], [float(x) for x in r[1]])

    res = {"scenarios": {}, "missing": [], "criterion": {"STAT_MAX": STAT_MAX, "SD_SCALE": SD_SCALE, "DDCD_SCALE": DDCD_SCALE}}; byname = {s["name"]: s for s in S}
    for s in S:
        if s["num"] not in M or s["den"] not in M: res["missing"].append(s["name"]); continue
        o = dict(kind=s["kind"], set=s["set"], grp=s["grp"], ref=s["ref"], desc=s["desc"], kND=kND(s)[0] if s["kind"] == "soft" else None, kND_from_xsec=kND(s)[1] if s["kind"] == "soft" else None)
        for sp in SPL:
            r = R_of(s, sp)
            if r is None: continue
            rp = R_pt(s, sp); r["R_pt"], r["R_pt_err"] = rp
            if s["kind"] == "soft":
                r["R_t1"] = (R_of(s, sp, name="pt_tag_t1") or {}).get("R")
            o[sp] = r
        if "D0" in o and o["D0"].get("R_pt"):
            o["N427929_D0_3_30"] = float(sum(n * r for n, r in zip(N_INCL_D0, o["D0"]["R_pt"][1:])) * AEPS * EXICAL)
        res["scenarios"][s["name"]] = o
    SC = res["scenarios"]

    # ---- niepewność: obwiednia w ln R po zbiorze A, osobno w grupach; M i P nie sumowane (większa z dwóch); + X, L, normalizacja σ w kwadraturze
    def env(grp, sp, sets=("A",)):
        up = dn = 0.0; used = []
        for n, o in SC.items():
            if o["grp"] != grp or o["set"] not in sets or sp not in o: continue
            ref = SC.get(o["ref"], {}).get(sp)
            if ref is None or not o[sp]["R"] > 0: continue
            if o[sp]["err"] / o[sp]["R"] > STAT_MAX: used.append(n + " (STAT, pominięty)"); continue
            x = math.log(o[sp]["R"] / ref["R"]); up = max(up, x); dn = min(dn, x); used.append(f"{n} {math.exp(x):.3f}")
        return up, dn, used

    unc = {}
    for sp in SPL:
        nom = next((s for s in S if s["set"] == "N" and s["kind"] == "soft" and s["name"] in SC and sp in SC[s["name"]]), None)
        if nom is None: continue
        g = {k: env(k, sp) for k in "MPXL"}
        r0 = SC[nom["name"]][sp]["R"]; nrm = [0.0, 0.0]
        for sc in ({"SD": SD_SCALE[0]}, {"DD": DDCD_SCALE[0]}, {"DD": DDCD_SCALE[1]}, {"CD": DDCD_SCALE[0]}, {"CD": DDCD_SCALE[1]}):
            x = math.log(R_of(nom, sp, sc)["R"] / r0); key = "+".join(f"{k}×{v}" for k, v in sc.items()); g.setdefault("norm_detail", {})[key] = math.exp(x)
            nrm[0] = math.hypot(nrm[0], max(x, 0)); nrm[1] = math.hypot(nrm[1], min(x, 0))
        shape = (max(g["M"][0], g["P"][0]), min(g["M"][1], g["P"][1]))
        up = math.sqrt(shape[0] ** 2 + g["X"][0] ** 2 + g["L"][0] ** 2 + nrm[0] ** 2); dn = math.sqrt(shape[1] ** 2 + g["X"][1] ** 2 + g["L"][1] ** 2 + nrm[1] ** 2)
        unc[sp] = dict(nominal=nom["name"], R=r0, R_stat=SC[nom["name"]][sp]["err"], factor_up=math.exp(up), factor_dn=math.exp(-dn), R_hi=r0 * math.exp(up), R_lo=r0 * math.exp(-dn),
                       groups={k: dict(up=math.exp(g[k][0]), dn=math.exp(g[k][1]), used=g[k][2]) for k in "MPXL"}, norm=dict(up=math.exp(nrm[0]), dn=math.exp(-nrm[1]), detail=g["norm_detail"]),
                       envelope_B={k: [math.exp(v) for v in env(k, sp, ("A", "B"))[:2]] for k in "MPX"})
        h = env("H", sp); hn = next((s["name"] for s in S if s["set"] == "N" and s["kind"] == "hard" and s["name"] in SC and sp in SC[s["name"]]), None)
        if hn: unc[sp]["hard"] = dict(nominal=hn, R=SC[hn][sp]["R"], factor_up=math.exp(h[0]), factor_dn=math.exp(h[1]), used=h[2])
    res["uncertainty"] = unc
    if "D0" in unc and "N427929_D0_3_30" in SC.get(unc["D0"]["nominal"], {}):
        n0 = SC[unc["D0"]["nominal"]]["N427929_D0_3_30"]; res["N427929_D0"] = dict(nominal=n0, lo=n0 * unc["D0"]["factor_dn"], hi=n0 * unc["D0"]["factor_up"], note="N_incl × R(pT) × 0,55 × 0,99; A·ε_AFP nie wariowane")

    L = ["# SD_var: niepewność modelowa R = σ(D + p)/σ_incl(D) — POZIOM CZĄSTEK (kryterium: PLAN.md §4)\n",
         f"Wejście: `{a.var}`, `{a.ref}`; scenariusze `{a.scen}`. Brak danych dla: {', '.join(res['missing']) or '—'}.\n",
         f"Nadpisane z `{a.over}` (domknięcie): {', '.join(f'{k} ({m.nj} części, {m.nok} zd.)' for k, m in over.items()) or '—'}.\n",
         "## 1. Scenariusze\n", "| scenariusz | zbiór | grupa | k_ND | R(D⁰) [%] | R/odniesienie | udział nie-SD [%] | R(D⁺)/R(D⁰) | R(D⁰) w binach pT 2–3 / 3–5 / 5–8 / 8–30 [%] | N(D⁰, 427929) | opis |", "|" + "---|" * 11]
    for n, o in SC.items():
        if "D0" not in o: continue
        d = o["D0"]; ref = SC.get(o["ref"], {}).get("D0"); dp = o.get("Dp")
        L.append(f"| {n} | {o['set']} | {o['grp']} | {'–' if o['kND'] is None else format(o['kND'], '.3f') + ('' if o['kND_from_xsec'] else ' (brak xsec)')} | {100 * d['R']:.3f} ± {100 * d['err']:.3f} | "
                 f"{'–' if ref is None else format(d['R'] / ref['R'], '.3f')} | {100 * d['non_main']:.1f} | {'–' if not dp else format(dp['R'] / d['R'], '.3f')} | "
                 f"{' / '.join(format(100 * x, '.3f') for x in d['R_pt']) if d.get('R_pt') else '–'} | {o.get('N427929_D0_3_30', float('nan')):.1f} | {o['desc']} |")
    L += ["\n## 2. Niepewność (obwiednia zbioru A w grupach; czynniki mnożące R nominalne)\n", "| gatunek | R nom [%] | M (pełne modele) | P (parametry kształtu) | X (układ X) | L (proton ND) | normalizacja σ | **razem** | R zakres [%] | obwiednia A+B (M / P / X) |", "|" + "---|" * 10]
    f2 = lambda u, d: f"×{d:.2f} … ×{u:.2f}"
    for sp, u in unc.items():
        g = u["groups"]
        L.append(f"| {SPL[sp]} | {100 * u['R']:.3f} ± {100 * u['R_stat']:.3f} | " + " | ".join(f2(g[k]["up"], g[k]["dn"]) for k in "MPXL") + f" | {f2(u['norm']['up'], u['norm']['dn'])} | **{f2(u['factor_up'], u['factor_dn'])}** | "
                 f"{100 * u['R_lo']:.3f} – {100 * u['R_hi']:.3f} | " + " / ".join(f2(u["envelope_B"][k][0], u["envelope_B"][k][1]) for k in "MPX") + " |")
    L.append("\nScenariusz twardy (osobno, nie wchodzi do sumy): " + "; ".join(f"{SPL[sp]} R = {100 * u['hard']['R']:.3f} % {f2(u['hard']['factor_up'], u['hard']['factor_dn'])}" for sp, u in unc.items() if "hard" in u))
    if "N427929_D0" in res: L.append(f"\nUzysk D⁰ w 427929 (pT 3–30, scenariusz M z TAGREG, tylko niepewność modelu R): {res['N427929_D0']['nominal']:.1f} ({res['N427929_D0']['lo']:.1f} – {res['N427929_D0']['hi']:.1f}).")
    open(os.path.join(a.out, "tables.md"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    json.dump(res, open(os.path.join(a.out, "results.json"), "w"), indent=1, ensure_ascii=False, default=float)
    print(f"OK: {len(SC)} scenariuszy ({len(res['missing'])} bez danych) → {a.out}/tables.md")


if __name__ == "__main__":
    main()
