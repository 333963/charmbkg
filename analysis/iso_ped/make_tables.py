#!/usr/bin/env python3
"""ISO_ped: tabele (out/tables.md) i werdykt (verdict.json) z out/iso_ped.json wg KRYTERIA.md."""
import json, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, "out", "iso_ped.json")))
CLS = ["9:MB", "8:all", "8:J12only", "8:MBTSbased", "8:other", "8:notJ12"]
BINS = ["3-5", "5-8", "8-30", "3-30"]
VARS = ["raw", "V1", "V2", "V3", "V1r"]
L = []
f = lambda v, e, n=3: "–" if v is None or not np.isfinite(v) else f"{v:.{n}f} ± {e:.{n}f}".replace(".", ",")  # noqa
g = lambda v, n=3: "–" if v is None or not np.isfinite(v) else f"{v:.{n}f}".replace(".", ",")  # noqa

L.append("## 0. Odtworzenie odniesienia (gałąź ptcone40 < 2 GeV, MC bez wag, pT 3–30)\n")
L.append("| klasa | ε fit | ε_MC | SF | det/valid | oczekiwane (TRKJET2) |\n|---|---|---|---|---|---|")
exp = {"9:MB": "0,98 ± 0,04", "8:all": "0,51 ± 0,02", "8:J12only": "0,46 ± 0,03"}
for n, r in R["ref"].items():
    L.append(f"| {n} | {f(r['fit']['eps'], r['fit']['err'])} | {g(r['eps_mc'])} | {f(r['SF'], r['SF_err'])} | {r['fit']['det']}/{r['fit']['valid']} | {exp[n]} |")

L.append("\n## 1. Progi przy równej ε_S^MC (MC bez wag)\n")
L.append("| bin | ε_S^MC | n_MC | raw | V1 | V2 | V3 | V1r |\n|---|---|---|---|---|---|---|---|")
for b in BINS:
    t = R["thr"][b]
    L.append(f"| {b} | {g(t['eps_mc'])} | {t['n_mc']} | 2 | " + " | ".join(g(t[v]) for v in VARS[1:]) + " |")
L.append("\nMapy P(n) [GeV] (biny ntrk05 " + str(R["NE"][:-1]) + "+):\n")
for k in ("427929", "428770", "MC"):
    L.append(f"- {k}: " + ", ".join(g(x, 2) for x in R["maps"][k]) + "  (n par: " + ", ".join(str(int(x)) for x in R["maps"][k + "_n"]) + ")")

L.append("\n## 2. SF per klasa / bin / wariant\n")
L.append("SF = ε_dane(best)/ε_MC(pT). Błąd SF: fit/SBS ⊕ MC. D = |SF_raw−1| − |SF_var−1|; σ_D z bootstrapu (SBS). „f” = fit, „S” = SBS (fit det/valid = False).\n")
L.append("| klasa | bin | N_sig | jakość | wariant | ε dane | met. | ε_MC | SF | SF (SBS) | D | D (SBS) | σ_D | D/σ_D |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for n in CLS:
    for b in BINS:
        r = R["cells"][f"{n}|{b}"]
        for v in VARS:
            x = r["var"][v]
            d = "" if v == "raw" else f"{g(x['D'])} | {g(x['D_sbs'])} | {g(x['D_err'])} | {g(x['D']/x['D_err'], 1) if x['D_err'] > 0 else '–'}"
            if v == "raw":
                d = "– | – | – | –"
            L.append(f"| {n} | {b} | {f(r['N_sig'], r['N_err'], 0)} | {r['quality']} | {v} | {f(x['eps'], x['eps_err'])} | {'f' if x['method']=='fit' else 'S'} | {g(x['eps_mc'])} | "
                     f"{f(x['SF'], x['SF_err'])} | {g(x['SF_sbs'])} | {d} |")

L.append("\nFity z det/valid = False (pokazane, w tabeli użyto SBS):\n")
bad = []
for n in CLS:
    for b in BINS:
        for v in VARS:
            x = R["cells"][f"{n}|{b}"]["var"][v]
            if not x["fit_ok"]:
                ft = x.get("fit", {})
                bad.append(f"- {n} {b} {v}: ε_fit {g(ft.get('eps'))} ± {g(ft.get('err'))}, det {ft.get('det')}, valid {ft.get('valid')}, {ft.get('reason','')}; SBS {f(x['sbs']['eps'], x['sbs']['err'])}")
L += bad if bad else ["- brak"]

L.append("\n## 3. Tło: ε_B na pasmach bocznych masy D⁰ (dane), MC mb_nd truth ≥ 8; zysk Z\n")
L.append("| klasa | bin | n_B dane | n_B MC | wariant | ε_B dane | ε_B MC | dane/MC | ε_B var/raw (dane) | var/raw (MC) | Z_var/Z_raw (dane) | Z_var/Z_raw (MC) |\n|---|---|---|---|---|---|---|---|---|---|---|---|")
for n in CLS:
    for b in BINS:
        r = R["cells"][f"{n}|{b}"]
        for v in VARS:
            x = r["var"][v]
            ra = x["epsB"] / x["epsB_mc"] if x["epsB_mc"] else np.nan
            rae = ra * np.hypot(x["epsB_err"] / max(x["epsB"], 1e-9), x["epsB_mc_err"] / max(x["epsB_mc"], 1e-9)) if np.isfinite(ra) else np.nan
            ex = "– | – | – | –" if v == "raw" else f"{f(x['rB'], x['rB_err'])} | {g(x['rB_mc'])} | {f(x['Zgain'], x['Zgain_err'])} | {g(x['Zgain_mc'])}"
            L.append(f"| {n} | {b} | {r['nb']} | {r['nb_mc']} | {v} | {f(x['epsB'], x['epsB_err'])} | {f(x['epsB_mc'], x['epsB_mc_err'])} | {f(ra, rae)} | {ex} |")

L.append("\n## 4. Rozkład piedestału sygnału D* (odjęcie pasm Δm), dane vs MC (pT)\n")
L.append("Średnie obcięte: perp, rnd ≤ 10 GeV; cone ≤ 30; near = cone − perp_avg ∈ [−10, 30].\n")
L.append("| klasa | bin | ⟨perp_avg⟩ dane | MC | dane/MC | ⟨near⟩ dane | MC | dane/MC | ⟨cone⟩ dane | MC | F_ped | f(perp<1) dane | MC | χ²/ndf rozkładu | ⟨ntrk05⟩ dane | MC |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
PED = {}
for n in CLS:
    for b in BINS:
        P = R["ped"][f"{n}|{b}"]
        if "perp" not in P:
            continue
        def rr(k):
            a, ea = P[k]["data"]; m, em = P[k]["mc"]
            return a, ea, m, em, a / m, np.hypot(ea / m, a * em / m ** 2)
        p, n_, c = rr("perp"), rr("near"), rr("cone")
        dc, edc = c[0] - c[2], np.hypot(c[1], c[3]); dp = p[0] - p[2]; edp = np.hypot(p[1], p[3])
        fped = f(dp / dc, abs(dp / dc) * np.hypot(edp / abs(dp) if dp else 0, edc / abs(dc))) if abs(dc) > 3 * edc else "– (mian. < 3σ)"
        f1 = P["f<1.0"]; h = P["hist"]; nt = P["ntrk05"]
        PED[f"{n}|{b}"] = {"perp_ratio": p[4], "perp_ratio_err": p[5], "near_ratio": n_[4], "near_ratio_err": n_[5], "d_cone": dc, "d_cone_err": edc, "d_perp": dp, "d_perp_err": edp,
                          "chi2": h["chi2"], "ndf": h["ndf"]}
        L.append(f"| {n} | {b} | {f(p[0], p[1], 2)} | {g(p[2], 2)} | {f(p[4], p[5], 2)} | {f(n_[0], n_[1], 2)} | {g(n_[2], 2)} | {f(n_[4], n_[5], 2)} | {f(c[0], c[1], 2)} | {g(c[2], 2)} | {fped} | "
                 f"{f(f1['data'][0], f1['data'][1])} | {g(f1['mc'][0])} | {g(h['chi2'], 1)}/{h['ndf']} | {f(nt['data'][0], nt['data'][1], 1)} | {g(nt['mc'][0], 1)} |")
L.append("\nRozkład perp_avg (ułamki w binach " + str(R["ped"]["9:MB|3-30"]["hist"]["edges"][:-1]) + "+ GeV), pT 3–30:\n")
L.append("| klasa | próbka | " + " | ".join(f"b{i}" for i in range(8)) + " |\n|---|---|" + "---|" * 8)
for n in CLS:
    h = R["ped"][f"{n}|3-30"]["hist"]
    L.append(f"| {n} | dane | " + " | ".join(f(a, e) for a, e in zip(h["data"], h["data_err"])) + " |")
    L.append(f"| {n} | MC | " + " | ".join(g(a) for a in h["mc"]) + " |")

L.append("\n## 5. Stałość SF po klasach rozłącznych (9:MB, J12only, MBTSbased, other)\n")
L.append("| bin | wariant | średnia | χ²/ndf | p |\n|---|---|---|---|---|")
for b in BINS:
    for v in VARS:
        c = R["const"][f"{b}|{v}"]
        L.append(f"| {b} | {v} | {f(c['mean'], c['err'])} | {g(c['chi2'], 1)}/{c['ndf']} | {g(c['p'], 4)} |")

# ---------------- werdykt wg KRYTERIA.md
V = {"variants": {}, "ped": PED}
for v in ("V1", "V2", "V3"):
    o = {"K1_cells": {}, "K2_cells": {}, "K3_cells": {}}
    for n in ("8:all", "9:MB"):
        for b in ("5-8", "8-30"):
            r = R["cells"][f"{n}|{b}"]; x = r["var"][v]
            ok = r["quality"] == "ok"
            z, zs = x["D"] / x["D_err"], x["D_sbs"] / x["D_err"]
            amb = (np.sign(x["D"]) != np.sign(x["D_sbs"])) or ((z >= 2) != (zs >= 2))
            cell = {"D": x["D"], "D_sbs": x["D_sbs"], "sigma": x["D_err"], "z": z, "z_sbs": zs, "quality": r["quality"], "ambiguous": bool(amb),
                    "SF_raw": r["var"]["raw"]["SF"], "SF_var": x["SF"], "SF_var_err": x["SF_err"]}
            if n == "8:all":
                cell["improve"] = bool(ok and not amb and z >= 2)
                o["K1_cells"][f"{n}|{b}"] = cell
            cell["not_worse"] = bool(z > -2 and zs > -2)
            o["K2_cells"][f"{n}|{b}"] = cell
            o["K3_cells"][f"{n}|{b}"] = {"rB": x["rB"], "rB_err": x["rB_err"], "pass": bool(x["rB"] <= 1.10)}
    k1 = any(c["improve"] for c in o["K1_cells"].values())
    k2 = all(c["not_worse"] for c in o["K2_cells"].values())
    k3 = all(c["pass"] for c in o["K3_cells"].values())
    closes = all(abs(c["SF_var"] - 1) < 2 * c["SF_var_err"] for c in o["K1_cells"].values())
    o.update({"K1": k1, "K2": k2, "K3": k3, "closes_8all": bool(closes),
              "verdict": "GORZEJ" if not k2 else ("BEZ POPRAWY" if not k1 else ("POPRAWA" if k3 else "POPRAWA SF, GORSZE TŁO"))})
    V["variants"][v] = o
V["R4"] = "potwierdzona" if any(o["verdict"] == "POPRAWA" for o in V["variants"].values()) else "niepotwierdzona"
V["R4_closes"] = any(o["verdict"] == "POPRAWA" and o["closes_8all"] for o in V["variants"].values())
json.dump(V, open(os.path.join(HERE, "verdict.json"), "w"), indent=1, ensure_ascii=False)
L.append("\n## 6. Werdykt wg kryteriów\n")
L.append("| wariant | K1 poprawa (8:all) | K2 nie gorzej | K3 tło ≤ 1,10 | domyka 8:all | werdykt |\n|---|---|---|---|---|---|")
for v, o in V["variants"].items():
    k1 = "; ".join(f"{k.split('|')[1]}: {g(c['z'],1)}σ (SBS {g(c['z_sbs'],1)}σ{', niejedn.' if c['ambiguous'] else ''}{'' if c['quality']=='ok' else ', niemiar.'})" for k, c in o["K1_cells"].items())
    k2 = "; ".join(f"{k}: {g(c['z'],1)}σ" for k, c in o["K2_cells"].items())
    k3 = "; ".join(f"{k}: {g(c['rB'],2)}" for k, c in o["K3_cells"].items())
    L.append(f"| {v} | {o['K1']} ({k1}) | {o['K2']} ({k2}) | {o['K3']} ({k3}) | {o['closes_8all']} | **{o['verdict']}** |")
L.append(f"\nR4 łącznie: **{V['R4']}**; domyka: {V['R4_closes']}")
open(os.path.join(HERE, "out", "tables.md"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
