"""D23/ANALIZA: tabele z out/wp.json + out/fit.json -> tables.md, out/verdict.json.  Użycie: python s03_tables.py"""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TOP = os.path.dirname(HERE)
WP = json.load(open(os.path.join(TOP, "out", "wp.json")))
FT = json.load(open(os.path.join(TOP, "out", "fit.json")))
Y35 = {"427929": (4331.0, 257.0), "428770": (5917.0, 429.0)}
grid = WP["wp"]["grid"]
eps = {w["wp"]: w["eps"] for w in grid}
C = {(c["run"], c["wp"]): c for c in FT["cells"]}
L = []
V = {}


def bo(c, var):
    b = c["boot"].get(var)
    if b is None:
        return None
    ok = np.array(b["ok"], bool)
    return dict(N=np.array(b["N"]), dmu=np.array(b["dmu"]), sc=np.array(b["sc"]), ok=ok)


L.append("# D23/ANALIZA — tabele (generowane przez scripts/s03_tables.py)\n")
L.append(f"Repliki bootstrapu (po zdarzeniach, WP ≤ 50): {FT['nrep']}. Info: {json.dumps(FT['info'])}\n")
L.append("## T1. Siatka WP (MC fd_nd 2–3 GeV) i pasma boczne danych (okno ±60 MeV zaślepione); B i S w ±45 MeV\n")
for run in ("427929", "428770"):
    L.append(f"\n**{run}** — wybór: WP{int(WP['choice'][run]['wp']*100)}\n")
    L.append("| WP | próg score | ε_S^MC | B(pol3) | B pol2 / pol4 | χ²/ndf pasm | S_oczek | S/B oczek | Z_oczek |\n|---|---|---|---|---|---|---|---|---|")
    for r in WP["wp"][run]:
        L.append(f"| {int(r['wp']*100)} | {r['thr']:.3f} | {r['eps']:.3f} | {r['B']:.0f} | {r['B_pol2']:.0f} / {r['B_pol4']:.0f} | {r['chi2']:.0f}/{r['ndf']} | {r['S_exp']:.0f} | {r['SB_exp']:.4f} | {r['Z_exp']:.2f} |")

L.append("\n## T2. Fity (pełna statystyka, 5 MeV, 1,70–2,20). N_S: kształt z MC (fix), tło exp(pol3); σ_boot = rms bootstrapu po zdarzeniach; Z = N/max(σ_HESSE, σ_boot)\n")
L.append("| run | WP | ε | N_S ± HESSE (± boot) | Z | S/B (±45 MeV) | χ²/ndf | N pol2 / pol4 | N (kształt wolny) | y = N/ε | udane repl. |\n|---|---|---|---|---|---|---|---|---|---|---|")
for run in ("427929", "428770"):
    V[run] = {"wp": {}}
    for w in grid:
        c = C[(run, w["wp"])]
        f = c["fit"]
        r3 = f["3|fix"]
        b = bo(c, "fix")
        if b is not None:
            sb = float(b["N"][b["ok"]].std())
            okf = float(b["ok"].mean())
            sig = max(r3["N_err"], sb)
            sbs = f"{sb:.0f}"
        else:
            sb, okf, sig, sbs = None, None, r3["N_err"], "—"
        Z = r3["N"] / sig
        fr = f["3|free"]
        V[run]["wp"][str(w["wp"])] = dict(eps=w["eps"], thr=w["thr"], N=r3["N"], hesse=r3["N_err"], boot=sb, Z=Z, SB=r3["SB"], chi2=r3["chi2"], ndf=r3["ndf"], valid=r3["valid"],
                                         N_pol2=f["2|fix"]["N"], N_pol4=f["4|fix"]["N"], N_free=fr["N"], N_free_err=fr["N_err"], y=r3["N"] / w["eps"], y_err=sig / w["eps"], ok_frac=okf)
        L.append(f"| {run} | {int(w['wp']*100)} | {w['eps']:.3f} | {r3['N']:.0f} ± {r3['N_err']:.0f} (± {sbs}) | {Z:.1f} | {r3['SB']:.4f} | {r3['chi2']:.0f}/{r3['ndf']} | {f['2|fix']['N']:.0f} / {f['4|fix']['N']:.0f} | "
                 f"{fr['N']:.0f} ± {fr['N_err']:.0f}{' ⚠gr.' if fr['at_limit'] else ''} | {r3['N']/w['eps']:.0f} ± {sig/w['eps']:.0f} | {'—' if okf is None else f'{100*okf:.0f} %'} |")

L.append("\n## T3. Pozycja i szerokość (fit z wolnym dmu, sc; tło pol3; błąd = rms bootstrapu; ⊕ błąd MC w nawiasie osobno)\n")
L.append("| run | WP | μ_MC [MeV] | σ_eff MC [MeV] | Δμ = dane − MC [MeV] (± boot; MC ±) | σ_dane/σ_MC (± boot; MC ±) | pol2 / pol4: Δμ | pol2 / pol4: sc | granica | udane repl. |\n|---|---|---|---|---|---|---|---|---|---|")
for run in ("427929", "428770"):
    for w in grid:
        c = C[(run, w["wp"])]
        b = bo(c, "free")
        if b is None:
            continue
        f = c["fit"]
        fr = f["3|free"]
        ok = b["ok"]
        e_mu, e_sc = float(b["dmu"][ok].std() * 1e3), float(b["sc"][ok].std())
        mcb = FT["mcshape_boot"][str(w["wp"])]
        V[run]["wp"][str(w["wp"])].update(dmu=fr["dmu"] * 1e3, dmu_err=e_mu, dmu_mc_err=mcb["mu_err"] * 1e3, sc=fr["sc"], sc_err=e_sc, sc_mc_err=mcb["rms_err"] / w["shape"]["rms_eff"],
                                         at_limit=fr["at_limit"], ok_frac_free=float(ok.mean()), dmu_pol2=f["2|free"]["dmu"] * 1e3, dmu_pol4=f["4|free"]["dmu"] * 1e3,
                                         sc_pol2=f["2|free"]["sc"], sc_pol4=f["4|free"]["sc"], mu_mc=w["shape"]["mu"] * 1e3, rms_mc=w["shape"]["rms_eff"] * 1e3)
        L.append(f"| {run} | {int(w['wp']*100)} | {w['shape']['mu']*1e3:.2f} | {w['shape']['rms_eff']*1e3:.2f} | {fr['dmu']*1e3:+.2f} ± {e_mu:.2f} (MC ± {mcb['mu_err']*1e3:.2f}) | "
                 f"{fr['sc']:.3f} ± {e_sc:.3f} (MC ± {mcb['rms_err']/w['shape']['rms_eff']:.3f}) | {f['2|free']['dmu']*1e3:+.2f} / {f['4|free']['dmu']*1e3:+.2f} | {f['2|free']['sc']:.3f} / {f['4|free']['sc']:.3f} | "
                 f"{'TAK' if fr['at_limit'] else 'nie'} | {100*ok.mean():.0f} % |")

L.append("\n## T4. Skan: y = N/ε_MC (fix, pol3), różnice sparowane względem WP wybranego i ilorazy ρ(a/b) = y_a/y_b (bootstrap sparowany)\n")
L.append("| run | para a/b | y_a | y_b | y_a − y_b ± σ_par | odchylenie [σ] | ρ(a/b) ± σ | czułość 2σ [%] | ρ pol4 |\n|---|---|---|---|---|---|---|---|---|")
for run in ("427929", "428770"):
    ch = WP["choice"][run]["wp"]
    bw = [w["wp"] for w in grid if (run, w["wp"]) in C and C[(run, w["wp"])]["boot"]]
    V[run]["pairs"] = {}
    pairs = [(ch, x) for x in bw if x != ch] + [(0.3, 0.5), (0.4, 0.5)]
    for a_, b_ in pairs:
        ba, bb = bo(C[(run, a_)], "fix"), bo(C[(run, b_)], "fix")
        ok = ba["ok"] & bb["ok"]
        ya, yb = V[run]["wp"][str(a_)]["y"], V[run]["wp"][str(b_)]["y"]
        dd = ba["N"][ok] / eps[a_] - bb["N"][ok] / eps[b_]
        rr = (ba["N"][ok] / eps[a_]) / (bb["N"][ok] / eps[b_])
        rho = ya / yb
        rho4 = (C[(run, a_)]["fit"]["4|fix"]["N"] / eps[a_]) / (C[(run, b_)]["fit"]["4|fix"]["N"] / eps[b_])
        V[run]["pairs"][f"{a_}/{b_}"] = dict(diff=ya - yb, diff_err=float(dd.std()), pull=float((ya - yb) / dd.std()), rho=rho, rho_err=float(rr.std()), sens2=float(200 * rr.std() / rho), rho_pol4=rho4)
        L.append(f"| {run} | {int(a_*100)}/{int(b_*100)} | {ya:.0f} | {yb:.0f} | {ya-yb:+.0f} ± {dd.std():.0f} | {(ya-yb)/dd.std():+.2f} | {rho:.3f} ± {rr.std():.3f} | {200*rr.std()/rho:.0f} | {rho4:.3f} |")

L.append("\n## T5. Werdykt K1–K3 i porównanie z 3–5 GeV\n")
L.append("| run | WP | N_S | Z | K1 | Δμ [MeV] | σ_d/σ_MC | K3 (sąsiedzi) | y(2–3) | y(3–5, WP70) | y(2–3)/y(3–5) (oczek. 0,6) |\n|---|---|---|---|---|---|---|---|---|---|---|")
for run in ("427929", "428770"):
    ch = WP["choice"][run]["wp"]
    v = V[run]["wp"][str(ch)]
    k1 = bool(v["Z"] >= 5 and v["valid"] and (v["ok_frac"] or 0) >= 0.8)
    nb = [x for x in (round(ch + 0.1, 1), round(ch - 0.1, 1)) if f"{ch}/{x}" in V[run]["pairs"]]
    pulls = {x: V[run]["pairs"][f"{ch}/{x}"]["pull"] for x in nb}
    k3 = all(abs(p) <= 2 for p in pulls.values())
    y35, e35 = Y35[run]
    rat = v["y"] / y35
    rat_e = rat * np.hypot(v["y_err"] / v["y"], e35 / y35)
    k2 = k1 and not v["at_limit"]
    V[run]["verdict"] = dict(wp=ch, K1=k1, K2_reportable=k2, K3=k3, pulls={str(k): p for k, p in pulls.items()}, y_ratio_35=rat, y_ratio_35_err=float(rat_e))
    pos = f"{v['dmu']:+.2f} ± {np.hypot(v['dmu_err'], v['dmu_mc_err']):.2f}" if k2 else ("nieokreślone (granica)" if k1 else "—")
    wid = f"{v['sc']:.3f} ± {np.hypot(v['sc_err'], v['sc_mc_err']):.3f}" if k2 else ("nieokreślone" if k1 else "—")
    L.append(f"| {run} | {int(ch*100)} | {v['N']:.0f} ± {v['N']/v['Z']:.0f} | {v['Z']:.1f} | {'TAK' if k1 else 'NIE'} | {pos} | {wid} | {'stabilny' if k3 else 'NIESTABILNY'} "
             f"({', '.join(f'WP{int(k*100)}: {p:+.2f}σ' for k, p in pulls.items())}) | {v['y']:.0f} ± {v['y_err']:.0f} | {y35:.0f} ± {e35:.0f} | {rat:.2f} ± {rat_e:.2f} |")

open(os.path.join(TOP, "tables.md"), "w").write("\n".join(L) + "\n")
json.dump(V, open(os.path.join(TOP, "out", "verdict.json"), "w"), indent=1)
print("\n".join(L))
