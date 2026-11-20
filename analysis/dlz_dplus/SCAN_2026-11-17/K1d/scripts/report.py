"""K1d: tabele i werdykt wg KRYTERIA_K1d.md z out/k1d.json, out/weights.json → tables.md, verdict.json."""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
K1D = os.path.dirname(HERE)
SCAN = os.path.dirname(K1D)
R = json.load(open(os.path.join(K1D, "out", "k1d.json")))
WI = json.load(open(os.path.join(K1D, "out", "weights.json")))
NOM = json.load(open(os.path.join(SCAN, "K1c", "out", "k1c_nominal.json")))
BN = ["3-5", "5-8", "8-30"]
CL = ["H", "DLZ"]
WPS = [0.5, 0.7, 0.8]
ORD = ["2", "3", "4"]
VARS = R["variants"]
REF = {"a": "427929", "b1_427929": "427929", "b0_427929": "427929", "b1_428770": "428770", "b0_428770": "428770"}
C = R["cells"]


def bs(var, bn, c, wp, o):
    q = C[f"{var}|{bn}|{c}|{wp}"]["fit"][o]
    B = q["boot"]
    N, e, ok, nt = (np.array(B[k], float) for k in ("N", "err", "ok", "nt"))
    return q, N, e, ok > 0, nt


def stats(var, bn, wp, o):
    """b, σ_b per klasyfikator; Δ, σ_Δ; zysk; pull width."""
    out = {}
    qs = {c: bs(var, bn, c, wp, o) for c in CL}
    for c in CL:
        q, N, e, ok, nt = qs[c]
        br = N[ok] / nt[ok]
        out[c] = dict(b=q["b"], sb=float(br.std(ddof=1)), b_boot_mean=float(br.mean()), f_ok=float(ok.mean()), valid=bool(q["valid"] and not q["at_limit"]),
                      pull_w_boot=float(((N[ok] - q["N"]) / e[ok]).std(ddof=1)), pull_w_toy=q["toy"].get("pull_w"), pull_mean_toy=q["toy"].get("pull_mean"),
                      toy_rms=q["toy"].get("rms"), rel_true=q["rel_true"], chi2=q["chi2"] / q["ndf"], N_true=C[f"{var}|{bn}|{c}|{wp}"]["N_true"],
                      rms_boot=float((N[ok] / q["N"]).std(ddof=1)))
    okb = qs["H"][3] & qs["DLZ"][3]
    dH, dD = qs["H"][1][okb] / qs["H"][4][okb], qs["DLZ"][1][okb] / qs["DLZ"][4][okb]
    out["d"] = float(out["DLZ"]["b"] - out["H"]["b"])
    out["sd"] = float((dD - dH).std(ddof=1))
    gr = 1 - (qs["DLZ"][2][okb] / qs["DLZ"][4][okb]) / (qs["H"][2][okb] / qs["H"][4][okb])
    out["g"] = float(1 - out["DLZ"]["rel_true"] / out["H"]["rel_true"])
    out["sg"] = float(gr.std(ddof=1))
    # zysk w definicji danych (σ/N_fit)
    gf = 1 - (qs["DLZ"][2][okb] / qs["DLZ"][1][okb]) / (qs["H"][2][okb] / qs["H"][1][okb])
    out["g_fit"] = float(1 - (qs["DLZ"][0]["N_err"] / qs["DLZ"][0]["N"]) / (qs["H"][0]["N_err"] / qs["H"][0]["N"]))
    out["sg_fit"] = float(gf.std(ddof=1))
    return out


def sbvalid(var, bn):
    run = REF[var]
    fv = R["fitvar"][bn]
    o = {}
    for c in CL:
        cell = C[f"{var}|{bn}|{c}|0.7"]
        sb_mc = cell["S_win"] / cell["B_win"]
        sb_d = NOM["nominal"][f"{run}|{bn}|{c}|0.7"][fv]["fits"]["3"]["SB"]
        o[c] = dict(sb_mc=sb_mc, sb_data=sb_d, ratio=sb_mc / sb_d)
    o["valid"] = bool(all(abs(o[c]["ratio"] - 1) < 0.30 for c in CL))
    return o


L = []
V = {"criteria": "KRYTERIA_K1d.md", "variants": {}}
L.append("# K1d — tabele (MC minbias z pełnym tłem; fit jak K1c)\n")
L.append(f"Repliki bootstrapu po zdarzeniach MC: {R['nrep']}; pseudo-eksperymenty Poissona: {R['ntoy']}; kształt: " + ", ".join(f"{bn} {R['fitvar'][bn]}" for bn in BN) + ".\n")
L.append("## 1. Ważność: S/B prawdy MC (okno ±45 MeV) vs S/B z fitu danych (K1c, pol3), WP70\n")
L.append("| wariant | odniesienie | pT | H: MC / dane (stosunek) | H+DLZ: MC / dane (stosunek) | ważny |\n|---|---|---|---|---|---|")
VAL = {}
for var in VARS:
    for bn in BN:
        s = sbvalid(var, bn)
        VAL[(var, bn)] = s
        L.append(f"| {var} | {REF[var]} | {bn} | " + " | ".join(f"{s[c]['sb_mc']:.3f} / {s[c]['sb_data']:.3f} ({s[c]['ratio']:.2f})" for c in CL) + f" | {'tak' if s['valid'] else 'NIE'} |")
L.append("\nDla (a) także względem 428770: " + "; ".join(
    f"{bn} H {C[f'a|{bn}|H|0.7']['S_win']/C[f'a|{bn}|H|0.7']['B_win']/NOM['nominal'][f'428770|{bn}|H|0.7'][R['fitvar'][bn]]['fits']['3']['SB']:.2f}"
    f", H+DLZ {C[f'a|{bn}|DLZ|0.7']['S_win']/C[f'a|{bn}|DLZ|0.7']['B_win']/NOM['nominal'][f'428770|{bn}|DLZ|0.7'][R['fitvar'][bn]]['fits']['3']['SB']:.2f}" for bn in BN) + " (nieważne).\n")
L.append("## 2. Przeważenie (b): wagi i domknięcie po cięciu (dane vs MC ważone)\n")
L.append("| wariant | pT | k sygnału | ⟨w_B⟩ | w p99 / max | n_eff/n | pasma: MC/dane H70, DLZ70 | χ²/n pasm H70, DLZ70 | okno: MC/dane H70, DLZ70 (±stat MC) | okno H50/DLZ50 | okno H80/DLZ80 |\n|---|---|---|---|---|---|---|---|---|---|---|")
for k, w in WI.items():
    var, bn = k.split("|")
    cl = w["closure"]
    e = lambda c, wp: cl[f"{c}|{wp}"]  # noqa: E731
    L.append(f"| {var} | {bn} | {w['k_sig']:.2f} | {w['w_mean']:.2f} | {w['w_p99']:.1f} / {w['w_max']:.0f} | {w['neff_frac']:.2f} | {e('H',0.7)['sb_ratio']:.3f}, {e('DLZ',0.7)['sb_ratio']:.3f} | "
             f"{e('H',0.7)['chi2_sb']:.2f}, {e('DLZ',0.7)['chi2_sb']:.2f} | {e('H',0.7)['win_ratio']:.3f} ± {e('H',0.7)['win_mc_err']/e('H',0.7)['win_mc']:.3f}, {e('DLZ',0.7)['win_ratio']:.3f} ± {e('DLZ',0.7)['win_mc_err']/e('DLZ',0.7)['win_mc']:.3f} | "
             f"{e('H',0.5)['win_ratio']:.3f} / {e('DLZ',0.5)['win_ratio']:.3f} | {e('H',0.8)['win_ratio']:.3f} / {e('DLZ',0.8)['win_ratio']:.3f} |")
L.append("\n„okno” = wszystkie kandydaci w ±45 MeV (S + B); dla H70 S jest dopasowane z konstrukcji (k), reszta to przewidywanie.\n")

ST = {}
for wp in WPS:
    L.append(f"## 3. Obciążenie b = N_fit/N_true, WP{int(wp*100)} (± bootstrap po zdarzeniach MC)\n")
    L.append("| wariant | pT | rząd | N_true H / DLZ | b H | b H+DLZ | Δ = DLZ − H | Δ/σ | pull width H / DLZ (boot; toy) | χ²/ndf H / DLZ | f_ok |\n|---|---|---|---|---|---|---|---|---|---|---|")
    for var in VARS:
        for bn in BN:
            for o in ORD:
                s = stats(var, bn, wp, o)
                ST[(var, bn, wp, o)] = s
                H, D = s["H"], s["DLZ"]
                pt_ = lambda x: "—" if x is None else f"{x:.2f}"  # noqa: E731
                L.append(f"| {var} | {bn} | pol{o} | {H['N_true']:.0f} / {D['N_true']:.0f} | {H['b']:.3f} ± {H['sb']:.3f} | {D['b']:.3f} ± {D['sb']:.3f} | {s['d']:+.3f} ± {s['sd']:.3f} | {s['d']/s['sd']:+.1f} | "
                         f"{H['pull_w_boot']:.2f} / {D['pull_w_boot']:.2f}; {pt_(H['pull_w_toy'])} / {pt_(D['pull_w_toy'])} | {H['chi2']:.2f} / {D['chi2']:.2f} | {min(H['f_ok'], D['f_ok']):.2f} |")
    L.append("")

L.append("## 4. Skład tego, co fit bierze za sygnał (leave-one-out wg prawdy), WP70, pol3; wartości w % N_true\n")
L.append("Δ_g = N_fit(wszystko) − N_fit(bez grupy g). „sygnał − 1” = Δ_sig/N_true − 1. „samo tło” = N_fit na histogramie bez sygnału (kształt ustalony); „sama kombinatoryka” = to samo tylko dla kodów 8–15.\n")
grp = R["groups"]
L.append("| wariant | pT | klas. | b − 1 | sygnał − 1 | 2–7 razem | " + " | ".join(grp[1:]) + " | samo tło | sama komb. | n(2–7) w oknie / B w oknie |\n|---|---|---|---|---|---|" + "---|" * (len(grp) - 1) + "---|---|---|")
LOO = {}
for var in VARS:
    for bn in BN:
        for c in CL:
            cell = C[f"{var}|{bn}|{c}|0.7"]
            q = cell["fit"]["3"]
            nt = cell["N_true"]
            lo = q["loo"]
            LOO[f"{var}|{bn}|{c}"] = {g: lo[g]["dN"] / nt for g in lo}
            L.append(f"| {var} | {bn} | {c} | {100*(q['b']-1):+.1f} | {100*(lo['sig']['dN']/nt-1):+.1f} | {100*lo['pc_2-7']['dN']/nt:+.1f} | " + " | ".join(f"{100*lo[g]['dN']/nt:+.1f}" for g in grp[1:]) +
                     f" | {100*lo['bkg_only_fixshape']['dN']/nt:+.1f} | {100*lo['comb_only_fixshape']['dN']/nt:+.1f} | {lo['pc_2-7']['n_win']:.0f} / {cell['B_win']:.0f} |")
L.append("")
L.append("## 5. Zysk H+DLZ nad H na MC (pol3/pol3), g = 1 − [σ/N_true]_DLZ/[σ/N_true]_H; obok definicja danych (σ/N_fit) i zysk z danych K1c [%]\n")
L.append("| wariant | WP | 3–5 | 5–8 | 8–30 | średnio | def. σ/N_fit: średnio | dane K1c pol3/pol3 (run odniesienia): 3–5 / 5–8 / 8–30 / średnio |\n|---|---|---|---|---|---|---|---|")
GA = {}
for var in VARS:
    for wp in WPS:
        g = [ST[(var, bn, wp, "3")] for bn in BN]
        run = REF[var]
        gd = []
        for bn in BN:
            fv = R["fitvar"][bn]
            f = {c: NOM["nominal"][f"{run}|{bn}|{c}|{wp}"][fv]["fits"]["3"] for c in CL}
            gd.append(1 - (f["DLZ"]["N_err"] / f["DLZ"]["N"]) / (f["H"]["N_err"] / f["H"]["N"]))
        av, sav = np.mean([x["g"] for x in g]), np.sqrt(sum(x["sg"] ** 2 for x in g)) / 3
        avf, savf = np.mean([x["g_fit"] for x in g]), np.sqrt(sum(x["sg_fit"] ** 2 for x in g)) / 3
        GA[f"{var}|{wp}"] = dict(g=[x["g"] for x in g], sg=[x["sg"] for x in g], mean=float(av), mean_err=float(sav), mean_fitdef=float(avf), mean_fitdef_err=float(savf), data=gd, data_mean=float(np.mean(gd)))
        L.append(f"| {var} | {int(wp*100)} | " + " | ".join(f"{100*x['g']:+.1f} ± {100*x['sg']:.1f}" for x in g) + f" | {100*av:+.1f} ± {100*sav:.1f} | {100*avf:+.1f} ± {100*savf:.1f} | " +
                 " / ".join(f"{100*x:+.1f}" for x in gd) + f" / {100*np.mean(gd):+.1f} |")
L.append("\nσ = błąd hesjanu fitu histogramu o oczekiwanej statystyce (dla b: statystyka danych). Błąd ze sparowanego bootstrapu po zdarzeniach MC.\n")

# --- werdykt
L.append("## 6. Kryteria (WP70, pol3)\n")
L.append("| wariant | pT | ważny | C1 H (|b−1|−σ) | C1 H+DLZ | C2 |Δ|/σ_Δ | moc (σ_b H / DLZ / σ_Δ) | b1 vs b0 (H / DLZ) |\n|---|---|---|---|---|---|---|---|")
for var in ["a", "b1_427929", "b1_428770"]:
    rows, okall, power, valid_all, robust = [], True, True, True, True
    nvalid = 0
    for bn in BN:
        s = ST[(var, bn, 0.7, "3")]
        v = VAL[(var, bn)]["valid"]
        c1 = {c: bool(abs(s[c]["b"] - 1) - s[c]["sb"] < 0.05 and s[c]["valid"] and s[c]["f_ok"] >= 0.8) for c in CL}
        c2 = bool(abs(s["d"]) < 2 * s["sd"])
        pw = bool(s["H"]["sb"] <= 0.05 and s["DLZ"]["sb"] <= 0.05 and s["sd"] <= 0.05)
        rb = None
        if var.startswith("b1"):
            s0 = ST[(var.replace("b1", "b0"), bn, 0.7, "3")]
            rb = {c: float(s[c]["b"] - s0[c]["b"]) for c in CL}
            if v and any(abs(x) >= 0.05 for x in rb.values()):
                robust = False
        if v:
            nvalid += 1
            okall &= c1["H"] and c1["DLZ"] and c2
            power &= pw
        else:
            valid_all = False
        rows.append(dict(bin=bn, valid=v, C1=c1, C2=c2, power=pw, b={c: s[c]["b"] for c in CL}, sb={c: s[c]["sb"] for c in CL}, d=s["d"], sd=s["sd"], b1_minus_b0=rb))
        L.append(f"| {var} | {bn} | {'tak' if v else 'NIE'} | " + " | ".join(f"{'PASS' if c1[c] else 'FAIL'} ({abs(s[c]['b']-1)-s[c]['sb']:+.3f})" for c in CL) +
                 f" | {'PASS' if c2 else 'FAIL'} ({abs(s['d'])/s['sd']:.1f}) | {'tak' if pw else 'nie'} ({s['H']['sb']:.3f} / {s['DLZ']['sb']:.3f} / {s['sd']:.3f}) | " +
                 ("—" if rb is None else f"{rb['H']:+.3f} / {rb['DLZ']:+.3f}") + " |")
    if nvalid == 0 or not robust:
        ver = "NIEROZSTRZYGNIĘTY"
    elif not okall:
        ver = "POZA NORMĄ"
    elif power and valid_all:
        ver = "W NORMIE (wykazane)"
    else:
        ver = "ZGODNE, MOC OGRANICZONA"
    V["variants"][var] = dict(verdict=ver, rows=rows, robust=robust)
    L.append(f"| **{var}** | | | | | | | **{ver}** |")
okv = ("W NORMIE (wykazane)", "ZGODNE, MOC OGRANICZONA")
V["overall_bias_ok"] = bool(V["variants"]["a"]["verdict"] in okv and V["variants"]["b1_427929"]["verdict"] in okv)
V["extrapolation_428770"] = V["variants"]["b1_428770"]["verdict"]
V["gain"] = GA
V["loo_wp70_pol3"] = LOO
V["bias"] = {f"{k[0]}|{k[1]}|{k[2]}|pol{k[3]}": dict(bH=s["H"]["b"], sbH=s["H"]["sb"], bD=s["DLZ"]["b"], sbD=s["DLZ"]["sb"], d=s["d"], sd=s["sd"],
                                                    pullw_boot=[s["H"]["pull_w_boot"], s["DLZ"]["pull_w_boot"]], pullw_toy=[s["H"]["pull_w_toy"], s["DLZ"]["pull_w_toy"]]) for k, s in ST.items()}
V["validity"] = {f"{k[0]}|{k[1]}": v for k, v in VAL.items()}
L.append(f"\nWerdykt łączny (a ∧ b1_427929): {'obciążenie w normie / zgodne' if V['overall_bias_ok'] else 'NIE wykazano'}; b1_428770 (ekstrapolacja): {V['extrapolation_428770']}.\n")
open(os.path.join(K1D, "tables.md"), "w").write("\n".join(L) + "\n")
json.dump(V, open(os.path.join(K1D, "verdict.json"), "w"), indent=1)
print("\n".join(L))
