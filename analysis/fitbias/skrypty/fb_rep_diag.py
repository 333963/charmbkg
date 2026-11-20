#!/usr/bin/env python3
"""FITBIAS — tabele diagnozy (diag_T.json, diag_W.json) → tabele/diag_T.md, tabele/diag_W.md."""
import json, os, math
HERE = os.path.dirname(os.path.abspath(__file__))
f = lambda x, n=3: "—" if x is None else f"{x:+.{n}f}".replace(".", ",")
u = lambda x, n=3: "—" if x is None else f"{x:.{n}f}".replace(".", ",")
R = json.load(open(os.path.join(HERE, "..", "out", "diag_T.json")))["res"]
L = ["# FITBIAS — diagnoza T (MC_templates, mb_nd z kategoriami prawdy; szablony fd_nd). b = N_fit/N_prawda − 1.", ""]
L += ["## T1. Skład w oknie piku i S/B", "", "| gat. | pT | kat. | N_prawda | S/B okno | f_refl fd / mb | okno: Sig | Refl | DMM | CM | PartComb | Comb |", "|---|---|---|---|---|---|---|---|---|---|---|---|"]
for r in R:
    w = r["win"]
    L.append(f"| {r['sp']} | {r['bin']} | {r['cat']} | {r['N_true']:.0f} | {u(r['SB_win'])} | {u(r['frefl_fd'])} / {u(r['frefl_mb'])} | " + " | ".join(f"{w[g]:.0f}" for g in ("Sig", "Refl", "DMM", "CM", "PC", "Comb")) + " |")
L += ["", "## T2. Dokładanie składowych (kolejność I: Sig → +Refl → +DMM → +CM → +PartComb → +Comb), model nominalny; b po każdym kroku. Błąd b(Sig) po odjęciu Poissona sygnału; b(wszystko) ± σ_fit/N.", "",
      "| gat. | pT | kat. | kształt | Sig | +Refl | +DMM | +CM | +PC | +Comb (= wszystko) | σ_fit/N | sc (wszystko) | kolejność II: Sig+Comb | χ²/ndf |", "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
for r in R:
    for sh in ("fix", "free"):
        rows = {x["add"]: x for x in r["seq"][f"I|{sh}"]}
        r2 = {x["add"]: x for x in r["seq"][f"II|{sh}"]}
        a = rows["Comb"]
        L.append(f"| {r['sp']} | {r['bin']} | {r['cat']} | {sh} | {f(rows['Sig']['b'])} ± {u(rows['Sig'].get('b_err'))} | " + " | ".join(f(rows[g]["b"]) if g in rows else "—" for g in ("Refl", "DMM", "CM", "PC", "Comb"))
                 + f" | {u((a['err'] or 0) / r['N_true'])} | {u(a['sc'], 2)} | {f(r2['Comb']['b'])} | {a['chi2']:.0f}/{a['ndf']} |")
L += ["", "## T3. Sygnał pozorny składowych tła osobno (szablon sygnału stały, tło exp nominalne): N_spur/N_prawda", "",
      "| gat. | pT | kat. | DMM | CM | PartComb | Comb | DMM_C | DMM_B | PC_421 | PC_411 | PC_431 | PC_413 | całe tło |", "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
for r in R:
    s = r["spur"]
    g_ = lambda k: f"{f(s[k]['rel'])} ± {u(s[k]['rel_err'])}" if k in s else "—"
    L.append(f"| {r['sp']} | {r['bin']} | {r['cat']} | " + " | ".join(g_(k) for k in ("DMM", "CM", "PC", "Comb", "DMM_C", "DMM_B", "PC_421", "PC_411", "PC_431", "PC_413", "bkg_all")) + " |")
L += ["", "## T4. Gładkość tła (całe tło prawdy bez Sig/Refl): p(χ²) fitu samego tła i N_spur/N_prawda vs rząd exp(Czebyszew k); zakres nominalny i alternatywny", "",
      "| gat. | pT | kat. | zakres | k = nom−1 | nom | nom+1 | nom+2 |", "|---|---|---|---|---|---|---|---|"]
for r in R:
    nb = 2 if r["sp"] == "D0" else 3
    for rn in ("nom", "alt"):
        L.append(f"| {r['sp']} | {r['bin']} | {r['cat']} | {rn} | " + " | ".join(f"p = {r['smooth'][f'bkg_all|{rn}|exp{k}']['p_0']:.2g}; {f(r['smooth'][f'bkg_all|{rn}|exp{k}']['rel'], 2)} ± {u(r['smooth'][f'bkg_all|{rn}|exp{k}']['rel_err'], 2)}" for k in range(nb - 1, nb + 3)) + " |")
L += ["", "## T5. Asimov bez szumu i rząd tła na pełnych pseudo-danych (b; kształt stały / wolny)", "",
      "| gat. | pT | kat. | model sygnału + gładkie tło(nom+2): fix / free | sygnał prawdy + gładkie tło: fix / free | wszystko, nom: fix / free | wszystko, nom+1 | wszystko, nom+2 |", "|---|---|---|---|---|---|---|---|"]
for r in R:
    nb = 2 if r["sp"] == "D0" else 3
    a = r["asimov"]
    L.append(f"| {r['sp']} | {r['bin']} | {r['cat']} | {f(a[f'model_sig+smooth{nb+2}|fix']['b'])} / {f(a[f'model_sig+smooth{nb+2}|free']['b'])} | {f(a[f'true_sig+smooth{nb+2}|fix']['b'])} / {f(a[f'true_sig+smooth{nb+2}|free']['b'])} | "
             f"{f(a['all|fix']['b'])} / {f(a['all|free']['b'])} | {f(a[f'all|fix|exp{nb+1}']['b'])} / {f(a[f'all|free|exp{nb+1}']['b'])} | {f(a[f'all|fix|exp{nb+2}']['b'])} / {f(a[f'all|free|exp{nb+2}']['b'])} |")
L += ["", "## T6. Obciążenie własne estymatora — zabawki z modelu (prawda = model, N = N_prawda, statystyka mb_nd), 300 zabawek", "",
      "| gat. | pT | kat. | S/B | fix: obciążenie | pull śr. / rms | free: obciążenie (średnia) | mediana | pull śr. / rms | pokrycie 68 % | rms N/N: free / fix | sc śr. ± rms | na granicy |", "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
for r in R:
    a, b = r["toys"]["fix"], r["toys"]["free"]
    L.append(f"| {r['sp']} | {r['bin']} | {r['cat']} | {u(r['SB_win'])} | {f(a['bias'])} ± {u(a['bias_err'])} | {f(a['pull_mean'], 2)} / {u(a['pull_rms'], 2)} | {f(b['bias'])} ± {u(b['bias_err'])} | {f(b['bias_med'])} | {f(b['pull_mean'], 2)} / {u(b['pull_rms'], 2)} | "
             f"{u(b['cov68'], 2)} | {u(b['rms_rel'])} / {u(a['rms_rel'])} | {u(b['sc_mean'], 2)} ± {u(b['sc_rms'], 2)} | {u(b['at_limit'], 2)} |")
open(os.path.join(HERE, "..", "tabele", "diag_T.md"), "w").write("\n".join(L) + "\n")

R = json.load(open(os.path.join(HERE, "..", "out", "diag_W.json")))["res"]
L = ["# FITBIAS — diagnoza W (histogramy WPSCAN, drzewo H; mb_nd: sygnał / odbicie / reszta). b = N_fit/N_prawda − 1.", ""]
L += ["## W1. Fit jednej kategorii: dokładanie składowych i warianty odbicia / rzędu tła", "",
      "| gat. | pT | kat. | N_prawda | S/B | f_refl fd / mb | kształt | Sig | Sig+Refl (f_fd) | Sig+Refl (f_mb) | Sig+reszta | wszystko (f_fd) = NOMINAŁ ± σ/N | wszystko (f_mb) | wszystko (odbicie mb: kształt+norm.) | nom+1 | nom+2 | sc |", "|" + "---|" * 17]
for r in R:
    nb = 2 if r["sp"] == "D0" else 3
    for n, o in r["single"].items():
        for sh in ("fix", "free"):
            g = lambda k: f(o[k]["b"]) if k in o else "—"
            a = o[f"all|{sh}|frefl_fd"]
            L.append(f"| {r['sp']} | {r['bin']} | {n} | {o['N_true']:.0f} | {u(o['SB_win'])} | {u(o['frefl_fd'])} / {u(o['frefl_mb'])} | {sh} | {g(f'Sig|{sh}')} ± {u(o[f'Sig|{sh}'].get('b_err'))} | {g(f'Sig+Refl|{sh}|frefl_fd')} | {g(f'Sig+Refl|{sh}|frefl_mb')} | {g(f'Sig+rest|{sh}')} | "
                     f"{f(a['b'])} ± {u((a['err'] or 0) / o['N_true'])} | {g(f'all|{sh}|frefl_mb')} | {g(f'all|{sh}|refl_mb_shape')} | {g(f'all|{sh}|frefl_fd|exp{nb+1}')} | {g(f'all|{sh}|frefl_fd|exp{nb+2}')} | {u(a['sc'], 2)} |")
L += ["", "## W2. Reszta tła (bez sygnału i odbicia): p(χ²) samego tła i N_spur/N_prawda vs rząd", "", "| gat. | pT | kat. | nom−1 | nom | nom+1 | nom+2 |", "|---|---|---|---|---|---|---|"]
for r in R:
    nb = 2 if r["sp"] == "D0" else 3
    for n, o in r["single"].items():
        L.append(f"| {r['sp']} | {r['bin']} | {n} | " + " | ".join(f"p = {o[f'rest|exp{k}']['p_0']:.2g}; {f(o[f'rest|exp{k}']['rel'], 2)} ± {u(o[f'rest|exp{k}']['rel_err'], 2)}" for k in range(nb - 1, nb + 3)) + " |")
L += ["", "## W3. Fit łączny A + B (pass/fail): SF_rel = f_fit/f_MC(fd), prawda = f_prawda(mb)/f_MC(fd); obciążenia N_A, N_B", "",
      "| gat. | pT | ε | SF prawda | wariant | SF_rel ± σ | b(N_A) | b(N_B) | sc |", "|---|---|---|---|---|---|---|---|---|"]
for r in R:
    for e, o in r["joint"].items():
        for k, v in o.items():
            if not isinstance(v, dict):
                continue
            L.append(f"| {r['sp']} | {r['bin']} | 0,{e} | {u(v['SF_true'])} | {k} | {u(v['SF'])} ± {u(v.get('SF_err'))} | {f(v['bA'])} | {f(v['bB'])} | {u(v['sc'], 2)} |")
open(os.path.join(HERE, "..", "tabele", "diag_W.md"), "w").write("\n".join(L) + "\n")
