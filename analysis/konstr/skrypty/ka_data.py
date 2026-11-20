#!/usr/bin/env python3
"""KONSTR A — dane: (1) pomiar A_R z bliźniaków poza skrzyżowaniem (bez fitu okna), (2) fity M0 / M8 / M8-R / M0-R (odbicie: kształt fd, normalizacja z danych przez f = A_R/N_S, iteracyjnie).
ka_data.py [RUN ...]  → out/ka_data.json, tabele/A_dane.md"""
import json, os, sys
import numpy as np
H = os.path.dirname(os.path.abspath(__file__)); T = os.path.dirname(H); W = os.path.dirname(T)
sys.path.insert(0, H); sys.path.insert(0, os.path.join(W, "FITBIAS_2026-11-19", "skrypty"))
import ka_lib as K
import fb_fit as F
import fb_var as V
V._init()
FD = np.load(f"{T}/in/kx_fd_nd.npz")
NREP = 60
out, L1, L2 = {}, [], []
runs = sys.argv[1:] or ["428770"]
for run in runs:
    D = np.load(f"{T}/in/kx_data_{run}.npz")
    for bn in ("5-8", "8-30"):
        for e in (0.7, 0.5):
            fd1 = K.cum(FD[f"G|fd_nd|{bn}|t1"], e)
            G2 = K.cum(D[f"G|data_{run}|{bn}|all"], e)
            a = K.measure2(G2, fd1)
            infl = max(1.0, np.sqrt(a["chi2"] / a["ndf"]))
            cat = V.wcat("D0", bn, e)
            H5 = V.wdata(cat, V.SRC[run])
            H5 = H5[0] if H5.ndim == 2 else H5
            chk = float(G2.sum() / H5.sum())
            res = {"A": a["A"], "eA": a["eA"] * infl, "chi2_shape": a["chi2"], "ndf_shape": a["ndf"], "frefl_mc": cat["frefl"], "frac_in_2D": chk, "y": a["y"].tolist(), "ey": a["ey"].tolist(), "f": (a["A"] * a["f"]).tolist()}
            rng = np.random.default_rng(5 + int(e * 100) + (0 if bn == "5-8" else 1000))
            reps = rng.poisson(H5, size=(NREP, len(H5))).astype(float)

            # UWAGA (oznaczone 20.11.2026): kscale=1.364 to skala κ WYCOFANA (docs/results/cver/RAPORT.md pkt 5; nominał = 1.0).
            # Wartość zostaje bez zmian, bo z nią policzono zapisane wyniki KONSTR A; nie przenosić do nowych rachunków.
            def run_fit(h, mname, frefl=None, kscale=1.364):
                c = dict(cat); c["kap"] = cat["kap"] * kscale
                if frefl is not None:
                    c["frefl"] = frefl
                return V.fit_set(h[None], [c], "D0", mname)

            def fitR(h, mname, A):
                f = cat["frefl"]
                for _ in range(6):
                    r = run_fit(h, mname, frefl=f)
                    fn = A / max(r["N"][0], 1.0)
                    if abs(fn - f) < 1e-4:
                        break
                    f = fn
                r["frefl_used"] = f
                return r
            for mname, isR in (("M0", False), ("M8", False), ("M8", True), ("M0", True)):
                r = fitR(H5, mname, a["A"]) if isR else run_fit(H5, mname)
                bs = []
                if bn == "8-30" and e == 0.7:                      # bootstrap tylko w komórce docelowej (czas); w -R: f = A_b/N z jednej iteracji od punktu zbieżności
                    for h in reps:
                        if isR:
                            Ab = a["A"] + rng.normal() * a["eA"] * infl
                            rb = run_fit(h, mname, frefl=Ab / r["N"][0]); rb = run_fit(h, mname, frefl=Ab / max(rb["N"][0], 1.0))
                        else:
                            rb = run_fit(h, mname)
                        if rb["valid"]:
                            bs.append(rb["N"][0])
                bs = bs or [float("nan")]
                print(run, bn, e, mname, isR, round(r["N"][0]), r["chi2"], flush=True)
                nm = mname + ("-R" if isR else "")
                res[nm] = {"N": r["N"][0], "hesse": (r.get("err") or [None])[0], "boot": float(np.std(bs)), "nboot": len(bs), "chi2": r["chi2"], "ndf": r["ndf"], "p": F.pval(r["chi2"], r["ndf"]), "sc": r["sc"], "dmu": r["dmu"],
                           "valid": r["valid"], "at_limit": r["at_limit"], "frefl": r.get("frefl_used", cat["frefl"])}
            out[f"{run}|{bn}|{e}"] = res
            n8 = res["M8"]["N"]
            L1.append(f"| {run} | {bn} | {int(e * 100)} | {a['A']:.0f} ± {a['eA'] * infl:.0f} | {a['chi2']:.1f}/{a['ndf']} | {cat['frefl']:.3f} | {a['A'] / res['M0']['N']:.3f} ± {a['eA'] * infl / res['M0']['N']:.3f} | {a['A'] / n8:.3f} ± {a['eA'] * infl / n8:.3f} | {a['A'] / n8 / cat['frefl']:.2f} ± {a['eA'] * infl / n8 / cat['frefl']:.2f} |")
            for nm in ("M0", "M8", "M8-R", "M0-R"):
                q = res[nm]
                he = q["hesse"] if q["hesse"] is not None else float("nan")
                L2.append(f"| {run} | {bn} | {int(e * 100)} | {nm} | {q['N']:.0f} | {he:.0f} | {q['boot']:.0f} | {q['chi2']:.1f}/{q['ndf']} | {q['p']:.2g} | {q['frefl']:.3f} | {q['sc']:.3f} | {'' if q['valid'] else 'nieważny '}{'granica' if q['at_limit'] else ''} |")
tag = "_".join(runs)
json.dump(out, open(f"{T}/out/ka_data_{tag}.json", "w"), indent=1)
txt = ("# A — dane " + tag + "\n\n## Pomiar odbicia z bliźniaków (E2, poza skrzyżowaniem; błąd × √(χ²/ndf) gdy > 1)\n\n| run | bin | WP | A_R | χ² kształtu fd/ndf | f_refl MC | A_R/N_S(M0) | A_R/N_S(M8) | (dane/MC, M8) |\n|---|---|---|---|---|---|---|---|---|\n" + "\n".join(L1) +
       "\n\n## Fity (σ: HESSE i bootstrap Poissona 60, tylko 8–30 WP70; w -R bootstrap zawiera błąd A_R)\n\n| run | bin | WP | model | N | HESSE | boot | χ²/ndf | p | f_refl | sc | uwagi |\n|---|---|---|---|---|---|---|---|---|---|---|---|\n" + "\n".join(L2) + "\n")
open(f"{T}/tabele/A_dane_{tag}.md", "w").write(txt)
print(txt)
