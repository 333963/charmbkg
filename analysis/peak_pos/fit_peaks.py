"""PEAK_pos: pozycja i szerokość piku D0 / D+ osobno dla D i D̄, dane (30 % lumi 428770) vs MC (prawda).
Wynik: out/peaks.json, out/peaks_rows.csv. Użycie: python fit_peaks.py"""
import json
import numpy as np, pandas as pd
from common import load, M0, PTB, OUT
from fitlib import fit_mc, fit_data

WIN = {"D0": {"nom": (1.76, 1.98), "narrow": (1.79, 1.94)}, "Dplus": {"nom": (1.74, 2.04), "narrow": (1.78, 1.96)}}
ETAB = [(-2.5, -1.0), (-1.0, 0.0), (0.0, 1.0), (1.0, 2.5)]
rows = []


def bins():
    for lo, hi in PTB:
        yield f"pT{lo}-{hi}", lambda L, lo=lo, hi=hi: (L.pt >= lo) & (L.pt < hi)
    yield "pT3-30", lambda L: np.ones(len(L), bool)
    for lo, hi in ETAB:
        yield f"eta{lo:+.0f}..{hi:+.0f}", lambda L, lo=lo, hi=hi: (L.eta >= lo) & (L.eta < hi)
    yield "abseta<1", lambda L: np.abs(L.eta) < 1
    yield "abseta>1", lambda L: np.abs(L.eta) >= 1


def run(sp, wp=0.7, veto=True, variants=True):
    D = load(sp, "data", wp=wp, veto=veto)
    S = load(sp, "mc", wp=wp, veto=veto)
    R = None
    if sp == "D0":
        R = S[S.hypo_swap == 1]          # refleksje (złe przypisanie prawdziwego D0)
        S = S[S.hypo_swap == 0]
    lo, hi = WIN[sp]["nom"]
    for bname, sel in bins():
        sb, db = np.asarray(sel(S)), np.asarray(sel(D))
        shape = fit_mc(S.m.to_numpy()[sb], lo, hi)          # kształt z MC (oba ładunki)
        for q in (0, 1, -1):
            qs = np.ones(len(S), bool) if q == 0 else (S.q == q).to_numpy()
            qd = np.ones(len(D), bool) if q == 0 else (D.q == q).to_numpy()
            mc = fit_mc(S.m.to_numpy()[sb & qs], lo, hi)
            base = dict(sp=sp, bin=bname, q=q, wp=wp, veto=veto, mu_mc=mc["mu"], e_mu_mc=mc["e_mu"], seff_mc=mc["seff"],
                        n_mc=int((sb & qs).sum()), shape_seff=shape["seff"])
            md = D.m.to_numpy()[db & qd]
            vlist = [("nom", "cheb2", "nom", None)]
            if variants and q != 0 or variants and bname.startswith("pT"):
                vlist += [("cheb3", "cheb3", "nom", None), ("exp1", "exp1", "nom", None), ("narrow", "cheb2", "narrow", None),
                          ("bw1", "cheb2", "nom", "bw1")]
                if sp == "D0" and veto is False:
                    pass
            if sp == "D0" and not veto:   # keep_both: szablon refleksji z MC, R = N_refl/N_sig w oknie
                rq = np.ones(len(R), bool) if q == 0 else (R.q == q).to_numpy()
                rb = np.asarray(sel(R)) & rq
                rm = R.m.to_numpy()[rb]
                ratio = ((rm >= lo) & (rm < hi)).sum() / max(((S.m.to_numpy()[sb & qs] >= lo) & (S.m.to_numpy()[sb & qs] < hi)).sum(), 1)
                # refleksja z q=−q' sygnału: przy podziale na q bierzemy refleksje oznaczone tym samym q (co widać w danych)
                vlist = [(v[0], v[1], v[2], v[3]) for v in vlist]
            for vname, bkg, wn, opt in vlist:
                wlo, whi = WIN[sp][wn]
                kw = {}
                if sp == "D0" and not veto:
                    kw["refl"] = (rm, ratio)
                r = fit_data(md, wlo, whi, shape, bkg=bkg, bw=0.001 if opt == "bw1" else 0.002, mu0=M0[sp], **kw)
                rows.append(dict(base, var=vname, mu=r["mu"], e_mu=r["e_mu"], k=r["k"], e_k=r["e_k"], seff=r["seff"],
                                 e_seff=r["e_seff"], Ns=r["Ns"], e_Ns=r["e_Ns"], ok=r["ok"], n=r["n"],
                                 refl_ratio=(ratio if sp == "D0" and not veto else np.nan)))
        print(sp, wp, veto, bname, "done", flush=True)


if __name__ == "__main__":
    run("D0"); run("Dplus")
    run("D0", veto=False, variants=False)
    run("D0", wp=0.5, variants=False); run("Dplus", wp=0.5, variants=False)
    T = pd.DataFrame(rows)
    T.to_csv(f"{OUT}/peaks_rows.csv", index=False)
    print(T[(T["var"] == "nom") & (T.wp == 0.7)].to_string(float_format=lambda x: f"{x:.5g}"))
