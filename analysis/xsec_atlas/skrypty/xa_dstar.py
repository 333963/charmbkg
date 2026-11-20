#!/usr/bin/env python3
"""XSEC_atlas A — sigma(D*+- -> D0 pi -> K pi pi) w 427929 z fitow Delta m (KRYTERIA.md §A; POST HOC wzgledem XSEC).
Wejscie: tabele par INCAL (D0__pair_*), prawda: out/truth_ngen.json (N_gen po prawdziwym trybie). Fit: charmtag.dstar_refit.fit_single."""
import json, os, sys
import numpy as np
H = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(H, "skrypty")); sys.path.insert(0, os.path.join(os.path.dirname(H), "AA_merge", "src"))
import xa_lib as X
from charmtag.dstar_refit import fit_single, shape_from_mc, FIT_SINGLE_ATTEMPTS, MPI

LUMI, LUMI_E = 0.4648, 0.0103           # nb^-1
BR, BR_E = 0.02665, 0.00028
BINS = [("3-5", 3, 5), ("3.5-5", 3.5, 5), ("5-8", 5, 8), ("8-30", 8, 30)]
NG = json.load(open(os.path.join(H, "out/truth_ngen.json")))
NEV = NG["nev"]
VAR = {"nom": dict(), "Lxy010": dict(lxy=0.10), "ptmin06": dict(ptmin=0.6), "okno55": dict(mw=0.055)}


def sel(P, lxy=0.25, mw=0.040, ptmin=1.0):
    return (P["Lxy"] >= lxy) & (np.abs(P["m"] - X.M0["D0"]) < mw) & (P["ptmin"] > ptmin) & (np.abs(P["eta"]) < 2.2)


def hist(dm, hi=156.0, step=0.25):
    e = np.arange(MPI + 0.9, hi + 1e-9, step)
    return np.histogram(dm * 1e3, e)[0].astype(float), e


def main():
    D = X.derived(X.load("data_427929", "D0__pair_"), "D0"); D = {k: v[(D["cls"] & 3) == 3] for k, v in D.items()}
    M = X.derived(X.load("mb_nd", "D0__pair_"), "D0"); F = X.derived(X.load("fd_nd", "D0__pair_"), "D0")
    out = {}
    for vn, kw in VAR.items():
        for bn, lo, hi in BINS:
            inb = lambda P: sel(P, **kw) & (P["ds_pt"] >= lo) & (P["ds_pt"] < hi)  # noqa: E731
            m, f, d = inb(M), inb(F), inb(D)
            sg = np.isin(M["ds_truth"], (0, 1)); pk = np.isin(M["ds_truth"], (2, 3, 6, 7)); sgf = np.isin(F["ds_truth"], (0, 1))
            rng = (M["ds_dm"] * 1e3 >= MPI + 0.9)
            n_sig = int((m & sg & rng).sum()); n_pk = int((m & pk & rng & (np.abs(M["ds_dm"] - 0.14543) < 0.003)).sum())
            ng = NG["bins"][f"Dstar|{bn}"]
            C = n_sig / ng["N_excl"]; C_e = np.sqrt(C * (1 - C) / ng["N_excl"]); C_dty = n_sig / ng["N_dty"]
            fpk = n_pk / n_sig
            shp = shape_from_mc(np.r_[M["ds_dm"][m & sg], F["ds_dm"][f & sgf]] * 1e3)
            shp_fd = shape_from_mc(F["ds_dm"][f & sgf] * 1e3)
            r = {"N_reco_truth": n_sig, "N_pk": n_pk, "f_pk": fpk, "C": C, "C_err": C_e, "C_dty": C_dty, "shape": shp, "N_gen_excl": ng["N_excl"], "f_b_gen": ng["f_b"],
                 "f_b_reco": float(np.isin(M["ds_truth"][m & sg & rng], (1,)).mean())}
            fits = {}
            hd, e = hist(D["ds_dm"][d]); fits["nom"] = fit_single(hd, e, shp)
            if vn == "nom":
                fits["ksztalt_fd"] = fit_single(hd, e, shp_fd)
                h2, e2 = hist(D["ds_dm"][d], hi=152.0); fits["zakres152"] = fit_single(h2, e2, shp)
                h3, e3 = hist(D["ds_dm"][d], step=0.5); fits["bin05"] = fit_single(h3, e3, shp)
                fits["tlo_expo"] = fit_single(hd, e, shp, attempts=FIT_SINGLE_ATTEMPTS[1:])
                r["hist"] = hd.tolist(); r["edges"] = e.tolist()
            hm, e = hist(M["ds_dm"][m]); fm = fit_single(hm, e, shp)          # domkniecie: caly mb_nd jako pseudodane
            r["closure"] = {"N_fit": fm["N"], "N_err": fm["N_err"], "N_true_sig": n_sig, "ratio_to_sig_pk": fm["N"] / (n_sig * (1 + fpk)),
                            "ratio_err": fm["N_err"] / (n_sig * (1 + fpk)), "valid": fm["valid"], "chi2": fm["chi2"], "ndf": fm["ndf"]}
            for k, ft in fits.items():
                ft.pop("attempts", None)
            N, Ne = fits["nom"]["N"], fits["nom"]["N_err"]
            Ns = N / (1 + fpk)
            sig = Ns / (C * BR * LUMI) / 1e3                                  # nb -> ub
            r.update({"fits": fits, "N": N, "N_err": Ne, "N_sig": Ns, "sigma_ub": sig, "stat_ub": sig * Ne / N, "n_cand": int(d.sum()),
                      "pythia_ub": ng["N_all"] / NEV * 55.587e3})             # sigma_ND(A3) = 55,587 mb (jak XSEC)
            out[f"{vn}|{bn}"] = r
            c = r["closure"]
            print(f"{vn:8s} {bn:6s} kand {r['n_cand']:6d} N = {N:7.1f} ± {Ne:5.1f} (s={fits['nom']['scale']:.2f}, chi2/ndf {fits['nom']['chi2']:.0f}/{fits['nom']['ndf']}, valid {fits['nom']['valid']}, {fits['nom']['bkg_model']})"
                  f"  C = {100 * C:.2f} ± {100 * C_e:.2f} % (C_dty {100 * C_dty:.2f})  f_pk = {fpk:.3f}  sigma = {sig:.1f} ± {sig * Ne / N:.1f} ub  Pythia {r['pythia_ub']:.1f}"
                  f"  domkn. {c['ratio_to_sig_pk']:.3f} ± {c['ratio_err']:.3f}" + ("  war: " + ", ".join(f"{k} {v['N']:.0f}±{v['N_err']:.0f}" for k, v in fits.items() if k != "nom") if vn == "nom" else ""), flush=True)
    json.dump(out, open(os.path.join(H, "out/dstar.json"), "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
