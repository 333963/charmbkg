#!/usr/bin/env python3
"""XSEC_atlas B — sigma(D0), sigma(D+) z selekcja cieciowa (bez klasyfikatora) vs WP70 drzewa H, te same tabele, ten sam kod fitu (FITBIAS fb_fit.fit;
D+ = M0, D0 = M8). KRYTERIA.md §B; POST HOC wzgledem XSEC. Wspolny bootstrap po zdarzeniach (uid) -> blad stosunku z korelacja.
Uzycie: xa_fit.py [--nboot 200] [--only nom]"""
import argparse, glob, json, math, os, sys
import numpy as np
H = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(H, "skrypty")); sys.path.insert(0, os.path.join(os.path.dirname(H), "FITBIAS_2026-11-19", "skrypty"))
import xa_lib as X
import fb_fit as F

LUMI = 0.4648
BR = {"D0": 0.03936, "Dplus": 0.0931}
THR = {"D0": {"3-5": 0.822, "5-8": 0.8719, "8-30": 0.8631}, "Dplus": {"3-5": 0.9415, "5-8": 0.9517, "8-30": 0.9346}}
LXY_NOM = {"D0": 0.25, "Dplus": 0.5}
GC = {"D0": (0.0008, 0.10), "Dplus": (0.0005, 0.10)}
MODEL = {"D0": dict(shape="gauss", peak=True, rng="alt"), "Dplus": dict(shape="free", peak=False, rng="fit")}   # M8 / M0
NG = json.load(open(os.path.join(H, "out/truth_ngen.json")))["bins"]


def load(sample, sp):
    out = {}
    for f in sorted(glob.glob(os.path.join(H, "in/parts", f"xa_{sample}_*.npz"))):
        z = np.load(f, allow_pickle=False)
        for k in z.files:
            if k.startswith(sp + "__t_"):
                out.setdefault(k[len(sp) + 4:], []).append(z[k])
    return X.derived({k: np.concatenate(v) for k, v in out.items()}, sp)


def sel_wp(T, sp, bn):
    return T["ok"].astype(bool) & (T["Lxy"] >= LXY_NOM[sp]) & (T["score"] >= THR[sp][bn])


def sel_cut(T, sp, lxy=None, isig=4.0, ptmin=None):
    ok = T["ok"].astype(bool)
    if sp == "Dplus":
        return ok & (T["ptmin"] > (0.8 if ptmin is None else ptmin)) & (T["chi2"] < 8) & (T["Lxy"] > (1.1 if lxy is None else lxy)) & (np.abs(T["impact"]) < 1.0) & \
            (T["isig"] < isig) & (T["cts"] > -0.8) & (T["dphi"] > 0.008)
    return ok & (T["ptmin"] > (1.0 if ptmin is None else ptmin)) & (T["chi2"] < 10) & (T["Lxy"] > (0.45 if lxy is None else lxy)) & (np.abs(T["impact"]) < 1.0) & \
        (T["isig"] < isig) & (T["cts"] > -0.8)


VARS = {"D0": {"ciecia": dict(), "Lxy035": dict(lxy=0.35), "Lxy060": dict(lxy=0.60), "bez_s3D": dict(isig=1e9), "pt08": dict(ptmin=0.8)},
        "Dplus": {"ciecia": dict(), "Lxy080": dict(lxy=0.80), "Lxy150": dict(lxy=1.50), "bez_s3D": dict(isig=1e9)}}


def components(sp, FD, MB, sf, sm, bn):
    """Szablony po danej selekcji (maski sf na fd, sm na mb) w binie; -> dict."""
    c = F.SP[sp]; lo, hi = c["mr"]
    inr = lambda T: (T["m"] >= lo) & (T["m"] < hi)  # noqa: E731
    tsig = lambda T: np.isin(T["truth"], (0, 1))  # noqa: E731
    sw_f = FD["swap"] == 1 if sp == "D0" else np.zeros(len(FD["m"]), bool)
    sw_m = MB["swap"] == 1 if sp == "D0" else np.zeros(len(MB["m"]), bool)
    fs = sf & tsig(FD) & ~sw_f & inr(FD)
    o = {"sigT": F.Tmpl(FD["m"][fs], np.ones(fs.sum()), lo, hi), "n_fd": int(fs.sum()), "refT": None, "frefl": 0.0}
    if sp == "D0":
        fr = sf & tsig(FD) & sw_f & inr(FD)
        o["refT"], o["frefl"] = F.Tmpl(FD["m"][fr], np.ones(fr.sum()), lo, hi), float(fr.sum() / fs.sum())
    ms = sm & tsig(MB) & ~sw_m & inr(MB)
    mp = sm & np.isin(MB["truth"], (2, 3, 6, 7)) & inr(MB)
    ed = np.arange(lo, hi + F.BW / 2, F.BW)
    o["peakT"] = F.TmplE(ed, np.histogram(MB["m"][mp], ed)[0]); o["kap"] = float(mp.sum() / max(ms.sum(), 1))
    o["N_truth"] = int(ms.sum()); o["mask_truth"] = ms
    return o


def dofit(h, comp, sp, kscale=1.0, start=None, hesse=True):
    c, S = F.SP[sp], MODEL[sp]
    rng = c[S["rng"]]
    edges = np.arange(rng[0], rng[1] + F.BW / 2, F.BW)
    t = [(comp["refT"], comp["frefl"])] if comp["refT"] is not None else []
    if S["peak"]:
        t.append((comp["peakT"], comp["kap"] * kscale))
    f = F.fit(np.asarray(h, float)[None, :], [comp["sigT"]], [t], edges, c["m0"], c["nbkg"], shape=S["shape"], gc=GC[sp], win=c["win"], start=start, hesse=hesse)
    f.pop("mu", None)
    return f


def mh(m, sp, w=None):
    rng = F.SP[sp][MODEL[sp]["rng"]]
    return np.histogram(m, np.arange(rng[0], rng[1] + F.BW / 2, F.BW), weights=w)[0].astype(float)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--nboot", type=int, default=200); ap.add_argument("--sp", default="D0,Dplus")
    a = ap.parse_args()
    out = {}
    for sp in a.sp.split(","):
        D = load("data_427929", sp); D = {k: v[(D["cls"] & 3) == 3] for k, v in D.items()}
        MB, FD = load("mb_nd", sp), load("fd_nd", sp)
        uq, inv = np.unique(D["uid"], return_inverse=True)
        print(sp, "wiersze: dane", len(D["m"]), "mb", len(MB["m"]), "fd", len(FD["m"]), "zdarzen danych", len(uq), flush=True)
        for bn, lo, hi in X.BINS:
            inb = lambda T: (T["pt"] >= lo) & (T["pt"] < hi) & (np.abs(T["eta"]) < 2.2)  # noqa: E731
            sels = {"WP70": (lambda T: sel_wp(T, sp, bn))}
            for vn, kw in VARS[sp].items():
                sels[vn] = (lambda T, kw=kw: sel_cut(T, sp, **kw))
            ng = NG[f"{sp}|{bn}"]
            res, keep = {}, {}
            for sn, fs in sels.items():
                d, m, f = inb(D) & fs(D), inb(MB) & fs(MB), inb(FD) & fs(FD)
                comp = components(sp, FD, MB, f, m, bn)
                ft = dofit(mh(D["m"][d], sp), comp, sp)
                fm = dofit(mh(MB["m"][m], sp), comp, sp)                       # domkniecie na mb_nd (kappa = 1 x MC z definicji)
                C = comp["N_truth"] / ng["N_excl"]
                r = {"n_cand": int(d.sum()), "N": ft["N"][0], "N_err": (ft.get("err") or [float("nan")])[0], "valid": ft["valid"], "chi2": ft["chi2"], "ndf": ft["ndf"], "p": F.pval(ft["chi2"], ft["ndf"]),
                     "sc": ft["sc"], "dmu": ft["dmu"], "at_limit": ft["at_limit"], "bkg_win": ft["bkg_win"][0], "kap": comp["kap"], "frefl": comp["frefl"], "n_fd": comp["n_fd"],
                     "N_truth_mb": comp["N_truth"], "N_gen_excl": ng["N_excl"], "N_gen_dty": ng["N_dty"], "C": C, "C_err": math.sqrt(C * (1 - C) / ng["N_excl"]), "C_dty": comp["N_truth"] / ng["N_dty"],
                     "closure": fm["N"][0] / comp["N_truth"], "closure_err": (fm.get("err") or [float("nan")])[0] / comp["N_truth"], "closure_valid": fm["valid"]}
                r["sigma_ub"] = r["N"] / (C * BR[sp] * LUMI) / 1e3; r["stat_ub"] = r["sigma_ub"] * r["N_err"] / r["N"] if r["N"] else float("nan")
                if sp == "D0":
                    f2 = dofit(mh(D["m"][d], sp), comp, sp, kscale=1.364); r["N_kap1364"] = f2["N"][0]
                res[sn] = r; keep[sn] = (d, comp, ft["x"], m)
                print(f"{sp} {bn} {sn:8s} kand {r['n_cand']:7d} N = {r['N']:7.1f} ± {r['N_err']:5.1f} p = {r['p']:.3f} sc {r['sc']:.2f} valid {r['valid']} C = {100 * C:.2f} % (dty {100 * r['C_dty']:.2f}) κ {r['kap']:.3f}"
                      f" σ = {r['sigma_ub']:.1f} ± {r['stat_ub']:.1f} µb  domkn. {r['closure']:.3f} ± {r['closure_err']:.3f}" + (f"  N(κ×1,364) {r['N_kap1364']:.1f}" if sp == "D0" else ""), flush=True)
            # stosunki sigma_ciecia / sigma_WP70
            k = res["WP70"]; mk = keep["WP70"][3] & keep["WP70"][1]["mask_truth"]
            for sn in VARS[sp]:
                c = res[sn]; mc_ = keep[sn][3] & keep[sn][1]["mask_truth"]
                na, nb_, nc = int((mc_ & mk).sum()), int((mc_ & ~mk).sum()), int((mk & ~mc_).sum())
                Nc, Nk = na + nb_, na + nc
                c["R"] = c["sigma_ub"] / k["sigma_ub"]
                c["R_mcstat_rel"] = math.sqrt(na * (1 / Nc - 1 / Nk) ** 2 + nb_ / Nc ** 2 + nc / Nk ** 2)
                c["R_closure_corr"] = c["R"] * k["closure"] / c["closure"]
                c["overlap_mc"] = {"oba": na, "tylko_ciecia": nb_, "tylko_WP70": nc}
                dd, dk = keep[sn][0], keep["WP70"][0]
                c["overlap_data_cand"] = {"oba": int((dd & dk).sum()), "tylko_ciecia": int((dd & ~dk).sum()), "tylko_WP70": int((dk & ~dd).sum())}
            # bootstrap nominalu
            if a.nboot:
                Rb, Nb = [], {"ciecia": [], "WP70": []}
                for rep in range(a.nboot):
                    w = np.random.default_rng(770000 + rep).poisson(1.0, len(uq))[inv].astype(float)
                    n = {}
                    for sn in ("ciecia", "WP70"):
                        d, comp, x0, _ = keep[sn]
                        ft = dofit(mh(D["m"][d], sp, w[d]), comp, sp, start=x0, hesse=False)
                        n[sn] = ft["N"][0] if ft["valid"] else float("nan")
                        Nb[sn].append(n[sn])
                    Rb.append(n["ciecia"] / n["WP70"] if n["WP70"] > 0 else float("nan"))
                Rb = np.array(Rb); g = np.isfinite(Rb)
                c = res["ciecia"]
                c["boot"] = {"n_ok": int(g.sum()), "R_rel_sd": float(np.std(Rb[g], ddof=1) / np.mean(Rb[g])), "R_mean_over_nom": float(np.mean(Rb[g]) * res["WP70"]["N"] / c["N"]),
                             "N_c_rel_sd": float(np.nanstd(Nb["ciecia"], ddof=1) / np.nanmean(Nb["ciecia"])), "N_k_rel_sd": float(np.nanstd(Nb["WP70"], ddof=1) / np.nanmean(Nb["WP70"])),
                             "rho": float(np.corrcoef(np.array(Nb["ciecia"])[g], np.array(Nb["WP70"])[g])[0, 1])}
                b = c["boot"]
                print(f"   R = σ_c/σ_k = {c['R']:.3f} ± {c['R'] * b['R_rel_sd']:.3f} (boot, ρ = {b['rho']:.2f}; sd N_c {100 * b['N_c_rel_sd']:.1f} %, N_k {100 * b['N_k_rel_sd']:.1f} %) ± {c['R'] * c['R_mcstat_rel']:.3f} (stat MC)   po poprawce domkn. {c['R_closure_corr']:.3f}", flush=True)
            for sn in VARS[sp]:
                if sn != "ciecia":
                    print(f"   wariant {sn}: R = {res[sn]['R']:.3f} (stat MC ± {res[sn]['R'] * res[sn]['R_mcstat_rel']:.3f}), N = {res[sn]['N']:.0f} ± {res[sn]['N_err']:.0f}", flush=True)
            out[f"{sp}|{bn}"] = res
            json.dump(out, open(os.path.join(H, f"out/ciecia_{a.sp.replace(',', '_')}.json"), "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
