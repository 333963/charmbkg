#!/usr/bin/env python3
"""ZADANIE 2 — testy skoku częstości protonów przypadkowych poniżej progu ξ_cal dawców (bez obszaru sygnału).

(i)   strona przeciwna / „lustro”: w zdarzeniach z ξ_cal,s > 0,14 (własny proton po stronie s wykluczony) częstość tagu po stronie s
      w funkcji ξ_cal,o strony przeciwnej, względem przewidywania z warstw nominalnych (R = obs/pred); oraz estymator „lustrzany”:
      nominał z dodatkową klasą warstw — dawca z ξ_cal,o w tym samym przedziale co ξ_cal,s biorcy.
(ii)  protony „niemożliwe”: dokładnie jeden proton po stronie, ξ_p w wycinku [a, b], ξ_cal,s > b + δ + m; obs/pred w funkcji ξ_cal,s
      (0,06–0,14), marginesu m i ξ_cal,o; wycinki także poniżej okna tagu (0,02–0,03).
(iii) zależność od nSelBL / μ / LB oraz wariant drobnych warstw nSelBL.
Żadna wielkość nie używa zdarzeń z tagiem przechodzącym cięcie ξ_cal (obszar sygnału) poza N_acc z mieszania.
Wymaga template.json (profil obszaru użycia). Wynik: jump_tests.json.
"""
from __future__ import annotations
import json, math, os, sys, time
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "4")
from dataclasses import replace
import numpy as np
from common import HERE, SAMPLES, M, load, jsonable

XB_O = np.array([0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.2, 0.3, 0.5, np.inf])
MIR = np.array([0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14])           # klasy lustrzane (wewnętrzne brzegi)
SLICES = ((0.02, 0.03), (0.03, 0.04), (0.04, 0.05), (0.05, 0.06), (0.06, 0.08), (0.08, 0.10))
XB_S = np.array([0.06, 0.07, 0.08, 0.09, 0.10, 0.12, 0.14])
MARG = np.array([0.01, 0.025, 0.04, np.inf])                         # klasy marginesu m_eff = ξ_cal − b − δ
OPP = np.array([0, 0.06, 0.14, np.inf])
FINE_ACT = (4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 25, 28, 30, 35, 40, 45, 50, 60, 80)
NB = 500


def rr(k, p, ep=0.0):
    """R = k/p z błędem: Poisson obs (max(k, p) — bez zaniżania przy fluktuacji w dół) ⊕ błąd przewidywania."""
    if not p > 0:
        return [float(k), float(p), None, None]
    return [float(k), float(p), k / p, math.sqrt(max(k, p) + ep ** 2) / p]


def run(s, usage, log):
    d = load(s, with_cands=False)
    Dn, T, cfg = d["Dn"], d["T"], d["cfg"]
    boot = M.Boot.for_events(Dn, cfg, NB, 0)
    out = {"sample": s, "run": d["run"]}
    ms, m0 = M.mix_region(T, Dn, cfg, "AorC")
    n0, rep0 = m0.total(), m0.boot(boot)
    side_n0 = {sd: (ms[sd].total(), ms[sd].boot(boot)) for sd in "AC"}
    out["nominal"] = {"n_acc": n0, "err": float(rep0.std(ddof=1)), "A": side_n0["A"][0], "C": side_n0["C"][0]}
    log(f"{s}: nominał N_acc {n0:.1f} ± {rep0.std(ddof=1):.1f}")

    BASE = {}

    def ratio(parts, base=None):
        mv = M.MixSum(parts.values())
        nv, rv = mv.total(), mv.boot(boot)
        res = {"n_acc": nv, "ratio": nv / n0, "ratio_err_corr": float(np.std(rv / rep0 - 1, ddof=1)), "n_acc_err": float(rv.std(ddof=1)),
               "level_fractions": {int(k): float(np.mean(mv.level == k)) for k in np.unique(mv.level)}, "n_tagged_donors": int(mv.ntd)}
        for sd, p in parts.items():
            res[f"ratio_{sd}"] = p.total() / side_n0[sd][0]
            res[f"ratio_{sd}_err_corr"] = float(np.std(p.boot(boot) / side_n0[sd][1] - 1, ddof=1))
        if base is not None:
            nb_, rb_, pb_ = BASE[base]
            res[f"ratio_to_{base}"] = nv / nb_
            res[f"ratio_to_{base}_err_corr"] = float(np.std(rv / rb_ - 1, ddof=1))
            for sd, p in parts.items():
                res[f"ratio_{sd}_to_{base}"] = p.total() / pb_[sd][0]
                res[f"ratio_{sd}_to_{base}_err_corr"] = float(np.std(p.boot(boot) / pb_[sd][1] - 1, ddof=1))
        res["_keep"] = (nv, rv, {sd: (p.total(), p.boot(boot)) for sd, p in parts.items()})
        return res

    # ================================================================ (i) strona przeciwna
    hasD = Dn["hasD"].to_numpy().astype(bool)
    t1 = {"edges_opp": XB_O.tolist(), "sides": {}}
    acc = {}
    for sd in "AC":
        o = "C" if sd == "A" else "A"
        pool = M.donor_pool(Dn, cfg, sd)
        P = Dn.loc[pool].reset_index(drop=True)
        mp = M.mix(P, Dn, cfg, sd, pairs=False)
        tag, p = mp.own_tag[sd], mp.p_tag
        ok = np.isfinite(p)
        xo, hd, npo = P[f"xical_{o}"].to_numpy(), P["hasD"].to_numpy().astype(bool), P[f"np_{o}"].to_numpy()
        b = np.digitize(xo, XB_O) - 1
        nb = len(XB_O) - 1
        tb = {}
        for lab, m in (("hasD", hd), ("noD", ~hd), ("all", np.ones(len(P), bool)), ("hasD_npo0", hd & (npo == 0)), ("hasD_npo1", hd & (npo >= 1)),
                       ("all_npo0", npo == 0), ("all_npo1", npo >= 1)):
            mm = m & ok
            tb[lab] = {"n": np.bincount(b[mm], minlength=nb).tolist(), "k": np.bincount(b[mm], tag[mm], nb).tolist(), "pred": np.bincount(b[mm], p[mm], nb).tolist()}
        XS2 = np.array([cfg.veto_thr, 0.2, 0.3, 0.5, np.inf])
        xs_ = P[f"xical_{sd}"].to_numpy()
        b2 = np.digitize(xs_, XS2) - 1
        for lab, m in (("hasD", hd), ("all", np.ones(len(P), bool))):
            mm = m & ok
            sh2 = (len(XS2) - 1, nb)
            f2 = lambda w: np.bincount(b2[mm] * nb + b[mm], w, sh2[0] * nb).reshape(sh2).tolist()
            tb[lab + "_2d"] = {"edges_same": XS2.tolist(), "k": f2(tag[mm].astype(float)), "pred": f2(p[mm]), "n": f2(np.ones(mm.sum())),
                               "sumln_same": f2(np.log(xs_[mm])), "sumln_opp": f2(np.log(np.maximum(xo[mm], 1e-6)))}
        # (iii) zależność R przy małym ξ_cal,o od nSelBL, μ, LB
        low = ok & (xo < 0.12)
        dep = {}
        for var, edges in (("nsel", np.array([0, 10, 15, 20, 25, 30, 40, 50, 80, 10 ** 6])), ("mu", np.quantile(P["mu"].to_numpy(), [0, 1 / 3, 2 / 3, 1]) + [0, 0, 0, 1]),
                           ("lb", np.quantile(P["lb"].to_numpy(), [0, 1 / 3, 2 / 3, 1]) + [0, 0, 0, 1])):
            bb = np.clip(np.digitize(P[var].to_numpy(), edges) - 1, 0, len(edges) - 2)
            dep[var] = {"edges": edges.tolist()}
            for lab, m in (("hasD", hd), ("all", np.ones(len(P), bool))):
                for reg, mk in (("low", low), ("high", ok & ~low)):
                    mm = m & mk
                    dep[var][f"{lab}_{reg}"] = {"k": np.bincount(bb[mm], tag[mm], len(edges) - 1).tolist(), "pred": np.bincount(bb[mm], p[mm], len(edges) - 1).tolist()}
        t1["sides"][sd] = {"n_pool": int(len(P)), "n_pool_hasD": int(hd.sum()), "tables": tb, "dependence_low_xical_opp": dep,
                           "mean_nsel_pool_hasD": float(P["nsel"].to_numpy()[hd].mean()), "mean_nsel_pool_hasD_low": float(P["nsel"].to_numpy()[hd & (xo < 0.12)].mean())}
        del mp, P
    out["test_i"] = t1
    # kontrola: nominał z dawcami bez protonu po stronie przeciwnej (izoluje efekt samego wymogu np_o == 0)
    parts = {sd: M.mix(T, Dn, cfg, sd, logic="AorC", donor_sel=(Dn[f"np_{'C' if sd == 'A' else 'A'}"].to_numpy() == 0)) for sd in "AC"}
    out["nominal_npo0"] = ratio(parts)
    BASE["nominal_npo0"] = out["nominal_npo0"].pop("_keep")
    log(f"{s}: nominal_npo0: ×{out['nominal_npo0']['ratio']:.3f} ± {out['nominal_npo0']['ratio_err_corr']:.3f}")
    del parts
    # estymator lustrzany: klasa = przedział ξ_cal,s biorcy ↔ przedział ξ_cal,o dawcy
    cfgm = replace(cfg, cls_extra=cfg.cls_extra + ("mir",))
    for lab, sel_np in (("mirror_all", False), ("mirror_npo0", True)):
        parts = {}
        for sd in "AC":
            o = "C" if sd == "A" else "A"
            Ts = T.assign(mir=np.digitize(T[f"xical_{sd}"].to_numpy(), MIR).astype(np.int8))
            Ds = Dn.assign(mir=np.digitize(Dn[f"xical_{o}"].to_numpy(), MIR).astype(np.int8))
            parts[sd] = M.mix(Ts, Ds, cfgm, sd, logic="AorC", donor_sel=(Dn[f"np_{o}"].to_numpy() == 0) if sel_np else None)
        out[lab] = ratio(parts, "nominal_npo0")
        out[lab].pop("_keep")
        log(f"{s}: {lab} / nominal_npo0: ×{out[lab]['ratio_to_nominal_npo0']:.3f} ± {out[lab]['ratio_to_nominal_npo0_err_corr']:.3f}")
        log(f"{s}: {lab}: ×{out[lab]['ratio']:.3f} ± {out[lab]['ratio_err_corr']:.3f} (A {out[lab]['ratio_A']:.3f} ± {out[lab]['ratio_A_err_corr']:.3f}, C {out[lab]['ratio_C']:.3f} ± {out[lab]['ratio_C_err_corr']:.3f}); poziomy {out[lab]['level_fractions']}")
        del parts
    # (iii) drobne warstwy nSelBL
    _, mv = M.mix_region(T, Dn, replace(cfg, act_edges=FINE_ACT), "AorC")
    rv = mv.boot(boot)
    out["act_fine"] = {"act_edges": list(FINE_ACT), "n_acc": mv.total(), "ratio": mv.total() / n0, "ratio_err_corr": float(np.std(rv / rep0 - 1, ddof=1)),
                       "level_fractions": {int(k): float(np.mean(mv.level == k)) for k in np.unique(mv.level)}}
    log(f"{s}: act_fine: ×{out['act_fine']['ratio']:.3f} ± {out['act_fine']['ratio_err_corr']:.3f}; poziomy {out['act_fine']['level_fractions']}")
    del mv, ms, m0

    # ================================================================ (ii) protony „niemożliwe”
    shp = (len(SLICES), len(XB_S) - 1, len(MARG) - 1, len(OPP) - 1, 2)
    OBS, PRED = np.zeros(shp), np.zeros(shp)
    KP = np.zeros((len(SLICES), 2))            # tagi w puli (wycinek × hasD) — błąd przewidywania
    for sd in "AC":
        o = "C" if sd == "A" else "A"
        xs_all = Dn[f"xical_{sd}"].to_numpy()
        N = Dn.loc[(xs_all >= XB_S[0]) & (xs_all <= cfg.veto_thr)].reset_index(drop=True)
        xs, xo, hd, nps = N[f"xical_{sd}"].to_numpy(), N[f"xical_{o}"].to_numpy(), N["hasD"].to_numpy().astype(int), N[f"np_{sd}"].to_numpy()
        ix = np.digitize(xs, XB_S) - 1
        io = np.digitize(xo, OPP) - 1
        Dp = Dn.loc[M.donor_pool(Dn, cfg, sd)]
        for i, (a, b_) in enumerate(SLICES):
            mq = M.mix(N, Dn, cfg, sd, window=(a, b_), pairs=False)
            tg, p = mq.own_tag[sd], mq.p_tag
            me = xs - b_ - cfg.delta
            im = np.digitize(me, MARG) - 1
            ok = np.isfinite(p) & (im >= 0) & (ix >= 0) & (ix < shp[1]) & ((nps == 0) | tg)
            np.add.at(OBS[i], (ix[ok & tg], im[ok & tg], io[ok & tg], hd[ok & tg]), 1.0)
            np.add.at(PRED[i], (ix[ok], im[ok], io[ok], hd[ok]), p[ok])
            tgd = M.side_tags(Dp, cfg, (a, b_))[sd][0]
            KP[i] += np.bincount(Dp["hasD"].to_numpy().astype(int)[tgd], minlength=2)
            del mq
        log(f"{s}: (ii) strona {sd}: obs {OBS.sum():.0f}, pred {PRED.sum():.1f}")
        del N
    out["test_ii"] = {"slices": SLICES, "edges_xical_same": XB_S.tolist(), "edges_margin": MARG.tolist(), "edges_xical_opp": OPP.tolist(),
                      "axes": "[wycinek ξ_p, ξ_cal strona tagu, margines, ξ_cal strona przeciwna, hasD]", "obs": OBS.tolist(), "pred": PRED.tolist(), "k_pool": KP.tolist()}
    return out


def main():
    t00 = time.time()
    log = lambda m: print(f"[{time.time() - t00:6.0f} s] {m}", flush=True)
    tj = {x["sample"]: x for x in json.load(open(os.path.join(HERE, "template.json")))["samples"]}
    res = []
    for s in (sys.argv[1].split(",") if len(sys.argv) > 1 else SAMPLES):
        res.append(run(s, tj[s]["usage_profile"], log))
        json.dump(jsonable({"status": "WSTEPNE", "n_boot": NB, "results": res}), open(os.path.join(HERE, "jump_tests.json"), "w"), indent=1)
    log("zapisano jump_tests.json")


if __name__ == "__main__":
    main()
