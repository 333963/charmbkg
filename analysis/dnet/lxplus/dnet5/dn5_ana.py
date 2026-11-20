#!/usr/bin/env python3
"""DNET etap 5a — porównania SPAROWANE między wariantami sieci i zysk kaskady H × E (kryteria: DNET/KRYTERIA_etap5.md). Tylko pojedyncze D,
bez tagu AFP, bez okna sygnału danych: sygnał = MC (fd_nd) w oknie, tło = MC truth ≥ 8 w oknie oraz dalekie pasma boczne danych obu runów.

  dn5_ana.py --data DIR1/data --sc DN_OUT,DN2_OUT,DN5_OUT [--sc4 DN4_OUT] --out dnet5_cmp.json [--nboot 200]

(1) pary (A, B): r = Z_A/Z_B przy WP70 (i WP80/50 opisowo), Z = ε_S/√ε_B, progi = kwantyle score sygnału MC per bin pT, wiersze wspólne,
    bootstrap Poissona po zdarzeniach — te same repliki dla A i B (błąd sparowany). „B = średnia ziaren” → średnia r po ziarnach odniesienia.
(2) kaskada: wiersze po cięciu drzewem H (odniesienie etapu 1) na WP70; G_E(wp) = ε_S(E|H)/√ε_B(E|H) = Z(H×E)/Z(H); punkt pracy E wybierany
    na tle MC (siatka 0,9…0,5), raportowany na pasmach danych; kontrola „samo zaostrzenie H do tej samej wydajności”; rzeźbienie E po H."""
import argparse
import datetime
import glob
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dn5_lib as L5  # noqa: E402
import dn_feat as DF  # noqa: E402
import dn_train as DT  # noqa: E402
import analyze as AN  # noqa: E402

WPS = [0.80, 0.70, 0.50]
EWPS = [0.90, 0.80, 0.70, 0.60, 0.50]
BN = [AN.bname(b) for b in AN.BINS]
BASE = [f"D0.tr.H+e04.mc.d{L5.LAM}.s{s}" for s in (1, 2, 3)]
LC50, LC25, PAIR = L5.BLOCKS["lc50"]["tag"], L5.BLOCKS["lc25"]["tag"], L5.BLOCKS["pair"]["tag"]
PAIRS = {  # nazwa: (A, [B…])
    "lc_100over50": (BASE, [LC50]), "lc_50over25": ([LC50], [LC25]), "lc_100over25": (BASE, [LC25]),
    "pair_over_base": ([PAIR], BASE), "seed2_over_seed1": ([BASE[1]], [BASE[0]]), "seed3_over_seed1": ([BASE[2]], [BASE[0]]),
}
ETAGS = {"D0": [L5.BLOCKS["eD0"]["tag"], L5.BLOCKS["eD0h"]["tag"]], "Dplus": [L5.BLOCKS["eDp"]["tag"], L5.BLOCKS["eDph"]["tag"]]}


def zrep(sc, th, bw):
    """ε = udział score ≥ th; → (ε, repliki ε)."""
    p = (sc >= th).astype(np.float32)
    return float(p.mean()) if len(p) else float("nan"), (bw @ p) / np.maximum(bw.sum(1), 1)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True); ap.add_argument("--sc", required=True); ap.add_argument("--sc4", default=""); ap.add_argument("--out", required=True)
    ap.add_argument("--nboot", type=int, default=200)
    a = ap.parse_args(argv)
    scf = {}
    for d in a.sc.split(","):
        for f in sorted(glob.glob(os.path.join(d, "sc", "*.npz"))):
            scf.setdefault(os.path.basename(f)[:-4], f)
    if a.sc4:
        for t in BASE:
            f = os.path.join(a.sc4, "sc", t + ".npz")
            if os.path.isfile(f) and os.path.isfile(os.path.join(a.sc4, "meta", t + ".json")):
                scf[t + "@4"] = f
                PAIRS[f"etap4_over_etap2_s{t[-1]}"] = ([t + "@4"], [t])
    res = {"generated": datetime.datetime.now().isoformat(timespec="seconds"), "nboot": a.nboot, "sources": {}, "pairs": {}, "cascade": {}, "missing": []}
    for sp in DF.SPECIES:
        scfg = DF.cfg["species"][sp]
        m0, win, (lo, hi), sbm, far = scfg["m0"], scfg["signal_window"], scfg["mass_range"], scfg["sideband_min"], DF.FAR[sp]
        S = DT.load_species(a.data, sp, {"m", "pt", "uid", "truth", "ok", "swap", "m_swap"})
        sig, bmc, _ = DT.masks(S, sp)
        fd, mb = S["fd_nd"], S["mb_nd"]
        sets = {"sig": ("fd_nd", sig & (np.abs(fd["a_m"] - m0) < win)), "mc": ("mb_nd", bmc & (np.abs(mb["a_m"] - m0) < win)), "bmc_all": ("mb_nd", bmc)}
        for run in (428770, 427929):
            D = S.get(f"data_{run}")
            if D is not None:
                okd = D["a_ok"] & ((np.abs(D["a_m_swap"] - m0) > DF.REFL_VETO) if sp == "D0" else True)
                sets[f"sb{run}"] = (f"data_{run}", okd & (np.abs(D["a_m"] - m0) > far))        # wyłącznie dalekie pasma boczne
        idx = {k: np.flatnonzero(mk) for k, (s, mk) in sets.items()}
        ib = {k: np.digitize(S[sets[k][0]]["a_pt"][idx[k]], [3, 5, 8, 30]) - 1 for k in sets}
        mm = {k: S[sets[k][0]]["a_m"][idx[k]].astype(np.float64) for k in sets}
        BW = {k: AN.boot_w(S[sets[k][0]]["a_uid"][idx[k]], a.nboot, 11 + i) for i, k in enumerate(sets) if k != "bmc_all"}
        bk = [k for k in sets if k not in ("sig", "bmc_all")]
        cache = {}

        def sco(t, k):
            if t not in cache:
                cache[t] = dict(np.load(scf[t]))
            key = sets[k][0] + "__a"
            return cache[t][key][idx[k]] if key in cache[t] else np.full(len(idx[k]), np.nan, np.float32)

        def zz(t, valid, bi, wp):
            """Z i repliki dla tagu t na wierszach valid, bin bi, punkt pracy wp → {tło: (Z, repliki)}."""
            v = {k: valid[k] & (ib[k] == bi) for k in ["sig"] + bk}
            ss = sco(t, "sig")[v["sig"]]
            if len(ss) < 50:
                return {}
            th = float(np.quantile(ss, 1 - wp))
            es, esr = zrep(ss, th, BW["sig"][:, v["sig"]])
            o = {}
            for k in bk:
                if v[k].sum() < 20:
                    continue
                eb, ebr = zrep(sco(t, k)[v[k]], th, BW[k][:, v[k]])
                if eb > 0:
                    o[k] = (es / math.sqrt(eb), esr / np.sqrt(np.maximum(ebr, 1e-9)))
            return o
        # ---- (1) pary
        if sp == "D0":
            for name, (A, B) in PAIRS.items():
                need = [t for t in A + B if t not in scf]
                if need:
                    res["missing"].append({"pair": name, "tags": need}); continue
                valid = {k: np.all([np.isfinite(sco(t, k)) for t in A + B], axis=0) for k in ["sig"] + bk}
                R = {"A": A, "B": B, "rows": {k: int(valid[k].sum()) for k in valid}}
                for bi, b in enumerate(BN):
                    for wp in WPS:
                        za, zb = [zz(t, valid, bi, wp) for t in A], [zz(t, valid, bi, wp) for t in B]
                        for k in bk:
                            if not all(k in z for z in za + zb):
                                continue
                            r0 = np.mean([z[k][0] for z in za]) / np.mean([z[k][0] for z in zb])
                            rr = np.mean([z[k][1] for z in za], axis=0) / np.mean([z[k][1] for z in zb], axis=0)
                            e = float(np.std(rr))
                            R[f"{b}|{wp}|{k}"] = {"r": float(r0), "err": e, "nsig": float((r0 - 1) / e) if e > 0 else None}
                res["pairs"][name] = R
                print(f"[para] {name}: " + "  ".join(f"{b} mc {R.get(f'{b}|0.7|mc', {}).get('r', float('nan')):.3f}±{R.get(f'{b}|0.7|mc', {}).get('err', float('nan')):.3f}"
                                                      for b in BN), flush=True)
            for t in set(sum([x + y for x, y in PAIRS.values()], [])):
                if t in scf:
                    res["sources"][t] = scf[t]
        # ---- (2) kaskada H × E
        ref = DF.REF[sp]
        full = BASE[0] if sp == "D0" else None
        for et in ETAGS[sp]:
            if et not in scf or ref not in scf:
                res["missing"].append({"cascade": et, "tags": [t for t in (et, ref) if t not in scf]}); continue
            res["sources"][et] = scf[et]
            keys = ["sig", "bmc_all"] + bk
            valid = {k: np.isfinite(sco(et, k)) & np.isfinite(sco(ref, k)) & (np.isfinite(sco(full, k)) if full and full in scf else True) for k in keys}
            C = {"rows": {k: int(valid[k].sum()) for k in valid}}
            hthr = []
            for bi, b in enumerate(BN):
                v = {k: valid[k] & (ib[k] == bi) for k in keys}
                hs = sco(ref, "sig")[v["sig"]]
                if len(hs) < 50:
                    hthr.append(float("inf")); continue
                th_h = float(np.quantile(hs, 0.30)); hthr.append(th_h)
                hp = {k: v[k] & (sco(ref, k) >= th_h) for k in keys}
                es_e = sco(et, "sig")[hp["sig"]]
                o = {"n_sig_H": int(hp["sig"].sum()), "eps_b_H": {k: float((sco(ref, k)[v[k]] >= th_h).mean()) for k in bk if v[k].sum()},
                     "n_b_H": {k: int(hp[k].sum()) for k in bk}}
                both = np.concatenate([es_e, sco(et, "mc")[hp["mc"]]])
                o["auc_E_afterH_mc"] = DF.auc(both, np.arange(len(both)) < len(es_e)) if hp["mc"].sum() > 20 else None
                for wp in EWPS:
                    th = float(np.quantile(es_e, 1 - wp))
                    es, esr = zrep(es_e, th, BW["sig"][:, hp["sig"]])
                    w = {"thr": th, "eps_s_E": es}
                    for k in bk:
                        if hp[k].sum() < 30:
                            continue
                        eb, ebr = zrep(sco(et, k)[hp[k]], th, BW[k][:, hp[k]])
                        if eb > 0:
                            g = es / math.sqrt(eb)
                            w["G_" + k], w["G_" + k + "_err"], w["eps_b_E_" + k] = g, float(np.std(esr / np.sqrt(np.maximum(ebr, 1e-9)))), eb
                        # kontrola: samo zaostrzenie drzewa H do tej samej wydajności całkowitej 0,7·wp
                        th2 = float(np.quantile(hs, 1 - 0.7 * wp))
                        e1, e2 = float((sco(ref, k)[v[k]] >= th_h).mean()), float((sco(ref, k)[v[k]] >= th2).mean())
                        if e1 > 0 and e2 > 0:
                            w["G_Htight_" + k] = wp * math.sqrt(e1 / e2)
                    o[f"wp{wp}"] = w
                cand = [(o[f"wp{wp}"].get("G_mc", -1), wp) for wp in EWPS]
                o["wp_nom"] = max(cand)[1] if hp["mc"].sum() >= 50 and max(cand)[0] > 0 else 0.80        # wybór na tle MC, nie na danych
                if full and full in scf:                                                                  # pełna sieć tr.H+e04 (etap 2, ziarno 1) na tych samych wierszach
                    zf, zh = zz(full, valid, bi, 0.7), zz(ref, valid, bi, 0.7)
                    o["G_fullnet"] = {k: zf[k][0] / zh[k][0] for k in bk if k in zf and k in zh}
                vb = hp["bmc_all"]
                if vb.sum() > 200:
                    s_, m_ = sco(et, "bmc_all")[vb].astype(np.float64), mm["bmc_all"][vb]
                    o["corr_E_m_afterH"] = float(np.corrcoef(s_, m_)[0, 1])
                    c = AN.closure(m_, s_ >= o[f"wp{o['wp_nom']}"]["thr"], m0, lo, hi, sbm, 0.0, win)
                    if c:
                        o["sculpt_afterH"] = c
                C[b] = o
            C["h_thr_wp70"] = hthr
            res["cascade"][et] = C
            gn = lambda o, k: o[f"wp{o['wp_nom']}"].get(k, float("nan"))  # noqa: E731
            print(f"[kaskada] {et}: " + "  ".join(f"{b} WP_E {C[b]['wp_nom']} G mc {gn(C[b], 'G_mc'):.3f} sb428770 {gn(C[b], 'G_sb428770'):.3f} "
                                                  f"sb427929 {gn(C[b], 'G_sb427929'):.3f}" for b in BN if b in C), flush=True)
        del S, cache
    json.dump(res, open(a.out, "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    print("zapisano", a.out, "; brak:", res["missing"], flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
