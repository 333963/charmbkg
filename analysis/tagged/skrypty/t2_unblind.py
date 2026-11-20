#!/usr/bin/env python3
"""TAGGED krok 4 — JEDEN odczyt regionu tagowanego wg ../KRYTERIA.md. Uruchamiany raz: zapisuje out/UNBLINDED i odmawia powtórki
(TG_FORCE=1 tylko do odtworzenia tabel z tych samych danych i tego samego kodu). Wyjście: out/unblind.json, out/hists.npz."""
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import tg_lib as L  # noqa: E402
import t1_prep as P  # noqa: E402


DRY = os.environ.get("TG_DRY") == "1"      # test mechaniki: okno sygnału pozostaje zamaskowane (puste)


def quad(*x):
    return math.sqrt(sum(t * t for t in x if t is not None and np.isfinite(t)))


def syst_shape(item, B, n, N0):
    """Źródła z KRYTERIA §2: max |ΔN| per źródło."""
    sp, bins = item[3], item[4]
    H, nb = B["H"], len(B["H"])
    out = {}

    def nfit(S, nn=n, Hh=None):
        Hh = H if Hh is None else Hh
        return L.TFit(S, len(Hh), nn).fit(Hh)[0][0]
    d = []
    for dsc in (1.1, 0.9):
        d.append(nfit(L.sig_shape(sp, B["fits"], bins=bins, dsc=dsc)) - N0)
    for dd in (+L.DMU_S[sp], -L.DMU_S[sp]):
        d.append(nfit(L.sig_shape(sp, B["fits"], bins=bins, ddmu=dd)) - N0)
    out["ksztalt_sc_dmu"] = float(max(abs(t) for t in d))
    if len(bins) > 1 and sum(B["wpt_sb"]) > 0:
        out["wagi_pT"] = float(abs(nfit(L.sig_shape(sp, B["fits"], bins=bins, w=np.array(B["wpt_sb"]) + 1e-9)) - N0))
    out["rzad_tla"] = float(abs(nfit(B["S"], nn=n + 1) - N0))
    Ba = P.build(item, blind=DRY, rng=L.ALT[sp])
    out["zakres"] = float(abs(nfit(Ba["S"], Hh=Ba["H"]) - N0))
    if sp == "D0":
        d = [nfit(L.sig_shape(sp, B["fits"], bins=bins, frs=f)) - N0 for f in (1.05, 0.95)] + [nfit(L.sig_shape(sp, B["fits"], bins=bins, ks=k)) - N0 for k in (1.5, 0.5)]
        out["odbicie_kappa"] = float(max(abs(t) for t in d))
    out["razem"] = quad(*out.values())
    return out


def k4(B, n, ntag, nrep=200, seed=4279):
    rng = np.random.default_rng(seed)
    Hi = B["Hincl"]
    p = Hi / Hi.sum()
    exp = B["N_incl"] * ntag / Hi.sum()
    t = L.TFit(B["S"], len(Hi), n)
    T = []
    for _ in range(nrep):
        r = t.full(rng.multinomial(int(ntag), p).astype(float), cl95=False)
        T.append([r["N"], r["lo68"], r["hi68"], r["hesse"], r["q0"]])
    s = L.toy_summary(np.array(T), exp)
    s["expected"] = float(exp)
    return s


def main():
    t0 = time.time()
    mark = os.path.join(L.OUT, "UNBLINDED")
    if os.path.exists(mark) and os.environ.get("TG_FORCE") != "1":
        sys.exit("STOP: region już odsłonięty (out/UNBLINDED) — jeden odczyt")
    prep = json.load(open(os.path.join(L.OUT, "prep.json")))
    if not os.path.exists(mark) and not DRY:
        open(mark, "w").write(time.strftime("%F %T %Z") + "\n")
    nf = os.path.join(L.OUT, "toys_null_q0.npz")
    null = np.load(nf) if os.path.exists(nf) else np.load(nf) if not DRY else type("E", (), {"files": []})()
    R, hs = {}, {}
    for item in P.ITEMS:
        key, run, reg, sp, bins, e, xsel, role = item
        pr = prep["items"][key]
        if role == "R2bin" and pr["nsb_tag"] < 2000:
            continue
        B = P.build(item, blind=DRY)
        n = pr["nbkg"]
        t = L.TFit(B["S"], len(B["H"]), n)
        f = t.full(B["H"])
        r = {k: pr[k] for k in ("run", "reg", "sp", "bins", "eps", "xsel", "role", "nbkg", "N_incl", "N_incl_err", "P_tag_SB", "nsb_tag", "nsb_incl")}
        r["toys"] = pr.get("toys")
        r["fit"] = f
        r["n_tot"] = float(B["H"].sum()); r["n_win"] = float(B["H"][B["wm"]].sum())
        mu0 = t.mu(np.array([0.0] + f["x"][1:]))
        r["bkg_win"] = float(mu0[B["wm"]].sum()); r["sig_win_frac"] = float(B["S"][B["wm"]].sum())
        # reguła pokrycia / obciążenia z zabawek (KRYTERIA §2)
        scale, bias = 1.0, 0.0
        if r["toys"]:
            ts = sorted(r["toys"], key=lambda q: q["ninj"])
            if min(q["cov68"] for q in ts) < 0.63:
                scale = max(q["pull_rms"] for q in ts)
            bias = float(np.interp(f["N"], [q["ninj"] for q in ts], [q["bias"] for q in ts]))
        r["err_scale"], r["toy_bias_at_N"] = scale, bias
        lo, hi = (f["N"] - f["lo68"]) * scale, (f["hi68"] - f["N"]) * scale
        r["stat_lo"], r["stat_hi"] = lo, hi
        kq = key.replace("|", "__")
        if kq in null.files:
            q = null[kq]
            r["p_toys"] = float((np.sum(q >= f["q0"]) + 1) / (len(q) + 1)); r["n_null"] = int(len(q))
        if role in ("main", "K1", "pt5", "reg", "K6"):
            r["syst"] = syst_shape(item, B, n, f["N"])
        sy = (r.get("syst") or {}).get("razem", 0.0)
        N, Ni, eNi = f["N"], pr["N_incl"], pr["N_incl_err"]
        # ułamki i tagi przypadkowe
        afp_trig = reg in ("R1", "R2", "MBAFP")
        if afp_trig:
            fa, efa = L.FACC["AFP"]
            nacc, enacc = N * fa, abs(N) * efa
        else:
            fa, efa = L.FACC[run]
            nacc, enacc = Ni * fa * pr["P_tag_SB"], Ni * efa * pr["P_tag_SB"]
        r["N_acc"], r["N_acc_err"], r["N_true"] = nacc, enacc, N - nacc
        if reg in ("MBTS", "J12", "MBTS_A", "MBTS_C") and Ni > 0:
            r["f_raw"] = N / Ni
            rel = eNi / Ni
            Nt = N - nacc
            r["f_true"] = Nt / Ni
            r["f_true_lo68"] = quad(lo, sy, enacc, Nt * rel) / Ni; r["f_true_hi68"] = quad(hi, sy, enacc, Nt * rel) / Ni
            r["f_stat_lo"], r["f_stat_hi"], r["f_syst"] = lo / Ni, hi / Ni, quad(sy, enacc, Nt * rel) / Ni
            l95, h95 = (N - f["lo95"]) * scale, (f["hi95"] - N) * scale
            r["f_true_lo95"] = quad(l95, 1.96 * sy, 1.96 * enacc, 1.96 * Nt * rel) / Ni; r["f_true_hi95"] = quad(h95, 1.96 * sy, 1.96 * enacc, 1.96 * Nt * rel) / Ni
            if xsel == "T2" and reg == "MBTS":
                r["pred"] = {}
                for nm, (c, a, b) in L.PRED[sp].items():
                    i68 = (r["f_true"] - r["f_true_lo68"], r["f_true"] + r["f_true_hi68"]); i95 = (r["f_true"] - r["f_true_lo95"], r["f_true"] + r["f_true_hi95"])
                    v = "zgodny" if (i68[0] <= b and i68[1] >= a) else ("napięcie" if (i95[0] <= b and i95[1] >= a) else "niezgodny")
                    sg = r["f_true_hi68"] if c > r["f_true"] else r["f_true_lo68"]
                    r["pred"][nm] = {"central": c, "lo": a, "hi": b, "werdykt": v, "pull_central": (r["f_true"] - c) / sg, "central_in_68": bool(i68[0] <= c <= i68[1])}
        if reg in L.LUMI or (run == "427929" and reg in ("MBTS", "MBAFP")):
            r["N_per_nb"] = N / L.LUMI.get(reg, L.LUMI["427929"])
        if role in ("main", "K1"):
            r["K4"] = k4(B, n, r["n_tot"])
        R[key] = r
        if role in ("main", "K1", "pt5", "reg", "K5", "K6"):
            hs[kq + "__H"] = B["H"]; hs[kq + "__mu"] = t.mu(np.array(f["x"])); hs[kq + "__bkg"] = mu0
        print(f"[{time.time() - t0:5.0f}s] {key}: N = {N:.1f} −{lo:.1f} +{hi:.1f} (syst {sy:.1f}), Z = {f['Z']:.2f}, n_okno {r['n_win']:.0f}, tło w oknie {r['bkg_win']:.1f}"
              + (f", f_raw = {100 * r['f_raw']:.3f} %" if "f_raw" in r else ""), flush=True)
    json.dump({"items": R, "time": time.strftime("%F %T"), "prep_time": prep["time"]}, open(os.path.join(L.OUT, "unblind.json"), "w"), indent=1, ensure_ascii=False)
    np.savez_compressed(os.path.join(L.OUT, "hists.npz"), **hs)
    print("zapisano out/unblind.json")


if __name__ == "__main__":
    main()
