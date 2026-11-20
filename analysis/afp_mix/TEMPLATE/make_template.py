#!/usr/bin/env python3
"""ZADANIE 1 — szablon tła przypadkowego tagu AFP (projekt B, AFP-MIX v2, WSTĘPNE).

Dla każdego runu × strony (A, C, AorC) × gatunku (Dplus, D0, Dstar) × przedziału pT: N_tag, N_acc, f_acc z rozbiciem błędu
oraz kształt masy składowej przypadkowej (Σ p_acc po kandydatach w binach masy) porównany z kształtem masy wszystkich
kandydatów bez wymogu tagu (χ² z kowariancją: bootstrap dawców ⊕ próbkowanie zdarzeń-biorców, bez podwójnego liczenia;
warianty „cons” i „samp”; kalibracja przez permutację wag).
Estymator = nominał v2 (``afp_mix.mix_region`` bez zmian); systematyki względne (pula, ekstrapolacja, ½Δ metod) z
``out_v2/result_v2.json``. NIE liczy uzysków sygnału, fitów masy ani przekrojów; nie histogramuje masy próbki tagowanej.

Użycie: make_template.py [--samples a,b] [--n-boot 500]
"""
from __future__ import annotations
import argparse, csv, json, math, os, time
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "4")
import numpy as np
from scipy import sparse, stats
from common import HERE, INDIR, SAMPLES, RD, M, load, jsonable

PT_BINS = (("3-5", 3.0, 5.0), ("5-8", 5.0, 8.0), ("8-30", 8.0, 30.0), ("3-30", 3.0, 30.0), ("3-150_v2check", 3.0, 150.0))
Q = 5e-5   # pół kwantu zapisu masy (0,0001 GeV) — brzegi binów między wartościami skwantowanymi
MASS_EDGES = {
    "Dplus": np.concatenate([1.72 + 0.01 * np.arange(30) - Q, [2.02 + Q]]),
    "D0": np.concatenate([1.75 + 0.01 * np.arange(25) - Q, [2.00 + Q]]),
    "Dstar": np.concatenate([0.139 + 0.001 * np.arange(21) - Q, [0.160 + Q]]),
}
REGIONS = ("A", "C", "AorC")
XB_S = np.array([0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.12, 0.14])
XB_O = np.array([0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.2, 0.3, 0.5, np.inf])


def run_sample(s, a, v2, log):
    d = load(s)
    Dn, T, ca, cfg = d["Dn"], d["T"], d["cands"], d["cfg"]
    nt = len(T)
    boot = M.Boot.for_events(Dn, cfg, a.n_boot, 0)
    t0 = time.time()
    ms, m0 = M.mix_region(T, Dn, cfg, "AorC")
    log(f"{s}: mieszanie {time.time() - t0:.0f} s")
    parts = {"AorC": list(ms.values()), "A": [ms["A"].relogic("A")], "C": [ms["C"].relogic("C")]}
    own = {"AorC": m0.own, "A": parts["A"][0].own, "C": parts["C"][0].own}
    # macierze par (otagowany dawca × biorca) i wagi bootstrapu
    PM = {}
    for reg, pl in parts.items():
        PM[reg] = []
        for mx in pl:
            k = mx.acc
            P = sparse.csr_matrix((mx.pw[k].astype(np.float64), (mx.pd[k], mx.pt[k])), shape=(mx.ntd, nt))
            PM[reg].append((P, boot.get(mx.TD["key"]).astype(np.float64)))
    p_acc = {reg: np.asarray(sum(P.sum(0) for P, _ in PM[reg])).ravel() for reg in REGIONS}
    # wariancja bootstrapowa p_acc per biorca (wagi Poissona(1) per dawca): v_t = Σ_d P[d,t]²
    v_acc = {reg: np.asarray(sum(P.multiply(P).sum(0) for P, _ in PM[reg])).ravel() for reg in REGIONS}
    # --- kontrola zgodności z out_v2
    chk = {}
    for reg in REGIONS:
        na = float(p_acc[reg].sum())
        rb = sum(W.T @ np.asarray(P.sum(1)).ravel() for P, W in PM[reg])
        e2 = v2["estimates_by_region"][reg]
        chk[reg] = {"n_tag": int(own[reg].sum()), "n_acc": na, "n_acc_err_stat": float(rb.std(ddof=1)), "v2_n_tag": e2["n_tag"], "v2_n_acc": e2["n_acc"],
                    "v2_n_acc_err_stat": e2["n_acc_err_stat"]}
        log(f"{s} [{reg}] zdarzenia: N_tag {chk[reg]['n_tag']} (v2 {e2['n_tag']:.0f}), N_acc {na:.2f} ± {chk[reg]['n_acc_err_stat']:.2f} (v2 {e2['n_acc']:.2f} ± {e2['n_acc_err_stat']:.2f})")
        assert abs(na - e2["n_acc"]) < 1e-3 * e2["n_acc"] and chk[reg]["n_tag"] == int(e2["n_tag"]), "niezgodność z out_v2"
    sy = v2["syst_rel"]
    rel = {"pool": sy["pool"], "extrap": sy["extrapolation"], "half": sy["half_diff_methods"]}
    rel_syst = math.sqrt(sum(x * x for x in rel.values()))
    # --- profil obszaru użycia (dla ZADANIA 2): rozkład wag par „acc” w (ξ_cal strona tagu, ξ_cal strona przeciwna), nSelBL
    usage = {}
    nsel_t = T["nsel"].to_numpy()
    for sd, mx in ms.items():
        o = "C" if sd == "A" else "A"
        k = mx.acc & mx.dtag[sd]
        xs, xo, w = mx.T[f"xical_{sd}"][mx.pt[k]], mx.T[f"xical_{o}"][mx.pt[k]], mx.pw[k].astype(float)
        H, _, _ = np.histogram2d(xs, xo, [XB_S, XB_O], weights=w)
        usage[sd] = {"edges_same": XB_S.tolist(), "edges_opp": XB_O.tolist(), "w2d": H.tolist(), "w_total": float(w.sum()),
                     "w_outside_same_edges": float(w[xs >= XB_S[-1]].sum()), "mean_ln_xical_same": float((w * np.log(np.maximum(xs, 1e-6))).sum() / w.sum()),
                     "mean_nsel": float((w * nsel_t[mx.pt[k]]).sum() / w.sum()),
                     "nsel_hist": np.histogram(nsel_t[mx.pt[k]], list(cfg.act_edges) and [0] + list(cfg.act_edges) + [10 ** 6], weights=w)[0].tolist(),
                     "frac_opp_proton": float(w[mx.T[f"np_{o}"][mx.pt[k]] >= 1].sum() / w.sum())}
    rows, shapes = [], []
    rng = np.random.default_rng(12345)
    tc, pt_c, sp_c = ca["t"].to_numpy(), ca["pt"].to_numpy(), ca["species"].to_numpy()
    for name, sp in RD.SPECIES.items():
        ms_ = sp_c == sp["code"]
        v = ca["dm" if name == "Dstar" else "m"].to_numpy()
        ed = MASS_EDGES[name]
        nb = len(ed) - 1
        for lab, lo, hi in PT_BINS:
            mb = ms_ & (pt_c >= lo) & (pt_c < hi)
            tt, vv = tc[mb], v[mb]
            pk = (vv >= sp["peak"][0]) & (vv <= sp["peak"][1])
            n_peak = np.bincount(tt[pk], minlength=nt).astype(np.float64)
            ib = np.digitize(vv, ed) - 1
            okb = (ib >= 0) & (ib < nb)
            N = sparse.csr_matrix((np.ones(okb.sum()), (tt[okb], ib[okb])), shape=(nt, nb))
            n_t = np.asarray(N.sum(1)).ravel()
            h_all = np.asarray(N.sum(0)).ravel()
            sh = h_all / h_all.sum()
            ipk = (ed[:-1] + Q >= sp["peak"][0] - 1e-9) & (ed[1:] - Q <= sp["peak"][1] + 1e-9)   # biny w całości w oknie piku
            for reg in REGIONS:
                # ---- liczności i f_acc
                n_tag = float(n_peak[own[reg]].sum())
                e_tag2 = float((n_peak[own[reg]] ** 2).sum())
                n_acc, reps = 0.0, np.zeros(a.n_boot)
                Cm = np.zeros((0, nb))
                h_acc, reps_h = np.zeros(nb), np.zeros((a.n_boot, nb))
                for P, W in PM[reg]:
                    c = P @ n_peak
                    n_acc += float(c.sum())
                    reps += W.T @ c
                    Cd = np.asarray((P @ N).todense())
                    h_acc += Cd.sum(0)
                    reps_h += W.T @ Cd
                e_acc = float(reps.std(ddof=1))
                f = n_acc / n_tag if n_tag > 0 else float("nan")
                r_don = e_acc / n_acc if n_acc > 0 else float("nan")
                r_tag = math.sqrt(e_tag2) / n_tag if n_tag > 0 else float("nan")
                row = {"run": d["run"], "region": reg, "species": name, "pt_bin": lab, "pt_lo": lo, "pt_hi": hi, "mass_var": sp["var"],
                       "win_lo": sp["peak"][0], "win_hi": sp["peak"][1], "n_cand_all_win": float(n_peak.sum()), "n_tag": n_tag, "n_acc": n_acc,
                       "n_acc_err_stat": e_acc, "n_acc_err_pool": n_acc * rel["pool"], "n_acc_err_extrap": n_acc * rel["extrap"], "n_acc_err_half": n_acc * rel["half"],
                       "n_acc_err_syst": n_acc * rel_syst, "f_acc": f, "f_acc_err_stat": f * math.hypot(r_don, r_tag) if f == f else float("nan"),
                       "f_acc_err_stat_donors": f * r_don if f == f else float("nan"), "f_acc_err_stat_ntag": f * r_tag if f == f else float("nan"),
                       "f_acc_err_pool": f * rel["pool"], "f_acc_err_extrap": f * rel["extrap"], "f_acc_err_half": f * rel["half"],
                       "f_acc_err_syst": f * rel_syst, "f_acc_err_tot": f * math.sqrt(r_don ** 2 + r_tag ** 2 + rel_syst ** 2) if f == f else float("nan"),
                       "mean_p_acc_all_win": n_acc / n_peak.sum() if n_peak.sum() > 0 else float("nan")}
                # ---- kształt masy: d = h_acc − Σh_acc · s_all ; kowariancja = bootstrap dawców ⊕ próbkowanie biorców
                tot = h_acc.sum()
                dvec = h_acc - tot * sh
                Db = reps_h - reps_h.sum(1)[:, None] * sh[None, :]
                cov_d = np.cov(Db, rowvar=False, ddof=1)
                p = p_acc[reg]
                pbar = float((p * n_t).sum() / n_t.sum())
                q = (p - pbar) ** 2
                A_ = np.asarray((N.T @ sparse.diags(q) @ N).todense())
                b_ = N.T @ (q * n_t)
                c_ = float((q * n_t * n_t).sum())
                cov_t = A_ - np.outer(b_, sh) - np.outer(sh, b_) + c_ * np.outer(sh, sh)
                # Kowariancja testu H0: „prawdziwa składowa przypadkowa (wagi π_t) ma kształt próbki nietagowanej”:
                #   C = C_boot + C_samp − C_overlap,
                #   C_boot    = bootstrap dawców przy ustalonych biorcach: Σ Var(p)·(zrealizowane odchylenie komórki)² — szum wag, także między
                #               komórkami o różnym kształcie masy (korelacja wag z jednej puli dawców),
                #   C_samp    = próbkowanie biorców przy zrealizowanych wagach: Σ_t (p_t − p̄)² g_t g_tᵀ,
                #   C_overlap = Σ_t v_t g_t g_tᵀ, v_t = Var_boot(p_t): część C_samp policzona już w C_boot (E[(p−p̄)²] = (π−π̄)² + v).
                # Warianty: „cons” = C_boot + C_samp (górna granica wariancji), „samp” = samo C_samp (H0 ścisłe: masa niezależna od wag;
                # kalibrowane permutacją wag).
                vq = v_acc[reg]
                bv = N.T @ (vq * n_t)
                cov_v = np.asarray((N.T @ sparse.diags(vq) @ N).todense()) - np.outer(bv, sh) - np.outer(sh, bv) + float((vq * n_t * n_t).sum()) * np.outer(sh, sh)
                cov = cov_d + cov_t - cov_v
                use = h_all > 0
                pinv = lambda c_: np.linalg.pinv(c_[np.ix_(use, use)], rcond=1e-10, hermitian=True)
                chi2 = float(dvec[use] @ pinv(cov) @ dvec[use])
                ndf = int(use.sum() - 1)
                chi2_cons = float(dvec[use] @ pinv(cov_t + cov_d) @ dvec[use])
                chi2_samp = float(dvec[use] @ pinv(cov_t) @ dvec[use])
                min_eig = float(np.linalg.eigvalsh(cov[np.ix_(use, use)])[1] / max(np.trace(cov) / use.sum(), 1e-300))
                # stosunek udziału okna piku: (Σ_pik h_acc / Σ h_acc) / (Σ_pik h_all / Σ h_all)
                gp = ipk.astype(float)
                fr_acc, fr_all = float(h_acc[ipk].sum() / tot) if tot > 0 else float("nan"), float(sh[ipk].sum())
                e_fr = math.sqrt(float(gp @ cov @ gp)) / tot if tot > 0 else float("nan")
                # nachylenie liniowe stosunku kształtów r_i = (h_acc_i/tot)/s_i względem masy (1/GeV), z kowariancją
                xm = 0.5 * (ed[:-1] + ed[1:])
                xw = xm - float((sh * xm).sum())
                den = float((sh * xw * xw).sum())
                slope = float((dvec * xw).sum() / tot / den) if tot > 0 else float("nan")
                e_slope = math.sqrt(float(xw @ cov @ xw)) / tot / den if tot > 0 else float("nan")
                # kontrola kalibracji χ²: p_acc przetasowane między biorcami kategorii (H0 z konstrukcji), kowariancja tylko próbkowania
                hasc = np.flatnonzero(n_t > 0)
                pp = p.copy()
                pp[hasc] = rng.permutation(p[hasc])
                qn = (pp - pbar) ** 2
                dn = np.asarray(N.T @ (pp - pbar)).ravel() - float(((pp - pbar) * n_t).sum()) * sh
                bn = N.T @ (qn * n_t)
                cn = np.asarray((N.T @ sparse.diags(qn) @ N).todense()) - np.outer(bn, sh) - np.outer(sh, bn) + float((qn * n_t * n_t).sum()) * np.outer(sh, sh)
                chi2_null = float(dn[use] @ np.linalg.pinv(cn[np.ix_(use, use)], rcond=1e-10, hermitian=True) @ dn[use])
                row.update({"shape_chi2_null_perm": chi2_null, "shape_chi2_cons": chi2_cons, "shape_p_cons": float(stats.chi2.sf(chi2_cons, ndf)),
                            "shape_chi2_samp": chi2_samp, "shape_p_samp": float(stats.chi2.sf(chi2_samp, ndf)), "shape_cov_min_eig_rel": min_eig,
                            "shape_overlap_frac": float(np.trace(cov_v) / np.trace(cov_t)) if np.trace(cov_t) > 0 else float("nan"), "shape_chi2": chi2, "shape_ndf": ndf, "shape_p": float(stats.chi2.sf(chi2, ndf)), "shape_n_acc_fullrange": float(tot),
                            "shape_n_all_fullrange": float(h_all.sum()), "peakfrac_ratio": fr_acc / fr_all if fr_all > 0 else float("nan"),
                            "peakfrac_ratio_err": e_fr / fr_all if fr_all > 0 else float("nan"), "shape_slope_per_GeV": slope, "shape_slope_err": e_slope,
                            "shape_boot_over_sampling_var": float(np.trace(cov_d) / np.trace(cov_t)) if np.trace(cov_t) > 0 else float("nan")})
                rows.append(row)
                eh = np.sqrt(np.diag(np.cov(reps_h, rowvar=False, ddof=1)))
                ed_sh = np.sqrt(np.clip(np.diag(cov), 0, None))
                shapes.append({"run": d["run"], "region": reg, "species": name, "pt_bin": lab, "mass_var": sp["var"], "edges": (ed[:-1] + Q).tolist() + [float(ed[-1] - Q)],
                               "w_acc": h_acc.tolist(), "w_acc_err_stat": eh.tolist(), "n_all": h_all.tolist(),
                               "ratio": (np.where(use, h_acc / np.maximum(tot * sh, 1e-300), np.nan)).tolist(),
                               "ratio_err": (np.where(use, ed_sh / np.maximum(tot * sh, 1e-300), np.nan)).tolist()})
            log(f"{s} {name} pT {lab}: " + ", ".join(f"{r['region']}: N_tag {r['n_tag']:.0f} N_acc {r['n_acc']:.2f} f {r['f_acc']:.4f} χ²/ndf {r['shape_chi2']:.1f}/{r['shape_ndf']}" for r in rows[-3:]))
    return {"sample": s, "run": d["run"], "check_vs_v2": chk, "syst_rel_run": dict(rel, syst=rel_syst, stat_incl=sy["stat"], total_incl=sy["total_rel"],
            kappa=sy["kappa"], kappa_err=sy["kappa_err"]), "usage_profile": usage, "n_targets": nt, "n_donors": len(Dn), "rows": rows, "shapes": shapes}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", default=",".join(SAMPLES))
    ap.add_argument("--n-boot", type=int, default=500)
    a = ap.parse_args()
    t00 = time.time()
    log = lambda m: print(f"[{time.time() - t00:6.0f} s] {m}", flush=True)
    V2 = {r["sample"]: r for r in json.load(open(os.path.join(INDIR, "out_v2", "result_v2.json")))["results"]}
    out = []
    for s in a.samples.split(","):
        out.append(run_sample(s, a, V2[s], log))
        json.dump(jsonable({"status": "WSTEPNE", "schema": "template v1 (README.md)", "method": "AFP-MIX v2 estymator A (charmbkg.afp_mix.mix_region), systematyki z out_v2/result_v2.json",
                            "n_boot": a.n_boot, "pt_bins": PT_BINS, "samples": out}), open(os.path.join(HERE, "template.json"), "w"), indent=1)
    rows = [r for o in out for r in o["rows"]]
    with open(os.path.join(HERE, "template.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})
    with open(os.path.join(HERE, "template_shapes.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["run", "region", "species", "pt_bin", "mass_var", "bin_lo", "bin_hi", "w_acc", "w_acc_err_stat", "n_all", "ratio_acc_over_all", "ratio_err"])
        for o in out:
            for sh in o["shapes"]:
                for i in range(len(sh["w_acc"])):
                    w.writerow([sh["run"], sh["region"], sh["species"], sh["pt_bin"], sh["mass_var"], f"{sh['edges'][i]:.4f}", f"{sh['edges'][i + 1]:.4f}",
                                f"{sh['w_acc'][i]:.6g}", f"{sh['w_acc_err_stat'][i]:.4g}", f"{sh['n_all'][i]:.0f}", f"{sh['ratio'][i]:.5g}", f"{sh['ratio_err'][i]:.4g}"])
    log("zapisano template.json, template.csv, template_shapes.csv")


if __name__ == "__main__":
    main()
