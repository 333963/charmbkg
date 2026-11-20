"""R3 — uruchomienie: (toys) test wstrzyknięcia, (data) fit nominalny + bootstrap sparowany po zdarzeniach + pass/fail tym samym silnikiem.
  python r3_run.py toys [NTOY]     -> out/toys.npz
  python r3_run.py data [NBOOT]    -> out/data.json, out/boot.npz
"""
import json, os, sys, time
import multiprocessing as _mp
Pool = _mp.get_context("fork").Pool
import numpy as np
import r3_common as R

OUT = os.path.join(R.HERE, "out"); os.makedirs(OUT, exist_ok=True)
PF_WPS = (0.95, 0.80, 0.70, 0.50)            # = granice kategorii 1..4 w binowaniu głównym
INP = None
WID = json.load(open(os.path.join(R.HERE, "mc_width.json")))


def setup(run, b, wps=R.WP_MAIN):
    d = INP["data"][run][b]; mc = INP["mc"][b]
    thr = R.edges_for(INP, b, wps)
    return d, mc, thr, R.cat_index(d["score"], thr), R.cat_index(mc["score"], thr)


def merge_pf(H, j):
    """pass/fail na granicy kategorii j: fail = kategorie < j, pass = >= j."""
    return np.stack([H[:j].sum(0), H[j:].sum(0)])


def r_for(b, wps=R.WP_MAIN):
    return np.array(WID[b]["r"]) if wps is R.WP_MAIN else np.ones(len(wps) + 1)


def pack(res, n):
    """wektor: S (n), good, valid, posdef, mu, s, nll, covS (n*n)"""
    cov = res["covS"] if res["covS"] is not None else np.full((n, n), np.nan)
    return np.r_[res["S"], float(res["good"]), float(res["valid"]), float(res["posdef"]), res["mu"], res["s"], res["nll"], cov.ravel()]


# ------------------------------------------------------------------------------------------------ toys
def toy_job(a):
    run, b, c_inj, ntoy, seed = a
    d, mc, thr, cd, cm = setup(run, b)
    sh = INP["shapes"][b]; nc = 6
    Hd = R.hist_cats(d["dm"], cd, nc)
    F0 = R.CatFit(Hd, sh)
    bexp = np.stack([(lambda v: v[0] * R.bkg_shape(*v[1:]))(F0.bkg_only(c, side=R.SIDE_GEN, mu0=145.5)) for c in range(nc)])
    rank = np.argsort(np.argsort(mc["score"])) / max(len(mc["score"]) - 1, 1)
    w = 1 + c_inj * (0.5 - rank)
    Hs = R.hist_cats(mc["dm"], cm, nc, w=w)
    inrange = Hs.sum() / w.sum()                                   # ułamek sygnału MC w zakresie fitu
    off = [r for r in R.official()["rows"] if r["run"] == run and r["bin"] == b and r["wp"] == "0.95" and r["cls"] in ("all", "all427929")][0]
    NS = off["N_sig"]
    sexp = Hs / Hs.sum() * NS
    f_mc = np.bincount(cm, minlength=nc) / len(cm)
    f_true = np.bincount(cm, weights=w, minlength=nc) / w.sum()
    rng = np.random.default_rng(seed)
    rows = []
    r = r_for(b)
    nco = R.ncoef_rule(sexp + bexp, sh["mu"]); nco_pf = R.ncoef_rule(merge_pf(sexp + bexp, 3), sh["mu"])
    for _ in range(ntoy):
        H = rng.poisson(sexp + bexp).astype(float)
        res = R.CatFit(H, sh, r, nco).fit()
        pf = R.CatFit(merge_pf(H, 3), sh, None, nco_pf).fit()
        rows.append(np.r_[pack(res, nc), pack(pf, 2)])
    return (run, b, c_inj), {"rows": np.array(rows), "f_mc": f_mc, "f_true": f_true, "NS": NS, "B_cat": bexp.sum(1),
                             "B_win": bexp[:, np.abs(R.CEN - 145.5) < 2.5].sum(1), "S_cat": sexp.sum(1), "inrange": inrange}


def run_toys(ntoy):
    jobs = [(run, b, c, ntoy, 1000 + 97 * i) for i, (run, b, c) in enumerate((run, b, c) for run in R.RUNS for b in R.BINS for c in (1.0, 0.0))]
    # podział na kawałki, by wykorzystać 8 rdzeni
    chunks = []
    for (run, b, c, n, seed) in jobs:
        for k in range(4):
            chunks.append((run, b, c, n // 4, seed + k))
    with Pool(8) as P:
        res = P.map(toy_job, chunks, chunksize=1)
    out = {}
    for key, v in res:
        k = f"{key[0]}_{key[1]}_{key[2]:g}"
        if k + "_rows" in out:
            out[k + "_rows"] = np.vstack([out[k + "_rows"], v["rows"]])
        else:
            for q, val in v.items():
                out[f"{k}_{q}"] = np.asarray(val)
    np.savez(os.path.join(OUT, "toys.npz"), **out)


# ------------------------------------------------------------------------------------------------ data
def data_job(a):
    run, b, nboot, seed = a
    d, mc, thr, cd, cm = setup(run, b)
    sh = INP["shapes"][b]; nc = 6
    H = R.hist_cats(d["dm"], cd, nc)
    r = r_for(b)
    nco = R.ncoef_rule(H, sh["mu"]); nco_pf = {j: R.ncoef_rule(merge_pf(H, j), sh["mu"]) for j in (1, 2, 3, 4)}
    nom = R.CatFit(H, sh, r).fit()
    nom_r1 = R.CatFit(H, sh).fit() if not np.allclose(r, 1) else nom           # kontrola: szerokość wspólna
    pfn = {j: R.CatFit(merge_pf(H, j), sh).fit() for j in (1, 2, 3, 4)}
    # wariant: równe ułamki sygnału
    thr_e = R.edges_for(INP, b, R.WP_EQ)
    ce, cme = R.cat_index(d["score"], thr_e), R.cat_index(mc["score"], thr_e)
    He = R.hist_cats(d["dm"], ce, nc)
    eq = R.CatFit(He, sh).fit()
    res = {"run": run, "bin": b, "thr": thr.tolist(), "n_cat": H.sum(1).tolist(), "r": r.tolist(), "ncoef": nco.tolist(),
           "f_mc": (np.bincount(cm, minlength=nc) / len(cm)).tolist(), "n_mc": len(cm),
           "nom": {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in nom.items()},
           "nom_r1": {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in nom_r1.items()},
           "pf": {str(j): {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in p.items()} for j, p in pfn.items()},
           "eq": {"thr": thr_e.tolist(), "f_mc": (np.bincount(cme, minlength=nc) / len(cme)).tolist(), "n_cat": He.sum(1).tolist(),
                  **{k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in eq.items()}},
           "H": H.tolist(), "lam": R.CatFit(H, sh, r).lam(nom["x"]).tolist(),
           "lam_bkg": R.CatFit(H, sh, r).lam(np.where(np.arange(len(nom["x"])) % 5 == 2, 0.0, nom["x"])).tolist()}
    rng = np.random.default_rng(seed)
    rows = []
    for _ in range(nboot):
        w = rng.poisson(1.0, d["nev"]).astype(float)[d["inv"]]
        Hb = R.hist_cats(d["dm"], cd, nc, w=w)
        fb = R.CatFit(Hb, sh, r, nco).fit(x0=nom["x"])
        if not fb["good"]:
            f2 = R.CatFit(Hb, sh, r, nco).fit()
            if f2["good"] or f2["nll"] < fb["nll"]:
                fb = f2
        row = [pack(fb, nc)]
        for j in (1, 2, 3, 4):
            Fp = R.CatFit(merge_pf(Hb, j), sh, None, nco_pf[j])
            p = Fp.fit(x0=np.array(pfn[j]["x"]))
            if not p["good"]:
                p2 = Fp.fit()
                if p2["good"] or p2["nll"] < p["nll"]:
                    p = p2
            row.append(pack(p, 2))
        rows.append(np.concatenate(row))
    return (run, b), res, np.array(rows)


def run_data(nboot):
    jobs = []
    for run in R.RUNS:
        for b in R.BINS:
            for k in range(2 if b == "8-30" else 3):                 # kawałki; pierwszy kawałek niesie wynik nominalny
                jobs.append((run, b, nboot // (2 if b == "8-30" else 3) + 1, 5000 + 31 * len(jobs)))
    with Pool(8) as P:
        res = P.map(data_job, jobs, chunksize=1)
    J, B = {}, {}
    for key, r, rows in res:
        k = f"{key[0]}_{key[1]}"
        J.setdefault(k, r)
        B[k] = np.vstack([B[k], rows]) if k in B else rows
    json.dump(J, open(os.path.join(OUT, "data.json"), "w"))
    np.savez(os.path.join(OUT, "boot.npz"), **B)


if __name__ == "__main__":
    t0 = time.time()
    INP = R.prep_inputs()
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    {"toys": run_toys, "data": run_data}[sys.argv[1]](n)
    print(sys.argv[1], "gotowe", round(time.time() - t0), "s")
