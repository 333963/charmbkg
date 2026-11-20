"""R3 — kontrola MC (bez danych): test gradientu silnika, ułamki f_c MC, zależność szerokości sygnału od binu score."""
import json, sys, time
import numpy as np
from iminuit import Minuit
from scipy.stats import norm
import r3_common as R

inp = R.prep_inputs()
out = {}
rng = np.random.default_rng(1)
for b in R.BINS:
    mc = inp["mc"][b]; sh = inp["shapes"][b]
    thr = R.edges_for(inp, b, R.WP_MAIN)
    cat = R.cat_index(mc["score"], thr)
    n = len(cat); f = np.bincount(cat, minlength=6) / n
    # szerokość rdzenia per kategoria przy (k, f) ustalonych z kształtu łącznego; fit nie-binowany w 143–148,5 MeV
    def fit_sc(x):
        x = x[(x > 143) & (x < 148.5)]
        def nll(mu, sc):
            c = lambda t: sh["f"] * norm.cdf(t, mu, sc) + (1 - sh["f"]) * norm.cdf(t, mu, sc * sh["k"])
            p = (sh["f"] * norm.pdf(x, mu, sc) + (1 - sh["f"]) * norm.pdf(x, mu, sc * sh["k"])) / (c(148.5) - c(143))
            return -np.sum(np.log(np.clip(p, 1e-300, None)))
        m = Minuit(nll, mu=sh["mu"], sc=sh["sc"]); m.errordef = 0.5; m.limits["sc"] = (0.15, 2); m.migrad(); m.hesse()
        return m.values["sc"], m.errors["sc"], m.values["mu"], m.errors["mu"]
    s_all = fit_sc(mc["dm"])
    rows = []
    for c in range(6):
        sc, e, mu, emu = fit_sc(mc["dm"][cat == c])
        r = sc / s_all[0]; er = e / s_all[0]
        rows.append({"cat": c, "n": int((cat == c).sum()), "f": float(f[c]), "sc": sc, "sc_err": e, "r": r, "r_err": er,
                     "nsig": (r - 1) / er, "mu": mu, "mu_err": emu})
    need = any(abs(x["r"] - 1) > 0.05 and abs(x["nsig"]) > 3 for x in rows)
    out[b] = {"n_mc": n, "thr": thr.tolist(), "sc_all": s_all[0], "sc_shape_json": sh["sc"], "rows": rows, "width_per_bin": bool(need),
              "r": [x["r"] if need else 1.0 for x in rows]}
    print(b, "n_mc", n, "thr", np.round(thr, 4), "f", np.round(f, 4), "width_per_bin:", need)
    for x in rows:
        print("   cat %d n %6d sc %.3f±%.3f r %.3f (%+.1fσ) mu %.3f±%.3f" % (x["cat"], x["n"], x["sc"], x["sc_err"], x["r"], x["nsig"], x["mu"], x["mu_err"]))
json.dump(out, open("mc_width.json", "w"), indent=1)

# test gradientu + szybkość na pseudo-próbce z MC (sygnał MC + płaskie tło)
b = "5-8"; mc = inp["mc"][b]; thr = R.edges_for(inp, b, R.WP_MAIN); cat = R.cat_index(mc["score"], thr)
H = R.hist_cats(mc["dm"], cat, 6) * 0.1
H = rng.poisson(H + 20 * R.bkg_shape(0.5, 1.0, 0.0)[None, :] * np.array([500, 100, 30, 20, 10, 5])[:, None]).astype(float)
F = R.CatFit(H, inp["shapes"][b]); x0 = F.start()
g = F.grad(x0); gn = np.array([(F.nll(x0 + d) - F.nll(x0 - d)) / (2 * h) for h, d in ((1e-5 * max(1, abs(v)), np.eye(len(x0))[i] * 1e-5 * max(1, abs(v))) for i, v in enumerate(x0))])
print("grad max rel diff", np.max(np.abs(g - gn) / (np.abs(gn) + 1e-3)))
t = time.time(); r = F.fit(x0); print("fit", time.time() - t, "s", r["valid"], r["posdef"], r["S"].round(1), "true", (np.bincount(cat, minlength=6) * 0.1).round(1), r["at_limit"])
