#!/usr/bin/env python3
"""Wielkości pochodne z jump_tests.json (+ out_v2/result_v2.json): tabele testów (i)–(iii), test symetrii stron, propozycja systematyki.
Zapis: jump_summary.json. Importowane przez make_report.py."""
from __future__ import annotations
import json, math, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
V2P = os.path.join(os.path.dirname(HERE), "runs", "lxplus", "out_v2", "result_v2.json")


def R(k, p, relp=0.0):
    k, p = float(k), float(p)
    if not p > 0:
        return {"k": k, "pred": p, "R": None, "err": None}
    return {"k": k, "pred": p, "R": k / p, "err": math.sqrt(max(k, p) / p ** 2 + (relp * k / p) ** 2)}


def wfit(x, y, e):
    x, y, e = map(np.asarray, (x, y, e))
    w = 1 / e ** 2
    xb = (w * x).sum() / w.sum()
    sxx = (w * (x - xb) ** 2).sum()
    b = (w * (x - xb) * y).sum() / sxx
    return float(b), float(1 / math.sqrt(sxx))


def derive():
    J = json.load(open(os.path.join(HERE, "jump_tests.json")))["results"]
    V2 = {r["sample"]: r for r in json.load(open(V2P))["results"]}
    out = {}
    for r in J:
        v2 = V2[r["sample"]]
        o = {"run": r["run"]}
        ed = r["test_i"]["edges_opp"]
        # ---------- (i) profil R(ξ_cal,o), strony razem
        for lab in ("hasD", "hasD_npo0", "all", "noD"):
            k = sum(np.array(r["test_i"]["sides"][s]["tables"][lab]["k"]) for s in "AC")
            p = sum(np.array(r["test_i"]["sides"][s]["tables"][lab]["pred"]) for s in "AC")
            o[f"i_profile_{lab}"] = {"edges": ed, "bins": [R(a, b) for a, b in zip(k, p)], "low": R(k[:7].sum(), p[:7].sum()), "high": R(k[7:].sum(), p[7:].sum())}
            lo, hi = o[f"i_profile_{lab}"]["low"], o[f"i_profile_{lab}"]["high"]
            o[f"i_profile_{lab}"]["low_over_high"] = [lo["R"] / hi["R"], lo["R"] / hi["R"] * math.sqrt(1 / lo["k"] + 1 / hi["k"])]
        for s in "AC":
            t = r["test_i"]["sides"][s]["tables"]["hasD"]
            k, p = np.array(t["k"]), np.array(t["pred"])
            o[f"i_side_{s}_hasD"] = {"low": R(k[:7].sum(), p[:7].sum()), "high": R(k[7:].sum(), p[7:].sum())}
        for key in ("mirror_all", "mirror_npo0", "nominal_npo0", "act_fine", "nominal"):
            o[key] = r[key]
        # ---------- (i-b) symetria stron w ćwiartce obu ξ_cal > 0,14 (hasD): nachylenia R względem ln ξ_cal własnej i przeciwnej strony
        for lab in ("hasD", "all"):
            K = sum(np.array(r["test_i"]["sides"][s]["tables"][lab + "_2d"]["k"]) for s in "AC")[:, 7:]
            Pd = sum(np.array(r["test_i"]["sides"][s]["tables"][lab + "_2d"]["pred"]) for s in "AC")[:, 7:]
            Nn = sum(np.array(r["test_i"]["sides"][s]["tables"][lab + "_2d"]["n"]) for s in "AC")[:, 7:]
            Ls = sum(np.array(r["test_i"]["sides"][s]["tables"][lab + "_2d"]["sumln_same"]) for s in "AC")[:, 7:]
            Lo = sum(np.array(r["test_i"]["sides"][s]["tables"][lab + "_2d"]["sumln_opp"]) for s in "AC")[:, 7:]
            res = {}
            for nm, ax, L in (("same", 1, Ls), ("opp", 0, Lo)):
                k, p, x = K.sum(ax)[:3], Pd.sum(ax)[:3], (L.sum(ax) / np.maximum(Nn.sum(ax), 1))[:3]
                b, eb = wfit(x, k / p, np.sqrt(np.maximum(k, p)) / p)
                res[nm] = {"k": k.tolist(), "pred": p.tolist(), "R": (k / p).tolist(), "mean_ln": x.tolist(), "slope": b, "slope_err": eb}
            res["diff"] = res["opp"]["slope"] - res["same"]["slope"]
            res["diff_err"] = math.hypot(res["opp"]["slope_err"], res["same"]["slope_err"])
            res["z"] = res["diff"] / res["diff_err"]
            o[f"ib_symmetry_{lab}"] = res
        # nachylenie względem strony przeciwnej w pełnym zakresie (wszystkie przedziały z k > 0, hasD) — przybliżone środki przedziałów
        # ---------- (iii) zależności przy małym / dużym ξ_cal,o
        dep = {}
        for var in ("nsel", "mu", "lb"):
            d0 = r["test_i"]["sides"]["A"]["dependence_low_xical_opp"][var]
            dep[var] = {"edges": d0["edges"]}
            for lab in ("hasD_low", "hasD_high", "all_low", "all_high"):
                k = sum(np.array(r["test_i"]["sides"][s]["dependence_low_xical_opp"][var][lab]["k"]) for s in "AC")
                p = sum(np.array(r["test_i"]["sides"][s]["dependence_low_xical_opp"][var][lab]["pred"]) for s in "AC")
                dep[var][lab] = [R(a, b) for a, b in zip(k, p)]
            dep[var]["hasD_low_over_high"] = [[a["R"] / b["R"], a["R"] / b["R"] * math.sqrt(1 / max(a["k"], 1) + 1 / max(b["k"], 1))] if a["R"] and b["R"] else [None, None]
                                              for a, b in zip(dep[var]["hasD_low"], dep[var]["hasD_high"])]
        o["iii_dep"] = dep
        o["iii_dep_note_edges_side_A"] = True
        # ---------- (ii)
        O, P, KP = np.array(r["test_ii"]["obs"]), np.array(r["test_ii"]["pred"]), np.array(r["test_ii"]["k_pool"])
        t2 = {"edges_xical_same": r["test_ii"]["edges_xical_same"], "edges_margin": r["test_ii"]["edges_margin"], "edges_xical_opp": r["test_ii"]["edges_xical_opp"], "slices": r["test_ii"]["slices"]}
        for lab, hs in (("hasD", [1]), ("all", [0, 1])):
            o_, p_ = O[..., hs], P[..., hs]
            # błąd przewidywania: wycinek × hasD niezależnie, względny 1/√K_puli
            def tot(mask_fn):
                ko = po = vp = 0.0
                for i in range(O.shape[0]):
                    for h in hs:
                        a, b = mask_fn(O[i, ..., h]), mask_fn(P[i, ..., h])
                        ko, po = ko + a, po + b
                        vp += b * b / max(KP[i, h], 1.0)
                d = R(ko, po)
                if d["R"] is not None:
                    d["err"] = math.sqrt(max(ko, po) + vp) / po
                return d
            t2[lab] = {"by_margin": [tot(lambda a, m=m: a[:, m, :].sum()) for m in range(3)],
                       "by_xical_margin": [[tot(lambda a, m=m, x=x: a[x, m, :].sum()) for m in range(3)] for x in range(O.shape[1])],
                       "by_opp_margin": [[tot(lambda a, m=m, q=q: a[:, m, q].sum()) for m in range(3)] for q in range(3)],
                       "by_slice_margin": [[R(O[i, :, m, :, hs].sum(), P[i, :, m, :, hs].sum()) for m in range(3)] for i in range(O.shape[0])],
                       "in_window_far": tot(lambda a: a[:, 2, :].sum()), "midfar": tot(lambda a: a[:, 1:, :].sum()),
                       "near_minus_far": None}
            n_, f_ = t2[lab]["by_margin"][0], t2[lab]["by_margin"][2]
            t2[lab]["near_minus_far"] = [n_["R"] - f_["R"], math.hypot(math.sqrt(max(n_["k"], n_["pred"])) / n_["pred"], math.sqrt(max(f_["k"], f_["pred"])) / f_["pred"])]
        o["ii"] = t2
        # ---------- propozycja
        sy = v2["syst_rel"]
        km, ekm = r["mirror_npo0"]["ratio_to_nominal_npo0"], r["mirror_npo0"]["ratio_to_nominal_npo0_err_corr"]
        far = t2["hasD"]["by_margin"][2]
        down = math.hypot(km - 1, ekm)
        up = math.hypot(far["R"] - 1, far["err"])
        rest = math.sqrt(sy["stat"] ** 2 + sy["pool"] ** 2 + sy["half_diff_methods"] ** 2)
        f = v2["estimate"]["f_acc"]
        o["proposal"] = {"f_acc": f, "v2": {"stat": sy["stat"], "pool": sy["pool"], "extrap": sy["extrapolation"], "half": sy["half_diff_methods"], "total": sy["total_rel"], "f_err": f * sy["total_rel"]},
                         "kappa_mirror": km, "kappa_mirror_err": ekm, "kappa_upper95": km + 1.645 * ekm, "R_far_hasD": far["R"], "R_far_hasD_err": far["err"],
                         "R_far_upper95": far["R"] + 1.645 * far["err"],
                         "xical_syst_down": down, "xical_syst_up": up, "xical_syst_sym": max(down, up), "rest": rest,
                         "total_down": math.hypot(rest, down), "total_up": math.hypot(rest, up), "total_sym": math.hypot(rest, max(down, up)),
                         "f_err_down": f * math.hypot(rest, down), "f_err_up": f * math.hypot(rest, up), "f_err_sym": f * math.hypot(rest, max(down, up)),
                         "npo0_shift": r["nominal_npo0"]["ratio"] - 1, "npo0_shift_err": r["nominal_npo0"]["ratio_err_corr"]}
        out[str(r["run"])] = o
    json.dump(out, open(os.path.join(HERE, "jump_summary.json"), "w"), indent=1)
    return out


if __name__ == "__main__":
    o = derive()
    for run, x in o.items():
        print("=====", run)
        for lab in ("hasD", "hasD_npo0", "all"):
            pr = x[f"i_profile_{lab}"]
            print(lab, "R bins", [None if b["R"] is None else round(b["R"], 3) for b in pr["bins"]], "\n   err", [None if b["err"] is None else round(b["err"], 3) for b in pr["bins"]],
                  "\n   low", pr["low"], "high", pr["high"], "low/high", np.round(pr["low_over_high"], 3))
        for lab in ("hasD", "all"):
            print("symmetry", lab, json.dumps(x[f"ib_symmetry_{lab}"], default=float)[:900])
        for var in ("nsel", "mu", "lb"):
            print("dep", var, "hasD low/high", [None if a[0] is None else (round(a[0], 2), round(a[1], 2)) for a in x["iii_dep"][var]["hasD_low_over_high"]])
        for lab in ("hasD", "all"):
            t = x["ii"][lab]
            print("ii", lab, "by margin", [(b["k"], round(b["pred"], 1), round(b["R"], 3), round(b["err"], 3)) for b in t["by_margin"]], "near-far", np.round(t["near_minus_far"], 3))
            print("   by opp×margin", [[(b["k"], round(b["R"], 2), round(b["err"], 2)) for b in row] for row in t["by_opp_margin"]])
        print("proposal", json.dumps(x["proposal"], indent=0, default=float).replace("\n", " "))
