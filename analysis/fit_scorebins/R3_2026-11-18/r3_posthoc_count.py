"""R3 — kontrola POST HOC: wydajki sygnału per kategoria ze zliczenia w oknie |dm - mu| < 3 MeV minus tło z fitu nominalnego
(niezależne od kształtu sygnału), vs wydajki z fitu. Wykrywa niedopasowanie wspólnego kształtu piku w kategorii."""
import json, numpy as np, r3_common as R
D = json.load(open("out/data.json")); out = {}
for k, d in D.items():
    H, lb, S = np.array(d["H"]), np.array(d["lam_bkg"]), np.array(d["nom"]["S"]); fm = np.array(d["f_mc"])
    w = np.abs(R.CEN - d["nom"]["mu"]) < 3.0
    cnt = (H[:, w] - lb[:, w]).sum(1); err = np.sqrt(H[:, w].sum(1))
    f = cnt / cnt.sum()
    out[k] = {"S_count": cnt.tolist(), "S_count_err_poisson": err.tolist(), "S_fit": S.tolist(), "R_count": (f / fm).tolist(), "R_fit": (S / S.sum() / fm).tolist()}
    print(k, "zliczenie:", " ".join(f"{c:.0f}±{e:.0f}" for c, e in zip(cnt, err)), "| fit:", " ".join(f"{s:.0f}" for s in S),
          "| R_zlicz:", np.round(f / fm, 2), "| R_fit:", np.round(S / S.sum() / fm, 2))
json.dump(out, open("posthoc_count.json", "w"), indent=1)
