#!/usr/bin/env python3
"""ISO_origin krok 3: ε_B tła pod pikiem — jak ważone jest MC per klasa w btrig_analyze epsb i ile z R = dane/MC to wagi.
MC tła (mb_nd, truth ≥ 8, |m−m0| < 55 MeV): ε_B z wagami (1) = 1, (2) emulacja A (class_prob, jak w nocy), (3) gęstość B
(p_dane,klasa/p_MC w (nSelBL, lead_pt) na kandydatach D0 danych tej klasy). Dane: pol2 z nocnego btrig_epsb.json."""
from iso_common import *
import btrig_emul as EM

O = os.path.join(R, "out")
J = json.load(open(os.path.join(O, "btrig_epsb.json")))
thr = J["thresholds"]
em = EM.load_map(os.path.join(O, "effmap_427929.json"))
z = np.load(os.path.join(R, "bt_mb_nd.npz"))
Mb = {k: z[k] for k in ("all_truth", "all_pt", "all_m", "all_score", "all_e_nSelBL", "all_e_lead_pt")}
dat = {}
for run in (428770, 427929):
    zd = np.load(os.path.join(R, f"bt_data_{run}.npz"))
    D = {k: zd[k] for k in ("all_pt", "all_m", "all_e_nSelBL", "all_e_lead_pt")}
    bits = {k[6:]: zd[k] for k in zd.files if k.startswith("all_e_tp_") or k == "all_e_grl"}
    dat[run] = (D, C.classify(bits, run))
rows = []
for b, (lo, hi) in C.BINS.items():
    mt = (Mb["all_truth"] >= 8) & (Mb["all_pt"] >= lo) & (Mb["all_pt"] < hi) & (np.abs(Mb["all_m"] - C.M0_D0) < 0.055)
    ns, ld, sc = Mb["all_e_nSelBL"][mt], Mb["all_e_lead_pt"][mt], Mb["all_score"][mt]
    for run, c in ((428770, "all"), (428770, "J12only"), (428770, "MBTSbased"), (428770, "other"), (427929, "MB427929")):
        D, cl = dat[run]
        sd = cl[c] & (D["all_pt"] >= lo) & (D["all_pt"] < hi) & (np.abs(D["all_m"] - C.M0_D0) > 0.075)      # pasma boczne danych (tło)
        wA = EM.class_prob(em, ns, ld, c)
        wB, _ = EM.density_weights(ns, ld, D["all_e_nSelBL"][sd], D["all_e_lead_pt"][sd])
        for wp in ("0.80", "0.70", "0.50"):
            p = sc >= thr[b][wp]
            d = [r for r in J["rows"] if r["run"] == run and r["cls"] == c and r["bin"] == b and r["wp"] == wp][0]
            e1, eA, eB = float(p.mean()), float((wA * p).sum() / wA.sum()), float((wB * p).sum() / wB.sum())
            rows.append({"run": run, "cls": c, "bin": b, "wp": wp, "eps_data": d["pol2"], "mc_w1": e1, "mc_emulA": eA, "mc_densB": eB,
                         "mc_night": d["mc_bkg"], "R_night": d["R_data_mc"], "R_w1": d["pol2"] / e1, "R_emulA": d["pol2"] / eA, "R_densB": d["pol2"] / eB,
                         "neff_A": EM.n_eff(wA), "neff_B": EM.n_eff(wB), "n_mc": int(mt.sum()), "wA_const": bool(np.allclose(wA, wA[0]))})
json.dump(rows, open(os.path.join(HERE, "out", "s3.json"), "w"), indent=1)
print("| run | klasa | bin | WP | ε_B dane (pol2) | MC w=1 | MC emul. A | MC gęstość B | R noc | R(w=1) | R(emul A) | R(gęst. B) | n_eff A/B (z n_MC) |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for r in rows:
    print(f"| {r['run']} | {r['cls']} | {r['bin']} | {r['wp']} | {r['eps_data']:.4f} | {r['mc_w1']:.4f} | {r['mc_emulA']:.4f} | {r['mc_densB']:.4f} | {r['R_night']:.3f} | "
          f"{r['R_w1']:.3f} | {r['R_emulA']:.3f} | {r['R_densB']:.3f} | {r['neff_A']:.0f}/{r['neff_B']:.0f} ({r['n_mc']}) |")
