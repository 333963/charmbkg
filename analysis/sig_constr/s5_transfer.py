#!/usr/bin/env python3
"""SIG_constr A(iii)/K(c3): przenoszenie D0←D* → inkluzywne D0 (tylko MC fd_nd, prawda): ε_iso(inkluzywne prawdziwe D0) vs ε_iso(D0←D*)
przeważone w pT; to samo dla score H (WP). Cięcie obszaru: |m − m0| < 40 MeV; izolacja w H06 (score > 0,6), score w PRE.
Uwaga: gałąź ptcone40 zawiera wolny pion D* (gdy leży w stożku) — wariant „−π_s” odejmuje pT pionu od stożka par (górne ograniczenie efektu)."""
from sc_common import *

z = np.load(os.path.join(CW, "B_trig", "runs", "lxplus", "bt_fd_nd.npz"))
A = {k[4:]: z[k] for k in z.files if k.startswith("all_")}
M = load("fd_nd", tj=False)
thr = json.load(open(os.path.join(CW, "B_trig", "runs", "lxplus", "out", "btrig_sf.json")))["thresholds"]
FINE = np.array([3, 3.25, 3.5, 3.75, 4, 4.5, 5, 5.5, 6, 7, 8, 9, 10, 12, 15, 20, 30.0])
rows = []


def cmp(name, passA, passP, selA, selP, b, lo, hi):
    a = selA & (A["pt"] >= lo) & (A["pt"] < hi); p = selP & (M["pt"] >= lo) & (M["pt"] < hi)
    ha = np.histogram(A["pt"][a], FINE)[0].astype(float); hp = np.histogram(M["pt"][p], FINE)[0].astype(float)
    w = (ha / np.maximum(hp, 1))[np.clip(np.digitize(M["pt"][p], FINE) - 1, 0, len(FINE) - 2)]
    eA = passA[a].mean(); sA = np.sqrt(eA * (1 - eA) / a.sum())
    eP0 = passP[p].mean()
    eP = np.average(passP[p], weights=w); sP = np.sqrt(np.sum(w ** 2 * (passP[p] - eP) ** 2)) / w.sum()
    d = eP - eA; s = np.hypot(sA, sP)
    rows.append({"test": name, "bin": b, "n_incl": int(a.sum()), "n_tag": int(p.sum()), "eps_incl": float(eA), "eps_tag_raw": float(eP0), "eps_tag_ptrw": float(eP),
                 "diff": float(d), "err": float(s), "rel": float(d / eA), "pull": float(d / s)})
    print(f"{name:26s} {b:5s} N_inkl {a.sum():6d} N_tag {p.sum():6d} ε_inkl {eA:.4f} ε_tag {eP0:.4f} → pT-rw {eP:.4f}  Δ {d:+.4f} ± {s:.4f} ({100*d/eA:+.1f} %, {d/s:+.1f}σ)")


sigA = np.isin(A["truth"], [0, 1]) & (np.abs(A["m"] - C.M0_D0) < MWIN)
sigP = np.isin(M["truth"], [0, 1]) & np.isin(M["truth_ds"], [0, 1]) & (np.abs(M["m"] - C.M0_D0) < MWIN) & (M["dm"] > 0.1435) & (M["dm"] < 0.1475)
print("inkluzywne prawdziwe D0:", int(sigA.sum()), " D0←D*:", int(sigP.sum()))
for b, (lo, hi) in list(PTB.items()) + [("3-30", (3.0, 30.0))]:
    for cut in (1.0, 2.0, 4.0):
        cmp(f"ptcone40<{cut:g} (H06)", A["ptcone40"] < cut, M["ptcone40"] < cut, sigA & (A["score"] > 0.6), sigP & (M["score"] > 0.6), b, lo, hi)
        cmp(f"ptcone40−π_s<{cut:g} (H06)", A["ptcone40"] < cut, (M["ptcone40"] - M["pis_pt"]) < cut, sigA & (A["score"] > 0.6), sigP & (M["score"] > 0.6), b, lo, hi)
    cmp("ptcone40/pT<0.5 (H06)", A["ptcone40"] / A["pt"] < 0.5, M["ptcone40"] / M["pt"] < 0.5, sigA & (A["score"] > 0.6), sigP & (M["score"] > 0.6), b, lo, hi)
    if b in thr:
        for wp in ("0.95", "0.80", "0.50"):
            cmp(f"score>WP{wp} (PRE)", A["score"] > thr[b][wp], M["score"] > thr[b][wp], sigA, sigP, b, lo, hi)
jdump({"rows": rows, "n_incl": int(sigA.sum()), "n_tag": int(sigP.sum())}, "s5_transfer.json")
