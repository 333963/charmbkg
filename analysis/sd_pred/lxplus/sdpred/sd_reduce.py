#!/usr/bin/env python3
"""sd_reduce.py — TSV z pythia_sd → histogramy (npz) per składowa × gatunek. Tylko numpy.
Użycie: sd_reduce.py part.tsv[.gz] part.tsv.json out.npz
Definicja (poziom cząstek): D = D0/D+/D*+ (z antycząstkami; D0 z D* liczone też jako D0; prompt + z b), 2 ≤ pT < 30 GeV, |η| < 2,2;
tag strony s ∈ {A (+z), C (−z)}: proton końcowy o najmniejszym ξ po stronie s z 0,030 ≤ ξ ≤ 0,100 (w oknie mieści się najwyżej jeden).
Każda strona osobno (D z tagiem A i C liczone dwa razy — jak region `tag` = A lub C; udział podwójnych: h `ntagAC`).
Zmienne „po stronie protonu": ηs = η·znak, ys = y·znak (znak = +1 dla A) → ηs > 0 = półkula protonu.
"""
import gzip, json, sys
import numpy as np

PTB = [2, 3, 5, 8, 30]
PTF = [2, 2.5, 3, 3.5, 4, 5, 6, 8, 10, 15, 30]
XIB = [0.030, 0.040, 0.050, 0.060, 0.075, 0.100]
EDGES = {
    "pt": PTB, "ptf": PTF, "xi": XIB,
    "ys": list(np.round(np.arange(-2.5, 2.51, 0.25), 2)), "etas": list(np.round(np.arange(-2.2, 2.21, 0.2), 2)),
    "t": list(np.round(np.arange(0, 2.001, 0.05), 3)) + [1e9],
    "nch": [0, 10, 20, 30, 45, 70, 1e9],
    "xiw": list(np.round(np.arange(0, 0.2001, 0.01), 3)),
    "d": [-1e9] + list(np.round(np.arange(-0.12, 0.3001, 0.005), 4)) + [1e9],      # ξ_central − ξ_p
    "r": list(np.round(np.arange(0, 3.001, 0.05), 3)) + [1e9],                      # ξ_central / ξ_p
    "xc": list(np.round(np.arange(0, 1.0001, 0.01), 3)) + [1e9],                    # ξ_central (bez wymogu protonu)
}
XILO, XIHI, DCUT = 0.030, 0.100, 0.02
SP = {421: "D0", 411: "Dp", 413: "Ds"}


def comp(code, hd):
    if code == 101: return "ND"
    if code in (103, 104): return "SD"
    if code == 105: return "DD"
    if code == 106: return "CD"
    return "HQ" + ("nd" if hd == 0 else "dif" if hd in (1, 2) else "dd")


def main(tsv, jsn, out):
    meta = json.load(open(jsn))
    ev = {}; rows = []
    op = gzip.open if tsv.endswith(".gz") else open
    with op(tsv, "rt") as f:
        for ln in f:
            c = ln.split("\t")
            if c[0] == "E":
                ev[c[1]] = (len(ev), comp(int(c[2]), int(c[3])), [float(x) for x in c[4:]])
            elif c[0] == "D":
                rows.append((c[1], int(c[2]), float(c[3]), float(c[4]), float(c[5]), int(c[7])))
    H = {}
    def fill(key, edges, x, w=None):
        h, _ = np.histogram(x, bins=edges, weights=w)
        H[key] = H.get(key, 0) + h
    def fill2(key, ex, ey, x, y):
        h, _, _ = np.histogram2d(x, y, bins=[ex, ey])
        H[key] = H.get(key, 0) + h.ravel()
    if rows:
        E = np.array([ev[r[0]][2] for r in rows])       # nch lead xiA tA xiC tC xcalA xcalC xtrkA xtrkC xtk1A xtk1C
        cmpn = np.array([ev[r[0]][1] for r in rows]); pdg = np.abs(np.array([r[1] for r in rows]))
        pt = np.array([r[2] for r in rows]); eta = np.array([r[3] for r in rows]); y = np.array([r[4] for r in rows]); fb = np.array([r[5] for r in rows])
        fid = (pt >= 2) & (pt < 30) & (np.abs(eta) < 2.2)
        for cn in np.unique(cmpn):
            for pid, sp in SP.items():
                m = fid & (cmpn == cn) & (pdg == pid)
                if not m.any(): continue
                k = f"{cn}|{sp}|"
                fill(k + "pt_incl", PTB, pt[m]); fill(k + "ptf_incl", PTF, pt[m]); fill(k + "ptb_incl", PTB, pt[m][fb[m] == 1])
                fill(k + "ys_incl", EDGES["ys"], y[m]); fill(k + "etas_incl", EDGES["etas"], eta[m]); fill(k + "nch_incl", EDGES["nch"], E[m, 0])
                wA = (E[:, 2] >= XILO) & (E[:, 2] <= XIHI); wC = (E[:, 4] >= XILO) & (E[:, 4] <= XIHI)
                fill(k + "pt_tagAC", PTB, pt[m & wA & wC])
                for sgn, ix, it, ic, itr, it1 in ((1, 2, 3, 6, 8, 10), (-1, 4, 5, 7, 9, 11)):
                    xi = E[:, ix]; tt = -E[:, it]; xc = E[:, ic]; xt = E[:, itr]; x1 = E[:, it1]
                    fill(k + "xc_incl", EDGES["xc"], xc[m]); fill(k + "xtrk_incl", EDGES["xc"], xt[m]); fill(k + "xtk1_incl", EDGES["xc"], x1[m])
                    fill(k + "xiw", EDGES["xiw"], xi[m])
                    fill2(k + "pt_xiw", PTB, EDGES["xiw"], pt[m], xi[m])
                    w = m & (xi >= XILO) & (xi <= XIHI)
                    if not w.any(): continue
                    fill(k + "pt_tag", PTB, pt[w]); fill(k + "ptf_tag", PTF, pt[w]); fill(k + "ptb_tag", PTB, pt[w][fb[w] == 1])
                    fill(k + "ys_tag", EDGES["ys"], sgn * y[w]); fill(k + "etas_tag", EDGES["etas"], sgn * eta[w])
                    fill(k + "xi_tag", XIB, xi[w]); fill(k + "t_tag", EDGES["t"], tt[w]); fill(k + "nch_tag", EDGES["nch"], E[w, 0])
                    fill(k + "pt_tag_t1", PTB, pt[w & (tt < 1.0)])
                    fill2(k + "pt_xi_tag", PTB, XIB, pt[w], xi[w])
                    fill(k + "pt_tag_cal", PTB, pt[w & (xc <= xi + DCUT)]); fill(k + "pt_tag_trk", PTB, pt[w & (xt <= xi + DCUT)])
                    fill(k + "xi_tag_cal", XIB, xi[w & (xc <= xi + DCUT)])
                    for nm, v in (("cal", xc), ("trk", xt), ("tk1", x1)):
                        fill(k + "d" + nm + "_tag", EDGES["d"], (v - xi)[w]); fill(k + "r" + nm + "_tag", EDGES["r"], (v / xi)[w])
                        fill2(k + "xi_d" + nm, XIB, EDGES["d"], xi[w], (v - xi)[w])
    np.savez_compressed(out, __meta__=np.array(json.dumps(meta)), __edges__=np.array(json.dumps(EDGES)),
                        **{k: np.asarray(v, dtype=np.float64) for k, v in H.items()})
    nt = sum(float(v.sum()) for k, v in H.items() if k.endswith("|pt_tag"))
    print(f"sd_reduce: {len(ev)} zdarzeń z D, {len(rows)} D, {nt:.0f} D fid. z tagiem (A+C) → {out}")


if __name__ == "__main__":
    main(*sys.argv[1:4])
