#!/usr/bin/env python3
"""DNET etap 6 — studium maski wolnego pionu (MC + pasma danych; samo numpy, bez modeli). Użycie: dn6_study.py --data DIR --env DIR --out JSON [--samples …]
(1) wzór: Δm tokenu prawdziwego π_s (tor o pT zgodnym z p_pis_pt w 1 %, jak dn4_pis) wobec p_dm tabeli; (2) dla siatki DM_MAX × znak:
udział usuniętych wolnych pionów (próbka tagu z prawdą D*), strata torów otoczenia w ΔR < 0,4 dla sygnału inkluzywnego / tła MC / pasm danych."""
import argparse, glob, json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
_Q = os.path.dirname(HERE)
for p in (os.environ.get("DN1") or os.path.join(_Q, "dnet1"), HERE):
    sys.path.insert(0, p)
import dn_feat as DF  # noqa: E402
import dn6_mask as MK  # noqa: E402

SP = "D0"
GRID = [0.148, 0.150, 0.152, 0.155, 0.160]


def load(ddir, edir, s, pre, cols, nmax=0):
    z = np.load(os.path.join(ddir, f"dn_{s}.npz"), allow_pickle=False)
    k = f"{SP}__{pre}"
    ix, E = [], []
    for f in sorted(glob.glob(os.path.join(edir, f"dn2_{s}_*.npz"))):
        e = np.load(f, allow_pickle=False)
        if k + "idx" in e.files:
            ix.append(e[k + "idx"]); E.append(e[k + "E"])
    ix, E = np.concatenate(ix), np.concatenate(E)
    o = np.argsort(ix); ix, E = ix[o], E[o]
    if nmax and len(ix) > nmax:
        st = len(ix) // nmax + 1; ix, E = ix[::st], E[::st]
    T = {c: z[k + c][ix] for c in cols if k + c in z.files}
    T["E"] = E
    return T


def loss(T, sel, tag):
    """strata torów e04 (ΔR < 0,4, pierwszych 12 tokenów) na wierszach sel, dla siatki."""
    E = T["E"][sel][:, :12]
    Dd, pt, eta, phi, m = (T[c][sel] for c in ("Dd", "pt", "eta", "phi", "m"))
    in04 = (E[:, :, 0] > 0) & (np.hypot(E[:, :, 1].astype(np.float64), E[:, :, 2].astype(np.float64)) < 0.4)
    o = {"n": int(sel.sum()), "tracks_per_cand_e04": float(in04.sum(1).mean())}
    for sg in (1, 0):
        for d in GRID:
            v = MK.veto(E, Dd, pt, eta, phi, m, d, sg) & in04
            o[f"dm<{d}|sign{sg}"] = {"frac_tracks_removed": float(v.sum() / max(in04.sum(), 1)), "frac_cands_touched": float(v.any(1).mean()),
                                     "removed_per_cand": float(v.sum(1).mean())}
    x = o[f"dm<{MK.DM_MAX}|sign{MK.SIGN}"] if f"dm<{MK.DM_MAX}|sign{MK.SIGN}" in o else None
    print(f"{tag}: n {o['n']}, torów e04/kandydat {o['tracks_per_cand_e04']:.2f}" + (f"; nominalnie usunięte {x['frac_tracks_removed']:.4f} torów, "
          f"{x['frac_cands_touched']:.4f} kandydatów" if x else ""), flush=True)
    return o


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True); ap.add_argument("--env", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--samples", default="fd_nd,mb_nd,data_427929,data_428770"); ap.add_argument("--nmax", type=int, default=400000)
    a = ap.parse_args()
    m0, win = DF.cfg["species"][SP]["m0"], DF.cfg["species"][SP]["signal_window"]
    far = DF.FAR[SP]
    res = {"DM_MAX": MK.DM_MAX, "SIGN": MK.SIGN, "grid": GRID}
    C = ["Dd", "m", "pt", "eta", "phi", "truth", "swap", "ok", "m_swap"]
    smp = a.samples.split(",")
    if "fd_nd" in smp:
        P = load(a.data, a.env, "fd_nd", "p_", C + ["pis_pt", "truth_ds", "dm"])
        sg = np.isin(P["truth"], DF.SIGNAL_CATS) & (P["swap"] == 0) & (np.abs(P["m"] - m0) < win) & np.isin(P["truth_ds"], DF.SIGNAL_CATS)
        Ept = P["E"][:, :, 0].astype(np.float64); pis = P["pis_pt"].astype(np.float64)
        rel = np.abs(Ept - pis[:, None]) / np.clip(pis[:, None], 1e-3, None); rel[Ept <= 0] = np.inf
        j = np.argmin(rel, 1); r = np.arange(len(j)); found = rel[r, j] < 0.01
        dmt = MK.token_dm(P["E"], P["pt"], P["eta"], P["phi"], P["m"])[r, j]
        qrel = P["E"][r, j, 3].astype(np.float64) * P["Dd"][:, 0, 7].astype(np.float64)
        dr = np.hypot(P["E"][r, j, 1].astype(np.float64), P["E"][r, j, 2].astype(np.float64))
        w = sg & found
        d = (dmt - P["dm"].astype(np.float64))[w] * 1e3
        q = lambda x, p: [float(v) for v in np.quantile(x, p)]  # noqa: E731
        res["formula"] = {"n_tag": int(sg.sum()), "found": float(found[sg].mean()), "dm_token_minus_table_MeV": {"mean": float(d.mean()), "rms": float(d.std()),
                          "q": q(d, [0.001, 0.01, 0.5, 0.99, 0.999])}, "dm_token_MeV_q": q(dmt[w] * 1e3, [0.001, 0.01, 0.5, 0.99, 0.995, 0.999]),
                          "dm_table_MeV_q": q(P["dm"][sg].astype(np.float64) * 1e3, [0.001, 0.01, 0.5, 0.99, 0.995, 0.999]),
                          "q_rel_minus1": float((qrel[w] < 0).mean()), "dr_lt04": float((dr[w] < 0.4).mean()), "slot_lt12": float((j[w] < 12).mean())}
        print("wzór:", json.dumps(res["formula"]), flush=True)
        seen = found & (dr < 0.4) & (j < 12)
        res["pis_removed"] = {}
        for sgn in (1, 0):
            for dmx in GRID:
                v = MK.veto(P["E"], P["Dd"], P["pt"], P["eta"], P["phi"], P["m"], dmx, sgn)
                rm = v[r, j]
                # po masce: czy w e04 został JAKIKOLWIEK tor o Δm < 0,160 dowolnego znaku (resztka)
                res["pis_removed"][f"dm<{dmx}|sign{sgn}"] = {"of_found": float(rm[sg & found].mean()), "of_seen_e04": float(rm[sg & seen].mean()),
                                                              "of_all_tag": float((rm & found)[sg].mean()), "seen_e04_left_of_all_tag": float((seen & ~rm)[sg].mean())}
                print(f"π_s usunięty dm<{dmx} znak{sgn}: {res['pis_removed'][f'dm<{dmx}|sign{sgn}']}", flush=True)
        res["pis_seen_e04"] = float(seen[sg].mean())
        for lo, hi in DF.BINS:
            b = sg & (P["pt"] >= lo) & (P["pt"] < hi)
            v = MK.veto(P["E"], P["Dd"], P["pt"], P["eta"], P["phi"], P["m"])[r, j]
            res.setdefault("pis_removed_bins", {})[f"{lo}-{hi}"] = {"n": int(b.sum()), "seen_e04": float(seen[b].mean()), "removed_of_seen": float(v[b & seen].mean()),
                                                                    "left_of_all": float((seen & ~v)[b].mean())}
        print("biny:", json.dumps(res["pis_removed_bins"]), flush=True)
        res["loss_tag"] = loss(P, sg, "próbka tagu (prawda D*)")
        del P
        A = load(a.data, a.env, "fd_nd", "a_", C)
        sa = np.isin(A["truth"], DF.SIGNAL_CATS) & (A["swap"] == 0) & (np.abs(A["m"] - m0) < win)
        res["loss_sig_incl"] = loss(A, sa, "sygnał inkluzywny fd_nd")
        del A
    if "mb_nd" in smp:
        A = load(a.data, a.env, "mb_nd", "a_", C, a.nmax)
        b = (A["truth"] >= DF.BKG_MIN) & (A["truth"] < 100) & A["ok"]
        res["loss_bkg_mc"] = loss(A, b, "tło MC mb_nd (cały zakres mas)")
        res["loss_bkg_mc_win"] = loss(A, b & (np.abs(A["m"] - m0) < win), "tło MC mb_nd (okno)")
        del A
    for run in (427929, 428770):
        if f"data_{run}" in smp:
            A = load(a.data, a.env, f"data_{run}", "a_", C, a.nmax)
            b = A["ok"] & (np.abs(A["m"] - m0) > far) & (np.abs(A["m_swap"] - m0) > DF.REFL_VETO)
            res[f"loss_sb{run}"] = loss(A, b, f"dalekie pasma {run}")
            del A
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    json.dump(res, open(a.out, "w"), indent=1)
    print("zapisano", a.out)


if __name__ == "__main__":
    sys.exit(main())
