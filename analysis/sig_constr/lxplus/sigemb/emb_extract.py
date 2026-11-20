#!/usr/bin/env python3
"""SIGEMB ekstrakcja (projekt B, SIG_constr, 17.11.2026): EMBEDDING — kinematyka (pT, η) prawdziwego D0←D* z MC (emb_donors.npz) wstawiana
w zdarzenia-nośniki (dane per klasa triggera; mb_nd jako kontrola) i otoczenie liczone z PRAWDZIWYCH torów nośnika (drzewo `Tracks`, ROOT v2.0).
Nośniki = wszystkie wpisy drzewa (zdarzenia z ≥ 1 luźnym kandydatem D — tak zbudowane są ntuple), co `--prescale`-ty wg hasza eventNumber.
Warianty kierunku (kolumna `var`), jeden dawca na nośnik, ten sam dla wszystkich wariantów:
  0 rnd    — η dawcy, φ losowe
  1 perpJ  — η dawcy, φ = φ(wiodący dżet torowy; gdy brak — wiodący tor) ± π/2 (znak wg parzystości eventNumber)
  2 nearJ  — przy dżecie torowym: dżet losowany jednostajnie spośród dżetów pT > 2 GeV (maks. 6 wiodących), oś D0 w ΔR = r·(0…0,8) od osi dżetu
             (r płasko), azymut wokół osi losowy; wymagane |η_D0| < 2,2. Rozkład (ΔR, pT dżetu) NIE jest fizyczny — analiza przeważa do rozkładu
             zmierzonego dla prawdziwych D0 z tagu D* (TRKJET: b_drn*, b_ptn*, o_*), osobno per klasa i bin pT.
Dla każdego wiersza jądro TRKJET2 (tj2_kernel.one_pair, bez zmian) liczy stożki, dżety A (D0 jako pseudocząstka) / B / O, półkule, stożki prostopadłe.
Córki D0 nie są wycinane (nie ma ich w nośniku) — wektor córek = 0. Wyjście npz: wiersz = (nośnik, wariant).
Użycie: emb_extract.py --sample S --rawdir DIR --donors emb_donors.npz --kernel-dir ../trkjet2 --part k --nparts N --prescale P --out OUT.npz [--max-files n]"""
import argparse
import glob
import json
import math
import os
import re
import sys
import time

import numpy as np

M0 = 1.86484
SEED = {"data_428770": 11, "data_427929": 12, "fd_nd": 13, "mb_nd": 14}
CHAINS = ["HLT_j20_L1MBTS_2", "HLT_mb_sptrk_pt4_L1MBTS_2", "HLT_mb_sptrk_pt2_L1MBTS_2", "HLT_mb_mbts_L1MBTS_1", "HLT_noalg_mb_L1MBTS_1",
          "HLT_noalg_mb_L1MBTS_2", "HLT_mb_sptrk_L1MBTS_1", "HLT_noalg_L1J12", "HLT_noalg_L1AFP_A_OR_C", "HLT_mb_sptrk_pt4_L1AFP_A_OR_C",
          "HLT_noalg_L1RD0_FILLED", "HLT_mb_sptrk_pt6_L1MBTS_2", "HLT_mb_sptrk_pt8_L1MBTS_2"]
EXTRA = {"nSelBL": ["nSelBL"], "lead_pt": ["leadTrackPtBL"], "actual_mu": ["actual_mu", "actualInteractionsPerCrossing"], "grl": ["passGRL_MinBiasAFP"],
         "npu": ["nPileupVertices"], "lb": ["lumiBlock"], "ntrk500": ["nTracksPt500Eta25"], "run": ["runNumber"]}
KEEPN = ["a_jD_pt", "a_dr_jD", "a_jD_n", "z", "a_j1_pt", "b_n5", "b_j1_pt", "b_dr_j1", "b_j2_pt", "b_drn2", "b_ptn2", "b_drn3", "b_ptn3", "b_drn5", "b_ptn5",
         "b_drn8", "b_ptn8", "b_ptmax04", "b_ptmax08", "c_n04", "c_s04", "c_s10", "v_s40_p2", "v_s40_p5", "v_s100_p2", "h_s_away", "h_n_away", "h_s_tow",
         "t1_pt", "t1_dr", "ntrk05", "o_n5", "o_j1_pt", "o_dr_j1", "o_drn3", "o_ptn3", "o_drn5", "o_ptn5", "o_ptmax08", "p_s2_p", "p_s2_m",
         "b_tow_pt", "b_away_pt", "o_tow_pt", "o_away_pt", "sumpt_p2", "ntrk02"]
NVAR = 3
PT_STRATA = np.array([3.0, 4.0, 5.0, 6.5, 8.0, 12.0, 30.0])


def col(keys, names):
    s = set(keys)
    for n in names:
        for pre in ("EventInfo_", ""):
            for suf in ("_NOSYS", ""):
                if pre + n + suf in s:
                    return pre + n + suf
    return None


def make_placer(K):
    from numba import njit
    antikt = K.antikt

    @njit(cache=False)
    def place(tpt, teta, tphi, off, dpt, deta, U, par, DD, AUX, pt_min, eta_max, R, mpi):
        """Dla nośnika e: dżety torowe (jak kol. B TRKJET), 3 wiersze DD (rnd, perpJ, nearJ). U[e] = 8 liczb losowych U(0,1).
        AUX[row] = (ok, pT dżetu wybranego, r, liczba dżetów pT > 2, pT dżetu wiodącego, pT toru wiodącego)."""
        ne = len(off) - 1
        for e in range(ne):
            a, b = off[e], off[e + 1]
            n = 0
            for i in range(a, b):
                if tpt[i] > pt_min and abs(teta[i]) < eta_max:
                    n += 1
            sp = np.empty(n); se = np.empty(n); sph = np.empty(n); sm = np.full(n, mpi); fl = np.zeros(n, np.int64)
            j = 0; t1 = 0.0; t1phi = 0.0
            for i in range(a, b):
                if tpt[i] > pt_min and abs(teta[i]) < eta_max:
                    sp[j] = tpt[i]; se[j] = teta[i]; sph[j] = tphi[i]; j += 1
                    if tpt[i] > t1:
                        t1 = tpt[i]; t1phi = tphi[i]
            nj2 = 0; j1pt = 0.0; j1phi = t1phi
            jpt = np.zeros(1); jeta = np.zeros(1); jphi = np.zeros(1)
            if n > 0:
                jpt, jeta, jphi, _jn, _jfl = antikt(sp, se, sph, sm, fl, R)
                j1pt = jpt[0]; j1phi = jphi[0]
                for q in range(len(jpt)):
                    if jpt[q] > 2.0 and nj2 < 6:
                        nj2 += 1
            for v in range(3):
                r = 3 * e + v
                DD[r, 0] = dpt[e]; DD[r, 3] = 1.86484
                AUX[r, 3] = nj2; AUX[r, 4] = j1pt; AUX[r, 5] = t1
                if v == 0:
                    DD[r, 1] = deta[e]; DD[r, 2] = -math.pi + 2.0 * math.pi * U[e, 0]; AUX[r, 0] = 1.0
                elif v == 1:
                    sgn = 1.0 if par[e] == 0 else -1.0
                    ph = j1phi + sgn * 0.5 * math.pi
                    if ph > math.pi:
                        ph -= 2.0 * math.pi
                    if ph < -math.pi:
                        ph += 2.0 * math.pi
                    DD[r, 1] = deta[e]; DD[r, 2] = ph; AUX[r, 0] = 1.0 if n > 0 else 0.0
                else:
                    AUX[r, 0] = 0.0; DD[r, 1] = deta[e]; DD[r, 2] = 0.0
                    if nj2 > 0:
                        q = min(int(U[e, 1] * nj2), nj2 - 1)
                        rr = 0.8 * U[e, 2]
                        for t in range(4):
                            psi = 2.0 * math.pi * U[e, 3 + t]
                            et = jeta[q] + rr * math.cos(psi)
                            if abs(et) < 2.2:
                                ph = jphi[q] + rr * math.sin(psi)
                                if ph > math.pi:
                                    ph -= 2.0 * math.pi
                                if ph < -math.pi:
                                    ph += 2.0 * math.pi
                                DD[r, 1] = et; DD[r, 2] = ph; AUX[r, 0] = 1.0; AUX[r, 1] = jpt[q]; AUX[r, 2] = rr
                                break
    return place


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True); ap.add_argument("--rawdir", required=True); ap.add_argument("--donors", required=True)
    ap.add_argument("--kernel-dir", required=True)
    ap.add_argument("--part", type=int, default=0); ap.add_argument("--nparts", type=int, default=1); ap.add_argument("--prescale", type=int, default=1)
    ap.add_argument("--out", required=True); ap.add_argument("--max-files", type=int, default=0); ap.add_argument("--step", type=int, default=20_000)
    a = ap.parse_args(argv)
    sys.path.insert(0, a.kernel_dir)
    import tj2_kernel as K
    import awkward as ak
    import uproot
    place = make_placer(K)
    t0 = time.time()
    dn = np.load(a.donors)
    dpt_all, deta_all = dn["pt"].astype(np.float64), dn["eta"].astype(np.float64)
    strata = [np.nonzero((dpt_all >= lo) & (dpt_all < hi))[0] for lo, hi in zip(PT_STRATA[:-1], PT_STRATA[1:])]
    files = sorted(f for f in glob.glob(os.path.join(a.rawdir, "*.root")) if "dpt3to150" in f)
    mine = files[a.part::a.nparts]
    if a.max_files:
        mine = mine[:a.max_files]
    print(f"[{a.sample}] część {a.part}/{a.nparts}: {len(mine)} z {len(files)} plików, prescale {a.prescale}, dawców {len(dpt_all)}, numba={K.HAVE_NUMBA}", flush=True)
    ix = [K.IX[n] for n in KEEPN]
    parts = []; stat = {"files": [], "n_entries": 0, "n_carriers": 0, "sec_per_file": {}}
    for path in mine:
        tf = time.time()
        fn = os.path.basename(path)
        m = re.search(r"\._(\d{6})\.", fn)
        tag = int(m.group(1)) if m else len(stat["files"])
        rng = np.random.default_rng([20261117, SEED.get(a.sample, 9), tag])
        with uproot.open(path) as R:
            T, C = R["Tracks"], R["CharmAnalysis"]
            if T.num_entries != C.num_entries:
                raise RuntimeError(f"{fn}: Tracks {T.num_entries} ≠ CharmAnalysis {C.num_entries}")
            keys = C.keys()
            cmap = {"event": col(keys, ["eventNumber"])}
            for k, v in EXTRA.items():
                c = col(keys, v)
                if c:
                    cmap["e_" + k] = c
            for ch in CHAINS:
                c = col(keys, ["trigPassed_" + ch])
                if c:
                    cmap["e_tp_" + ch] = c
            for l1 in ("L1_J12", "L1_MBTS_1", "L1_MBTS_2"):
                c = col(keys, ["trigTBP_" + l1])
                if c:
                    cmap["x_tbp_" + l1] = c
            if cmap["event"] is None:
                raise RuntimeError(f"{fn}: brak eventNumber")
            stat["n_entries"] += int(T.num_entries)
            for s in range(0, T.num_entries, a.step):
                e = min(T.num_entries, s + a.step)
                cv = C.arrays(list(cmap.values()), entry_start=s, entry_stop=e, library="np")
                evn = cv[cmap["event"]].astype(np.int64)
                h = (evn * 2654435761 + tag * 40503) & 0xFFFFFFFF
                keep = ((h >> 7) % a.prescale) == 0
                le = np.nonzero(keep)[0]
                if len(le) == 0:
                    continue
                tr = T.arrays(["Tracks_pt_NOSYS", "Tracks_eta_NOSYS", "Tracks_phi_NOSYS"], entry_start=s, entry_stop=e)[le]
                cnt = ak.to_numpy(ak.num(tr["Tracks_pt_NOSYS"]))
                off = np.concatenate([[0], np.cumsum(cnt)]).astype(np.int64)
                tpt = ak.to_numpy(ak.flatten(tr["Tracks_pt_NOSYS"])).astype(np.float64) * 1e-3
                teta = ak.to_numpy(ak.flatten(tr["Tracks_eta_NOSYS"])).astype(np.float64)
                tphi = ak.to_numpy(ak.flatten(tr["Tracks_phi_NOSYS"])).astype(np.float64)
                del tr
                ne = len(le)
                st = rng.integers(0, len(strata), ne)
                di = np.array([strata[k][rng.integers(0, len(strata[k]))] for k in st], np.int64)
                U = rng.random((ne, 8))
                DD = np.zeros((NVAR * ne, K.ND)); AUX = np.zeros((NVAR * ne, 6))
                DD[:, 13:] = rng.uniform(-np.pi, np.pi, (NVAR * ne, K.NRND))
                place(tpt, teta, tphi, off, dpt_all[di], deta_all[di], U, (evn[le] % 2).astype(np.int64), DD, AUX, K.PT_MIN, K.ETA_MAX, K.R_JET, K.M_PI)
                pev = np.repeat(np.arange(ne), NVAR)
                O = K.run(tpt, teta, tphi, off, pev, DD)
                d = {"var": np.tile(np.arange(NVAR, dtype=np.int8), ne), "ok": AUX[:, 0].astype(np.int8), "pj_pt": AUX[:, 1].astype(np.float32),
                     "pj_r": AUX[:, 2].astype(np.float32), "nj2": AUX[:, 3].astype(np.int8), "j1_pt": AUX[:, 4].astype(np.float32), "trk1_pt": AUX[:, 5].astype(np.float32),
                     "pt": DD[:, 0].astype(np.float32), "eta": DD[:, 1].astype(np.float32), "phi": DD[:, 2].astype(np.float32),
                     "donor": np.repeat(di, NVAR).astype(np.int32), "file_tag": np.full(NVAR * ne, tag, np.int32), "entry": np.repeat(s + le, NVAR).astype(np.int64)}
                for k, c in cmap.items():
                    x = cv[c][le]
                    if k == "e_lead_pt":
                        x = np.where(np.nanmedian(x) > 100, x * 1e-3, x) if len(x) else x
                    d[k] = np.repeat(x, NVAR).astype(np.int64 if k == "event" else (np.int8 if k.startswith(("e_tp_", "x_tbp_")) or k == "e_grl" else np.float32))
                for n_, i_ in zip(KEEPN, ix):
                    d[n_] = O[:, i_]
                parts.append(d); stat["n_carriers"] += ne
        stat["files"].append(fn); stat["sec_per_file"][fn] = round(time.time() - tf, 1)
        print(f"[{a.sample}] {fn}: nośników łącznie {stat['n_carriers']}, {time.time() - tf:.1f}s", flush=True)
    keys = set(parts[0]) if parts else set()
    for p in parts:
        keys &= set(p)
    out = {k: np.concatenate([p[k] for p in parts]) for k in sorted(keys)}
    stat.update({"sample": a.sample, "part": a.part, "nparts": a.nparts, "prescale": a.prescale, "n_files_total": len(files), "n_rows": int(len(out.get("var", []))),
                 "names": KEEPN, "variants": {"0": "rnd", "1": "perpJ", "2": "nearJ"}, "seconds": round(time.time() - t0, 1),
                 "frac_ok": {str(v): float(np.mean(out["ok"][out["var"] == v])) for v in range(NVAR)} if parts else {}})
    out["meta"] = np.array(json.dumps(stat))
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    np.savez_compressed(a.out, **out)
    print(f"[{a.sample}] zapisano {a.out}: wierszy {stat['n_rows']}, nośników {stat['n_carriers']} z {stat['n_entries']} wpisów, {os.path.getsize(a.out) / 1e6:.1f} MB, "
          f"{time.time() - t0:.0f}s; ok per wariant {stat['frac_ok']}", flush=True)


if __name__ == "__main__":
    sys.exit(main())
