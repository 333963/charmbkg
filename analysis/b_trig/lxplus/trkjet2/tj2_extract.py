#!/usr/bin/env python3
"""TRKJET2 ekstrakcja (projekt B, 17.11.2026): dżety torowe i otoczenie D0 dla par D0←D* z B_trig — rozszerzenie TRKJET o trzy testy
z literatury (RECON_lit_B2 T2, T3, T8): stożki prostopadłe/losowy (piedestał), dżety w obszarach Δφ (podpróbka „tryger po przeciwnej
stronie”), ξ z torów + ξ_cal z ROOT (zastępnik energii neutralnej), wszystkie bity trigTBP i trigPassed MBTS/mb/RD0/dżetowe.
Kolumny TRKJET są zachowane bez zmian (TRKJET2 = nadzbiór). φ stożka losowego: generator z ziarnem (numer pliku, próbka) → powtarzalne.

Pary NIE są budowane od nowa: bierzemy bt_<próbka>.npz z etapu b_trig (ta sama selekcja i ten sam score H co „ext”, co do wiersza),
a dla każdej pary: wiersz h5 v2.0 (row → event_entry, córki D0; row_ds → wolny pion) → wpis drzewa `Tracks` pliku ROOT v2.0 o tym samym
numerze _NNNNNN_ → tj2_kernel (anty-kT R = 0,4 własna implementacja + stożki + półkule). Kontrola złączenia: eventNumber z CharmAnalysis
== event pary (każda niezgodność = błąd), pT wolnego pionu z h5 == pis_pt z bt.
Wyjście npz: kolumny pary z bt (KEEP + e_tp_*/e_ps wybrane) + bt_idx (indeks w bt npz) + tj2_kernel.NAMES + x_tp_/x_tbp_ (bity dżetowe z ROOT,
dane) + entry, uid + meta (json). Pamięć: h5 i ROOT czytane porcjami (--h5-chunk wierszy, --step wpisów).
Użycie: tj2_extract.py --sample S --bt bt_S.npz --h5dir DIR --rawdir DIR --part k --nparts N --out OUT.npz [--max-files n]
"""
import argparse
import glob
import json
import os
import re
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tj2_kernel as K  # noqa: E402

KEEP = ["run", "event", "file_idx", "row", "row_ds", "m", "pt", "eta", "phi", "score", "fold", "truth", "truth_ds", "weight", "m_swap",
        "dm", "pt_ds", "eta_ds", "pis_pt", "ptcone40", "ptcone100", "Lxy",
        "e_nSelBL", "e_lead_pt", "e_actual_mu", "e_grl", "e_npu", "e_lb", "e_ntrk500"]
KEEP_PS = ["e_ps_HLT_noalg_L1J12", "e_ps_HLT_j20_L1MBTS_2", "e_ps_HLT_mb_sptrk_pt4_L1MBTS_2", "e_ps_HLT_mb_mbts_L1MBTS_1"]
XTRIG = re.compile(r"^EventInfo_trig(Passed|TBP)_(.*?)(_NOSYS)?$")
XSEL = re.compile(r"(J\d+|_j\d+|MBTS|_mb_|RD0)")      # trigPassed: dżetowe + MBTS/mb/RD0; trigTBP: wszystkie
XI = re.compile(r"^EventInfo_(xi_cal_.*?)(_NOSYS)?$")
SEED = {"data_428770": 1, "data_427929": 2, "fd_nd": 3, "mb_nd": 4}
H5F = ["event_entry", "event_uid", "event_number", "Dmeson_pt", "Dmeson_ptcone40"] + \
      [f"Dmeson_daughter{i}_{v}" for i in (0, 1, 2) for v in ("pt", "eta", "phi")]


def h5_rows(path, rows, fields, chunk):
    """Pola `fields` dla posortowanych unikalnych wierszy `rows` (czytane porcjami)."""
    import h5py
    out = {}
    with h5py.File(path, "r") as f:
        ds = f["events"]
        fields = [c for c in fields if c in ds.dtype.names]
        n = ds.shape[0]
        if len(rows) and rows[-1] >= n:
            raise RuntimeError(f"{path}: wiersz {rows[-1]} ≥ {n}")
        parts = {c: [] for c in fields}
        for a in range(0, n, chunk):
            b = min(n, a + chunk)
            i0, i1 = np.searchsorted(rows, [a, b])
            if i1 == i0:
                continue
            x = ds.fields(fields)[a:b]
            loc = rows[i0:i1] - a
            for c in fields:
                parts[c].append(x[c][loc])
        for c in fields:
            out[c] = np.concatenate(parts[c]) if parts[c] else np.zeros(0)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True); ap.add_argument("--bt", required=True)
    ap.add_argument("--h5dir", required=True); ap.add_argument("--rawdir", required=True)
    ap.add_argument("--part", type=int, default=0); ap.add_argument("--nparts", type=int, default=1)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-files", type=int, default=0, help="test: tylko pierwsze n plików części")
    ap.add_argument("--h5-chunk", type=int, default=200_000); ap.add_argument("--step", type=int, default=20_000)
    a = ap.parse_args(argv)
    import awkward as ak
    import uproot

    t0 = time.time()
    z = np.load(a.bt, allow_pickle=False)
    meta = json.loads(str(z["meta"]))
    fmap = {}
    for g in meta["groups"]:
        for i, fn in enumerate(g["files"]):
            fmap[g["group"] * 1000 + i] = fn
    files_all = sorted(fmap.items(), key=lambda kv: kv[1])
    mine = files_all[a.part::a.nparts]
    if a.max_files:
        mine = mine[:a.max_files]
    file_idx = z["file_idx"]
    cols = [c for c in KEEP + KEEP_PS if c in z.files] + [c for c in z.files if c.startswith("e_tp_")]
    sel_all = np.nonzero(np.isin(file_idx, [k for k, _ in mine]))[0]
    P = {c: z[c][sel_all] for c in cols}
    del z
    npair = len(sel_all)
    print(f"[{a.sample}] część {a.part}/{a.nparts}: {len(mine)} z {len(files_all)} plików, par {npair}, numba={K.HAVE_NUMBA}, "
          f"wczytanie bt {time.time() - t0:.0f}s", flush=True)
    OUT = np.zeros((npair, K.NOUT), np.float32)
    entry = np.full(npair, -1, np.int64); uid = np.zeros(npair, np.uint64)
    h5cone = np.full(npair, np.nan, np.float32); same_entry = np.zeros(npair, np.int8)
    X = {}; XIV = {}
    stat = {"files": [], "pairs_per_file": {}, "sec_per_file": {}, "n_event_mismatch": 0, "n_pis_mismatch": 0, "n_d0pt_mismatch": 0}
    for fidx, fn in mine:
        tf = time.time()
        loc = np.nonzero(P["file_idx"] == fidx)[0]
        stat["files"].append(fn); stat["pairs_per_file"][fn] = int(len(loc))
        m = re.match(r"(.+)_(\d{6})_none\.h5$", fn)
        if not m:
            raise RuntimeError(f"nazwa h5 {fn}")
        tag = m.group(2)
        if len(loc) == 0:
            stat["sec_per_file"][fn] = 0.0
            continue
        rf = [f for f in glob.glob(os.path.join(a.rawdir, f"*._{tag}.*.root")) if "dpt3to150" in f]
        if len(rf) != 1:
            raise RuntimeError(f"{a.sample} {tag}: plików ROOT dpt3to150 {len(rf)} (oczekiwany 1)")
        r0 = P["row"][loc].astype(np.int64); r1 = P["row_ds"][loc].astype(np.int64)
        rows = np.unique(np.concatenate([r0, r1]))
        H = h5_rows(os.path.join(a.h5dir, fn), rows, H5F, a.h5_chunk)
        i0 = np.searchsorted(rows, r0); i1 = np.searchsorted(rows, r1)
        u = 1e-3 if np.nanmedian(H["Dmeson_daughter0_pt"]) > 100 else 1.0          # MeV → GeV, gdyby h5 miał MeV
        up = 1e-3 if np.nanmedian(H["Dmeson_pt"]) > 100 else 1.0
        ent = H["event_entry"][i0].astype(np.int64)
        same = H["event_entry"][i1].astype(np.int64) == ent        # MC: (run, event) nie jest kluczem → b_trig mógł skleić D0 i D* z różnych zdarzeń
        same_entry[loc] = same
        if a.sample.startswith("data_") and not same.all():
            raise RuntimeError(f"{fn}: D0 i D* pary w różnych wpisach drzewa (dane)")
        if np.any(H["event_number"][i0].astype(np.int64) != P["event"][loc].astype(np.int64)):
            raise RuntimeError(f"{fn}: event_number h5 ≠ event pary (zły plik/wiersz)")
        entry[loc] = ent; uid[loc] = H["event_uid"][i0].astype(np.uint64)
        h5cone[loc] = H["Dmeson_ptcone40"][i0] * up
        stat["n_d0pt_mismatch"] += int(np.sum(np.abs(H["Dmeson_pt"][i0] * up / P["pt"][loc] - 1) > 1e-3))
        pis = H["Dmeson_daughter2_pt"][i1] * u
        stat["n_pis_mismatch"] += int(np.sum(np.abs(pis / np.maximum(P["pis_pt"][loc], 1e-6) - 1) > 1e-3))
        DD = np.column_stack([P["pt"][loc], P["eta"][loc], P["phi"][loc], P["m"][loc],
                              H["Dmeson_daughter0_pt"][i0] * u, H["Dmeson_daughter0_eta"][i0], H["Dmeson_daughter0_phi"][i0],
                              H["Dmeson_daughter1_pt"][i0] * u, H["Dmeson_daughter1_eta"][i0], H["Dmeson_daughter1_phi"][i0],
                              pis, H["Dmeson_daughter2_eta"][i1], H["Dmeson_daughter2_phi"][i1]]).astype(np.float64)
        rng = np.random.default_rng([20261117, SEED.get(a.sample, 9), int(tag)])           # powtarzalne: zależy tylko od pliku i kolejności par w bt
        DD = np.hstack([DD, rng.uniform(-np.pi, np.pi, (len(loc), K.NRND))])
        with uproot.open(rf[0]) as R:
            T, C = R["Tracks"], R["CharmAnalysis"]
            if T.num_entries != C.num_entries:
                raise RuntimeError(f"{rf[0]}: Tracks {T.num_entries} ≠ CharmAnalysis {C.num_entries}")
            if ent.max() >= T.num_entries:
                raise RuntimeError(f"{rf[0]}: wpis {ent.max()} ≥ {T.num_entries}")
            xb = {}
            for k in C.keys():
                mm = XTRIG.match(k)
                if mm and (mm.group(1) == "TBP" or XSEL.search(mm.group(2))):
                    xb[("x_tp_" if mm.group(1) == "Passed" else "x_tbp_") + mm.group(2)] = k
            for k in xb:
                X.setdefault(k, np.full(npair, -1, np.int8))
            xib = {XI.match(k).group(1): k for k in C.keys() if XI.match(k)}
            for k in xib:
                XIV.setdefault(k, np.full(npair, np.nan, np.float32))
            order = np.argsort(ent, kind="stable")
            ents = ent[order]
            ue = np.unique(ents)
            for s in range(0, T.num_entries, a.step):
                e = min(T.num_entries, s + a.step)
                j0, j1 = np.searchsorted(ue, [s, e])
                if j1 == j0:
                    continue
                le = ue[j0:j1] - s
                tr = T.arrays(["Tracks_pt_NOSYS", "Tracks_eta_NOSYS", "Tracks_phi_NOSYS"], entry_start=s, entry_stop=e)[le]
                cnt = ak.to_numpy(ak.num(tr["Tracks_pt_NOSYS"]))
                off = np.concatenate([[0], np.cumsum(cnt)]).astype(np.int64)
                tpt = ak.to_numpy(ak.flatten(tr["Tracks_pt_NOSYS"])).astype(np.float64) * 1e-3
                teta = ak.to_numpy(ak.flatten(tr["Tracks_eta_NOSYS"])).astype(np.float64)
                tphi = ak.to_numpy(ak.flatten(tr["Tracks_phi_NOSYS"])).astype(np.float64)
                del tr
                p0, p1 = np.searchsorted(ents, [s, e])
                pi = order[p0:p1]                                   # pary (indeksy w loc) z tej porcji
                pev = np.searchsorted(ue[j0:j1], ent[pi])
                OUT[loc[pi]] = K.run(tpt, teta, tphi, off, pev, DD[pi])
                cv = C.arrays(["EventInfo_eventNumber"] + list(xb.values()) + list(xib.values()), entry_start=s, entry_stop=e, library="np")
                evn = cv["EventInfo_eventNumber"][le][pev].astype(np.int64)
                stat["n_event_mismatch"] += int(np.sum(evn != P["event"][loc[pi]].astype(np.int64)))
                for k, br in xb.items():
                    X[k][loc[pi]] = cv[br][le][pev].astype(np.int8)
                for k, br in xib.items():
                    XIV[k][loc[pi]] = cv[br][le][pev].astype(np.float32)
        stat["sec_per_file"][fn] = round(time.time() - tf, 1)
        print(f"[{a.sample}] {fn}: par {len(loc)}, zdarzeń {len(ue)}, {time.time() - tf:.1f}s", flush=True)
    if stat["n_event_mismatch"]:
        raise RuntimeError(f"eventNumber ROOT ≠ para: {stat['n_event_mismatch']} (złe wyrównanie Tracks/CharmAnalysis/h5)")
    out = dict(P)
    out["bt_idx"] = sel_all.astype(np.int64); out["entry"] = entry; out["uid"] = uid; out["h5_ptcone40"] = h5cone
    out["same_entry"] = same_entry          # 0 = D0 i D* z różnych wpisów drzewa (fałszywa para MC; otoczenie liczone we wpisie D0)
    for i, nme in enumerate(K.NAMES):
        out[nme] = OUT[:, i]
    out.update(X); out.update(XIV)
    # zgodność własnego ptcone z gałęzią (pary z obiema córkami dopasowanymi)
    val = {}
    if npair:
        ok = (OUT[:, K.IX["nmatch"]] + OUT[:, K.IX["nd_lo"]] == 2) & (same_entry == 1)      # komplet córek wycięty
        ref = P["ptcone40"]
        for v in ("v_s40_p2", "v_s40_p2x", "v_s40_p5", "v_s40_p5x", "c_s04"):
            d = OUT[ok, K.IX[v]] - ref[ok]
            val[v] = {"frac_abs_lt_1MeV": float(np.mean(np.abs(d) < 1e-3)), "frac_abs_lt_50MeV": float(np.mean(np.abs(d) < 0.05)),
                      "median": float(np.median(d)), "rms": float(np.sqrt(np.mean(d ** 2)))}
        d = OUT[ok, K.IX["v_s100_p2"]] - P["ptcone100"][ok] if "ptcone100" in P else np.zeros(1)
        val["v_s100_p2_vs_ptcone100"] = {"frac_abs_lt_1MeV": float(np.mean(np.abs(d) < 1e-3)), "median": float(np.median(d))}
        d = OUT[ok, K.IX["v_s100_p5"]] - P["ptcone100"][ok] if "ptcone100" in P else np.zeros(1)
        val["v_s100_p5_vs_ptcone100"] = {"frac_abs_lt_1MeV": float(np.mean(np.abs(d) < 1e-3)), "median": float(np.median(d))}
        val["frac_same_entry"] = float(np.mean(same_entry == 1))
        val["frac_daughters_ok"] = float(np.mean(OUT[:, K.IX["nmatch"]] + OUT[:, K.IX["nd_lo"]] == 2)); val["frac_pis_match"] = float(np.mean(OUT[:, K.IX["pis_match"]] == 1))
        val["frac_pis_pt_gt_0p2"] = float(np.mean(P["pis_pt"] > 0.2))
    stat.update({"sample": a.sample, "part": a.part, "nparts": a.nparts, "n_files_total": len(files_all), "n_pairs": int(npair),
                 "names": K.NAMES, "xtrig": sorted(X), "xi": sorted(XIV), "nrnd": K.NRND, "rnd_min": K.RND_MIN,
                 "frac_rnd_ok": float(np.mean(OUT[:, K.IX["r_ok"]] == 1)) if npair else None, "validation": val, "numba": K.HAVE_NUMBA, "seconds": round(time.time() - t0, 1),
                 "jet": {"algo": "anti-kT (własna implementacja, schemat E)", "R": K.R_JET, "pt_min": K.PT_MIN, "eta_max": K.ETA_MAX,
                         "match_dr": K.MATCH_DR, "match_rel": K.MATCH_REL},
                 "thresholds": meta.get("thresholds"), "m0": meta.get("m0")})
    out["meta"] = np.array(json.dumps(stat))
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    np.savez_compressed(a.out, **out)
    print(f"[{a.sample}] zapisano {a.out}: par {npair}, {os.path.getsize(a.out) / 1e6:.1f} MB, {time.time() - t0:.0f}s; walidacja: "
          f"{json.dumps(val)}; niezgodności pis/pt: {stat['n_pis_mismatch']}/{stat['n_d0pt_mismatch']}", flush=True)


if __name__ == "__main__":
    sys.exit(main())
