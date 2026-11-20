#!/usr/bin/env python3
"""Test jednostkowy TAGTRIG na syntetycznym drzewie CharmAnalysis (uproot): pakowanie bitów, tag AFP, liczniki kandydatów, preskale per LB, scalanie.
Uruchomienie (wymaga uproot, awkward, pyarrow): python analysis/tagtrig/test_tt.py   (katalog roboczy tymczasowy, kasowany)."""
import json, os, shutil, subprocess, sys, tempfile
import numpy as np, awkward as ak, uproot, pandas as pd, pyarrow.parquet as pq

HERE = os.path.dirname(os.path.abspath(__file__)); PK = HERE
sys.path.insert(0, PK)
import tt_bits, tt_merge  # noqa: E402

NCH, NL1 = 70, 5       # 70 łańcuchów → 2 słowa uint64
CH = tt_merge.TRIG11 + [f"HLT_fake{j:02d}_L1X" for j in range(NCH - 11)]
L1 = ["L1_AFP_A_OR_C", "L1_J12", "L1_MBTS_1", "L1_MBTS_2", "L1_RD0_FILLED"]


def synth(fn, seed, n, run=428770, drop=()):
    r = np.random.default_rng(seed)
    d = {"EventInfo_runNumber": np.full(n, run, np.uint32), "EventInfo_eventNumber": (np.arange(n) + seed * 1_000_000).astype(np.uint64),
         "EventInfo_lumiBlock": r.integers(4, 9, n).astype(np.uint32),
         "EventInfo_actualInteractionsPerCrossing": r.uniform(0.04, 0.06, n).astype(np.float32),
         "EventInfo_averageInteractionsPerCrossing": r.uniform(0.04, 0.06, n).astype(np.float32),
         "EventInfo_passGRL_MinBiasAFP": (r.random(n) < 0.8).astype(np.int8),
         "EventInfo_nPileupVertices_NOSYS": (r.random(n) < 0.1).astype(np.int32), "EventInfo_nSelBL_NOSYS": r.integers(0, 60, n).astype(np.int32),
         "EventInfo_xi_cal_A_NOSYS": r.uniform(0, 0.15, n).astype(np.float32), "EventInfo_xi_cal_C_NOSYS": r.uniform(0, 0.15, n).astype(np.float32)}
    npr = r.integers(0, 4, n)
    d["AfpProtons_xi_NOSYS"] = ak.unflatten(r.uniform(0.01, 0.13, npr.sum()).astype(np.float32), npr)
    d["AfpProtons_side_NOSYS"] = ak.unflatten(r.integers(0, 2, npr.sum()).astype(np.int32), npr)
    nd = r.integers(0, 6, n); k = nd.sum()
    d["DMesons_decayType_NOSYS"] = ak.unflatten(r.choice([4, 1, 10, 11], k).astype(np.int32), nd)
    d["DMesons_m_NOSYS"] = ak.unflatten(r.uniform(1650, 2100, k).astype(np.float32), nd)
    d["DMesons_DeltaMass_NOSYS"] = ak.unflatten(r.uniform(130, 170, k).astype(np.float32), nd)
    d["DMesons_fitOutput__Lxy_NOSYS"] = ak.unflatten(r.uniform(-0.2, 1.0, k).astype(np.float32), nd)
    for j, c in enumerate(CH):
        if c in drop:
            continue
        d["EventInfo_trigPassed_" + c + "_NOSYS"] = (r.random(n) < (0.5 if j < 3 else 0.03)).astype(np.int8)
        if j != 20:                                    # jeden łańcuch bez trigPS
            d["EventInfo_trigPS_" + c] = np.where(d["EventInfo_lumiBlock"] < 6, 120.0, 60.0).astype(np.float32) if j % 2 else np.ones(n, np.float32)
    for c in L1:
        d["EventInfo_trigTBP_" + c] = (r.random(n) < 0.4).astype(np.int8)
    with uproot.recreate(fn) as f:
        f["CharmAnalysis"] = d
    return d


def brute(d, chains):
    """Referencja pętlami po zdarzeniach: wiersze tagowane (GRL) i liczniki."""
    rows = {}
    xi, sd = ak.to_list(d["AfpProtons_xi_NOSYS"]), ak.to_list(d["AfpProtons_side_NOSYS"])
    dt, m, dm, lx = (ak.to_list(d["DMesons_" + k + "_NOSYS"]) for k in ("decayType", "m", "DeltaMass", "fitOutput__Lxy"))
    for i in range(len(xi)):
        if d["EventInfo_passGRL_MinBiasAFP"][i] <= 0:
            continue
        tag, x = {}, {}
        for s, S in ((0, "A"), (1, "C")):
            v = [np.float32(a) for a, b in zip(xi[i], sd[i]) if b == s]
            tag[S] = len(v) == 1 and np.float32(0.030) <= v[0] <= np.float32(0.100)
            x[S] = v[0] if tag[S] else -1.0
        if not (tag["A"] or tag["C"]):
            continue
        n0 = sum(1 for t, a, l in zip(dt[i], m[i], lx[i]) if t == 10 and l >= 0.25 and 1.75 <= np.float32(a) * 0.001 <= 2.00)
        n1 = sum(1 for t, a, l in zip(dt[i], m[i], lx[i]) if t == 4 and l >= 0 and 1.72 <= np.float32(a) * 0.001 <= 2.02)
        n2 = sum(1 for t, a, l in zip(dt[i], dm[i], lx[i]) if t == 1 and l >= 0 and 0.139 <= np.float32(a) * 0.001 <= 0.160)
        bits = [c for c in chains if ("EventInfo_trigPassed_" + c + "_NOSYS") in d and d["EventInfo_trigPassed_" + c + "_NOSYS"][i] > 0]
        tb = [c for c in L1 if d["EventInfo_trigTBP_" + c][i] > 0]
        rows[int(d["EventInfo_eventNumber"][i])] = (tag["A"], tag["C"], x["A"], x["C"], n0, n1, n2, bits, tb)
    return rows


def main():
    # 1. pakowanie bitów
    r = np.random.default_rng(1)
    for k in (1, 31, 64, 65, 106, 128, 129):
        M = r.random((500, k)) < 0.3
        W = tt_bits.pack_bits(M)
        assert W.shape == (500, (k + 63) // 64) and W.dtype == np.uint64
        assert np.array_equal(tt_bits.unpack_bits(W, k), M)
        assert all(int(W[i, j // 64]) >> (j % 64) & 1 == int(M[i, j]) for i in range(0, 500, 97) for j in range(k))
    src = ["b", "a", "d"]; dst = ["a", "b", "c", "d"]
    M = r.random((50, 3)) < 0.5; X = tt_merge.remap(M, src, dst)
    assert np.array_equal(X[:, 1], M[:, 0]) and np.array_equal(X[:, 0], M[:, 1]) and np.array_equal(X[:, 3], M[:, 2]) and not X[:, 2].any()
    print("OK pack/unpack/remap")
    wd = tempfile.mkdtemp(prefix="tt_test_")
    try:
        raw = os.path.join(wd, "raw", "data_428770", "sub"); os.makedirs(raw); os.makedirs(os.path.join(wd, "out", "parts"))
        D = []
        for i in range(3):   # plik 2 (część 0 przy nparts = 2: pliki 0 i 2; część 1: plik 1) — część 1 bez jednego łańcucha → remap w merge
            D.append(synth(os.path.join(raw, f"user.x.123456.EXT0._{i + 1:06d}.CharmAnalysis.dpt3to150.root"), i + 1, 4000,
                           drop=("HLT_fake05_L1X",) if i == 1 else ()))
        synth(os.path.join(raw, "user.x.123456.EXT0._000009.CharmAnalysis.dpt2to3.root"), 9, 100)      # musi być pominięty
        for k in range(2):
            subprocess.run([sys.executable, os.path.join(PK, "tt_bits.py"), "--sample", "data_428770", "--rootdir", os.path.dirname(raw), "--part", str(k),
                            "--nparts", "2", "--no-copy", "--step", "1MB" if k else "100MB", "--out-prefix", os.path.join(wd, "out", "parts", f"data_428770.p{k}")],
                           check=True, stdout=subprocess.DEVNULL)
        m1 = json.load(open(os.path.join(wd, "out", "parts", "data_428770.p1.meta.json")))
        assert len(m1["chains"]) == NCH - 1 and m1["n_files_total"] == 3 and len(m1["files"]) == 1
        # merge: zła liczba plików → odmowa; poprawna → wynik
        cmd = [sys.executable, os.path.join(PK, "tt_merge.py"), "--out", os.path.join(wd, "out"), "--wd", os.path.join(wd, "fin"), "--samples", "data_428770", "--nparts", "2"]
        assert subprocess.run(cmd + ["--nfiles", "data_428770:4"], capture_output=True).returncode != 0
        assert subprocess.run(cmd + ["--nfiles", "data_428770:3", "--max-mb", "0.001"], capture_output=True).returncode == 3
        subprocess.run(cmd + ["--nfiles", "data_428770:3"], check=True)
        fin = os.path.join(wd, "fin")
        S = json.load(open(os.path.join(fin, "summary.json")))["samples"]["data_428770"]
        C = json.load(open(os.path.join(fin, "chains_428770.json")))
        assert C["chains"] == sorted(CH) and C["l1tbp"] == sorted(L1) and C["ps_missing"] == [CH[20]] and not S["chain_lists_identical"]
        ev = pq.read_table(os.path.join(fin, "bits_428770.parquet")).to_pandas()
        assert not any("mass" in c.lower() or c.lower().endswith("_m") for c in ev.columns), "kolumna masy w wyjściu!"
        ref = {}
        for d in D:
            ref.update(brute(d, CH))
        assert len(ev) == len(ref) == S["n_T0"] and set(ev["event_number"]) == set(ref) and S["n_dup_event_number"] == 0
        Mb = tt_bits.unpack_bits(ev[["trig_w0", "trig_w1"]].to_numpy(), NCH); Tb = tt_bits.unpack_bits(ev[["tbp_w0"]].to_numpy(), NL1)
        for i, e in enumerate(ev.itertuples()):
            a, c, xa, xc, n0, n1, n2, bits, tb = ref[e.event_number]
            assert (bool(e.tag_A), bool(e.tag_C)) == (a, c) and e.n_D0 == n0 and e.n_Dp == n1 and e.n_Dst == n2, e
            assert abs(e.xi_A - xa) < 1e-7 and abs(e.xi_C - xc) < 1e-7
            assert sorted(np.array(C["chains"])[Mb[i]]) == sorted(bits) and sorted(np.array(C["l1tbp"])[Tb[i]]) == sorted(tb), e.event_number
        print(f"OK zdarzenia tagowane: {len(ev)} wierszy, tag/liczniki/bity zgodne z referencją (remap części bez łańcucha)")
        # preskale per LB × łańcuch
        ps = pd.read_csv(os.path.join(fin, "ps_perLB_428770.csv"))
        lb = np.concatenate([d["EventInfo_lumiBlock"] for d in D]); g = np.concatenate([d["EventInfo_passGRL_MinBiasAFP"] for d in D]) > 0
        for c in (CH[0], CH[1], CH[20], CH[33]):
            p = np.concatenate([d["EventInfo_trigPassed_" + c + "_NOSYS"] for d in D]) > 0
            for u in np.unique(lb):
                row = ps[(ps.chain == c) & (ps.lumiBlock == u) & (ps.kind == "HLT")]
                assert len(row) == 1 and row.n_events.iloc[0] == (g & (lb == u)).sum() and row.n_pass.iloc[0] == (g & (lb == u) & p).sum()
                exp = -1.0 if c == CH[20] else ((120.0 if u < 6 else 60.0) if CH.index(c) % 2 else 1.0)
                assert row.ps_min.iloc[0] == row.ps_max.iloc[0] == row.ps_mode.iloc[0] == exp and row.ps_nvals.iloc[0] == 1, (c, u, row)
        row = ps[(ps.chain == "HLT_fake05_L1X")]                # łańcuch nieobecny w pliku 2: liczony tylko z plików 1 i 3
        l2 = np.concatenate([D[0]["EventInfo_lumiBlock"], D[2]["EventInfo_lumiBlock"]]); g2 = np.concatenate([D[0]["EventInfo_passGRL_MinBiasAFP"], D[2]["EventInfo_passGRL_MinBiasAFP"]]) > 0
        assert row.n_events.sum() == g2.sum() and len(row) == len(np.unique(l2))
        t = np.concatenate([d["EventInfo_trigTBP_L1_J12"] for d in D]) > 0
        assert ps[(ps.chain == "L1_J12") & (ps.kind == "L1TBP")].n_pass.sum() == (g & t).sum()
        ls = pd.read_csv(os.path.join(fin, "lb_summary_428770.csv"))
        assert ls.n_all.sum() == 12000 and ls.n_grl.sum() == g.sum() == S["n_grl"] and ls.n_tag.sum() == len(ev) and ls.n_T2.sum() == S["n_T2"]
        # T1/T2 i klasa none z referencji
        T1 = ev.npu == 0
        T2 = T1 & ((ev.tag_A == 0) | (ev.xi_cal_A <= ev.xi_A + np.float32(0.02))) & ((ev.tag_C == 0) | (ev.xi_cal_C <= ev.xi_C + np.float32(0.02)))
        i11 = [C["chains"].index(c) for c in tt_merge.TRIG11]
        none11 = ~Mb[:, i11].any(1)
        assert S["n_T1"] == T1.sum() and S["n_T2"] == T2.sum() and S["none"]["T2_none11"] == (T2 & none11).sum()
        assert S["none"]["T2_noneAll"] == (T2 & ~Mb.any(1)).sum() and 0 < S["none"]["T2_none11"] < S["n_T2"]
        top = dict(S["none"]["chains_T2_none11"]); j = C["chains"].index(CH[40])
        assert top[CH[40]] == int(Mb[(T2 & none11).to_numpy(), j].sum()) and not any(c in top for c in tt_merge.TRIG11)
        assert sum(S["none"]["n_other_bits_hist_T2_none11"]) == S["none"]["T2_none11"]
        print(f"OK preskale per LB, lb_summary, summary: T0 {S['n_T0']} T1 {S['n_T1']} T2 {S['n_T2']} none11 {S['none']['T2_none11']} noneAll {S['none']['T2_noneAll']}")
        # --max-events
        subprocess.run([sys.executable, os.path.join(PK, "tt_bits.py"), "--sample", "data_428770", "--rootdir", os.path.dirname(raw), "--no-copy", "--max-files", "1",
                        "--max-events", "1000", "--out-prefix", os.path.join(wd, "x")], check=True, stdout=subprocess.DEVNULL)
        assert json.load(open(os.path.join(wd, "x.meta.json")))["n_events"] == 1000
        print("OK --max-events")
    finally:
        shutil.rmtree(wd, ignore_errors=True)
    print("WSZYSTKO OK")


if __name__ == "__main__":
    main()
