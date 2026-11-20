"""B_trig — histogramy Δm D* (|m(Kπ)−m_D0|<55 MeV) per łańcuch triggera z ROOT v2.0 (tylko odczyt EOS). Użycie (lxplus, LCG):
python3 dshist_root.py N_428770 N_427929 > dshist.json   (N = liczba plików, równomiernie; 97 458 = całość, ~35 min)"""
import sys, glob, json
import numpy as np, uproot, awkward as ak
TR = ["HLT_mb_sptrk_L1MBTS_1", "HLT_noalg_L1J12", "HLT_j20_L1MBTS_2", "HLT_mb_sptrk_pt4_L1MBTS_2", "HLT_mb_mbts_L1MBTS_1", "HLT_noalg_mb_L1MBTS_1", "HLT_noalg_mb_L1MBTS_2"]
def res(keys, n):
    for x in (n + "_NOSYS", n):
        if x in keys: return x
out = {}
E = np.arange(139.6, 170.0001, 0.1)
import time; T0 = time.time()
for samp, nf in (("data_428770", int(sys.argv[1])), ("data_427929", int(sys.argv[2]))):
    fl = sorted(glob.glob(f"/eos/user/k/klyczak/charm_raw_v20/{samp}/*dpt3to150*.root"))
    sel = [fl[i] for i in np.linspace(0, len(fl) - 1, nf).astype(int)]
    acc = {"files": len(sel), "nfiles_total": len(fl), "entries": 0, "ds": {}}
    for f in sel:
        t = uproot.open(f)["CharmAnalysis"]; k = t.keys()
        br = {"dt": res(k, "DMesons_decayType"), "pt": res(k, "DMesons_pt"), "dm": res(k, "DMesons_DeltaMass"), "d0i": res(k, "DMesons_D0Index"), "m": res(k, "DMesons_m")}
        tb = {c: res(k, f"EventInfo_trigPassed_{c}") for c in TR}; tb = {c: b for c, b in tb.items() if b}
        ps = {c: res(k, f"EventInfo_trigPS_{c}") for c in TR}; ps = {c: b for c, b in ps.items() if b}
        g = res(k, "EventInfo_passGRL_MinBiasAFP")
        a = t.arrays(list(br.values()) + list(tb.values()) + list(ps.values()) + ([g] if g else []), library="ak")
        acc["entries"] += len(a)
        dt = a[br["dt"]]; pt = a[br["pt"]] / 1e3; dm = a[br["dm"]] / 1e3; d0i = a[br["d0i"]]
        md0 = (a[br["m"]] - a[br["dm"]]) / 1e3
        isd = (dt == 1) & (d0i >= 0) & (abs(md0 - 1.86484) < 0.055)
        if not acc.get("mchk"): acc["mchk"] = [float(x) for x in np.percentile(ak.to_numpy(ak.flatten(md0[dt == 1])), [1, 50, 99])]
        cls = {"all": np.ones(len(a), bool)}
        for c, b in tb.items(): cls[c] = ak.to_numpy(a[b]).astype(bool)
        if g: cls["grl"] = ak.to_numpy(a[g]).astype(bool)
        for c, b in ps.items():
            v = ak.to_numpy(a[b]).astype(float); acc.setdefault("ps", {}).setdefault(c, {})
            for u, n in zip(*np.unique(np.round(v, 1), return_counts=True)): acc["ps"][c][str(u)] = acc["ps"][c].get(str(u), 0) + int(n)
        for lo, hi in ((3, 5), (5, 8), (8, 30)):
            m = isd & (pt >= lo) & (pt < hi)
            for c, ev in cls.items():
                x = ak.to_numpy(ak.flatten(dm[m][ev])) * 1e3
                h, _ = np.histogram(x, E)
                d = acc["ds"].setdefault(f"{c}|{lo}-{hi}", np.zeros(len(E) - 1, np.int64).tolist())
                acc["ds"][f"{c}|{lo}-{hi}"] = (np.asarray(d) + h).tolist()
        for c, ev in cls.items():
            acc.setdefault("nev", {}); acc["nev"][c] = acc["nev"].get(c, 0) + int(ev.sum())
    acc["sec"] = time.time() - T0
    out[samp] = acc
print(json.dumps(out))
