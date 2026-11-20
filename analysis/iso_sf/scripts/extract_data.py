"""ISO_sf: D0 (wiersz Kπ, decayType 10) z tagiem D* w danych 428770 — tylko pliki lokalne z kolumnami izolacji (h5 v1.2).
Złączenie: HYP_d0/runs/lxplus_dstar/data_428770.parquet (file,row → wiersz h5; row globalny = offset pliku + lokalny).
Krotność zdarzenia (brak n_trk w h5): n_cand = liczba kandydatów (wszystkie decayType) w bloku zdarzenia; n_trkid = liczba
różnych trackId córek wszystkich kandydatów w zdarzeniu (proxy liczby torów). Wyjście: out/data_tag.parquet
"""
import json, os, sys
import h5py, numpy as np, pandas as pd, pyarrow.parquet as pq
from evtools import event_mult

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HYP = "/home/k/Documents/charm/charm_work/HYP_d0/runs/lxplus_dstar"
H5 = "/home/k/Documents/charm/charm_work/FEAT_charmnet/data"


def run(sample, h5dir, out):
    meta = json.load(open(f"{HYP}/{sample}.meta.json"))
    files, offs = meta["files"], meta["file_offsets"]
    local = sorted(os.listdir(f"{H5}/{h5dir}"))
    fidx = [i for i, f in enumerate(files) if f in local]
    cols = ["file", "row", "evblk", "m", "pt", "dm", "score", "Lxy", "pis_pt", "dstar_conflict", "dstar_n", "truth", "hypo_swap",
            "hypo_from_tag", "d_trk0", "d_trk1"]
    t = pq.read_table(f"{HYP}/{sample}.parquet", columns=cols, filters=[("file", "in", fidx)]).to_pandas()
    parts = []
    for fi in fidx:
        s = t[t.file == fi].copy()
        lr = s.row.to_numpy() - (offs[fi] if s.row.max() >= offs[fi] and offs[fi] > 0 else 0)
        with h5py.File(f"{H5}/{h5dir}/{files[fi]}", "r") as f:
            ds = f["events"]
            ev = ds.fields(["event_number"])[:]["event_number"]
            mult = event_mult(ds, ev)       # per-row: n_cand, n_trkid
            a = ds.fields(["Dmeson_ptcone40", "Dmeson_ptcone100", "Dmeson_mass", "Dmeson_pt", "Dmeson_decayType", "actual_mu",
                           "Dmeson_truth_category"])[:]
        sel = lambda c: a[c][lr]
        assert np.allclose(sel("Dmeson_mass"), s.m.to_numpy(), atol=1e-4), "mass mismatch"
        assert (sel("Dmeson_decayType") == 10).all()
        s["ptcone40"], s["ptcone100"] = sel("Dmeson_ptcone40"), sel("Dmeson_ptcone100")
        s["mu"], s["truth_cat"] = sel("actual_mu"), sel("Dmeson_truth_category")
        s["n_cand"], s["n_trkid"] = mult["n_cand"][lr], mult["n_trkid"][lr]
        parts.append(s)
        print(sample, files[fi], len(s), flush=True)
    d = pd.concat(parts, ignore_index=True)
    d.to_parquet(out)
    print(d.describe().T.to_string())


if __name__ == "__main__":
    smp = sys.argv[1]
    run(smp, {"data_428770": "data_428770", "fd_nd": "fd_nd", "mb_nd": "mb_nd"}[smp], os.path.join(BASE, "out", f"{smp}_tag.parquet"))
