#!/usr/bin/env python3
"""TAGTRIG / finish (18.11.2026, projekt B): scalenie części tt_bits.py → pliki końcowe + summary.json. numpy/pandas/pyarrow.
Wejście: <out>/parts/<próbka>.p<k>.{events,ps,lb}.parquet + .meta.json (k = 0..nparts−1). Wyjście w --wd (na EOS kładzie tt_job.sh):
  bits_<run>.parquet        — zdarzenia tagowane (kolumny jak w tt_bits.py; bity wg chains_<run>.json)
  ps_perLB_<run>.csv        — run, lumiBlock, kind (HLT | L1TBP), chain, n_events (GRL w LB), n_pass, ps_min, ps_max, ps_mode, ps_nvals
  lb_summary_<run>.csv      — lumiBlock, n_all, n_grl, mean_mu_grl, n_tag (T0), n_T1, n_T2, n_T2_D0
  chains_<run>.json         — chains (HLT, kolejność bitów trig_w*), l1tbp (kolejność bitów tbp_w*), ps_missing
  summary.json              — liczby per run, rozkład klasy `none` (żaden z 11 łańcuchów zapisanych w h5) na łańcuchy — tylko liczniki
Kontrole: komplet części, każdy plik ROOT dokładnie raz, liczba plików = --nfiles; łączny rozmiar ≤ --max-mb (inaczej kod 3, bez DONE).
T0 = GRL ∧ tag(A ∨ C); T1 = T0 ∧ npu = 0; T2 = T1 ∧ ξ_cal,s ≤ ξ_s + 0,02 na każdej tagowanej stronie.
"""
import argparse, json, os, sys
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tt_bits import pack_bits, unpack_bits, codec  # noqa: E402

TRIG11 = ["HLT_j20_L1MBTS_2", "HLT_mb_mbts_L1MBTS_1", "HLT_mb_sptrk_L1MBTS_1", "HLT_mb_sptrk_pt2_L1MBTS_2", "HLT_mb_sptrk_pt4_L1AFP_A_OR_C",
          "HLT_mb_sptrk_pt4_L1MBTS_2", "HLT_noalg_L1AFP_A_OR_C", "HLT_noalg_L1J12", "HLT_noalg_L1RD0_FILLED", "HLT_noalg_mb_L1MBTS_1",
          "HLT_noalg_mb_L1MBTS_2"]          # = charmpp_h5_v20.TRIG_STORE (łańcuchy zapisane w h5 v2.0)
DELTA = 0.02


def words(tb, pre):
    cols = sorted((c for c in tb.column_names if c.startswith(pre)), key=lambda c: int(c[len(pre):]))
    return np.stack([tb.column(c).to_numpy() for c in cols], 1).astype(np.uint64) if cols else np.zeros((tb.num_rows, 0), np.uint64)


def remap(M, src, dst):
    """Macierz bool w kolejności src → kolejność dst (brakujące = 0)."""
    if src == dst:
        return M
    out = np.zeros((M.shape[0], len(dst)), bool)
    pos = {c: j for j, c in enumerate(dst)}
    for j, c in enumerate(src):
        out[:, pos[c]] = M[:, j]
    return out


def regions(d):
    tagA, tagC = d["tag_A"] > 0, d["tag_C"] > 0
    ok = (~tagA | (d["xi_cal_A"] <= d["xi_A"] + DELTA)) & (~tagC | (d["xi_cal_C"] <= d["xi_C"] + DELTA))
    T1 = d["npu"] == 0
    return T1, T1 & ok


def top(cnt, names, k=40):
    o = np.argsort(-cnt, kind="stable")[:k]
    return [[names[j], int(cnt[j])] for j in o if cnt[j] > 0]


def merge_sample(s, out, wd, nparts, nexp, allow, cd):
    P = os.path.join(out, "parts")
    metas = []
    for k in range(nparts):
        fn = os.path.join(P, f"{s}.p{k}.meta.json")
        if not os.path.exists(fn):
            sys.exit(f"{s}: brak części {k} ({fn})")
        metas.append(json.load(open(fn)))
    files = [f for m in metas for f in m["files"]]
    if len(set(files)) != len(files):
        sys.exit(f"{s}: plik ROOT w więcej niż jednej części")
    if not allow and len(files) != nexp:
        sys.exit(f"{s}: {len(files)} plików ROOT, oczekiwane {nexp}")
    runs = {m["run_number"] for m in metas}
    if len(runs) != 1:
        sys.exit(f"{s}: różne run_number w częściach {runs}")
    run = runs.pop()
    chains = sorted(set(c for m in metas for c in m["chains"]))
    l1 = sorted(set(c for m in metas for c in m["l1tbp"]))
    same = all(m["chains"] == chains and m["l1tbp"] == l1 for m in metas)
    i11 = [chains.index(c) for c in TRIG11 if c in chains]
    R = {"run": run, "sample": s, "n_files": len(files), "n_events": sum(m["n_events"] for m in metas), "n_grl": sum(m["n_grl"] for m in metas),
         "n_chains": len(chains), "n_l1tbp": len(l1), "chain_lists_identical": same, "trig11_missing": [c for c in TRIG11 if c not in chains],
         "ps_missing": sorted(set(c for m in metas for c in m["ps_missing"])), "seconds_parts": [m["seconds"] for m in metas]}
    nm = ["T0", "T1", "T2", "T2_D0"]
    N = {k: 0 for k in nm + ["T2_A_only", "T2_C_only", "T2_AC", "T2_Dp", "T2_Dst"]}
    cnt = {k: np.zeros(len(chains), np.int64) for k in ("T2", "T2_D0", "T2_none11", "T2_none11_D0")}
    tbp = {k: np.zeros(len(l1), np.int64) for k in ("T2", "T2_none11", "T2_noneAll")}
    nn = {"T2_none11": 0, "T2_none11_D0": 0, "T2_noneAll": 0, "T2_noneAll_D0": 0, "T2_none11_AC": 0}
    nbits = np.zeros(8, np.int64)        # T2 ∧ none11: liczba ustawionych bitów pozostałych łańcuchów (0..6, ≥7)
    lbt, evn, wr = [], [], None
    ofn = os.path.join(wd, f"bits_{run}.parquet")
    for k, m in enumerate(metas):
        tb = pq.read_table(os.path.join(P, f"{s}.p{k}.events.parquet"))
        M = remap(unpack_bits(words(tb, "trig_w"), len(m["chains"])), m["chains"], chains)
        T = remap(unpack_bits(words(tb, "tbp_w"), len(m["l1tbp"])), m["l1tbp"], l1)
        if not same:
            cols = {c: tb.column(c) for c in tb.column_names if not (c.startswith("trig_w") or c.startswith("tbp_w"))}
            for pre, X in (("trig_w", M), ("tbp_w", T)):
                W = pack_bits(X)
                for i in range(W.shape[1]):
                    cols[f"{pre}{i}"] = pa.array(W[:, i])
            tb = pa.table(cols)
        if wr is None:
            wr = pq.ParquetWriter(ofn, tb.schema, compression=cd)
        wr.write_table(tb)
        d = {c: tb.column(c).to_numpy() for c in ("lumiBlock", "npu", "tag_A", "tag_C", "xi_A", "xi_C", "xi_cal_A", "xi_cal_C", "n_D0", "n_Dp", "n_Dst")}
        evn.append(tb.column("event_number").to_numpy())
        T1, T2 = regions(d)
        hD0 = d["n_D0"] > 0
        a, c = d["tag_A"] > 0, d["tag_C"] > 0
        N["T0"] += len(T1); N["T1"] += int(T1.sum()); N["T2"] += int(T2.sum()); N["T2_D0"] += int((T2 & hD0).sum())
        N["T2_Dp"] += int((T2 & (d["n_Dp"] > 0)).sum()); N["T2_Dst"] += int((T2 & (d["n_Dst"] > 0)).sum())
        N["T2_A_only"] += int((T2 & a & ~c).sum()); N["T2_C_only"] += int((T2 & c & ~a).sum()); N["T2_AC"] += int((T2 & a & c).sum())
        none11 = ~M[:, i11].any(1) if i11 else np.ones(len(T1), bool)
        noneall = ~M.any(1)
        sets = {"T2": T2, "T2_D0": T2 & hD0, "T2_none11": T2 & none11, "T2_none11_D0": T2 & none11 & hD0}
        for key, msk in sets.items():
            cnt[key] += M[msk].sum(0)
        for key, msk in (("T2", T2), ("T2_none11", T2 & none11), ("T2_noneAll", T2 & noneall)):
            tbp[key] += T[msk].sum(0)
        nn["T2_none11"] += int((T2 & none11).sum()); nn["T2_none11_D0"] += int((T2 & none11 & hD0).sum())
        nn["T2_noneAll"] += int((T2 & noneall).sum()); nn["T2_noneAll_D0"] += int((T2 & noneall & hD0).sum())
        nn["T2_none11_AC"] += int((T2 & none11 & a & c).sum())
        nbits += np.bincount(np.minimum(M[T2 & none11].sum(1), 7), minlength=8)
        lbt.append(pd.DataFrame({"lumiBlock": d["lumiBlock"], "n_tag": 1, "n_T1": T1.astype(int), "n_T2": T2.astype(int),
                                 "n_T2_D0": (T2 & hD0).astype(int)}).groupby("lumiBlock", as_index=False).sum())
    wr.close()
    evn = np.concatenate(evn)
    R.update({"n_" + k: v for k, v in N.items()})
    R["n_dup_event_number"] = int(len(evn) - len(np.unique(evn)))
    R["none"] = dict(nn, def_none11="żaden bit z 11 łańcuchów zapisanych w h5 v2.0 (TRIG11)", def_noneAll="żaden bit z wszystkich łańcuchów drzewa",
                     n_other_bits_hist_T2_none11=nbits.tolist(),
                     chains_T2_none11=top(cnt["T2_none11"], chains), chains_T2_none11_D0=top(cnt["T2_none11_D0"], chains),
                     l1tbp_T2_none11=top(tbp["T2_none11"], l1), l1tbp_T2_noneAll=top(tbp["T2_noneAll"], l1))
    R["chains_T2"] = top(cnt["T2"], chains, 200); R["chains_T2_D0"] = top(cnt["T2_D0"], chains, 200); R["l1tbp_T2"] = top(tbp["T2"], l1, 200)
    # preskale per LB × łańcuch
    ps = []
    for k, m in enumerate(metas):
        p = pq.read_table(os.path.join(P, f"{s}.p{k}.ps.parquet")).to_pandas()
        if not same:
            mc = np.array([chains.index(c) for c in m["chains"]] or [0]); ml = np.array([l1.index(c) for c in m["l1tbp"]] or [0])
            idx = p["idx"].to_numpy()
            p["idx"] = np.where(p["kind"].to_numpy() == 0, mc[np.minimum(idx, len(mc) - 1)], ml[np.minimum(idx, len(ml) - 1)])
        ps.append(p)
    ps = pd.concat(ps, ignore_index=True).groupby(["kind", "idx", "lumiBlock", "ps"], as_index=False)[["n_ev", "n_pass"]].sum()
    ps = ps.sort_values(["kind", "idx", "lumiBlock", "n_ev", "ps"], ascending=[True, True, True, False, True])
    g = ps.groupby(["kind", "idx", "lumiBlock"], sort=True)
    t = g.agg(n_events=("n_ev", "sum"), n_pass=("n_pass", "sum"), ps_min=("ps", "min"), ps_max=("ps", "max"), ps_mode=("ps", "first"),
              ps_nvals=("ps", "nunique")).reset_index()
    names = np.array(chains + [""], object); l1n = np.array(l1 + [""], object)
    t["chain"] = np.where(t["kind"] == 0, names[np.minimum(t["idx"], len(chains))], l1n[np.minimum(t["idx"], len(l1))])
    t["kind"] = np.where(t["kind"] == 0, "HLT", "L1TBP")
    t.insert(0, "run", run)
    t[["run", "lumiBlock", "kind", "chain", "n_events", "n_pass", "ps_min", "ps_max", "ps_mode", "ps_nvals"]].to_csv(
        os.path.join(wd, f"ps_perLB_{run}.csv"), index=False)
    R["n_lb_grl"] = int(t["lumiBlock"].nunique()); R["n_ps_rows"] = int(len(t))
    R["n_lb_chain_ps_not_constant"] = int(((t["kind"] == "HLT") & (t["ps_nvals"] > 1)).sum())
    h = t[t["kind"] == "HLT"]
    R["ps_sets_TRIG11"] = {c: {str(v): int(n) for v, n in h[h["chain"] == c]["ps_mode"].value_counts().sort_index().items()} for c in TRIG11 if c in chains}
    # LB
    lb = pd.concat([pq.read_table(os.path.join(P, f"{s}.p{k}.lb.parquet")).to_pandas() for k in range(nparts)]).groupby("lumiBlock", as_index=False).sum()
    lb["mean_mu_grl"] = lb["sum_mu_grl"] / lb["n_grl"].where(lb["n_grl"] > 0)
    tg = pd.concat(lbt).groupby("lumiBlock", as_index=False).sum() if lbt else pd.DataFrame(columns=["lumiBlock", "n_tag", "n_T1", "n_T2", "n_T2_D0"])
    lb = lb.drop(columns="sum_mu_grl").merge(tg, on="lumiBlock", how="left").fillna({"n_tag": 0, "n_T1": 0, "n_T2": 0, "n_T2_D0": 0})
    for c in ("n_tag", "n_T1", "n_T2", "n_T2_D0"):
        lb[c] = lb[c].astype(int)
    lb.to_csv(os.path.join(wd, f"lb_summary_{run}.csv"), index=False)
    if int(lb["n_grl"].sum()) != R["n_grl"] or int(lb["n_all"].sum()) != R["n_events"] or int(lb["n_tag"].sum()) != R["n_T0"]:
        sys.exit(f"{s}: niezgodne sumy LB vs meta")
    json.dump({"run": run, "chains": chains, "l1tbp": l1, "ps_missing": R["ps_missing"], "trig11": TRIG11,
               "bit_layout": "łańcuch chains[j] → kolumna trig_w{j // 64}, bit j % 64 (LSB = 0); L1 l1tbp[j] → tbp_w{j // 64}"},
              open(os.path.join(wd, f"chains_{run}.json"), "w"), indent=1)
    R["files_out"] = {f: os.path.getsize(os.path.join(wd, f)) for f in (f"bits_{run}.parquet", f"ps_perLB_{run}.csv", f"lb_summary_{run}.csv", f"chains_{run}.json")}
    return R


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True); ap.add_argument("--wd", required=True)
    ap.add_argument("--samples", required=True); ap.add_argument("--nparts", type=int, required=True)
    ap.add_argument("--nfiles", default=""); ap.add_argument("--max-mb", type=float, default=300.0)
    ap.add_argument("--allow-incomplete", action="store_true", help="test: bez kontroli liczby plików")
    a = ap.parse_args()
    nexp = {kv.split(":")[0]: int(kv.split(":")[1]) for kv in a.nfiles.replace(",", " ").split()}
    os.makedirs(a.wd, exist_ok=True)
    S = {"samples": {}, "max_mb": a.max_mb, "allow_incomplete": a.allow_incomplete}
    for s in a.samples.replace(",", " ").split():
        S["samples"][s] = merge_sample(s, a.out, a.wd, a.nparts, nexp.get(s, -1), a.allow_incomplete, codec())
    # węzeł lumi (niekrytyczny): status i lista plików, jeśli coś zostawił w <out>
    st = os.path.join(a.out, "lumi_STATUS.txt")
    S["lumi"] = {"status_file": os.path.exists(st), "status": open(st, errors="replace").read()[-4000:] if os.path.exists(st) else
                 "brak lumi_STATUS.txt (węzeł lumi wyłączony, nie wystartował albo padł przed zapisem) — świetlności i preskal z COOL NIE MA",
                 "files": sorted(f for f in (os.listdir(a.out) if os.path.isdir(a.out) else []) if f.startswith(("cool_", "coolps_", "ilumicalc_", "lumi_")))}
    tot = sum(sum(r["files_out"].values()) for r in S["samples"].values())
    tot += sum(os.path.getsize(os.path.join(a.out, f)) for f in S["lumi"]["files"])
    S["total_bytes"] = tot; S["total_mb"] = round(tot / 1e6, 2); S["size_ok"] = tot <= a.max_mb * 1e6
    json.dump(S, open(os.path.join(a.wd, "summary.json"), "w"), indent=1, ensure_ascii=False)
    for s, r in S["samples"].items():
        print(f"{s}: run {r['run']} plików {r['n_files']} zdarzeń {r['n_events']} GRL {r['n_grl']} T0 {r['n_T0']} T1 {r['n_T1']} T2 {r['n_T2']} "
              f"T2∧D0 {r['n_T2_D0']} łańcuchów {r['n_chains']}+{r['n_l1tbp']} none11(T2) {r['none']['T2_none11']} noneAll(T2) {r['none']['T2_noneAll']} "
              f"dup {r['n_dup_event_number']}")
    print(f"razem {S['total_mb']} MB (limit {a.max_mb}) lumi_STATUS: {'jest' if S['lumi']['status_file'] else 'BRAK'}")
    if not S["size_ok"]:
        sys.exit(3)


if __name__ == "__main__":
    main()
