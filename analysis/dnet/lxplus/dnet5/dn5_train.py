#!/usr/bin/env python3
"""DNET etap 5a — trening ZRÓWNOLEGLONY PO FOLDACH: jeden fold = jeden job (fold), osobny węzeł łączący (merge) liczy score OOF i meta
w formacie etapu 4 (sc/<tag>.npz, meta/<tag>.json, models/<tag>/{fold<k>.pt, meta.json}) → analyze.py / dn4_apply bez zmian.

  dn5_train.py fold  --data DIR1/data --env DIR2/parts --env-dp DIR4/parts_dp --hsc DIR1 --out DIR --block B --fold k [--device cpu|auto]
  dn5_train.py merge --data … --env … --env-dp … --hsc DIR1 --models DIR/models --folds-dir DIR/folds --out DIR --block B
  dn5_train.py list
  dn5_train.py equiv --dir KATALOG [--folds 3 --epochs 2 --n 4000]    test: stara ścieżka (dn4_train.py train, foldy sekwencyjnie) ≡ nowa (fold × K + merge)
                                                                      na danych syntetycznych, D0 i D+; potem wszystkie bloki 1 epoka (mechanika)
  dn5_train.py bench [--threads 1,2,4,8] [--n 40960]                   czas epoki (przeliczony na pulę etapu 4) vs liczba wątków, tr i trp

Równoważność ze starą ścieżką (dn4_train.run_tr): pula, foldy, wagi, standaryzacja, DM.fit i ziarno 1000·s + k są te same; inicjalizacja wag sieci
foldu k w starej ścieżce zależy od stanu generatora torch po foldzie k − 1 (DM.fit kończy epokę wywołaniem torch.manual_seed(1000·s + (k − 1) + epok)),
a foldu 0 — od stanu generatora świeżego procesu, który NIE jest ustalony (test 19.11: dwa uruchomienia starej ścieżki dają różne wagi foldu 0).
Job foldu ustawia ten stan jawnie (init_rng; dla foldu 0 ta sama reguła z k − 1 = −1), więc: foldy k ≥ 1 są bitowo równe starej ścieżce przy tej
samej liczbie wątków, a fold 0 jest w nowej ścieżce powtarzalny (w starej nie był) — równoważny statystycznie, nie bitowo."""
import argparse
import json
import os
import socket
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dn5_lib as L5  # noqa: E402
import dn_feat as DF  # noqa: E402
import dn_train as DT  # noqa: E402

MC = ("fd_nd", "mb_nd")


def init_rng(c, k, epochs):
    import torch
    torch.manual_seed(1000 * c["seed"] + (k - 1) + epochs)


def prepare(c, a, samples, log):
    edir = a.env if c["sp"] == "D0" else a.env_dp
    S, full = L5.load(a.data, edir, c["sp"], samples)
    if "fd_nd" not in S or "mb_nd" not in S:
        raise SystemExit(f"brak torów otoczenia {c['sp']} dla fd_nd/mb_nd w {edir}")
    log(f"wczytano {c['sp']} ({edir}): " + ", ".join(f"{s} a={len(T.get('a_m', []))}/{full.get((s, 'a_'))}" for s, T in S.items()))
    hp, thr = None, None
    if c["hcut"]:
        hp, thr, _ = L5.hpass(S, c["sp"], a.hsc, c["hcut"])
        log(f"cięcie H WP{int(100 * c['hcut'])}: progi {['%.4f' % t for t in thr]}; przechodzi " + ", ".join(f"{s} {v.mean():.3f}" for s, v in hp.items()))
    P = L5.pool(S, c["sp"], c, a.cap_sig, a.cap_bkg, hp)
    return S, full, P, thr


def pool_inputs(c, S, P):
    g = lambda k: np.concatenate([S[s]["a_" + k][i] for s, i, _ in P["src"]])  # noqa: E731
    return L5.inputs(c, g("X"), g("Dd"), g("E"), g("eta"), g("phi"), g("pt"))


def cmd_fold(a, log):
    import torch
    import dn_models as DM
    import dn2_models as D2
    import dn5_models as D5
    c, k, K = L5.BLOCKS[a.block], a.fold, a.folds
    torch.set_num_threads(max(1, a.threads))
    DM.HP["epochs"] = a.epochs
    dev, devname = a.device, "cpu"
    if dev == "auto":
        dev = "cuda" if torch.cuda.is_available() else "cpu"
    if dev == "cuda":
        try:                                    # karta bez jądra CUDA dla tej wersji torch (np. V100) → CPU zamiast padać (jak dn2_train)
            _t = torch.randn(64, 8, device="cuda")
            float((torch.cat([_t, _t], -1) @ torch.randn(16, 4, device="cuda")).sum())
            devname = torch.cuda.get_device_name(0)
        except Exception as e:  # noqa: BLE001
            log(f"UWAGA: CUDA nie działa ({type(e).__name__}: {e}) → CPU")
            dev = "cpu"
    t0 = time.time()
    S, full, P, thr = prepare(c, a, MC, log)
    tok, mask, G = pool_inputs(c, S, P)
    fold = DF.fold_of(P["uid"], K)
    val = DF.unit_hash(P["uid"], salt=23) >= 0.9
    tload = time.time() - t0
    log(f"[{c['tag']}] fold {k}/{K}: pula {len(P['y'])} (S {int(P['y'].sum())}, B {int((P['role'] == 'B').sum())}), tokenów {tok.shape[1]}, "
        f"średnio ważnych {mask.sum(1).mean():.2f}, glob {G.shape[1]}, urządzenie {dev} ({devname}), wątków {a.threads}, torch {torch.__version__}")
    trn = fold != k
    u = trn & (P["w"] > 0)
    med, sc = DM.prep_stats(G[u])
    mv = mask[u] > 0
    tv = tok[u][:, :, :D2.NCONT][mv]
    tmean, tstd = tv.mean(0), tv.std(0)
    init_rng(c, k, a.epochs)
    net = D5.build(c["arch"], G.shape[1], tmean, np.where(tstd > 0, tstd, 1.0), med, sc)
    net, h = DM.fit(net, (tok[trn], mask[trn], G[trn]), P["y"][trn], P["w"][trn], P["m"][trn], P["dmask"][trn], P["dw"][trn], val[trn],
                    c["dec"][0], float(c["lam"]), 1000 * c["seed"] + k, log, device=dev)
    L5.save_fold(os.path.join(a.out, "models", c["tag"]), k, net, G.shape[1], c["arch"])
    info = {**h, "fold": k, "device": dev, "device_name": devname, "threads": a.threads, "epochs": a.epochs, "sec_epoch": round(h["sec"] / a.epochs, 2),
            "sec_load": round(tload, 1), "sec_total": round(time.time() - t0, 1), "n_pool": int(len(P["y"])), "n_train": int(trn.sum()),
            "nparams": DM.nparams(net), "host": socket.gethostname(), "t_start": int(t0), "torch": torch.__version__}
    os.makedirs(os.path.join(a.out, "folds", c["tag"]), exist_ok=True)
    json.dump(info, open(os.path.join(a.out, "folds", c["tag"], f"fold{k}.json"), "w"), indent=1)
    log(f"[{c['tag']}] fold {k} gotowy: {h['sec']:.0f}s trening ({info['sec_epoch']:.1f}s/epokę), {info['sec_total']:.0f}s razem")
    return 0


def cmd_merge(a, log):
    import torch
    import dn_models as DM
    import dn2_models as D2
    c, K = L5.BLOCKS[a.block], a.folds
    torch.set_num_threads(max(1, a.threads))
    dev = a.device
    if dev == "auto":
        dev = "cuda" if torch.cuda.is_available() else "cpu"
    t0 = time.time()
    S, full, P, thr = prepare(c, a, None, log)
    _, models = L5.load_models(a.models, c["tag"], K)
    info = [json.load(open(os.path.join(a.folds_dir, c["tag"], f"fold{k}.json"))) for k in range(K)]
    out = {}
    for s, T in S.items():
        for pre in ("a_", "p_"):
            if pre + "uid" not in T:
                continue
            sc_ = np.full(full[(s, pre)], np.nan, np.float32)
            sc_[T[pre + "idx"]] = L5.score_rows(c, models, T, pre, device=dev)
            out[f"{s}__{pre[0]}"] = sc_
        log(f"  score {s}: " + ", ".join(f"{k} {int(np.isfinite(v).sum())}/{len(v)}" for k, v in out.items() if k.startswith(s + "__")))
    nd = int(DF.cfg["species"][c["sp"]]["n_daughters"])
    ng = 2 if c["inputs"] == "env" else len(DF.H)
    meta = {**c, "features": (["log_pt", "abs_eta"] if c["inputs"] == "env" else DF.H), "token": D2.TOK,
            "n_feat": f"{ng} + {(0 if c['inputs'] == 'env' else nd) + 12}×{len(D2.TOK)}", "n_pool": int(len(P["y"])), "n_sig": int(P["y"].sum()),
            "n_bkg": int((P["role"] == "B").sum()), "n_src": {f"{s}:{r}": int(len(i)) for s, i, r in P["src"]}, "nparams": DM.nparams(models[0]),
            "npar_unit": "wagi", "folds": K, "fit": info, "threads": info[0]["threads"], "device": info[0]["device"], "hp": {**DM.HP, "epochs": info[0]["epochs"]},
            "h_thr": thr, "parallel_folds": True, "device_merge": dev, "device_name": info[0].get("device_name"), "sec": round(sum(i["sec_total"] for i in info) + time.time() - t0, 1),
            "sec_wall_fold_max": max(i["sec_total"] for i in info), "sec_merge": round(time.time() - t0, 1)}
    for d in ("sc", "meta", os.path.join("models", c["tag"])):
        os.makedirs(os.path.join(a.out, d), exist_ok=True)
    np.savez_compressed(os.path.join(a.out, "sc", c["tag"] + ".npz"), **out)
    json.dump(meta, open(os.path.join(a.out, "models", c["tag"], "meta.json"), "w"), indent=1)
    json.dump(meta, open(os.path.join(a.out, "meta", c["tag"] + ".json"), "w"), indent=1)
    log(f"[{c['tag']}] złączone: {K} foldów, najdłuższy fold {meta['sec_wall_fold_max']:.0f}s, łączenie {meta['sec_merge']:.0f}s")
    return 0


# ----------------------------------------------------------------------------------------------------------- testy
def cmd_equiv(a, log):
    """Stara ścieżka (dn4_train.py train: foldy sekwencyjnie w jednym procesie) wobec nowej (fold × K w osobnych procesach + merge)."""
    import torch
    import dn4_train as T4
    import dn5_models as D5
    d = a.dir
    ddir, edirs = os.path.join(d, "data"), {"D0": os.path.join(d, "parts"), "Dplus": os.path.join(d, "parts_dp")}
    T4.synth(ddir, edirs, n=a.n)
    env = {**os.environ, "DN_REDO": "1", "DN_DONE_DIR": "", "OMP_NUM_THREADS": str(a.threads), "MKL_NUM_THREADS": str(a.threads)}
    com = ["--data", ddir, "--env", edirs["D0"], "--env-dp", edirs["Dplus"], "--folds", str(a.folds), "--epochs", str(a.epochs), "--threads", str(a.threads)]
    py = [sys.executable, "-u"]
    run = lambda cmd: subprocess.run(cmd, check=True, env=env, stdout=subprocess.DEVNULL)  # noqa: E731
    old, new = os.path.join(d, "old"), os.path.join(d, "new")
    res = {"folds": a.folds, "epochs": a.epochs, "threads": a.threads, "n": a.n, "torch": torch.__version__, "host": socket.gethostname(), "cmp": {}}
    ok = True
    for b4, b5 in (("d0s1", "base"), ("dps1", "baseDp")):
        tag = L5.BLOCKS[b5]["tag"]
        t = time.time()
        run(py + [os.path.join(L5.DN4, "dn4_train.py"), "train", "--out", old, "--block", b4] + com)
        t_old = time.time() - t
        t = time.time()
        for k in range(a.folds):
            run(py + [os.path.abspath(__file__), "fold", "--out", new, "--block", b5, "--fold", str(k), "--hsc", old] + com)
        run(py + [os.path.abspath(__file__), "merge", "--out", new, "--block", b5, "--hsc", old, "--models", os.path.join(new, "models"),
                  "--folds-dir", os.path.join(new, "folds")] + com)
        sp = L5.BLOCKS[b5]["sp"]
        S, full = L5.load(ddir, edirs[sp], sp)
        zo, zn = np.load(os.path.join(old, "sc", tag + ".npz")), np.load(os.path.join(new, "sc", tag + ".npz"))
        assert sorted(zo.files) == sorted(zn.files), (zo.files, zn.files)
        d1, d0, n1 = 0.0, 0.0, 0
        for kf in zo.files:
            assert np.array_equal(np.isfinite(zo[kf]), np.isfinite(zn[kf])), f"{tag} {kf}: wzór NaN"
            s_, pre = kf.split("__")
            T = S[s_]
            fo = DF.fold_of(T[pre + "_uid"], a.folds)
            dd = np.abs(zo[kf][T[pre + "_idx"]] - zn[kf][T[pre + "_idx"]])
            d1 = max(d1, float(dd[fo > 0].max())); d0 = max(d0, float(dd[fo == 0].max())); n1 += int((fo > 0).sum())
        wd = []
        for k in range(a.folds):
            so = torch.load(os.path.join(old, "models", tag, f"fold{k}.pt"), map_location="cpu")["state"]
            sn = torch.load(os.path.join(new, "models", tag, f"fold{k}.pt"), map_location="cpu")["state"]
            wd.append(max(float((so[x].float() - sn[x].float()).abs().max()) for x in so))
        # powtarzalność: fold 0 nowej ścieżki drugi raz; stara ścieżka drugi raz (jej fold 0 nie jest ustalony ziarnem)
        run(py + [os.path.abspath(__file__), "fold", "--out", new + "2", "--block", b5, "--fold", "0", "--hsc", old] + com)
        sn = torch.load(os.path.join(new, "models", tag, "fold0.pt"), map_location="cpu")["state"]
        s2 = torch.load(os.path.join(new + "2", "models", tag, "fold0.pt"), map_location="cpu")["state"]
        rep_new = max(float((sn[x].float() - s2[x].float()).abs().max()) for x in sn)
        run(py + [os.path.join(L5.DN4, "dn4_train.py"), "train", "--out", old + "2", "--block", b4] + com)
        rep_old = []
        for k in range(a.folds):
            so = torch.load(os.path.join(old, "models", tag, f"fold{k}.pt"), map_location="cpu")["state"]
            s2 = torch.load(os.path.join(old + "2", "models", tag, f"fold{k}.pt"), map_location="cpu")["state"]
            rep_old.append(max(float((so[x].float() - s2[x].float()).abs().max()) for x in so))
        mo, mn = json.load(open(os.path.join(old, "meta", tag + ".json"))), json.load(open(os.path.join(new, "meta", tag + ".json")))
        same_ep = [f["best_epoch"] for f in mo["fit"]][1:] == [f["best_epoch"] for f in mn["fit"]][1:]
        miss = sorted(set(mo) - set(mn))
        res["cmp"][tag] = {"max_abs_dweight_per_fold": wd, "max_abs_dscore_folds_ge1": d1, "rows_folds_ge1": n1, "max_abs_dscore_fold0": d0,
                           "new_fold0_rerun_max_abs_dweight": rep_new, "old_rerun_max_abs_dweight_per_fold": rep_old, "best_epoch_same_folds_ge1": same_ep,
                           "meta_keys_missing": miss, "sec_old": round(t_old, 1), "sec_new": round(time.time() - t, 1)}
        log(f"równoważność {tag}: max|Δwaga| per fold (stara vs nowa) = {wd}; max|Δscore| foldy ≥ 1 = {d1:.3g} na {n1} wierszach (fold 0: {d0:.3g}); "
            f"powtórka foldu 0 nowej ścieżki: {rep_new:.3g}; powtórka starej ścieżki per fold: {rep_old}; epoki wyboru (k ≥ 1) zgodne: {same_ep}; brak kluczy meta: {miss}")
        ok &= max(wd[1:]) <= a.tol and d1 <= a.tol and rep_new <= a.tol and same_ep and not miss
        # modele nowej ścieżki czytane kodem etapu 4 (dn4_lib.load_models) → te same score
        import dn4_lib as L4
        _, ms = L4.load_models(os.path.join(new, "models"), tag)
        T = S["data_428770"]
        v = L4.score_rows("tr", ms, T["a_uid"], T["a_X"], T["a_Dd"], T["a_E"], T["a_eta"], T["a_phi"])
        d4 = float(np.max(np.abs(v - zn["data_428770__a"][T["a_idx"]])))
        res["cmp"][tag]["dn4_lib_read_max_abs"] = d4
        ok &= d4 < 1e-6
    res["equivalent"] = bool(ok)
    # mechanika pozostałych bloków: drzewo H (score do cięcia H) + fold × 2 + merge, 1 epoka; maska par
    for b in ("xgbD0", "xgbDp"):
        run(py + [os.path.join(L5.DN4, "dn4_train.py"), "train", "--out", old, "--block", b] + com)
    com1 = com[:6] + ["--folds", "2", "--epochs", "1", "--threads", str(a.threads)]
    mech = os.path.join(d, "mech")
    for b, c in L5.BLOCKS.items():
        if b.startswith("base"):
            continue
        for k in range(2):
            run(py + [os.path.abspath(__file__), "fold", "--out", mech, "--block", b, "--fold", str(k), "--hsc", old] + com1)
        run(py + [os.path.abspath(__file__), "merge", "--out", mech, "--block", b, "--hsc", old, "--models", os.path.join(mech, "models"),
                  "--folds-dir", os.path.join(mech, "folds")] + com1)
        z = np.load(os.path.join(mech, "sc", c["tag"] + ".npz"))
        assert all(np.isfinite(z[k]).any() and np.nanmin(z[k]) >= 0 and np.nanmax(z[k]) <= 1 for k in z.files), c["tag"]
        m = json.load(open(os.path.join(mech, "meta", c["tag"] + ".json")))
        log(f"mechanika {b} ({c['tag']}): pula {m['n_pool']} (S {m['n_sig']}, B {m['n_bkg']}), wag {m['nparams']}, klucze score {len(z.files)}")
    tok = torch.randn(5, 6, 11); tok[..., 8:] = 0; tok[:, 0, 8] = 1; tok[:, 1, 9] = 1; tok[:, 2:, 10] = 1
    mask = torch.ones(5, 6); mask[:, 5] = 0
    f, pm = D5.pair_feats(tok, mask)
    assert float(pm[:, :2, :2].abs().max()) == 0 and float(pm[:, 5].abs().max()) == 0 and float(pm[:, :, 5].abs().max()) == 0, "maska par: zera"
    assert bool((pm[:, 0, 2:5] == 1).all()) and bool((pm[:, 2, 3] == 1).all()) and float(pm[:, 2, 2].abs().max()) == 0, "maska par: jedynki"
    assert torch.isfinite(f).all()
    log("maska par OK: para córka–córka i tokeny nieważne wyzerowane, cechy skończone")
    res["mechanics_ok"] = True
    print("EQUIV " + json.dumps(res), flush=True)
    if a.json:
        json.dump(res, open(a.json, "w"), indent=1)
    return 0 if ok else 3


def cmd_bench(a, log):
    """Czas jednej epoki DM.fit (wsad 4096, DisCo λ = 5) na losowych tokenach o obsadzeniu jak w danych, przeliczony na pulę etapu 4
    (716 430 wierszy, 4/5 w treningu, z tego 90 % poza walidacją) — w funkcji liczby wątków."""
    import torch
    import dn_models as DM
    import dn5_models as D5
    rng = np.random.default_rng(1)
    n, T = a.n, 14
    tok = rng.normal(size=(n, T, 11)).astype(np.float32); tok[..., 8:] = 0; tok[:, 0, 8] = 1; tok[:, 1, 9] = 1; tok[:, 2:, 10] = 1
    mask = np.zeros((n, T), np.float32); mask[:, :2] = 1
    mask[:, 2:] = (np.arange(12)[None, :] < rng.poisson(2.1, n)[:, None])
    G = rng.normal(size=(n, 18)).astype(np.float32)
    y = rng.random(n) < 0.3
    w = np.ones(n, np.float32); m = rng.uniform(1.75, 2.0, n).astype(np.float32)
    val = rng.random(n) < 0.1
    scale = 716430 * 0.8 / n
    hp0 = dict(DM.HP)
    DM.HP["epochs"] = 1
    res = {"n": n, "host": socket.gethostname(), "ncpu": os.cpu_count(), "load1": os.getloadavg()[0], "torch": torch.__version__, "sec_epoch_full": {}}
    try:
        for arch in ("tr", "trp"):
            for nt in [int(x) for x in a.threads_list.split(",")]:
                torch.set_num_threads(nt)
                best = 1e99
                for _ in range(a.rep):
                    net = D5.build(arch, 18, np.zeros(8, np.float32), np.ones(8, np.float32), np.zeros(18, np.float32), np.ones(18, np.float32))
                    _, h = DM.fit(net, (tok, mask, G), y, w, m, ~y, w, val, "d", 5.0, 1, lambda s: None)
                    best = min(best, h["sec"])
                res["sec_epoch_full"][f"{arch}|{nt}"] = round(best * scale, 1)
                log(f"bench {arch} wątków {nt}: {best:.2f}s na {n} wierszy → {best * scale:.0f}s/epokę pełnej puli")
    finally:
        DM.HP.clear(); DM.HP.update(hp0)
    print("BENCH " + json.dumps(res), flush=True)
    if a.json:
        json.dump(res, open(a.json, "w"), indent=1)
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["fold", "merge", "list", "equiv", "bench"])
    ap.add_argument("--data"); ap.add_argument("--env"); ap.add_argument("--env-dp"); ap.add_argument("--hsc"); ap.add_argument("--out")
    ap.add_argument("--models"); ap.add_argument("--folds-dir"); ap.add_argument("--block"); ap.add_argument("--fold", type=int, default=0)
    ap.add_argument("--folds", type=int, default=int(os.environ.get("DN_FOLDS", "5")))
    ap.add_argument("--cap-sig", type=int, default=int(os.environ.get("DN5_CAP_SIG", "300000")))
    ap.add_argument("--cap-bkg", type=int, default=int(os.environ.get("DN5_CAP_BKG", "500000")))
    ap.add_argument("--epochs", type=int, default=int(os.environ.get("DN5_EPOCHS", "20")))
    ap.add_argument("--threads", type=int, default=int(os.environ.get("DN_CPUS", "4")))
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--dir", default=""); ap.add_argument("--n", type=int, default=0); ap.add_argument("--tol", type=float, default=0.0)
    ap.add_argument("--eq-folds", type=int, default=3); ap.add_argument("--eq-epochs", type=int, default=2)
    ap.add_argument("--threads-list", default="1,2,4,8"); ap.add_argument("--rep", type=int, default=2); ap.add_argument("--json", default="")
    a = ap.parse_args(argv)
    if a.cmd == "list":
        for b, c in L5.BLOCKS.items():
            print(b, c["tag"])
        return 0
    t0 = time.time()
    log = lambda s: print(f"[{time.time() - t0:6.0f}s] {s}", flush=True)  # noqa: E731
    if a.cmd == "equiv":
        a.n = a.n or 4000
        if not a.dir:
            raise SystemExit("--dir")
        a.folds, a.epochs = a.eq_folds, a.eq_epochs
        return cmd_equiv(a, log)
    if a.cmd == "bench":
        a.n = a.n or 40960
        return cmd_bench(a, log)
    if a.block not in L5.BLOCKS:
        raise SystemExit(f"nieznany blok {a.block}; są: {' '.join(L5.BLOCKS)}")
    return cmd_fold(a, log) if a.cmd == "fold" else cmd_merge(a, log)


if __name__ == "__main__":
    sys.exit(main())
