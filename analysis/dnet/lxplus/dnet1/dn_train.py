#!/usr/bin/env python3
"""DNET etap 1 — sklejanie części i trening bloku konfiguracji (dn_feat.BLOCKS) na TYCH SAMYCH zdarzeniach, foldach i wagach.

  dn_train.py merge --parts DIR --out DIR                     części dn_<próbka>_<g>.npz → dn_<próbka>.npz
  dn_train.py train --data DIR --out DIR --block NAZWA        dla każdej konfiguracji bloku: K foldów po zdarzeniach (hash event_uid),
      score OOF dla wszystkich tabel (a_ inkluzywne, p_ pary D0←D*) wszystkich próbek → sc/<tag>.npz, meta/<tag>.json, onnx/<tag>.onnx

Próbka treningowa (wspólna dla drzew i sieci): sygnał = fd truth 0/1 (D0: hypo_swap ≠ 1; D+: weto D*), cały zakres mas;
tło „mc” = mb truth ≥ 8; tło „sb” = dalekie SB danych 428770 (D0: > 100 MeV + weto refleksji; D+: > 110 MeV + weto D*);
„sbmc” = BCE na SB danych, człon dekorelacji na tle MC w całym zakresie mas (SB nie ma zdarzeń pod pikiem — sieć z cechami
rekonstruującymi masę nauczyłaby się okna masy bez kary). Wagi: tło przeważone do widma (pT, |η|) sygnału (charmbkg reweight_bins),
klasy zrównoważone. Limity po zdarzeniach (hash): DN_CAP_SIG, DN_CAP_BKG."""
import argparse
import glob
import json
import os
import re
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dn_feat as DF  # noqa: E402

SAMPLES = ["fd_nd", "mb_nd", "data_428770", "data_427929"]


def merge(pdir, out):
    by = {}
    for f in sorted(glob.glob(os.path.join(pdir, "dn_*_*.npz"))):
        m = re.match(r"dn_(\w+?)_(\d+)\.npz$", os.path.basename(f))
        if m:
            by.setdefault(m.group(1), []).append((int(m.group(2)), f))
    os.makedirs(out, exist_ok=True)
    for s, lst in by.items():
        cols, metas = {}, []
        for g, f in sorted(lst):
            z = np.load(f, allow_pickle=False)
            metas.append({"group": g, **json.loads(str(z["meta"]))})
            for k in z.files:
                if k != "meta":
                    cols.setdefault(k, []).append(z[k])
        o = {k: np.concatenate(v) for k, v in cols.items() if len(v) == len(lst)}
        o["meta"] = np.array(json.dumps({"sample": s, "groups": metas, "feats": metas[0]["feats"], "dau": metas[0]["dau"],
                                         "frac": metas[0]["frac"], "lxy_min": metas[0]["lxy_min"]}))
        np.savez(os.path.join(out, f"dn_{s}.npz"), **o)
        print(f"{s}: {len(lst)} części; " + ", ".join(f"{k}={len(v)}" for k, v in o.items() if k.endswith("_m")), flush=True)
    return sorted(by)


def load_species(ddir, sp, keys=None):
    """{próbka: {kolumna: tablica}} dla gatunku (klucze bez prefiksu gatunku; a_*/p_*)."""
    S = {}
    for s in SAMPLES:
        p = os.path.join(ddir, f"dn_{s}.npz")
        if not os.path.isfile(p):
            continue
        z = np.load(p, allow_pickle=False)
        S[s] = {k.split("__", 1)[1]: z[k] for k in z.files if k.startswith(sp + "__") and (keys is None or k.split("__", 1)[1][2:] in keys)}
    return S


def masks(S, sp):
    """Maski wierszy tabel a_: sygnał (fd), tło MC (mb), tło SB (dane 428770)."""
    m0 = DF.cfg["species"][sp]["m0"]
    fd, mb, da = S["fd_nd"], S["mb_nd"], S.get("data_428770")
    sig = np.isin(fd["a_truth"], DF.SIGNAL_CATS) & fd["a_ok"]
    if sp == "D0":
        sig &= fd["a_swap"] != 1
    bmc = (mb["a_truth"] >= DF.BKG_MIN) & (mb["a_truth"] < 100) & mb["a_ok"]
    bsb = None
    if da is not None:
        bsb = da["a_ok"] & (np.abs(da["a_m"] - m0) > DF.FAR[sp])
        if sp == "D0":
            bsb &= np.abs(da["a_m_swap"] - m0) > DF.REFL_VETO
    return sig, bmc, bsb


def reweight(ptb, etab, wb, pts, etas, ws):
    bins = DF.cfg["training"]["reweight_bins"]
    bp, be = np.asarray(bins["pt"]), np.asarray(bins["abs_eta"])
    hs, _, _ = np.histogram2d(pts, np.abs(etas), [bp, be], weights=ws)
    hb, _, _ = np.histogram2d(ptb, np.abs(etab), [bp, be], weights=wb)
    hs /= max(hs.sum(), 1e-12); hb /= max(hb.sum(), 1e-12)
    r = np.divide(hs, hb, out=np.zeros_like(hs), where=hb > 0)
    ip = np.clip(np.digitize(ptb, bp) - 1, 0, len(bp) - 2)
    ie = np.clip(np.digitize(np.abs(etab), be) - 1, 0, len(be) - 2)
    w = wb * r[ip, ie]
    return w * (ws.sum() / max(w.sum(), 1e-12))


def pool(S, sp, bkg, cap_sig, cap_bkg):
    """Pula treningowa: dict tablic (src = (próbka, indeksy a_)), y, w (BCE), dmask, dw, mass, uid — wspólna dla wszystkich modeli."""
    sig, bmc, bsb = masks(S, sp)
    if bkg in ("sb", "sbmc") and bsb is None:
        raise SystemExit("brak data_428770 dla tła sb")

    def take(s, mk, cap, salt):
        i = np.flatnonzero(mk)
        if cap and len(i) > cap:
            i = i[DF.unit_hash(S[s]["a_uid"][i], salt) < cap / len(i)]
        return i
    src = [("fd_nd", take("fd_nd", sig, cap_sig, 11), "S")]
    if bkg == "mc":
        src.append(("mb_nd", take("mb_nd", bmc, cap_bkg, 12), "B"))
    else:
        src.append(("data_428770", take("data_428770", bsb, cap_bkg, 13), "B"))
        if bkg == "sbmc":
            src.append(("mb_nd", take("mb_nd", bmc, cap_bkg // 2, 12), "D"))
    g = lambda k: np.concatenate([S[s]["a_" + k][i] for s, i, _ in src])  # noqa: E731
    role = np.concatenate([np.full(len(i), r) for _, i, r in src])
    P = {"src": src, "role": role, "y": role == "S", "m": g("m").astype(np.float32), "uid": g("uid"), "pt": g("pt"), "eta": g("eta")}
    w0 = np.where(role == "S", g("w"), 1.0).astype(np.float64)
    w = np.zeros(len(role))
    s_ = role == "S"
    w[s_] = w0[s_]
    for r in ("B", "D"):
        b_ = role == r
        if b_.any():
            w[b_] = reweight(P["pt"][b_], P["eta"][b_], w0[b_], P["pt"][s_], P["eta"][s_], w0[s_])
    w /= w[s_ | (role == "B")].mean()
    drole = "D" if bkg == "sbmc" else "B"
    P["dmask"] = role == drole
    P["dw"] = np.where(P["dmask"], w, 0.0).astype(np.float32)
    P["w"] = np.where(role == "D", 0.0, w).astype(np.float32)
    return P


def pool_x(S, P, key="X", cols=None):
    X = np.concatenate([S[s]["a_" + key][i] for s, i, _ in P["src"]])
    return X[:, cols] if cols is not None else X


def train_xgb(X, y, w, val, seed, threads):
    import xgboost as xgb
    p = dict(DF.cfg["training"]["xgb"])
    es = p.pop("early_stopping_rounds"); p.pop("n_jobs")
    if os.environ.get("DN_FAST") == "1":
        p.update(n_estimators=40, learning_rate=0.3)
    X = np.where(np.isfinite(X), X, np.nan)
    m = xgb.XGBClassifier(objective="binary:logistic", eval_metric="logloss", early_stopping_rounds=es, n_jobs=threads, random_state=seed, **p)
    m.fit(X[~val], y[~val], sample_weight=w[~val], eval_set=[(X[val], y[val])], sample_weight_eval_set=[w[val]], verbose=False)
    return m


def run_config(c, S, a, log):
    import dn_models as DM
    sp, K = c["sp"], a.folds
    feats = DF.SETS[sp][c["feats"]]
    cols = [DF.ALLF[sp].index(f) for f in feats]
    P = pool(S, sp, c["bkg"], a.cap_sig, a.cap_bkg)
    X = pool_x(S, P, "X", cols)
    fold = DF.fold_of(P["uid"], K)
    val = DF.unit_hash(P["uid"], salt=23) >= 0.9
    dec, lam = c["dec"][0], float(c["dec"][1:] or 0)
    t0 = time.time()
    log(f"[{c['tag']}] pula {len(P['y'])} (S {int(P['y'].sum())}, B {int((P['role'] == 'B').sum())}, D {int((P['role'] == 'D').sum())}), "
        f"{len(feats)} cech, dek. {dec} λ={lam}")
    models, info = [], []
    for k in list(range(K)) + (["full"] if c["full"] else []):
        trn = (fold != k) if k != "full" else np.ones(len(fold), bool)
        seed = 1000 * c["seed"] + (K if k == "full" else k)
        if c["model"] == "xgb":
            ub = trn & (P["w"] > 0)
            tm = time.time()
            m = train_xgb(X[ub], P["y"][ub], P["w"][ub], val[ub], seed, a.threads)
            models.append(("xgb", m)); info.append({"best_it": int(m.best_iteration), "sec": round(time.time() - tm, 1)})
            log(f"  fold {k}: best_it {m.best_iteration}, {time.time() - tm:.0f}s")
        else:
            med, sc = DM.prep_stats(X[trn & (P["w"] > 0)])
            net = DM.MLP(len(feats), c["model"], med, sc)
            net, h = DM.fit(net, (X[trn],), P["y"][trn], P["w"][trn], P["m"][trn], P["dmask"][trn], P["dw"][trn], val[trn], dec, lam, seed, log,
                            device=a.device)
            models.append(("nn", net)); info.append(h)
    npar = DM.nparams(models[0][1]) if c["model"] != "xgb" else int(sum(len(m.get_booster().get_dump()) for _, m in models[:1]))

    def pred(mod, Xv):
        kind, m = mod
        if kind == "xgb":
            return m.predict_proba(np.where(np.isfinite(Xv), Xv, np.nan))[:, 1].astype(np.float32)
        return DM.predict(m, (Xv,))
    out = {}
    for s, T in S.items():
        for pre in ("a_", "p_"):
            if pre + "uid" not in T:
                continue
            f = DF.fold_of(T[pre + "uid"], K)
            Xv = T[pre + "X"][:, cols]
            sc_ = np.full(len(f), np.nan, np.float32)
            for k in range(K):
                mk = f == k
                if mk.any():
                    sc_[mk] = pred(models[k], Xv[mk])
            out[f"{s}__{pre[0]}"] = sc_
    os.makedirs(os.path.join(a.out, "sc"), exist_ok=True); os.makedirs(os.path.join(a.out, "meta"), exist_ok=True)
    np.savez_compressed(os.path.join(a.out, "sc", c["tag"] + ".npz"), **out)
    meta = {**c, "features": feats, "n_feat": len(feats), "n_pool": int(len(P["y"])), "n_sig": int(P["y"].sum()), "n_bkg": int((P["role"] == "B").sum()),
            "n_dec": int(P["dmask"].sum()), "nparams": npar, "npar_unit": "wagi" if c["model"] != "xgb" else "drzewa (fold 0)", "folds": K,
            "fit": info, "sec": round(time.time() - t0, 1), "threads": a.threads, "device": a.device,
            "hp": (DM.HP if c["model"] != "xgb" else {k: v for k, v in DF.cfg["training"]["xgb"].items()})}
    if c["full"] and c["model"] != "xgb":
        os.makedirs(os.path.join(a.out, "onnx"), exist_ok=True)
        meta["onnx"] = DM.export_onnx(models[K][1], (X[:256],), os.path.join(a.out, "onnx", c["tag"] + ".onnx"), ["features"])
        meta["onnx"]["inputs"] = feats
        log(f"  ONNX: {meta['onnx']}")
    json.dump(meta, open(os.path.join(a.out, "meta", c["tag"] + ".json"), "w"), indent=1)
    log(f"[{c['tag']}] gotowe, {time.time() - t0:.0f}s, parametrów {npar}")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["merge", "train", "list"])
    ap.add_argument("--parts"); ap.add_argument("--data"); ap.add_argument("--out")
    ap.add_argument("--block"); ap.add_argument("--only", default="")
    ap.add_argument("--folds", type=int, default=int(os.environ.get("DN_FOLDS", "5")))
    ap.add_argument("--cap-sig", type=int, default=int(os.environ.get("DN_CAP_SIG", "400000")))
    ap.add_argument("--cap-bkg", type=int, default=int(os.environ.get("DN_CAP_BKG", "800000")))
    ap.add_argument("--epochs", type=int, default=int(os.environ.get("DN_EPOCHS", "30")))
    ap.add_argument("--threads", type=int, default=int(os.environ.get("DN_CPUS", "4")))
    ap.add_argument("--device", default="cpu")
    a = ap.parse_args(argv)
    if a.cmd == "list":
        for b, cs in DF.BLOCKS.items():
            print(b, len(cs), " ".join(c["tag"] for c in cs))
        return 0
    if a.cmd == "merge":
        merge(a.parts, a.out)
        return 0
    t0 = time.time()
    log = lambda s: print(f"[{time.time() - t0:6.0f}s] {s}", flush=True)  # noqa: E731
    cfgs = [c for c in DF.BLOCKS[a.block] if not a.only or c["tag"] in a.only.split(",")]
    if any(c["model"] != "xgb" for c in cfgs):
        import torch
        import dn_models as DM
        torch.set_num_threads(max(1, a.threads))
        DM.HP["epochs"] = a.epochs
    cache = {}
    for c in cfgs:
        if os.path.isfile(os.path.join(os.environ.get("DN_DONE_DIR") or a.out, "meta", c["tag"] + ".json")) and os.environ.get("DN_REDO") != "1":
            log(f"[{c['tag']}] już jest — pomijam")
            continue
        if c["sp"] not in cache:
            cache.clear()
            cache[c["sp"]] = load_species(a.data, c["sp"])
            log(f"wczytano {c['sp']}: " + ", ".join(f"{s} a={len(T.get('a_m', []))} p={len(T.get('p_m', []))}" for s, T in cache[c["sp"]].items()))
        run_config(c, cache[c["sp"]], a, log)
    return 0


if __name__ == "__main__":
    sys.exit(main())
