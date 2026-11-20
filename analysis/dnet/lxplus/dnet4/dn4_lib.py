"""DNET etap 4 — kod wspólny treningu (dn4_train.py) i zastosowania (dn4_apply.py): zapis/odczyt modeli foldów i liczenie score.
Reguła score (wszędzie ta sama): wiersz dostaje wynik modelu foldu fold_of(event_uid, K) — dla danych użytych jako tło treningowe
(pasma boczne; drzewo H „sb”, sieć „sb2mc”) jest to wynik spoza treningu (podział po zdarzeniach), dla reszty — dowolny model foldu.
Kod etapów 1–3 (../dnet1, ../dnet2, ../dnet3) NIE jest zmieniany."""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
_Q = os.path.dirname(HERE)
DN2 = os.environ.get("DN2") or os.path.join(_Q, "dnet2")
DN1 = os.environ.get("DN1") or os.path.join(_Q, "dnet1")
for _p in (DN1, DN2, HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import dn_feat as DF  # noqa: E402
import dn2_train as D2T  # noqa: E402

ENV = "e04"                                    # otoczenie: tory ΔR < 0,4 spośród 12 o największym pT w ΔR < 1,0 (jak etap 2)


def to_lg16(s):
    """score ∈ (0, 1) → logit float16 (zapis zwarty; rozdzielczość względna ~5·10⁻⁴ także przy score → 1)."""
    s = np.clip(np.asarray(s, np.float64), 1e-7, 1 - 1e-7)
    return np.log(s / (1 - s)).astype(np.float16)


def save_tr(mdir, k, net, nglob, arch):
    import torch
    os.makedirs(mdir, exist_ok=True)
    torch.save({"state": net.state_dict(), "nglob": int(nglob), "arch": arch}, os.path.join(mdir, f"fold{k}.pt"))


def load_models(mroot, tag):
    """→ (meta, [model foldu 0..K−1]); meta = models/<tag>/meta.json (kopia meta/<tag>.json)."""
    mdir = os.path.join(mroot, tag)
    meta = json.load(open(os.path.join(mdir, "meta.json")))
    K = int(meta["folds"])
    if meta["model"] == "xgb":
        import xgboost as xgb
        ms = []
        for k in range(K):
            m = xgb.XGBClassifier()
            m.load_model(os.path.join(mdir, f"fold{k}.json"))
            ms.append(m)
        return meta, ms
    import torch
    import dn2_models as D2
    ms = []
    for k in range(K):
        d = torch.load(os.path.join(mdir, f"fold{k}.pt"), map_location="cpu")
        ng = d["nglob"]
        net = D2.SetNet(ng, d["arch"], np.zeros(D2.NCONT, np.float32), np.ones(D2.NCONT, np.float32), np.zeros(ng, np.float32), np.ones(ng, np.float32))
        net.load_state_dict(d["state"])
        ms.append(net.eval())
    return meta, ms


def score_rows(kind, models, uid, X, Dd=None, E=None, eta=None, phi=None, chunk=200_000, device="cpu"):
    """Score float32 [n] wg reguły foldów. kind: 'xgb' (X = cechy modelu) | 'tr' (X = 18 cech H, Dd, E, eta, phi → tokeny e04)."""
    K = len(models)
    f = DF.fold_of(uid, K)
    out = np.full(len(f), np.nan, np.float32)
    if kind == "xgb":
        for k in range(K):
            w = f == k
            if w.any():
                Xv = X[w]
                out[w] = models[k].predict_proba(np.where(np.isfinite(Xv), Xv, np.nan))[:, 1].astype(np.float32)
        return out
    import dn_models as DM
    for a0 in range(0, len(f), chunk):
        sl = slice(a0, a0 + chunk)
        tk, mk = D2T.tokens(Dd[sl], E[sl], eta[sl], phi[sl], ENV)
        fs, Xs = f[sl], X[sl]
        sub = np.full(len(fs), np.nan, np.float32)
        for k in range(K):
            w = fs == k
            if w.any():
                sub[w] = DM.predict(models[k], (tk[w], mk[w], Xs[w]), device=device)
        out[sl] = sub
    return out
