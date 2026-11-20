"""DNET etap 5a — kod wspólny: konfiguracje wariantów, pula treningowa, wejścia modeli, zapis/odczyt modeli foldów, score OOF.
Reguła score jak w etapie 4 (dn4_lib): wiersz dostaje wynik modelu foldu fold_of(event_uid, K). Format models/<tag>/fold<k>.pt jak w etapie 4
({"state", "nglob", "arch"}) → warianty „tr” na pełnych wejściach (krzywa uczenia) czyta bez zmian dn4_lib.load_models / dn4_apply.py.
Kod etapów 1–4 NIE jest zmieniany (import)."""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
_Q = os.path.dirname(HERE)
DN4 = os.environ.get("DN4") or os.path.join(_Q, "dnet4")
DN2 = os.environ.get("DN2") or os.path.join(_Q, "dnet2")
DN1 = os.environ.get("DN1") or os.path.join(_Q, "dnet1")
for _p in (DN1, DN2, DN4, HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import dn_feat as DF  # noqa: E402
import dn_train as DT  # noqa: E402
import dn2_train as D2T  # noqa: E402

ENV = "e04"
LAM = DF.LAM[0]
BINS = [3, 5, 8, 30]
SALT_SFRAC = 31                                 # hash podpróbki sygnału (krzywa uczenia); 25 % ⊂ 50 % ⊂ 100 %


def _v(sp, model="tr", inputs="full", sfrac=1.0, hcut=0.0, seed=1):
    """tag = <gatunek>.<model>.<cechy>.<tło>.<dek>.s<ziarno>. inputs: full = 18 cech H + córki + otoczenie e04 (jak etap 4);
    env = WYŁĄCZNIE tokeny otoczenia e04 + (log pT, |η|) kandydata (tagger E; bez córek, bez cech wierzchołka, bez DisCo)."""
    env = inputs == "env"
    feats = (f"E04h{int(round(100 * hcut))}" if hcut else "E04") if env else f"H+{ENV}"
    bkg = "mc" if sfrac >= 1 else f"mcS{int(round(100 * sfrac))}"
    dec, lam = ("n", 0.0) if env else (f"d{LAM}", float(LAM))
    return {"sp": sp, "model": model, "feats": feats, "env": ENV, "bkg": bkg, "dec": dec, "seed": seed, "full": False, "lam": lam, "inputs": inputs,
            "sfrac": float(sfrac), "hcut": float(hcut), "arch": "trp" if model == "trp" else "tr", "tag": f"{sp}.{model}.{feats}.{bkg}.{dec}.s{seed}"}


BLOCKS = {
    "base": _v("D0"),                              # ≡ blok d0s1 etapu 4 (D0.tr.H+e04.mc.d5.s1): test równoważności i sonda GPU; nie w biegu domyślnym
    "baseDp": _v("Dplus"),                         # ≡ dps1 etapu 4: tylko test równoważności
    "lc25": _v("D0", sfrac=0.25), "lc50": _v("D0", sfrac=0.50),          # krzywa uczenia (pozycja #5 planu)
    "pair": _v("D0", model="trp"),                                         # cechy par w uwadze (#4)
    "eD0": _v("D0", model="trE", inputs="env"), "eD0h": _v("D0", model="trE", inputs="env", hcut=0.70),      # tagger otoczenia E (#3)
    "eDp": _v("Dplus", model="trE", inputs="env"), "eDph": _v("Dplus", model="trE", inputs="env", hcut=0.70),
}


def load(ddir, edir, sp, samples=None):
    """dn2_train.load z ograniczeniem listy próbek (job foldu potrzebuje tylko fd_nd i mb_nd)."""
    old = list(DT.SAMPLES)
    try:
        if samples:
            DT.SAMPLES[:] = [s for s in old if s in samples]
        return D2T.load(ddir, edir, sp)
    finally:
        DT.SAMPLES[:] = old


def hpass(S, sp, hdir, wp):
    """Cięcie drzewem H (odniesienie etapu 1, score OOF) na punkcie pracy wp: progi = kwantyle score sygnału MC w oknie sygnału per bin pT
    (jak analyze.py; tylko MC). → ({próbka: maska wierszy a_}, progi)."""
    z = np.load(os.path.join(hdir, "sc", DF.REF[sp] + ".npz"))
    scfg = DF.cfg["species"][sp]
    sig, _, _ = DT.masks(S, sp)
    fd = S["fd_nd"]
    hs = {s: z[f"{s}__a"][T["a_idx"]] for s, T in S.items() if "a_idx" in T and f"{s}__a" in z.files}
    w = sig & (np.abs(fd["a_m"] - scfg["m0"]) < scfg["signal_window"])
    ib = np.digitize(fd["a_pt"], BINS) - 1
    thr = [float(np.quantile(hs["fd_nd"][w & (ib == b)], 1 - wp)) for b in range(3)]
    out = {}
    for s, v in hs.items():
        b = np.clip(np.digitize(S[s]["a_pt"], BINS) - 1, 0, 2)
        out[s] = v >= np.asarray(thr)[b]
    return out, thr, hs


def pool(S, sp, c, cap_sig, cap_bkg, hp=None):
    """Pula treningowa. Bez podpróbki i bez cięcia H → dokładnie dn_train.pool(…, „mc”) (ta sama co w etapach 2–4). Inaczej ten sam przepis
    (limity po zdarzeniach, tło przeważone do widma (pT, |η|) sygnału, klasy zrównoważone) na zawężonych maskach."""
    if c["sfrac"] >= 1 and not c["hcut"]:
        return DT.pool(S, sp, "mc", cap_sig, cap_bkg)
    sig, bmc, _ = DT.masks(S, sp)
    if c["hcut"]:
        sig, bmc = sig & hp["fd_nd"], bmc & hp["mb_nd"]

    def take(s, mk, cap, salt):
        i = np.flatnonzero(mk)
        if cap and len(i) > cap:
            i = i[DF.unit_hash(S[s]["a_uid"][i], salt) < cap / len(i)]
        return i
    i_s = take("fd_nd", sig, cap_sig, 11)
    if c["sfrac"] < 1:
        i_s = i_s[DF.unit_hash(S["fd_nd"]["a_uid"][i_s], SALT_SFRAC) < c["sfrac"]]
    src = [("fd_nd", i_s, "S"), ("mb_nd", take("mb_nd", bmc, cap_bkg, 12), "B")]
    g = lambda k: np.concatenate([S[s]["a_" + k][i] for s, i, _ in src])  # noqa: E731
    role = np.concatenate([np.full(len(i), r) for _, i, r in src])
    P = {"src": src, "role": role, "y": role == "S", "m": g("m").astype(np.float32), "uid": g("uid"), "pt": g("pt"), "eta": g("eta")}
    w0 = np.where(role == "S", g("w"), 1.0).astype(np.float64)
    w = np.zeros(len(role))
    s_, b_ = role == "S", role == "B"
    w[s_] = w0[s_]
    w[b_] = DT.reweight(P["pt"][b_], P["eta"][b_], w0[b_], P["pt"][s_], P["eta"][s_], w0[s_])
    w /= w.mean()
    P["dmask"] = b_
    P["dw"] = np.where(b_, w, 0.0).astype(np.float32)
    P["w"] = w.astype(np.float32)
    return P


def inputs(c, X, Dd, E, eta, phi, pt):
    """→ (tok, mask, glob) dla modelu wariantu c."""
    tok, mask = D2T.tokens(Dd, E, eta, phi, ENV)
    if c["inputs"] == "env":
        nd = Dd.shape[1]
        G = np.stack([np.log(np.clip(pt, 0.5, None)), np.abs(eta)], 1).astype(np.float32)
        return tok[:, nd:], mask[:, nd:], G
    return tok, mask, X


def save_fold(mdir, k, net, nglob, arch):
    import torch
    os.makedirs(mdir, exist_ok=True)
    torch.save({"state": net.state_dict(), "nglob": int(nglob), "arch": arch}, os.path.join(mdir, f"fold{k}.pt"))


def load_fold(path):
    import torch
    import dn2_models as D2
    import dn5_models as D5
    d = torch.load(path, map_location="cpu")
    ng = d["nglob"]
    net = D5.build(d["arch"], ng, np.zeros(D2.NCONT, np.float32), np.ones(D2.NCONT, np.float32), np.zeros(ng, np.float32), np.ones(ng, np.float32))
    net.load_state_dict(d["state"])
    return net.eval()


def load_models(mroot, tag, K=None):
    mdir = os.path.join(mroot, tag)
    mp = os.path.join(mdir, "meta.json")
    meta = json.load(open(mp)) if os.path.isfile(mp) else {}
    K = int(K or meta["folds"])
    return meta, [load_fold(os.path.join(mdir, f"fold{k}.pt")) for k in range(K)]


def score_rows(c, models, T, pre, chunk=200_000, device="cpu"):
    """Score float32 [n] wierszy tabeli T (prefiks a_/p_) wg reguły foldów."""
    import dn_models as DM
    K = len(models)
    f = DF.fold_of(T[pre + "uid"], K)
    out = np.full(len(f), np.nan, np.float32)
    pch = 20_000 if c["arch"] == "trp" else 200_000                  # trp: tensor par [B, T, T, 32] — mniejsze porcje
    for a0 in range(0, len(f), chunk):
        sl = slice(a0, a0 + chunk)
        tk, mk, G = inputs(c, T[pre + "X"][sl], T[pre + "Dd"][sl], T[pre + "E"][sl], T[pre + "eta"][sl], T[pre + "phi"][sl], T[pre + "pt"][sl])
        fs = f[sl]
        sub = np.full(len(fs), np.nan, np.float32)
        for k in range(K):
            w = fs == k
            if w.any():
                sub[w] = DM.predict(models[k], (tk[w], mk[w], G[w]), chunk=pch, device=device)
        out[sl] = sub
    return out
