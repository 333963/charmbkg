#!/usr/bin/env python3
"""DNET etap 3 — skąd bierze się niedoszacowanie efektywności sygnału (SF 0,90 w runie jetowym) transformera `tr.H+e04` z etapu 2:
ablacje wejść otoczenia i dekorelacja od „aktywności”. Te same wiersze, pula, foldy, wagi, model (dn2_models.SetNet „tr”) i człon DisCo
masy (λ = DN_LAM[0]) co w etapie 2; dane = tabele etapu 1 + tory otoczenia etapu 2 (bez ponownej ekstrakcji).

  dn3_train.py train --data DIR1/data --env DIR2/parts --out DIR --block NAZWA [--device auto]
  dn3_train.py list
  dn3_train.py selftest [--need-torch] [--dir KATALOG]     dane syntetyczne: logika tokenów (numpy) + pełna mechanika (torch, 1 epoka, CPU)

Co JEST zapisane dla toru otoczenia (dn2_extract.E): pT, Δη, Δφ względem D, ładunek, dE/dx, npts — NIE MA d0/z0 (drzewo Tracks ich nie ma),
więc podziału „tor z wierzchołka pierwotnego / odsunięty” nie da się zrobić. Aktywność: nE[:, 0] = liczba torów otoczenia w ΔR < 0,4
(wszystkie, bez obcięcia), pc40 = Σ pT torów otoczenia (pT > 0,5 GeV, ΔR < 0,4).

Warianty tokenów otoczenia (VAR; córki zawsze bez zmian, 8 cech + flagi):
  e04     jak etap 2: tory otoczenia z ΔR < 0,4 spośród 12 o największym pT w ΔR < 1,0 (tak tnie dn2_train.tokens)
  e04g    A1: tokeny otoczenia bez cech pędowych — log pT, dE/dx (∝ pęd przez Bethego–Blocha) i npts wyzerowane; zostaje Δη, Δφ, ładunek względem K
  e04np   (poza domyślną listą) tylko log pT wyzerowany, dE/dx i npts zostają
  e04pt1  A4: tory otoczenia z pT > 1 GeV
  e04c    (poza domyślną listą) ΔR < 0,4 wybierane ze wszystkich 20 zapisanych torów, potem 12 o największym pT (bez zależności od ΔR < 1,0)
Dekorelacja od aktywności (dn3_models.fit3): + λ_act · ½[dCor²(score, log(1 + n)) + dCor²(score, log(1 + Σ pT))] na sygnale MC (aS)
albo na tle MC (aB). Tag: D0.tr.H+<wariant>.mc.d<λ>[a<S|B><λ_act>].s<ziarno>."""
import argparse
import glob
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
_Q = os.path.dirname(HERE)
DN2 = os.environ.get("DN2") or os.path.join(_Q, "dnet2")
DN1 = os.environ.get("DN1") or os.path.join(_Q, "dnet1")
for _p in (DN1, DN2, HERE):
    sys.path.insert(0, _p)
import dn_feat as DF  # noqa: E402
import dn_train as DT  # noqa: E402
import dn2_train as D2T  # noqa: E402

LAM = DF.LAM[0]
TOK = ["log_pt", "deta", "dphi", "q_rel", "dedx", "npts", "d0", "z0", "is_K", "is_pi", "is_env"]   # = dn2_models.TOK (sprawdzane przy treningu)
NCONT = 8                                                                                          # = dn2_models.NCONT
VAR = {
    "e04": {},
    "e04g": {"drop": ("log_pt", "dedx", "npts")},
    "e04np": {"drop": ("log_pt",)},
    "e04pt1": {"ptmin": 1.0},
    "e04c": {"compact": True},
}
ACT = ["log1p_n_env04", "log1p_sumpt_env04"]


def _c(env="e04", act="", lam_act=0.0, rep=False, seed=1):
    dec = f"d{LAM}" + (f"a{act}{lam_act:g}" if act else "") + ("rep" if rep else "")
    return {"sp": "D0", "model": "tr", "feats": f"H+{env}", "env": env, "bkg": "mc", "dec": dec, "seed": seed, "full": False, "lam": float(LAM),
            "act": act, "lam_act": float(lam_act), "fit3": bool(act or rep), "tag": f"D0.tr.H+{env}.mc.{dec}.s{seed}"}


BLOCKS = {
    "tr04g": [_c("e04g")],                       # A1
    "tr04pt1": [_c("e04pt1")],                   # A4
    "actS2": [_c(act="S", lam_act=2)],           # B: DisCo aktywności na sygnale MC
    "actS10": [_c(act="S", lam_act=10)],
    "actB10": [_c(act="B", lam_act=10)],         # B: DisCo aktywności na tle MC
    # poza domyślną listą DN3_BLOCKS:
    "tr04np": [_c("e04np")],
    "tr04c": [_c("e04c")],
    "ref": [_c(rep=True)],                       # C: powtórka tr.H+e04.mc.d5.s1 przez fit3 (λ_act = 0), tag …d5rep.s1
}


# ----------------------------------------------------------------------------------------------------------- tokeny (numpy)
def compact(E, rmax=0.4):
    """Tory z pT > 0 i ΔR < rmax na początek (kolejność po pT zachowana) — wybór ΔR < 0,4 ze wszystkich zapisanych torów."""
    Ef = np.nan_to_num(E.astype(np.float32), nan=0.0, posinf=0.0, neginf=0.0)
    ok = (Ef[:, :, 0] > 0) & (np.hypot(Ef[:, :, 1], Ef[:, :, 2]) < rmax)
    o = np.argsort(~ok, axis=1, kind="stable")
    return np.take_along_axis(E, o[:, :, None], 1)


def tokens3(Dd, E, eta, phi, env):
    """→ tok [n, nd + 12, 11], mask [n, nd + 12] jak dn2_train.tokens(…, "e04"), z filtrem/ablacją wariantu `env` (VAR)."""
    v = VAR[env]
    if v.get("compact"):
        E = compact(E)
    tok, mask = D2T.tokens(Dd, E, eta, phi, "e04")
    nd = Dd.shape[1]
    ne = tok.shape[1] - nd
    if "ptmin" in v:
        pt = np.nan_to_num(E[:, :ne, 0].astype(np.float32), nan=0.0, posinf=0.0, neginf=0.0)
        keep = (pt > v["ptmin"]).astype(np.float32)
        mask[:, nd:] *= keep
        tok[:, nd:] *= keep[:, :, None]
    for c in v.get("drop", ()):
        tok[:, nd:, TOK.index(c)] = 0.0
    return tok, mask


def tok_stats(tok, mask, nd, env):
    """Stałe standaryzacji tokenów jak w etapie 2 (po ważnych tokenach); dla cech wyzerowanych w otoczeniu — tylko z córek."""
    mv = mask > 0
    tv = tok[:, :, :NCONT][mv]
    tmean, tstd = tv.mean(0), tv.std(0)
    drop = VAR[env].get("drop", ())
    if drop:
        td = tok[:, :nd, :NCONT].reshape(-1, NCONT)
        for c in drop:
            j = TOK.index(c)
            tmean[j], tstd[j] = td[:, j].mean(), td[:, j].std()
    return tmean, np.where(tstd > 0, tstd, 1.0)


def act_vars(nE, pc40):
    return np.stack([np.log1p(np.clip(nE[:, 0].astype(np.float32), 0, None)), np.log1p(np.clip(pc40.astype(np.float32), 0, None))], 1)


def dcor2(x, y, nmax=3000, seed=1):
    """dCor² (numpy; diagnostyka w meta) — jak analyze.dcor2."""
    n = len(x)
    if n < 10:
        return float("nan")
    if n > nmax:
        i = np.random.default_rng(seed).choice(n, nmax, replace=False)
        x, y = x[i], y[i]
    x, y = x.astype(np.float64), y.astype(np.float64)
    a = np.abs(x[:, None] - x[None, :]); b = np.abs(y[:, None] - y[None, :])
    A = a - a.mean(0) - a.mean(1)[:, None] + a.mean(); B = b - b.mean(0) - b.mean(1)[:, None] + b.mean()
    return float((A * B).mean() / np.sqrt(max((A * A).mean() * (B * B).mean(), 1e-30)))


def load3(ddir, edir, sp="D0"):
    """dn2_train.load + pc40 (Σ pT otoczenia) z części dn2_<próbka>_*.npz, w tej samej kolejności wierszy (kontrola po idx)."""
    S, full = D2T.load(ddir, edir, sp)
    for s, T in S.items():
        parts = sorted(glob.glob(os.path.join(edir, f"dn2_{s}_*.npz")))
        for pre in ("a_", "p_"):
            if pre + "idx" not in T:
                continue
            k = f"{sp}__{pre}"
            ix, pc = [], []
            for f in parts:
                e = np.load(f, allow_pickle=False)
                if k + "idx" in e.files:
                    if k + "pc40" not in e.files:
                        pc = None
                        break
                    ix.append(e[k + "idx"]); pc.append(e[k + "pc40"])
            if pc is None:
                continue
            ix, pc = np.concatenate(ix), np.concatenate(pc)
            o = np.argsort(ix)
            if not np.array_equal(ix[o], T[pre + "idx"]):
                raise RuntimeError(f"{s} {pre}: kolejność wierszy pc40 niezgodna z dn2_train.load")
            T[pre + "pc40"] = pc[o]
    return S, full


# ----------------------------------------------------------------------------------------------------------- trening
def run_config(c, S, full, a, log):
    import torch  # noqa: F401
    import dn_models as DM
    import dn2_models as D2
    import dn3_models as D3
    assert D2.TOK == TOK and D2.NCONT == NCONT, "dn2_models.TOK/NCONT zmienione — popraw dn3_train"
    K = a.folds
    P = DT.pool(S, "D0", c["bkg"], a.cap_sig, a.cap_bkg)
    g = lambda k: np.concatenate([S[s]["a_" + k][i] for s, i, _ in P["src"]])  # noqa: E731
    G = g("X")
    Dd = g("Dd")
    tok, mask = tokens3(Dd, g("E"), g("eta"), g("phi"), c["env"])
    nd = Dd.shape[1]
    del Dd
    fold = DF.fold_of(P["uid"], K)
    val = DF.unit_hash(P["uid"], salt=23) >= 0.9
    dec, lam = "d", float(c["lam"])
    have_act = all("a_pc40" in S[s] for s, _, _ in P["src"])
    A = act_vars(g("nE"), g("pc40")) if have_act else None
    am = aw = None
    if c["act"]:
        if A is None:
            raise SystemExit(f"[{c['tag']}] brak pc40 w częściach dn2 — wariant z dekorelacją aktywności niemożliwy")
        am = (P["y"] & (P["w"] > 0)) if c["act"] == "S" else P["dmask"]
        aw = np.where(am, np.where(P["dmask"], P["dw"], P["w"]), 0.0).astype(np.float32)
    t0 = time.time()
    log(f"[{c['tag']}] pula {len(P['y'])} (S {int(P['y'].sum())}, B {int((P['role'] == 'B').sum())}), tokenów {tok.shape[1]}, średnio ważnych "
        f"{mask.sum(1).mean():.2f} (otoczenie S {mask[P['y'], nd:].sum(1).mean():.2f}, B {mask[~P['y'], nd:].sum(1).mean():.2f}), wariant {VAR[c['env']]}, "
        f"aktywność {c['act'] or '-'} λ_act {c['lam_act']:g}" + (f" (wierszy {int(am.sum())})" if am is not None else "") + f", urządzenie {a.device}")
    models, info = [], []
    for k in range(K):
        trn = fold != k
        u = trn & (P["w"] > 0)
        med, sc = DM.prep_stats(G[u])
        tmean, tstd = tok_stats(tok[u], mask[u], nd, c["env"])
        net = D2.SetNet(G.shape[1], c["model"], tmean, tstd, med, sc)
        args = (net, (tok[trn], mask[trn], G[trn]), P["y"][trn], P["w"][trn], P["m"][trn], P["dmask"][trn], P["dw"][trn], val[trn], dec, lam,
                1000 * c["seed"] + k, log)
        if c["fit3"]:
            kw = dict(act=A[trn], amask=am[trn], aw=aw[trn], lam_act=c["lam_act"]) if c["act"] else {}
            net, h = D3.fit3(*args, device=a.device, **kw)
        else:
            net, h = DM.fit(*args, device=a.device)            # ablacje wejść: ta sama pętla co w etapie 2
        models.append(net); info.append(h)
    out, subs = {}, {}
    for s, T in S.items():
        for pre in ("a_", "p_"):
            if pre + "uid" not in T:
                continue
            f = DF.fold_of(T[pre + "uid"], K)
            sc_ = np.full(full[(s, pre)], np.nan, np.float32)
            sub = np.full(len(f), np.nan, np.float32)
            for a0 in range(0, len(f), 400_000):
                sl = slice(a0, a0 + 400_000)
                tk, mk = tokens3(T[pre + "Dd"][sl], T[pre + "E"][sl], T[pre + "eta"][sl], T[pre + "phi"][sl], c["env"])
                for k in range(K):
                    w = f[sl] == k
                    if w.any():
                        sub[sl][w] = DM.predict(models[k], (tk[w], mk[w], T[pre + "X"][sl][w]), device=a.device)
            sc_[T[pre + "idx"]] = sub
            out[f"{s}__{pre[0]}"] = sc_
            subs[(s, pre)] = sub
    os.makedirs(os.path.join(a.out, "sc"), exist_ok=True); os.makedirs(os.path.join(a.out, "meta"), exist_ok=True)
    np.savez_compressed(os.path.join(a.out, "sc", c["tag"] + ".npz"), **out)
    meta = {**c, "features": DF.H, "token": D2.TOK, "n_feat": f"{len(DF.H)} + {tok.shape[1]}×{len(D2.TOK)}", "n_pool": int(len(P["y"])), "n_sig": int(P["y"].sum()),
            "n_bkg": int((P["role"] == "B").sum()), "nparams": DM.nparams(models[0]), "npar_unit": "wagi", "folds": K, "fit": info,
            "sec": round(time.time() - t0, 1), "threads": a.threads, "device": a.device, "hp": DM.HP, "arch": c["model"], "variant": VAR[c["env"]],
            "act_vars": ACT, "mean_env_tokens": {"S": float(mask[P["y"], nd:].sum(1).mean()), "B": float(mask[~P["y"], nd:].sum(1).mean())}}
    if A is not None:                    # diagnostyka: dCor²(score OOF, aktywność) na sygnale i tle puli (dla każdego wariantu, także bez członu)
        so = np.concatenate([subs[(s, "a_")][i] for s, i, _ in P["src"]])
        meta["act_dcor2"] = {r: [dcor2(so[mk], A[mk, j]) for j in range(A.shape[1])] for r, mk in (("sig", P["y"]), ("bkg", P["role"] == "B"))}
        log(f"  dCor²(score, aktywność): {meta['act_dcor2']}")
    json.dump(meta, open(os.path.join(a.out, "meta", c["tag"] + ".json"), "w"), indent=1)
    log(f"[{c['tag']}] gotowe, {time.time() - t0:.0f}s, parametrów {meta['nparams']}")
    return meta


# ----------------------------------------------------------------------------------------------------------- test mechaniki
def _synth_table(rng, n, kind):
    nf = len(DF.ALLF["D0"])
    truth = {"fd": np.where(rng.random(n) < 0.85, 0, 8), "mb": np.full(n, 8), "data": np.full(n, -1)}[kind].astype(np.int16)
    sig = (truth >= 0) & (truth < 8)
    X = rng.normal(size=(n, nf)).astype(np.float32)
    X[:, :3] += sig[:, None] * 1.0
    pt = np.clip(3.0 + rng.exponential(3.0, n), 3.0, 29.0)
    eta, phi = rng.uniform(-2.4, 2.4, n), rng.uniform(-np.pi, np.pi, n)
    m = np.where(sig, 1.86484 + rng.normal(0, 0.02, n), rng.uniform(1.75, 2.0, n))
    Dd = np.zeros((n, 2, 8), np.float32)
    fr = rng.uniform(0.3, 0.7, n)
    Dd[:, 0, 0], Dd[:, 1, 0] = pt * fr, pt * (1 - fr)
    Dd[:, :, 1] = eta[:, None] + rng.normal(0, 0.2, (n, 2)); Dd[:, :, 2] = phi[:, None] + rng.normal(0, 0.2, (n, 2))
    Dd[:, :, 3] = rng.normal(0, 0.1, (n, 2)); Dd[:, :, 4] = rng.normal(0, 0.2, (n, 2)); Dd[:, :, 5] = rng.uniform(0.5, 3, (n, 2))
    Dd[:, :, 6] = rng.integers(5, 20, (n, 2)); q = rng.choice([-1.0, 1.0], n); Dd[:, 0, 7], Dd[:, 1, 7] = q, -q
    T = {"X": X, "Dd": Dd, "m": m.astype(np.float32), "pt": pt.astype(np.float32), "eta": eta.astype(np.float32), "phi": phi.astype(np.float32),
         "uid": rng.integers(1, 2 ** 62, n).astype(np.uint64), "truth": truth, "ok": np.ones(n, bool), "swap": np.zeros(n, np.int8),
         "m_swap": (m + 0.1).astype(np.float32), "w": np.ones(n, np.float32)}
    # otoczenie: tło ma więcej torów (→ sieć uczy się aktywności); do 20 torów w ΔR < 1, malejąco w pT
    ntrk = np.minimum(rng.poisson(np.where(sig, 4.0, 10.0)), D2T_NT)
    pe = (0.1 + rng.exponential(0.8, (n, D2T_NT))) * (np.arange(D2T_NT)[None, :] < ntrk[:, None])
    pe = -np.sort(-pe, axis=1)
    r, an = np.sqrt(rng.random((n, D2T_NT))), rng.uniform(0, 2 * np.pi, (n, D2T_NT))
    E = np.stack([pe, r * np.cos(an), r * np.sin(an), rng.choice([-1.0, 1.0], (n, D2T_NT)), rng.uniform(0.5, 3, (n, D2T_NT)),
                  rng.integers(3, 20, (n, D2T_NT)).astype(float)], 2)
    E *= (pe > 0)[:, :, None]
    in04 = (pe > 0) & (r < 0.4)
    T["E"], T["nE"] = E.astype(np.float16), np.stack([in04.sum(1), ntrk], 1).astype(np.int16)
    T["pc40"] = (pe * (in04 & (pe > 0.5))).sum(1).astype(np.float32)
    return T


D2T_NT = 20          # dn2_extract.NT (liczba zapisanych torów otoczenia na kandydata)


def synth(ddir, edir, seed=5, n=6000):
    """Tabele dn_<próbka>.npz (jak etap 1) i części dn2_<próbka>_<k>.npz (jak etap 2) z liczb losowych."""
    rng = np.random.default_rng(seed)
    os.makedirs(ddir, exist_ok=True); os.makedirs(edir, exist_ok=True)
    for s, kind in (("fd_nd", "fd"), ("mb_nd", "mb"), ("data_428770", "data")):
        z, parts = {}, [{}, {}]
        for pre, nn_, frac in (("a_", n, 0.8), ("p_", n // 4, 1.0)):
            if pre == "p_" and kind == "mb":
                continue
            T = _synth_table(rng, nn_, kind)
            k = f"D0__{pre}"
            for c in D2T.KEYS:
                z[k + c] = T[c]
            idx = rng.permutation(np.flatnonzero(rng.random(nn_) < frac))          # części w dowolnej kolejności wierszy
            for j, ii in enumerate((idx[: len(idx) // 2], idx[len(idx) // 2:])):
                parts[j].update({k + "idx": ii, k + "E": T["E"][ii], k + "nE": T["nE"][ii], k + "pc40": T["pc40"][ii]})
        np.savez(os.path.join(ddir, f"dn_{s}.npz"), **z)
        for j, p in enumerate(parts):
            np.savez(os.path.join(edir, f"dn2_{s}_{j}.npz"), **p)


def test_numpy(S, log):
    """Logika wariantów tokenów, standaryzacji, aktywności i wyrównania wierszy — bez torch."""
    T = S["mb_nd"]
    Dd, E, eta, phi = T["a_Dd"], T["a_E"], T["a_eta"], T["a_phi"]
    nd, n = Dd.shape[1], len(Dd)
    t0, m0 = D2T.tokens(Dd, E, eta, phi, "e04")
    ne = t0.shape[1] - nd
    Ef = E.astype(np.float32)
    ok_all = (Ef[:, :, 0] > 0) & (np.hypot(Ef[:, :, 1], Ef[:, :, 2]) < 0.4)
    t, m = tokens3(Dd, E, eta, phi, "e04")
    assert np.array_equal(t, t0) and np.array_equal(m, m0), "e04 ≠ dn2_train.tokens"
    assert np.array_equal(m0[:, nd:] > 0, ok_all[:, :ne]) and m0[:, :nd].all()
    for env in ("e04g", "e04np"):
        t, m = tokens3(Dd, E, eta, phi, env)
        dj = [TOK.index(c) for c in VAR[env]["drop"]]
        kj = [j for j in range(len(TOK)) if j not in dj]
        assert np.array_equal(m, m0) and np.array_equal(t[:, :nd], t0[:, :nd]), f"{env}: maska/córki zmienione"
        assert not t[:, nd:, dj].any(), f"{env}: cechy pędowe otoczenia niewyzerowane"
        assert np.array_equal(t[:, nd:, kj], t0[:, nd:, kj]) and t0[:, nd:, dj].any(), f"{env}: pozostałe cechy zmienione"
        mu, sd = tok_stats(t, m, nd, env)
        td = t[:, :nd, :NCONT].reshape(-1, NCONT)
        assert np.allclose(mu[dj], td[:, dj].mean(0)) and np.allclose(sd[dj], td[:, dj].std(0)) and np.isfinite(mu).all() and (sd > 0).all()
    mu0, sd0 = tok_stats(t0, m0, nd, "e04")
    tv = t0[:, :, :NCONT][m0 > 0]
    assert np.array_equal(mu0, tv.mean(0)) and np.array_equal(sd0, np.where(tv.std(0) > 0, tv.std(0), 1.0)), "standaryzacja e04 ≠ etap 2"
    t, m = tokens3(Dd, E, eta, phi, "e04pt1")
    exp = ok_all[:, :ne] & (Ef[:, :ne, 0] > 1.0)
    assert np.array_equal(m[:, nd:] > 0, exp) and np.array_equal(t[:, :nd], t0[:, :nd]) and m[:, :nd].all(), "e04pt1: maska"
    assert not t[:, nd:][~exp].any() and np.array_equal(t[:, nd:][exp], t0[:, nd:][exp]), "e04pt1: tokeny"
    assert (t[:, nd:, 0][exp] > 0).all() and exp.sum() < (m0[:, nd:] > 0).sum(), "e04pt1: próg pT"
    t, m = tokens3(Dd, E, eta, phi, "e04c")
    nv = (m[:, nd:] > 0).sum(1)
    assert np.array_equal(nv, np.minimum(ok_all.sum(1), ne)) and (nv >= (m0[:, nd:] > 0).sum(1)).all(), "e04c: liczba torów"
    assert (m[:, nd:][:, :-1] >= m[:, nd:][:, 1:]).all(), "e04c: ważne tokeny nie na początku"
    assert (np.diff(np.where(m[:, nd:] > 0, t[:, nd:, 0], -1e9), axis=1) <= 1e-6).all(), "e04c: kolejność pT"
    assert np.array_equal(compact(E)[ok_all.sum(1) == 0], E[ok_all.sum(1) == 0]), "compact: zmienia wiersze bez torów w ΔR < 0,4"
    # aktywność i wyrównanie pc40 (load3) — pc40 odtworzone z E (w danych syntetycznych ≤ 20 torów, więc dokładnie)
    for s, T in S.items():
        for pre in ("a_", "p_"):
            if pre + "idx" not in T:
                continue
            e = T[pre + "E"].astype(np.float32)
            i4 = (e[:, :, 0] > 0) & (np.hypot(e[:, :, 1], e[:, :, 2]) < 0.4)
            # zaokrąglenie float16 przesuwa pojedyncze tory przez próg ΔR = 0,4 / pT = 0,5 → zgodność w ≥ 97 % wierszy (przesunięcie dałoby ~0 %)
            assert pre + "pc40" in T, f"{s} {pre}: brak pc40"
            fp = np.mean(np.abs(T[pre + "pc40"] - (e[:, :, 0] * (i4 & (e[:, :, 0] > 0.5))).sum(1)) < 0.02)
            fn = np.mean(T[pre + "nE"][:, 0] == i4.sum(1))
            assert fp > 0.97 and fn > 0.97, f"{s} {pre}: pc40/nE przesunięte względem E (zgodnych {fp:.3f}/{fn:.3f})"
    A = act_vars(T["a_nE"], T["a_pc40"])
    assert A.shape == (len(T["a_m"]), 2) and np.isfinite(A).all() and (A >= 0).all()
    x = np.random.default_rng(1).normal(size=2000)
    assert dcor2(x, x) > 0.99 and dcor2(x, np.random.default_rng(2).normal(size=2000)) < 0.02
    tags = [c["tag"] for cs in BLOCKS.values() for c in cs]
    assert len(set(tags)) == len(tags) and f"D0.tr.H+e04.mc.d{LAM}.s1" not in tags, "tag koliduje z etapem 2"
    log(f"numpy OK: warianty tokenów {list(VAR)}, standaryzacja, pc40/nE wyrównane, tagi: {' '.join(tags)}")


def test_torch(S, full, odir, log):
    import argparse as _ap
    import torch
    import dn_models as DM
    import dn2_models as D2
    import dn3_models as D3
    torch.set_num_threads(2)
    hp0 = dict(DM.HP)
    DM.HP.update(epochs=1, batch=256)
    a = _ap.Namespace(folds=2, cap_sig=0, cap_bkg=0, device="cpu", threads=2, out=odir)
    try:
        for b, cs in BLOCKS.items():
            for c in cs:
                meta = run_config(c, S, full, a, log)
                z = np.load(os.path.join(odir, "sc", c["tag"] + ".npz"))
                for s, T in S.items():
                    for pre in ("a_", "p_"):
                        if pre + "idx" not in T:
                            continue
                        v = z[f"{s}__{pre[0]}"]
                        fin = np.isfinite(v)
                        assert len(v) == full[(s, pre)] and fin.sum() == len(T[pre + "idx"]) and fin[T[pre + "idx"]].all(), f"{c['tag']} {s} {pre}: wzór NaN"
                        assert (v[fin] >= 0).all() and (v[fin] <= 1).all()
                assert len(meta["fit"]) == 2 and "act_dcor2" in meta
                if c["act"]:
                    h = meta["fit"][0]["hist"][0]
                    assert len(h) == 4 and np.isfinite(h).all() and h[3] > 0, f"{c['tag']}: człon aktywności {h}"
        # fit3(λ_act = 0) ≡ dn_models.fit (człon masy bez zmian)
        P = DT.pool(S, "D0", "mc", 0, 0)
        g = lambda k: np.concatenate([S[s]["a_" + k][i] for s, i, _ in P["src"]])  # noqa: E731
        G, Dd = g("X"), g("Dd")
        tok, mask = tokens3(Dd, g("E"), g("eta"), g("phi"), "e04")
        nd = Dd.shape[1]
        val = DF.unit_hash(P["uid"], salt=23) >= 0.9
        med, sc = DM.prep_stats(G); tmean, tstd = tok_stats(tok, mask, nd, "e04")
        A = act_vars(g("nE"), g("pc40"))
        am = P["y"] & (P["w"] > 0)
        aw = np.where(am, P["w"], 0.0).astype(np.float32)
        DM.HP.update(epochs=2)

        def mk():
            torch.manual_seed(0)
            return D2.SetNet(G.shape[1], "tr", tmean, tstd, med, sc)
        args = lambda net: (net, (tok, mask, G), P["y"], P["w"], P["m"], P["dmask"], P["dw"], val, "d", 5.0, 1000, log)  # noqa: E731
        n1, h1 = DM.fit(*args(mk()))
        n2, h2 = D3.fit3(*args(mk()))
        p1, p2 = DM.predict(n1, (tok, mask, G)), DM.predict(n2, (tok, mask, G))
        d = float(np.max(np.abs(p1 - p2)))
        log(f"fit3(λ_act=0) vs dn_models.fit: max|Δscore| = {d:.2e}, hist {h1['hist'][-1]} / {h2['hist'][-1][:3]}")
        assert d < 1e-4 and h1["best_epoch"] == h2["best_epoch"], "fit3 przy λ_act = 0 różni się od dn_models.fit — człon masy NIE jest bez zmian"
        # działanie członu aktywności (informacyjnie): dCor²(score, aktywność) na sygnale z λ_act = 50 vs 0
        DM.HP.update(epochs=4)
        n3, _ = D3.fit3(*args(mk()), act=A, amask=am, aw=aw, lam_act=50.0)
        n4, _ = D3.fit3(*args(mk()))
        d3 = [dcor2(DM.predict(n3, (tok, mask, G))[am], A[am, j]) for j in range(2)]
        d4 = [dcor2(DM.predict(n4, (tok, mask, G))[am], A[am, j]) for j in range(2)]
        log(f"dCor²(score, aktywność | sygnał): λ_act=0 {d4}, λ_act=50 {d3}" + ("" if sum(d3) < sum(d4) else "  UWAGA: człon aktywności nie obniżył dCor² (test informacyjny)"))
    finally:
        DM.HP.clear(); DM.HP.update(hp0)
    log("torch OK: wszystkie bloki (1 epoka, 2 foldy, CPU), wzór NaN score, fit3 ≡ fit przy λ_act = 0")


def selftest(a):
    import tempfile
    t0 = time.time()
    log = lambda s: print(f"[{time.time() - t0:6.0f}s] {s}", flush=True)  # noqa: E731
    ctx = None
    if a.dir:
        d = a.dir
        os.makedirs(d, exist_ok=True)
    else:
        ctx = tempfile.TemporaryDirectory(prefix="dn3_selftest_")
        d = ctx.name
    try:
        synth(os.path.join(d, "data"), os.path.join(d, "parts"))
        S, full = load3(os.path.join(d, "data"), os.path.join(d, "parts"))
        log("dane syntetyczne: " + ", ".join(f"{s} a={len(T.get('a_m', []))}/{full.get((s, 'a_'))} p={len(T.get('p_m', []))}/{full.get((s, 'p_'))}" for s, T in S.items()))
        test_numpy(S, log)
        try:
            import torch  # noqa: F401
        except ImportError:
            if a.need_torch:
                raise
            log("torch niedostępny — część torch (trening, fit3, score) NIE przetestowana")
            return 0
        test_torch(S, full, os.path.join(d, "out"), log)
    finally:
        if ctx is not None:
            ctx.cleanup()
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["train", "list", "selftest"])
    ap.add_argument("--data"); ap.add_argument("--env"); ap.add_argument("--out"); ap.add_argument("--block"); ap.add_argument("--only", default="")
    ap.add_argument("--folds", type=int, default=int(os.environ.get("DN_FOLDS", "5")))
    ap.add_argument("--cap-sig", type=int, default=int(os.environ.get("DN3_CAP_SIG", "300000")))
    ap.add_argument("--cap-bkg", type=int, default=int(os.environ.get("DN3_CAP_BKG", "500000")))
    ap.add_argument("--epochs", type=int, default=int(os.environ.get("DN3_EPOCHS", "20")))
    ap.add_argument("--threads", type=int, default=int(os.environ.get("DN_CPUS", "8")))
    ap.add_argument("--device", default="auto")
    ap.add_argument("--need-torch", action="store_true"); ap.add_argument("--dir", default="")
    a = ap.parse_args(argv)
    if a.cmd == "list":
        for b, cs in BLOCKS.items():
            print(b, len(cs), " ".join(c["tag"] for c in cs))
        return 0
    if a.cmd == "selftest":
        return selftest(a)
    if a.block not in BLOCKS:
        raise SystemExit(f"nieznany blok {a.block}; są: {' '.join(BLOCKS)}")
    import torch
    import dn_models as DM
    torch.set_num_threads(max(1, a.threads))
    DM.HP["epochs"] = a.epochs
    if a.device == "auto":
        a.device = "cuda" if torch.cuda.is_available() else "cpu"
    if a.device == "cuda":  # test dymny jak w dn2_train: karta bez jądra CUDA dla tej wersji torch (np. V100) → CPU zamiast padać
        try:
            _t = torch.randn(64, 8, device="cuda")
            float((torch.cat([_t, _t], -1) @ torch.randn(16, 4, device="cuda")).sum())
            float(torch.cdist(_t, _t, p=1).sum())
        except Exception as e:  # noqa: BLE001
            print(f"UWAGA: CUDA nie działa na {torch.cuda.get_device_name(0)} ({type(e).__name__}) → CPU", flush=True)
            a.device = "cpu"
    t0 = time.time()
    log = lambda s: print(f"[{time.time() - t0:6.0f}s] {s}", flush=True)  # noqa: E731
    log(f"torch {torch.__version__}, urządzenie {a.device}" + (f" ({torch.cuda.get_device_name(0)})" if a.device == "cuda" else ""))
    S, full = load3(a.data, a.env)
    log("wczytano: " + ", ".join(f"{s} a={len(T.get('a_m', []))}/{full.get((s, 'a_'))} p={len(T.get('p_m', []))}/{full.get((s, 'p_'))}" for s, T in S.items()))
    for c in BLOCKS[a.block]:
        if a.only and c["tag"] not in a.only.split(","):
            continue
        if os.path.isfile(os.path.join(os.environ.get("DN_DONE_DIR") or a.out, "meta", c["tag"] + ".json")) and os.environ.get("DN_REDO") != "1":
            log(f"[{c['tag']}] już jest — pomijam")
            continue
        run_config(c, S, full, a, log)
    return 0


if __name__ == "__main__":
    sys.exit(main())
