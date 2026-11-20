#!/usr/bin/env python3
"""ISO_pis (19.11.2026, projekt B) — czy wolny pion z D*→D0π wchodzi do ptcone kandydata D0 i czy przez to ε_iso / SF izolacji z próbki tagu D*
różni się od inkluzywnej. TYLKO MC (fd_nd, mb_nd), samo numpy, bez modeli i bez fitów.
Wejście (EOS, tylko odczyt): tabele DNET etapu 1 `dn_<próbka>.npz` (wiersze a_ = D0 inkluzywne, p_ = pary D0←D*; X zawiera gałęzie
Dmeson_ptcone40/100 = kolumny ptcone40_iso/ptcone100_iso) + tory otoczenia etapu 2 `dn2_<próbka>_*.npz` (E: 20 torów o największym pT w ΔR < 1,0,
bez córek D0; drzewo Tracks ma pT > 0,2 GeV).
Grupy (sygnał z prawdy: truth ∈ {0,1}, swap = 0, |m − m0| < 55 MeV; biny pT D0):
  a   próbka tagu: wiersze p_ z prawdą D* (truth_ds ∈ {0,1}), jeden wiersz na D0
  b   pozostałe D0 sygnałowe: wiersze a_ bez zrekonstruowanej pary D* z prawdą (zawiera D0←D* bez tagu → b1/b2)
  b1  b z torem zgodnym z wolnym pionem w stożku 0,4 (maska etapu 6: Δm < 0,155 GeV, ładunek = −q_K);  b2 = b bez takiego toru
  c   a po odjęciu pT wolnego pionu (prawda tabeli: p_pis_pt, tor znaleziony w tokenach) z ptcone, gdy pion leży w stożku
  am / bm / im   a / b / inkluzywne po odjęciu torów maski (definicja bez prawdy, stosowalna w danych)
  i   inkluzywne (wszystkie a_)
Użycie: iso_pis.py --data DIR --env DIR --out OUT.json"""
import argparse, glob, json, os, sys, time
import numpy as np

M0, WIN = 1.86484, 0.055
MPI = 0.13957
DM_MASK = 0.155
SIG = (0, 1)
BINS = [(3, 5), (5, 8), (8, 30), (3, 30)]
CUTS = [("ptcone40 < 2 GeV", 40, "abs", 2.0), ("ptcone40 < 4 GeV", 40, "abs", 4.0), ("ptcone40/pT < 0.5", 40, "rel", 0.5),
        ("ptcone40/pT < 1.0", 40, "rel", 1.0), ("ptcone100 < 8 GeV", 100, "abs", 8.0), ("ptcone100/pT < 1.0", 100, "rel", 1.0)]
TARGETS = [0.98, 0.95, 0.90, 0.76, 0.51]


def load(ddir, edir, s, pre, extra=()):
    z = np.load(os.path.join(ddir, "dn_%s.npz" % s), allow_pickle=False)
    feats = json.loads(str(z["meta"]))
    feats = (feats.get("feats") or feats["groups"][0]["feats"])["D0"]
    k = "D0__" + pre
    ix, E, pc = [], [], []
    for f in sorted(glob.glob(os.path.join(edir, "dn2_%s_*.npz" % s))):
        e = np.load(f, allow_pickle=False)
        if k + "idx" in e.files:
            ix.append(e[k + "idx"]); E.append(e[k + "E"]); pc.append(e[k + "pc40"])
    ix, E, pc = np.concatenate(ix), np.concatenate(E), np.concatenate(pc)
    o = np.argsort(ix)
    ix, E, pc = ix[o], E[o], pc[o]
    n_full = len(z[k + "m"])
    T = {c: z[k + c][ix] for c in ["m", "pt", "eta", "phi", "w", "truth", "swap", "row", "ftag", "Dd"] + list(extra)}
    X = z[k + "X"]
    T["c40"] = X[ix, feats.index("ptcone40_iso")].astype(np.float64)
    T["c100"] = X[ix, feats.index("ptcone100_iso")].astype(np.float64)
    del X
    for c in ("m", "pt", "eta", "phi", "w"):
        T[c] = T[c].astype(np.float64)
    T["E"], T["pc40"], T["n_full"] = E, pc.astype(np.float64), n_full
    if pre == "p_":      # klucze WSZYSTKICH par z prawdą D* i D0 (także bez torów otoczenia) → flaga „D0 ma tag” w tabeli a_
        f = np.isin(z[k + "truth"], SIG) & (z[k + "swap"] == 0) & (np.abs(z[k + "m"].astype(np.float64) - M0) < WIN) & np.isin(z[k + "truth_ds"], SIG)
        T["key_all"] = (z[k + "ftag"].astype(np.int64) * (1 << 40) + z[k + "row"].astype(np.int64))[f]
    return T


def token_dm(E, pt, eta, m):
    E = E.astype(np.float64)
    pt, eta, m = (np.asarray(x, np.float64)[:, None] for x in (pt, eta, m))
    tpt = E[:, :, 0]
    teta = eta + E[:, :, 1]
    ED = np.sqrt(m ** 2 + (pt * np.cosh(eta)) ** 2)
    Et = np.sqrt(MPI ** 2 + (tpt * np.cosh(teta)) ** 2)
    m2 = m ** 2 + MPI ** 2 + 2 * (ED * Et - pt * tpt * np.cos(E[:, :, 2]) - pt * np.sinh(eta) * tpt * np.sinh(teta))
    dm = np.sqrt(np.clip(m2, 0, None)) - m
    return np.where(tpt > 0, dm, np.nan)


def env(T):
    """Dopisuje: dr tokenów, maskę (bez prawdy), sumy pT torów maski w 0,4 / 1,0, własny stożek z tokenów."""
    E = T["E"]
    tpt = E[:, :, 0].astype(np.float64)
    dr = np.hypot(E[:, :, 1].astype(np.float64), E[:, :, 2].astype(np.float64))
    dm = token_dm(E, T["pt"], T["eta"], T["m"])
    qK = np.nan_to_num(T["Dd"][:, 0, 7].astype(np.float64))
    v = np.isfinite(dm) & (dm < DM_MASK) & ((E[:, :, 3].astype(np.float64) * qK[:, None]) < 0)
    T["tpt"], T["dr"], T["dmtok"], T["veto"] = tpt, dr, dm, v
    T["mk40"] = (tpt * (v & (dr < 0.4))).sum(1)
    T["mk100"] = (tpt * (v & (dr < 1.0))).sum(1)
    T["own40"] = (tpt * ((dr < 0.4) & (tpt > 0))).sum(1)
    T["n40"] = ((dr < 0.4) & (tpt > 0)).sum(1)
    return T


def eff(pas, w):
    sw = float(w.sum())
    if sw <= 0 or len(w) < 10:
        return None
    e = float((w * pas).sum() / sw)
    neff = sw ** 2 / float((w ** 2).sum())
    return [e, float(np.sqrt(max(e * (1 - e), 1e-9) / neff)), int(len(w))]


def ratio(x, y):
    if not x or not y or y[0] <= 0 or x[0] <= 0:
        return None
    r = x[0] / y[0]
    return [r, r * float(np.hypot(x[1] / x[0], y[1] / y[0]))]


def wq(x, w, qs):
    if len(x) == 0:
        return [None] * len(qs)
    o = np.argsort(x)
    cw = np.cumsum(w[o]); cw /= cw[-1]
    return [float(x[o][min(np.searchsorted(cw, q), len(x) - 1)]) for q in qs]


def passcut(c40, c100, pt, cut):
    _, R, kind, thr = cut
    v = c40 if R == 40 else c100
    return (v < thr) if kind == "abs" else (v / pt < thr)


def bisect(f, lo, hi, target, n=40):
    """f malejąca; zwraca parametr, przy którym f = target, albo None (cel nieosiągalny w [lo, hi])."""
    if f(hi) > target or f(lo) < target:
        return None
    for _ in range(n):
        mid = 0.5 * (lo + hi)
        if f(mid) > target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def sample(ddir, edir, s):
    t0 = time.time()
    A = env(load(ddir, edir, s, "a_"))
    P = env(load(ddir, edir, s, "p_", extra=("pis_pt", "truth_ds", "dm")))
    sig = lambda T: np.isin(T["truth"], SIG) & (T["swap"] == 0) & (np.abs(T["m"] - M0) < WIN)  # noqa: E731
    rows = np.arange(len(P["m"]))
    pis = P["pis_pt"].astype(np.float64)
    rel = np.abs(P["tpt"] - pis[:, None]) / np.clip(pis[:, None], 1e-3, None)
    rel[P["tpt"] <= 0] = np.inf
    j = np.argmin(rel, 1)
    found = rel[rows, j] < 0.01
    drp = np.where(found, P["dr"][rows, j], np.nan)
    dmp = np.where(found, P["dmtok"][rows, j], np.nan)
    in40, in100 = found & (drp < 0.4), found & (drp < 1.0)
    vetop = found & P["veto"][rows, j]
    sp_ = sig(P) & np.isin(P["truth_ds"], SIG)
    key = lambda T: T["ftag"].astype(np.int64) * (1 << 40) + T["row"].astype(np.int64)  # noqa: E731
    kp = key(P)
    _, first = np.unique(np.where(sp_, kp, -1 - rows), return_index=True)
    one = np.zeros(len(kp), bool); one[first] = True
    tag = sp_ & one
    sa = sig(A)
    in_tag = np.isin(key(A), P["key_all"])
    # stożki po odjęciu wolnego pionu
    P["c40c"] = np.clip(P["c40"] - pis * in40, 0, None)
    P["c100c"] = np.clip(P["c100"] - pis * in100, 0, None)
    # wariant górny: pion nieznaleziony w tokenach (pT < 0,2 GeV / poza 20 torami) też odjęty, jeśli gałąź ≥ pT pionu
    um = (~found) & (P["c40"] >= pis - 1e-3)
    P["c40x"] = np.clip(P["c40c"] - pis * um, 0, None)
    for T in (A, P):
        T["c40m"] = np.clip(T["c40"] - T["mk40"], 0, None)
        T["c100m"] = np.clip(T["c100"] - T["mk100"], 0, None)
    b = sa & ~in_tag
    b1 = b & (A["mk40"] > 0)
    res = {"sample": s, "rows": {"a_env": len(A["m"]), "a_full": A["n_full"], "p_env": len(P["m"]), "p_full": P["n_full"]},
           "n_sig_incl": int(sa.sum()), "n_tag": int(tag.sum()), "n_incl_tagged": int((sa & in_tag).sum()), "n_tag_keys_all": int(len(np.unique(P["key_all"]))), "bins": {}}
    # zgodność własnego stożka (tokeny, pT > 0,2) z gałęzią
    d = A["own40"][sa] - A["c40"][sa]
    res["own_vs_branch_sig_incl"] = {"abs_lt_1MeV": float(np.mean(np.abs(d) < 2e-3)), "abs_lt_200MeV": float(np.mean(np.abs(d) < 0.2)),
                                     "own_le_branch": float(np.mean(d < 2e-3)), "pc40_le_branch": float(np.mean(A["pc40"][sa] <= A["c40"][sa] + 2e-3))}
    d = (P["c40"] - pis)[tag & in40]
    res["branch_minus_pis_in40"] = {"frac_negative_gt_5MeV": float(np.mean(d < -5e-3)), "frac_equal_pis_within_2MeV": float(np.mean(np.abs(d) < 2e-3))}
    for lo, hi in BINS:
        bp = tag & (P["pt"] >= lo) & (P["pt"] < hi)
        ba = (A["pt"] >= lo) & (A["pt"] < hi)
        if bp.sum() < 50:
            continue
        w = P["w"][bp]
        fr = lambda m_: float((w * m_[bp]).sum() / w.sum())  # noqa: E731
        o = {"n_tag": int(bp.sum()), "n_b": int((b & ba).sum()), "n_b1": int((b1 & ba).sum()), "n_incl": int((sa & ba).sum()),
             "tag_frac_of_incl": float((A["w"][sa & ba & in_tag]).sum() / A["w"][sa & ba].sum()),
             "fb_tag": float((w * (P["truth"][bp] == 1)).sum() / w.sum()), "fb_b": float((A["w"][b & ba] * (A["truth"][b & ba] == 1)).sum() / A["w"][b & ba].sum()),
             "b1_frac_of_b": float((A["w"][b1 & ba]).sum() / max(A["w"][b & ba].sum(), 1e-9))}
        # --- pkt 1
        o["pis"] = {
            "pt_q10_25_50_75_90": wq(pis[bp], w, [0.1, 0.25, 0.5, 0.75, 0.9]),
            "frac_pt_lt_0.2": fr(pis < 0.2), "frac_pt_gt_0.5": fr(pis > 0.5), "frac_pt_gt_1.0": fr(pis > 1.0),
            "found_token": fr(found), "in40": fr(in40), "in100": fr(in100), "in40_pt_gt_0.5": fr(in40 & (pis > 0.5)), "in40_pt_gt_1.0": fr(in40 & (pis > 1.0)),
            "found_but_dr_gt_0.4": fr(found & ~in40), "notfound_pt_gt_0.2": fr(~found & (pis > 0.2)),
            "veto_mask_removes_of_in40": float((w * (vetop & in40)[bp]).sum() / max((w * in40[bp]).sum(), 1e-9)),
            "dr_q50_90_99": wq(drp[bp & found], P["w"][bp & found], [0.5, 0.9, 0.99]),
            "contrib40_GeV_q10_25_50_75_90": wq(pis[bp & in40], P["w"][bp & in40], [0.1, 0.25, 0.5, 0.75, 0.9]),
            "contrib40_mean_GeV_all_tag": float((w * (pis * in40)[bp]).sum() / w.sum()),
            "ptcone40_mean_tag": float((w * P["c40"][bp]).sum() / w.sum()), "ptcone40_mean_tag_nopis": float((w * P["c40c"][bp]).sum() / w.sum()),
            "ptcone40_mean_b": float((A["w"][b & ba] * A["c40"][b & ba]).sum() / A["w"][b & ba].sum()),
            "ptcone40_median_tag_c_b": [wq(P["c40"][bp], w, [0.5])[0], wq(P["c40c"][bp], w, [0.5])[0], wq(A["c40"][b & ba], A["w"][b & ba], [0.5])[0]],
            "share_of_ptcone40_q25_50_75": wq((pis / np.clip(P["c40"], 1e-6, None))[bp & in40 & (P["c40"] > 0)], P["w"][bp & in40 & (P["c40"] > 0)], [0.25, 0.5, 0.75]),
            "cone_only_pis": float((w * (in40 & (P["c40c"] < 2e-3))[bp]).sum() / w.sum()),
            "cone_empty_tag": fr(P["c40"] < 2e-3), "cone_empty_tag_nopis": fr(P["c40c"] < 2e-3),
            "cone_empty_b": float((A["w"][b & ba] * (A["c40"][b & ba] < 2e-3)).sum() / A["w"][b & ba].sum()),
            "pis_crosses_2GeV": fr((P["c40"] >= 2.0) & (P["c40c"] < 2.0)),
        }
        # przeważenie próbki tagu do widma pT(D0) grupy b (20 binów log) — kontrola kinematyki
        ed = np.exp(np.linspace(np.log(lo), np.log(hi), 21))
        hb, _ = np.histogram(A["pt"][b & ba], ed, weights=A["w"][b & ba])
        ht, _ = np.histogram(P["pt"][bp], ed, weights=w)
        rw = (hb / np.clip(ht, 1e-9, None))[np.clip(np.digitize(P["pt"][bp], ed) - 1, 0, 19)]
        wrw = w * rw
        # --- pkt 2
        o["cuts"] = {}
        for cut in CUTS:
            pa = passcut(A["c40"], A["c100"], A["pt"], cut)
            pam = passcut(A["c40m"], A["c100m"], A["pt"], cut)
            pp = passcut(P["c40"], P["c100"], P["pt"], cut)
            pc = passcut(P["c40c"], P["c100c"], P["pt"], cut)
            px = passcut(P["c40x"], P["c100c"], P["pt"], cut)
            pm = passcut(P["c40m"], P["c100m"], P["pt"], cut)
            g = {"a": eff(pp[bp], w), "a_in_incl": eff(pa[sa & ba & in_tag], A["w"][sa & ba & in_tag]), "b": eff(pa[b & ba], A["w"][b & ba]),
                 "b1": eff(pa[b1 & ba], A["w"][b1 & ba]), "b2": eff(pa[b & ~b1 & ba], A["w"][b & ~b1 & ba]),
                 "c": eff(pc[bp], w), "c_max": eff(px[bp], w), "i": eff(pa[sa & ba], A["w"][sa & ba]),
                 "am": eff(pm[bp], w), "bm": eff(pam[b & ba], A["w"][b & ba]), "im": eff(pam[sa & ba], A["w"][sa & ba]),
                 "a_rw": eff(pp[bp], wrw), "c_rw": eff(pc[bp], wrw)}
            g["a/b"], g["c/b"], g["c/a"], g["a/i"], g["am/bm"], g["am/im"] = (ratio(g["a"], g["b"]), ratio(g["c"], g["b"]), ratio(g["c"], g["a"]),
                                                                             ratio(g["a"], g["i"]), ratio(g["am"], g["bm"]), ratio(g["am"], g["im"]))
            g["a_rw/b"], g["c_rw/b"], g["c/b2"], g["c_max/b"] = ratio(g["a_rw"], g["b"]), ratio(g["c_rw"], g["b"]), ratio(g["c"], g["b2"]), ratio(g["c_max"], g["b"])
            for tc in (0, 1):       # podział wg prawdy: 0 = c (prompt), 1 = b (z rozpadów hadronów b)
                tp_, ta_ = bp & (P["truth"] == tc), b & ba & (A["truth"] == tc)
                g["a_t%d" % tc], g["c_t%d" % tc], g["b_t%d" % tc] = eff(pp[tp_], P["w"][tp_]), eff(pc[tp_], P["w"][tp_]), eff(pa[ta_], A["w"][ta_])
                g["a/b_t%d" % tc], g["c/b_t%d" % tc] = ratio(g["a_t%d" % tc], g["b_t%d" % tc]), ratio(g["c_t%d" % tc], g["b_t%d" % tc])
            # --- pkt 3: efekt drugiego rzędu. „Dane” = MC z przesuniętą aktywnością spoza wolnego pionu; pion bez zmian (tak samo modelowany).
            #   tag:  u = ptcone − π_s(prawda);  inkluzywne: u = ptcone − tory maski (π_s bez prawdy).  Modele: (S) u → k·u;  (T) waga exp(α·u).
            _, R, kind, thr = cut
            uP, piP = (P["c40c"], pis * in40) if R == 40 else (P["c100c"], pis * in100)
            uA, piA = (A["c40m"], A["mk40"]) if R == 40 else (A["c100m"], A["mk100"])
            uP, piP, ptP = uP[bp], piP[bp], P["pt"][bp]
            sel = sa & ba
            uA, piA, ptA, wA, tgA = uA[sel], piA[sel], A["pt"][sel], A["w"][sel], in_tag[sel]
            dP, dA = (1.0, 1.0) if kind == "abs" else (ptP, ptA)

            def e_scale(k, u, pi, d, ww):
                return float((ww * ((k * u + pi) / d < thr)).sum() / ww.sum())

            def e_shift(dl, u, pi, d, ww):       # piedestał: u → u + δ [GeV] dla wszystkich kandydatów (także pustych stożków)
                return float((ww * ((u + dl + pi) / d < thr)).sum() / ww.sum())
            e0P, e0A = e_scale(1.0, uP, piP, dP, w), e_scale(1.0, uA, piA, dA, wA)
            e0B = e_scale(1.0, uA[~tgA], piA[~tgA], dA if kind == "abs" else dA[~tgA], wA[~tgA])
            so = {"eps0_tag": e0P, "eps0_incl": e0A, "eps0_b": e0B, "scale": {}, "shift": {}}
            for name, fe, lo_, hi_ in (("scale", e_scale, 1.0, 30.0), ("shift", e_shift, 0.0, 30.0)):
                for tg in TARGETS:
                    par = bisect(lambda x: fe(x, uP, piP, dP, w) / e0P, lo_, hi_, tg)
                    if par is None:
                        so[name]["%.2f" % tg] = None
                        continue
                    sft = fe(par, uP, piP, dP, w) / e0P
                    sfi = fe(par, uA, piA, dA, wA) / e0A
                    sfb = fe(par, uA[~tgA], piA[~tgA], dA if kind == "abs" else dA[~tgA], wA[~tgA]) / e0B
                    so[name]["%.2f" % tg] = {"par": par, "SF_tag": sft, "SF_incl": sfi, "SF_b": sfb, "SF_incl/SF_tag": sfi / sft, "SF_b/SF_tag": sfb / sft}
            g["second_order"] = so
            o["cuts"][cut[0]] = g
        # wrażliwość ε inkluzywnej na wydajność rekonstrukcji wolnego pionu (jedyny wkład pierwszego rzędu w próbce inkluzywnej):
        # Δε_i/ε_i przy utracie 10 % torów maski w b1 (cięcie ptcone40 < 2 GeV)
        cut = CUTS[0]
        e_i = eff(passcut(A["c40"], A["c100"], A["pt"], cut)[sa & ba], A["w"][sa & ba])
        e_im = eff(passcut(A["c40m"], A["c100m"], A["pt"], cut)[sa & ba], A["w"][sa & ba])
        o["pis_eff_sens_ptcone40lt2"] = {"eps_i": e_i[0], "eps_i_all_mask_removed": e_im[0], "d_eps_rel_per_10pct_pis_loss": 0.1 * (e_im[0] - e_i[0]) / e_i[0]}
        res["bins"]["%d-%d" % (lo, hi)] = o
        c = o["cuts"][CUTS[0][0]]
        print("[%s %d-%d] tag %d, b %d; π_s w 0,4: %.3f (pT>0,5: %.3f); ε(<2 GeV) a %.4f±%.4f  b %.4f±%.4f  c %.4f  i %.4f;  a/b %.4f±%.4f  c/b %.4f±%.4f"
              % (s, lo, hi, o["n_tag"], o["n_b"], o["pis"]["in40"], o["pis"]["in40_pt_gt_0.5"], c["a"][0], c["a"][1], c["b"][0], c["b"][1], c["c"][0], c["i"][0],
                 c["a/b"][0], c["a/b"][1], c["c/b"][0], c["c/b"][1]), flush=True)
    res["sec"] = round(time.time() - t0, 1)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True); ap.add_argument("--env", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--samples", default="fd_nd,mb_nd")
    a = ap.parse_args()
    out = {"m0": M0, "win": WIN, "dm_mask": DM_MASK, "cuts": [c[0] for c in CUTS], "samples": {}}
    for s in a.samples.split(","):
        out["samples"][s] = sample(a.data, a.env, s)
        json.dump(out, open(a.out, "w"), indent=1)
    print("zapisano", a.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
