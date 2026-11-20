#!/usr/bin/env python3
"""BKGCAT (20.11.2026, projekt B; POST HOC względem FITBIAS) — mb_nd: histogramy masa × przedział score drzewa H × kategoria prawdy kandydata
oraz tablica rozpadów z prawdy generatora (TruthDmesons). Jeden plik wyjściowy na plik wejściowy (NNNNNN); suma: bkg_merge.py.
Łańcuch: app etapu 4 (uid, m, pt, score) ↔ h5 v2.0 (ten sam uid, m, pt co do bitu; gałęzie prawdy kandydata) ↔ ROOT DMesons (entry, m) ↔ TruthDmesons (entry, truthPt).
Progi WP: z wps_hist.npz (WPSCAN), przedział s = liczba progów ≤ score (jak wps_hist.py); zakres mas i bin 1 MeV jak WPSCAN.
Kategoria (klucz tekstowy):
  Sig{C,B}|<tryb>  Refl{C,B}|<tryb>   truth 0/1, swap ≠ 1 / = 1 (D0)        <tryb> = topologia prawdziwego D (lista TruthDmesons_daughterPdgId:
  CM{C,B}|<|pdg|>|<tryb albo d:córki>  truth 2/3 (CharmMismatch)             tylko π±, K±, π0, D0) z sufiksem: '' (m_D − m_wid < 30 MeV), '+X0' (m_brak < 100 MeV),
  DMM{C,B}|<tryb>|d:<córki>            truth 6/7 (DecayModeMismatch)         '+K0', '+eta' (m_brak w ±3 MeV), '+X' (inne); 'nolink' gdy brak dopasowania
  PC|<rodzic>|{NoB,SomeB,AllB}         truth 8–11 (PartComb)                 d:<córki> = |truthPdgId| torów kandydata (posortowane)
  MB|<|pdg|>, ML|<|pdg|>               truth 4/5 (wszystkie tory z jednego hadronu b / lekkiego)
  OtherNoB, OtherSomeB, OtherAllB      truth 12, 13, 14–15
Tablica rozpadów: zliczenia TruthDmesons per (|pdg|, kod topologii, klasa braku, fromBottom, region) + histogram m_brak (5 MeV) dla |η| < 2,5.
Bez AFP, bez danych. Użycie: bkg_hist.py --fid NNNNNN --app DIR --h5 DIR --raw DIR --thr wps_hist.npz --out OUT.npz"""
import argparse
import glob
import json
import os
import re
import sys
import time

import numpy as np

BINS = [(3, 5), (5, 8), (8, 30)]
SP = {"D0": {"mr": (1.70, 2.05), "dt": 10}, "Dplus": {"mr": (1.70, 2.20), "dt": 4}}
TAG = {"D0": "D0.xgb.H.sb.n.s1", "Dplus": "Dplus.xgb.H.sb.n.s1"}
MBW = 0.001
MPI, MK, MPI0 = 139.57039, 493.677, 134.9768
PDGS = (421, 411, 413)
MCN = ["", "+X0", "+K0", "+eta", "+X"]


def mode_str(code):
    code = int(code)
    n = [(code // 10 ** i) % 10 for i in range(6)]            # K-, K+, pi+, pi-, pi0, D0
    s = "".join(x * k for x, k in zip(("D0 ", "K- ", "K+ ", "pi+ ", "pi- ", "pi0 "), (n[5], n[0], n[1], n[2], n[3], n[4]))).strip()
    return s if s else "(nic)"


def res(keys, name):
    for n in (name + "_NOSYS", name):
        if n in keys:
            return n
    raise KeyError(name)


def truth_tables(rootfile):
    """→ (T: prawdziwe D [entry, pdg, pt, eta, code, mc, fromB, mmiss], R: kandydaci reco ROOT [entry, dt, m_GeV, truthPt, truthPdg])."""
    import awkward as ak
    import uproot
    t = uproot.open(rootfile)["CharmAnalysis"]
    keys = set(t.keys())
    tb = ["pdgId", "decayType", "pt", "eta", "phi", "m", "fromBottom", "daughterPdgId", "daughterPt", "daughterEta", "daughterPhi"]
    rb = {k: res(keys, "DMesons_" + k) for k in ("m", "decayType", "truthPt", "truthPdgId")}
    A = t.arrays(["TruthDmesons_" + b for b in tb] + list(rb.values()))
    g = lambda b: A["TruthDmesons_" + b]  # noqa: E731
    cnt = ak.to_numpy(ak.num(g("pdgId")))
    entry = np.repeat(np.arange(len(cnt), dtype=np.int64), cnt)
    fl = lambda b: ak.to_numpy(ak.flatten(g(b))).astype(np.float64)  # noqa: E731
    pdg = ak.to_numpy(ak.flatten(g("pdgId"))).astype(np.int64)
    pt, eta, phi, m = fl("pt"), fl("eta"), fl("phi"), fl("m")
    did, dp, de, dph = (ak.flatten(g(b), axis=1) for b in ("daughterPdgId", "daughterPt", "daughterEta", "daughterPhi"))
    sg = ak.Array(np.sign(pdg))
    d = did * sg
    c = lambda mask: np.minimum(ak.to_numpy(ak.sum(mask, axis=1)), 9).astype(np.int64)  # noqa: E731
    nself = (np.abs(pdg) == 421).astype(np.int64)
    nd0 = np.clip(ak.to_numpy(ak.sum(abs(did) == 421, axis=1)) - nself, 0, 9)
    code = c(d == -321) + 10 * c(d == 321) + 100 * c(d == 211) + 1000 * c(d == -211) + 10000 * c(abs(did) == 111) + 100000 * nd0
    n_other = int(ak.sum(~((abs(did) == 321) | (abs(did) == 211) | (abs(did) == 111) | (abs(did) == 421))))
    vis = abs(did) != 421
    dm = ak.where(abs(did) == 321, MK, ak.where(abs(did) == 211, MPI, MPI0))
    dp, de, dph, dm = dp[vis], de[vis], dph[vis], dm[vis]
    px, py, pz = dp * np.cos(dph), dp * np.sin(dph), dp * np.sinh(de)
    E = np.sqrt(px * px + py * py + pz * pz + dm * dm)
    s = lambda x: ak.to_numpy(ak.sum(x, axis=1)).astype(np.float64)  # noqa: E731
    Ev, pxv, pyv, pzv = s(E), s(px), s(py), s(pz)
    mvis = np.sqrt(np.maximum(Ev ** 2 - pxv ** 2 - pyv ** 2 - pzv ** 2, 0.0))
    Px, Py, Pz = pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(np.clip(eta, -20, 20))
    ED = np.sqrt(Px ** 2 + Py ** 2 + Pz ** 2 + m ** 2)
    mm2 = (ED - Ev) ** 2 - (Px - pxv) ** 2 - (Py - pyv) ** 2 - (Pz - pzv) ** 2
    mmiss = np.sign(mm2) * np.sqrt(np.abs(mm2))
    mc = np.full(len(pdg), 4, np.int64)
    mc[mmiss < 100] = 1
    mc[np.abs(mmiss - 497.611) < 3] = 2
    mc[np.abs(mmiss - 547.862) < 3] = 3
    mc[(m - mvis) < 30] = 0
    T = {"entry": entry, "pdg": pdg, "pt": pt, "eta": eta, "code": code, "mc": mc, "fromB": ak.to_numpy(ak.flatten(g("fromBottom"))).astype(np.int64),
         "mmiss": mmiss, "dty": ak.to_numpy(ak.flatten(g("decayType"))).astype(np.int64), "n_other": n_other, "nev": len(cnt)}
    rc = ak.to_numpy(ak.num(A[rb["m"]]))
    R = {"entry": np.repeat(np.arange(len(rc), dtype=np.int64), rc), "m": ak.to_numpy(ak.flatten(A[rb["m"]])).astype(np.float64) / 1000.0,
         "dt": ak.to_numpy(ak.flatten(A[rb["decayType"]])).astype(np.int64), "tpt": ak.to_numpy(ak.flatten(A[rb["truthPt"]])).astype(np.float64),
         "tpdg": ak.to_numpy(ak.flatten(A[rb["truthPdgId"]])).astype(np.int64)}
    return T, R


def br_tables(T):
    """Zliczenia per (pdgidx, code, mc, fromB) w 3 regionach + histogram m_brak (|η| < 2,5, tylko mc ≠ 0)."""
    pi = np.select([np.abs(T["pdg"]) == p for p in PDGS], [0, 1, 2], default=3)
    comp = ((T["code"] * 8 + T["mc"]) * 2 + (T["fromB"] > 0)) * 4 + pi
    regs = [np.ones(len(pi), bool), np.abs(T["eta"]) < 2.5, (np.abs(T["eta"]) < 2.2) & (T["pt"] >= 3000.0)]
    o = {}
    for r, msk in enumerate(regs):
        u, n = np.unique(comp[msk], return_counts=True)
        o[f"br{r}_key"], o[f"br{r}_n"] = u, n.astype(np.int64)
    msk = regs[1] & (T["mc"] != 0)
    mb = np.clip(np.floor(T["mmiss"][msk] / 5.0), 0, 399).astype(np.int64)
    u, n = np.unique((T["code"][msk] * 4 + pi[msk]) * 400 + mb, return_counts=True)
    o["mm_key"], o["mm_n"] = u, n.astype(np.int64)
    # kontrola: decayType sygnału vs klasa braku (D0: 10, D+: 4) w regionie 1
    for p, dty in ((421, 10), (411, 4)):
        s = regs[1] & (np.abs(T["pdg"]) == p) & (T["dty"] == dty)
        o[f"dty{dty}_mc"] = np.bincount(T["mc"][s], minlength=5).astype(np.int64)
        o[f"dty{dty}_codes"] = np.unique(T["code"][s], return_counts=True)[0]
    o["pdg_other"] = np.array([int((pi == 3).sum())])
    return o


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--fid", required=True); ap.add_argument("--app", required=True); ap.add_argument("--h5", required=True)
    ap.add_argument("--raw", required=True); ap.add_argument("--thr", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    import h5py
    t0 = time.time()
    fa = glob.glob(os.path.join(a.app, f"app_mb_nd_*_{a.fid}.npz"))
    fr = glob.glob(os.path.join(a.raw, f"*._{a.fid}.*.root"))
    fh = os.path.join(a.h5, f"mb_nd_{a.fid}_none.h5")
    if len(fa) != 1 or len(fr) != 1 or not os.path.exists(fh):
        raise SystemExit(f"{a.fid}: app {len(fa)}, root {len(fr)}, h5 {os.path.exists(fh)}")
    Z = np.load(a.thr, allow_pickle=False)
    T, R = truth_tables(fr[0])
    out = br_tables(T)
    st = {"fid": a.fid, "nev": T["nev"], "n_truthD": int(len(T["pdg"])), "n_other_daughters": T["n_other"]}
    # klucze dopasowania
    tk = T["entry"] * 1e7 + T["pt"]
    to = np.argsort(tk); tks = tk[to]
    z = np.load(fa[0], allow_pickle=False)
    meta = json.loads(str(z["meta"]))
    with h5py.File(fh, "r") as f:
        ds = f["events"]
        cols = ["event_uid", "event_entry", "Dmeson_mass", "Dmeson_pt", "Dmeson_decayType", "Dmeson_truth_category", "Dmeson_hypo_swap", "Dmeson_truthPdgId",
                "Dmeson_truth_channel"] + [f"Dmeson_daughter{i}_truthPdgId" for i in range(3)]
        H = ds.fields(cols)[:]
    keys, kid = [], {}
    tab = []
    for isp, (sp, c) in enumerate(SP.items()):
        tags = meta["tags"].get(sp) or []
        if TAG[sp] not in tags:
            raise SystemExit(f"{a.fid}: brak modelu {TAG[sp]}")
        k = sp + "__"
        A = {x: z[k + x] for x in ("uid", "m", "pt", "ok", "truth", "swap", "w")}
        lg = z[k + "lg"][:, tags.index(TAG[sp])].astype(np.float64)
        h = H[H["Dmeson_decayType"] == c["dt"]]
        # app ↔ h5: (uid, bity m, bity pt)
        def key3(uid, m, pt):
            o = np.empty(len(uid), dtype=[("u", "u8"), ("m", "u4"), ("p", "u4")])
            o["u"], o["m"], o["p"] = uid, np.asarray(m, np.float32).view(np.uint32), np.asarray(pt, np.float32).view(np.uint32)
            return o
        kh = key3(h["event_uid"], h["Dmeson_mass"], h["Dmeson_pt"])
        ka = key3(A["uid"], A["m"], A["pt"])
        _, inv = np.unique(np.concatenate([kh, ka]), return_inverse=True)
        inv_h, inv_a = inv[:len(kh)], inv[len(kh):]
        first_h = np.full(int(inv.max()) + 1, -1, np.int64)
        first_h[inv_h] = np.arange(len(kh))
        ih = first_h[inv_a]
        hit = ih >= 0
        dup = np.bincount(inv_h, minlength=len(first_h))[inv_a] > 1
        n = len(ka)
        st[f"{sp}|app_rows"], st[f"{sp}|h5_hit"], st[f"{sp}|h5_dup"] = int(n), int(hit.sum()), int(dup.sum())
        hh = h[ih]
        tcat = hh["Dmeson_truth_category"].astype(np.int64)
        st[f"{sp}|truth_mismatch"] = int((hit & (tcat != A["truth"])).sum())
        # h5 ↔ ROOT reco (entry, m) → truthPt → TruthDmesons
        r = R["dt"] == c["dt"]
        rk = R["entry"][r] * 10.0 + R["m"][r]
        ro = np.argsort(rk); rks = rk[ro]
        hk = hh["event_entry"].astype(np.float64) * 10.0 + hh["Dmeson_mass"].astype(np.float64)
        p = np.clip(np.searchsorted(rks, hk), 1, max(len(rks) - 1, 1))
        p = np.where(np.abs(rks[p - 1] - hk) < np.abs(rks[np.minimum(p, len(rks) - 1)] - hk), p - 1, np.minimum(p, len(rks) - 1))
        rhit = hit & (np.abs(rks[p] - hk) < 3e-6)
        ir = np.flatnonzero(r)[ro][p]
        st[f"{sp}|root_hit"] = int(rhit.sum())
        qk = R["entry"][ir] * 1e7 + R["tpt"][ir]
        q = np.clip(np.searchsorted(tks, qk), 1, max(len(tks) - 1, 1))
        q = np.where(np.abs(tks[q - 1] - qk) < np.abs(tks[np.minimum(q, len(tks) - 1)] - qk), q - 1, np.minimum(q, len(tks) - 1))
        it = to[q]
        thit = rhit & (np.abs(tks[q] - qk) < 1e-2) & (np.abs(T["pdg"][it]) == np.abs(R["tpdg"][ir])) & (R["tpt"][ir] > 0)
        # tryb prawdziwego D
        code, mc = T["code"][it], T["mc"][it]
        um, inv = np.unique(code * 8 + mc, return_inverse=True)
        ms = np.array([mode_str(x // 8) + MCN[x % 8] for x in um], dtype=object)[inv]
        ms[~thit] = "nolink"
        dpdg = np.stack([np.abs(hh[f"Dmeson_daughter{i}_truthPdgId"].astype(np.int64)) for i in range(2 if sp == "D0" else 3)], axis=1)
        dpdg.sort(axis=1)
        ud, invd = np.unique(dpdg, axis=0, return_inverse=True)
        dstr = np.array(["d:" + ",".join(str(x) for x in row) for row in ud], dtype=object)[np.ravel(invd)]
        Dt = np.abs(hh["Dmeson_truthPdgId"].astype(np.int64))
        ch = hh["Dmeson_truth_channel"]
        uc, invc = np.unique(ch, return_inverse=True)
        pcl = []
        for s_ in uc:
            s_ = s_.decode() if isinstance(s_, bytes) else str(s_)
            mm = re.search(r"MisMatch(NoB|SomeB|AllB)(\d+)", s_)
            pcl.append(f"PC|{mm.group(2)}|{mm.group(1)}" if mm else "PC|?|?")
        pcs = np.array(pcl, dtype=object)[np.ravel(invc)]
        tr = A["truth"].astype(np.int64)
        cb = np.where(np.isin(tr, (0, 2, 6)), "C", "B").astype(object)
        key = np.full(n, "", dtype=object)
        sig = np.isin(tr, (0, 1))
        sw = sig & (A["swap"] == 1) if sp == "D0" else np.zeros(n, bool)
        S = lambda x: np.array([str(v) for v in x], dtype=object)  # noqa: E731
        key[sig & ~sw] = ("Sig" + cb + "|" + ms)[sig & ~sw]
        key[sw] = ("Refl" + cb + "|" + ms)[sw]
        cm = np.isin(tr, (2, 3))
        key[cm] = ("CM" + cb + "|" + S(Dt) + "|" + np.where(thit, ms, dstr))[cm]
        dmm = np.isin(tr, (6, 7))
        key[dmm] = ("DMM" + cb + "|" + ms + "|" + dstr)[dmm]
        pc = (tr >= 8) & (tr <= 11)
        key[pc] = pcs[pc]
        for tt, nm_ in ((4, "MB|"), (5, "ML|")):
            key[tr == tt] = (nm_ + S(Dt))[tr == tt]
        key[tr == 12] = "OtherNoB"; key[tr == 13] = "OtherSomeB"; key[(tr == 14) | (tr == 15)] = "OtherAllB"
        rest = key == ""
        key[rest] = ("t" + S(tr))[rest]
        key[~hit] = ("nomatch|t" + S(tr))[~hit]
        # histogramy
        lo, hi = c["mr"]
        nm = int(round((hi - lo) / MBW))
        ib = np.floor((A["m"].astype(np.float64) - lo) / MBW + 1e-9).astype(np.int64)
        base = A["ok"] & (ib >= 0) & (ib < nm) & (A["truth"] < 100)
        st[f"{sp}|truth_ge100"] = int((A["truth"] >= 100).sum()); st[f"{sp}|w_ne1"] = int((A["w"] != 1).sum())
        st[f"{sp}|link"] = {nm_: [int((m_ & base).sum()), int((m_ & base & thit).sum())] for nm_, m_ in (("sig", sig), ("cm", cm), ("dmm", dmm))}
        for jb, b in enumerate(BINS):
            thr = Z[f"{sp}|{TAG[sp]}|{b[0]}-{b[1]}|thr"]
            inb = base & (A["pt"] >= b[0]) & (A["pt"] < b[1])
            if not inb.any():
                continue
            s = np.searchsorted(thr, lg[inb], side="right")
            uk, ki = np.unique(key[inb], return_inverse=True)
            gk = np.array([kid.setdefault(x, len(kid)) for x in uk], np.int64)[ki]
            flat = (gk * 32 + s) * 1024 + ib[inb]
            u, iv = np.unique(flat, return_inverse=True)
            cn = np.bincount(iv, minlength=len(u)); cw = np.bincount(iv, weights=A["w"][inb].astype(np.float64), minlength=len(u))
            tab.append(np.stack([np.full(len(u), isp), np.full(len(u), jb), u // (32 * 1024), (u // 1024) % 32, u % 1024, cn], axis=1).astype(np.int64))
            out.setdefault("sw", []).append(cw)
    out["tab"] = np.concatenate(tab) if tab else np.zeros((0, 6), np.int64)        # kolumny: gatunek, bin, klucz, s, bin masy, n
    out["sw"] = np.concatenate(out["sw"]) if tab else np.zeros(0)
    out["keys"] = np.array(sorted(kid, key=kid.get))
    st["sec"] = round(time.time() - t0, 1)
    out["meta"] = np.array(json.dumps(st))
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    np.savez_compressed(a.out, **out)
    print(json.dumps(st), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
