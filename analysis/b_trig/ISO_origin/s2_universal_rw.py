#!/usr/bin/env python3
"""ISO_origin krok 2+4: (a) uniwersalność ε_iso(zmienna) między 427929-MB / 428770-J12only / 428770-MBTSbased (χ²),
(b) ten sam trigger w obu runach, (c) przeważenie MC do rozkładu SYGNAŁU w danych (wagi odejmowania pasm Δm) i SF per klasa,
(d) sprawdzenia: połowy LB, okno masy, próg score_H. Wynik: out/s2.json."""
from iso_common import *
from s1_data_eff import samples

INF = 1e9
UNI = {"e_nSelBL": [0, 25, 35, 45, 60, 80, INF], "e_ntrk500": [0, 30, 45, 60, 80, 110, INF], "e_lead_pt": [0, 4, 5, 6.5, 8, 12, INF],
       "pt": [3, 4, 5, 6.5, 8, 12, 30], "z": [0, 0.8, 1.0, 1.2, 1.5, INF]}
SCHEMES = {
    "brak": [],
    "nSelBL": [("e_nSelBL", [0, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 100, INF])],
    "ntrk500": [("e_ntrk500", [0, 20, 30, 40, 50, 60, 70, 80, 100, 120, INF])],
    "lead_pt": [("e_lead_pt", [0, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, INF])],
    "pT(D0)": [("pt", [3, 3.5, 4, 4.5, 5, 6, 7, 8, 10, 12, 15, 20, 30])],
    "pT×nSelBL": [("pt", [3, 4, 5, 6.5, 8, 12, 30]), ("e_nSelBL", [0, 25, 35, 45, 60, 80, INF])],
    "pT×lead_pt": [("pt", [3, 4, 5, 6.5, 8, 12, 30]), ("e_lead_pt", [0, 4, 6, 8, 12, 20, INF])],
    "pT×lead_pt×nSelBL": [("pt", [3, 5, 8, 12, 30]), ("e_lead_pt", [0, 5, 8, 12, INF]), ("e_nSelBL", [0, 35, 55, INF])],
    "pT×pierścień(ptcone100−40)": [("pt", [3, 4, 5, 6.5, 8, 12, 30]), ("ann", [0, 2, 4, 6, 9, 13, 20, INF])],
}
CUTS = [("ptcone40", 1.0), ("ptcone40", 2.0), ("ptcone40", 4.0), ("rel40", 0.5)]
CLS = ["9:MB", "9:MB&j20&!J12", "9:MB&!j20&!J12", "8:all", "8:J12only", "8:MBTSbased", "8:other"]
DISJ = ["9:MB", "8:J12only", "8:MBTSbased", "8:other"]


def base(D, mwin, sthr, lo=3, hi=30):
    return (np.abs(D["m"] - C.M0_D0) < mwin) & (D["score"] > sthr) & (D["pt"] >= lo) & (D["pt"] < hi)


def mc_sel(M, mwin, sthr, lo=3, hi=30):
    return np.isin(M["truth"], [0, 1]) & (np.abs(M["m"] - C.M0_D0) < mwin) & (M["dm"] > 0.1435) & (M["dm"] < 0.1475) & (M["score"] > sthr) & (M["pt"] >= lo) & (M["pt"] < hi)


def cell_index(X, scheme, upto=None):
    idx = np.zeros(len(X["pt"]), np.int64)
    for var, ed in scheme[:upto]:
        ed = np.asarray(ed, float)
        idx = idx * (len(ed) - 1) + np.clip(np.digitize(X[var], ed) - 1, 0, len(ed) - 2)
    return idx


def mc_eps_for(D, md, M, mm, pasM, scheme, min_mc=25):
    """ε_MC komórki dla każdego kandydata danych (z hierarchicznym cofaniem do mniejszej liczby wymiarów gdy n_MC < min_mc).
    Zwraca (ε, var_ε, maska cofnięcia)."""
    n = int(md.sum())
    eps = np.full(n, np.nan); var = np.zeros(n); fb = np.zeros(n, bool)
    Dm = {v: D[v][md] for v, _ in scheme} | {"pt": D["pt"][md]}
    Mm = {v: M[v][mm] for v, _ in scheme} | {"pt": M["pt"][mm]}
    pm = pasM[mm].astype(float)
    todo = np.ones(n, bool)
    for upto in range(len(scheme), -1, -1):
        if upto == 0:
            e0 = pm.mean(); eps[todo] = e0; var[todo] = e0 * (1 - e0) / len(pm); fb[todo] = len(scheme) > 0; break
        nc = int(np.prod([len(ed) - 1 for _, ed in scheme[:upto]]))
        im = cell_index(Mm, scheme, upto); idd = cell_index(Dm, scheme, upto)
        cnt = np.bincount(im, minlength=nc).astype(float); sp = np.bincount(im, weights=pm, minlength=nc)
        e = sp / np.maximum(cnt, 1)
        good = cnt[idd] >= min_mc
        take = todo & good
        eps[take] = e[idd[take]]; var[take] = (e * (1 - e) / np.maximum(cnt, 1))[idd[take]]
        if upto < len(scheme):
            fb[take] = True
        todo &= ~good
        if not todo.any():
            break
    return eps, var, fb


def sf_rw(D, md, pas, M, mm, pasM, scheme):
    """SF = Σw·pass / Σw·ε_MC(komórka) — MC przeważone do rozkładu sygnału w danych. Błąd: delta (dane) ⊕ stat. MC komórek."""
    w = sbs_w(D["dm"][md], D["pt"][md], pas[md])
    p = pas[md].astype(float)
    eps, var, fb = mc_eps_for(D, md, M, mm, pasM, scheme)
    A, B = (w * p).sum(), (w * eps).sum()
    if B <= 0 or w.sum() <= 0:
        return None
    sf = A / B
    v_data = (w ** 2 * (p - sf * eps) ** 2).sum() / B ** 2
    # stat. MC: komórki wspólne → sumuj wagi per unikalna (ε,var)
    key = np.round(eps, 12)
    u, inv = np.unique(key, return_inverse=True)
    wc = np.bincount(inv, weights=w); vc = np.zeros(len(u)); vc[inv] = var
    v_mc = (sf ** 2) * ((wc ** 2) * vc).sum() / B ** 2
    return {"SF": float(sf), "err": float(np.sqrt(v_data + v_mc)), "err_data": float(np.sqrt(v_data)), "err_mc": float(np.sqrt(v_mc)),
            "eps_data": float(A / w.sum()), "eps_mc_rw": float(B / w.sum()), "S": float(w.sum()), "frac_fallback": float((w * fb).sum() / w.sum())}


def chi2_pair(a, b):
    v = [(x["eps"] - y["eps"]) ** 2 / (x["err"] ** 2 + y["err"] ** 2) for x, y in zip(a, b)
         if x and y and np.isfinite(x["eps"]) and np.isfinite(y["eps"]) and x["S"] > 25 and y["S"] > 25 and x["err"] > 0 and y["err"] > 0]
    from scipy.stats import chi2 as C2
    return {"chi2": float(np.sum(v)), "ndf": len(v), "p": float(C2.sf(np.sum(v), len(v))) if v else None}


def run_all(M, S, mwin=0.040, sthr=0.6, lbpar=None, full=True):
    out = {"sf": {}, "uni": {}, "same_trig": {}, "eps": {}}
    def sel(name, lo=3, hi=30):
        D, cm = S[name]
        m = base(D, mwin, sthr, lo, hi) & cm
        if lbpar is not None:
            m &= (D["e_lb"].astype(int) % 2) == lbpar
        return D, m
    # (c) SF per klasa × cięcie × schemat
    for v, cut in CUTS:
        pasM = M[v] < cut
        for name in CLS:
            for pb, (lo, hi) in ({"3-30": (3, 30), "3-5": (3, 5), "5-8": (5, 8), "8-30": (8, 30)} if full else {"3-30": (3, 30)}).items():
                D, md = sel(name, lo, hi); mm = mc_sel(M, mwin, sthr, lo, hi)
                if md.sum() < 50:
                    continue
                pas = D[v] < cut
                for sn, sch in SCHEMES.items():
                    if not full and sn not in ("brak", "nSelBL", "pT×lead_pt", "pT×lead_pt×nSelBL"):
                        continue
                    r = sf_rw(D, md, pas, M, mm, pasM, sch)
                    if r:
                        out["sf"][f"{v}<{cut:g}|{name}|{pb}|{sn}"] = r
    # χ² domknięcia (SF vs 1) na klasach rozłącznych
    out["closure"] = {}
    for v, cut in CUTS:
        for sn in SCHEMES:
            rr = [out["sf"].get(f"{v}<{cut:g}|{c}|3-30|{sn}") for c in DISJ]
            if all(rr):
                x2 = sum(((r["SF"] - 1) / r["err"]) ** 2 for r in rr)
                cc = chi2_const([r["SF"] for r in rr], [r["err"] for r in rr])
                out["closure"][f"{v}<{cut:g}|{sn}"] = {"chi2_vs1": float(x2), "ndf": len(rr), "chi2_spojnosc": cc["chi2"], "p_spojnosc": cc["p"], "srednia": cc["mean"]}
    # (a) uniwersalność ε(zmienna): sygnał (wagi SBS) i tło (Δm>150, bez score)
    v, cut = "ptcone40", 2.0
    for var, ed in UNI.items():
        tab = {}
        for name in ["9:MB", "8:J12only", "8:MBTSbased"]:
            D, m = sel(name); pas = D[v] < cut
            w = sbs_w(D["dm"][m], D["pt"][m], pas[m]); xv = D[var][m]; pp = pas[m].astype(float)
            cm = S[name][1]; bb = cm & (np.abs(D["m"] - C.M0_D0) < mwin) & (D["pt"] >= 3) & (D["pt"] < 30) & (D["dm"] > 0.150)
            rows = []
            for lo, hi in zip(ed[:-1], ed[1:]):
                k = (xv >= lo) & (xv < hi)
                e, er, T, eT = eff_w(w[k], pp[k])
                kb = bb & (D[var] >= lo) & (D[var] < hi); nb = int(kb.sum()); fbk = float(pas[kb].mean()) if nb else np.nan
                rows.append({"lo": lo, "hi": hi, "eps": e, "err": er, "S": T, "S_err": eT, "bkg_f": fbk, "bkg_n": nb,
                             "bkg_err": float(np.sqrt(max(fbk * (1 - fbk), 1e-9) / nb)) if nb else np.nan})
            tab[name] = rows
        mm = mc_sel(M, mwin, sthr); rows = []
        for lo, hi in zip(ed[:-1], ed[1:]):
            k = mm & (M[var] >= lo) & (M[var] < hi); e = float((M[v][k] < cut).mean()) if k.sum() else np.nan
            rows.append({"lo": lo, "hi": hi, "eps": e, "err": float(np.sqrt(max(e * (1 - e), 0) / max(k.sum(), 1))), "S": int(k.sum())})
        tab["MC"] = rows
        ch = {}
        for a, b in (("9:MB", "8:J12only"), ("9:MB", "8:MBTSbased"), ("8:J12only", "8:MBTSbased"), ("9:MB", "MC"), ("8:J12only", "MC"), ("8:MBTSbased", "MC")):
            ch[f"{a} vs {b}"] = chi2_pair(tab[a], tab[b])
        for a, b in (("9:MB", "8:J12only"), ("9:MB", "8:MBTSbased"), ("8:J12only", "8:MBTSbased")):
            vv = [(x["bkg_f"] - y["bkg_f"]) ** 2 / (x["bkg_err"] ** 2 + y["bkg_err"] ** 2) for x, y in zip(tab[a], tab[b]) if x["bkg_n"] > 100 and y["bkg_n"] > 100]
            ch[f"TŁO {a} vs {b}"] = {"chi2": float(np.sum(vv)), "ndf": len(vv)}
        out["uni"][var] = {"table": tab, "chi2": ch}
    # (b) ten sam trigger w obu runach, w binach pT
    for a, b in (("9:MB&j20&!J12", "8:j20&!J12"), ("9:MB&J12", "8:J12any"), ("9:MB&sptrk4&!j20&!J12", "8:sptrk4&!j20&!J12")):
        rows = []
        for lo, hi in ((3, 30), (3, 5), (5, 8), (8, 30)):
            r = {"pt": f"{lo}-{hi}"}
            for nm in (a, b):
                D, m = sel(nm, lo, hi); pas = D[v] < cut
                if m.sum() > 20:
                    w = sbs_w(D["dm"][m], D["pt"][m], pas[m]); e, er, T, eT = eff_w(w, pas[m].astype(float))
                else:
                    e, er, T, eT = np.nan, np.nan, 0.0, 0.0
                cm = S[nm][1]; bb = cm & (np.abs(D["m"] - C.M0_D0) < mwin) & (D["pt"] >= lo) & (D["pt"] < hi) & (D["dm"] > 0.150)
                if lbpar is not None:
                    bb &= (D["e_lb"].astype(int) % 2) == lbpar
                nb = int(bb.sum()); f = float(pas[bb].mean()) if nb else np.nan
                r[nm] = {"eps": e, "err": er, "S": T, "S_err": eT, "bkg_f": f, "bkg_n": nb, "bkg_err": float(np.sqrt(max(f * (1 - f), 1e-9) / nb)) if nb else np.nan}
            rows.append(r)
        out["same_trig"][f"{a} || {b}"] = rows
    return out


def main():
    M = load_pairs("bt_fd_nd")
    D9, D8, S = samples()
    res = {"nominal": run_all(M, S)}
    for tag, kw in (("LB_parzyste", dict(lbpar=0)), ("LB_nieparzyste", dict(lbpar=1)), ("okno_m_25MeV", dict(mwin=0.025)), ("okno_m_55MeV", dict(mwin=0.055)),
                    ("score>0.3", dict(sthr=0.3)), ("score>0.8", dict(sthr=0.8)), ("score>0.0", dict(sthr=0.0))):
        res[tag] = run_all(M, S, full=False, **kw)
        print(tag, "ok", flush=True)
    # rozkłady sygnału (dane, wagi SBS) vs MC — do rysunku i tabeli
    dist = {}
    for var, ed in (("pt", [3, 4, 5, 6.5, 8, 12, 30]), ("e_nSelBL", [0, 25, 35, 45, 60, 80, INF]), ("e_lead_pt", [0, 4, 6, 8, 12, 20, INF])):
        dist[var] = {"edges": ed}
        for name in ["9:MB", "8:J12only", "8:MBTSbased", "8:other"]:
            D, cm = S[name]; m = base(D, 0.040, 0.6) & cm; pas = D["ptcone40"] < 2
            w = sbs_w(D["dm"][m], D["pt"][m], pas[m]); h = np.histogram(D[var][m], ed, weights=w)[0]; h2 = np.histogram(D[var][m], ed, weights=w ** 2)[0]
            dist[var][name] = {"frac": (h / h.sum()).tolist(), "err": (np.sqrt(h2) / h.sum()).tolist()}
        mm = mc_sel(M, 0.040, 0.6); h = np.histogram(M[var][mm], ed)[0]
        dist[var]["MC"] = {"frac": (h / h.sum()).tolist()}
    res["sig_dist"] = dist
    json.dump(res, open(os.path.join(HERE, "out", "s2.json"), "w"), indent=1, default=lambda o: float(o))
    print("zapisano out/s2.json")


if __name__ == "__main__":
    main()
