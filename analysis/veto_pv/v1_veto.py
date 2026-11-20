#!/usr/bin/env python3
"""Wydajność weta wierzchołków `nPileupVertices == 0` na próbce INKLUZYWNEJ (bez tagu protonu AFP) i na MC.
Wejście (tylko odczyt): B_trig/runs/lxplus/{bt,trg}_*.npz (h5 v2.0 → etap ORCH btrig, 17.11). Wyjście: out/veto.json, out/tables.md.
Nie używa żadnej zmiennej AFP (protony, ξ). Klasa trygera „other” (ani J12, ani MBTS; głównie łańcuchy L1AFP) — tylko ułamek zdarzeń-tła,
bez uzysku D*. Zmienne zdarzeniowe brane z wiersza kandydata (e_*, all_e_*): w fd_nd klucz (file_idx, event) tabeli ev_ nie zgadza się
z tabelami kandydatów (v0d_explore.py), więc tabela ev_ służy tylko do zliczeń łącznych."""
import sys, os, json, time
sys.dont_write_bytecode = True
import numpy as np
from scipy.optimize import minimize
HERE = os.path.dirname(os.path.abspath(__file__))
CW = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(CW, "B_trig", "ISO_origin"))
import iso_common as I  # noqa  (sbs_w: wagi odejmowania tła z pasm Δm; C.classify: klasy trygera)
C = I.C
R = os.path.join(CW, "B_trig", "runs", "lxplus")
OUT = os.path.join(HERE, "out"); os.makedirs(OUT, exist_ok=True)
RNG = np.random.default_rng(20261118)
NB = 300
MWIN, STHR = 0.040, 0.6                      # |m(Kπ) − m_D0| < 40 MeV, H > 0,6 (jak B_trig/ISO_origin „base”)
PT_BINS = [(3, 5), (5, 8), (8, 30)]
NSEL_E = [0, 10, 20, 30, 40, 50, 65, 80, 1000]
LXY_E = [0.25, 0.4, 0.6, 1.0, 2.0, 1e3]
MU_E = {428770: [0, 0.025, 0.035, 0.045, 0.055, 0.065, 0.08, 0.12, 0.3], 427929: [0, 0.003, 0.004, 0.005, 0.007, 0.012, 0.02, 0.03]}
LB_E = {428770: [0, 100, 250, 290, 400], 427929: [119, 170, 220, 270, 325]}
res = {}
md = []


def P(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); md.append(s)


def binom(k, n):
    if n <= 0:
        return (np.nan, np.nan, 0)
    p = k / n
    v = p * (1 - p) if 0 < k < n else 1.0 / (n + 1)        # k = 0 lub n: ~1/n zamiast zera
    return (float(p), float(np.sqrt(v / n)), int(n))


def ev_frac(passed, mask=None):
    """Ułamek zdarzeń przechodzących weto; błąd dwumianowy (równoważny bootstrapowi po zdarzeniach dla wag jednostkowych — test w logu)."""
    if mask is not None:
        passed = passed[mask]
    return binom(int(passed.sum()), len(passed))


def boot_eps(uid, w, pas, raw=False):
    """ε = Σ w·pass / Σ w po kandydatach; bootstrap Poissona po ZDARZENIACH (wagi SBS stałe)."""
    if len(w) == 0 or w.sum() < 25:
        return (np.nan, np.nan, float(w.sum()) if len(w) else 0.0)
    u, inv = np.unique(uid, return_inverse=True)
    W = np.bincount(inv, w, len(u)); Wp = np.bincount(inv, w * pas, len(u))
    e0 = Wp.sum() / W.sum()
    g = RNG.poisson(1.0, (NB, len(u)))
    eb = (g @ Wp) / np.maximum(g @ W, 1e-9)
    if raw:
        return e0, eb
    return (float(e0), float(eb.std(ddof=1)), float(W.sum()))


def fmt(t):
    e, s, n = t
    return "–" if not np.isfinite(e) else f"{100*e:.2f} ± {100*s:.2f}"


def ff(f):
    return [f"{100*f['f']:.2f} ± {100*f['f_err']:.2f}", f"{f['eps_v']:.3f} ± {f['eps_v_err']:.3f}", f"{f['rho']:.2f}", f"{f['chi2']:.1f}/{f['ndf']}"]


def uid_of(fi, ev, ns, lp):
    """Identyfikator zdarzenia odporny na powtarzające się numery zdarzeń w MC (fd_nd: 1…10000 w każdym pliku)."""
    return (fi.astype(np.int64) << 44) ^ (ev.astype(np.int64) << 12) ^ (ns.astype(np.int64) << 2) ^ (np.round(lp * 1000).astype(np.int64) << 24)


def tab(z, pre, extra=()):
    d = {"uid": uid_of(z[pre + "file_idx"], z[pre + "event"], z[pre + "e_nSelBL"], z[pre + "e_lead_pt"])}
    for k in ("npu", "actual_mu", "lb", "grl", "nSelBL"):
        d[k] = z[pre + "e_" + k]
    for k in z.files:
        if k.startswith(pre + "e_tp_"):
            d[k[len(pre) + 2:]] = z[k]
    for k in extra:
        d[k] = z[pre + k]
    d["pas"] = d["npu"] == 0
    return d


def pois_fit(mu, k, n):
    """k_i ~ Bin(n_i, (1−f)·exp(−ε_v·μ_i)); zwraca f, ε_v, błędy (hesjan numeryczny), χ²/ndf."""
    mu, k, n = map(lambda a: np.asarray(a, float), (mu, k, n))

    def nll(p):
        q = np.clip((1 - p[0]) * np.exp(-p[1] * mu), 1e-9, 1 - 1e-9)
        return -(k * np.log(q) + (n - k) * np.log(1 - q)).sum()
    r = minimize(nll, [0.01, 0.8], method="Nelder-Mead", options={"xatol": 1e-7, "fatol": 1e-7, "maxiter": 4000})
    p = r.x; h = np.array([1e-4, 1e-3]); H = np.zeros((2, 2))
    for i in range(2):
        for j in range(2):
            ei, ej = np.eye(2)[i] * h[i], np.eye(2)[j] * h[j]
            H[i, j] = (nll(p + ei + ej) - nll(p + ei - ej) - nll(p - ei + ej) + nll(p - ei - ej)) / (4 * h[i] * h[j])
    try:
        cov = np.linalg.inv(H); err = np.sqrt(np.maximum(np.diag(cov), 0)); rho = cov[0, 1] / max(err[0] * err[1], 1e-30)
    except np.linalg.LinAlgError:
        err = np.array([np.nan, np.nan]); rho = np.nan
    q = (1 - p[0]) * np.exp(-p[1] * mu); chi2 = (((k - n * q) ** 2) / np.maximum(n * q * (1 - q), 1e-9)).sum()
    return {"f": float(p[0]), "f_err": float(err[0]), "eps_v": float(p[1]), "eps_v_err": float(err[1]), "rho": float(rho),
            "chi2": float(chi2), "ndf": int(len(mu) - 2)}


def mu_bins_fit(mu, pas, nb=16):
    ed = np.unique(np.quantile(mu, np.linspace(0, 1, nb + 1))); ed[-1] += 1e-9
    ib = np.clip(np.digitize(mu, ed) - 1, 0, len(ed) - 2)
    n = np.bincount(ib, minlength=len(ed) - 1); k = np.bincount(ib, pas, len(ed) - 1); m = np.bincount(ib, mu, len(ed) - 1) / np.maximum(n, 1)
    ok = n > 50
    return pois_fit(m[ok], k[ok], n[ok])


def table(title, rows, header):
    P(f"\n**{title}**\n"); P("| " + " | ".join(header) + " |"); P("|" + "---|" * len(header))
    for r in rows:
        P("| " + " | ".join(str(x) for x in r) + " |")


t0 = time.time()
for s in ["data_428770", "data_427929", "fd_nd", "mb_nd"]:
    z = np.load(os.path.join(R, f"bt_{s}.npz"), allow_pickle=False)
    isdata = s.startswith("data"); run = int(s.split("_")[1]) if isdata else 0
    rs = res[s] = {}
    P(f"\n## {s}\n")
    Pp = tab(z, "", ("m", "pt", "score", "dm", "Lxy", "truth", "truth_ds"))
    A = tab(z, "all_", ("pt",))
    # zdarzenia z kandydatem D0 (tabela all_: okno 1,75–2,00 GeV, pT ≥ 3 GeV, bez cięcia H) — jeden wiersz na zdarzenie
    _, first = np.unique(A["uid"], return_index=True)
    E = {k: v[first] for k, v in A.items()}

    def good_of(T):
        return ((T["grl"] == 1) & (T["actual_mu"] > 0)) if isdata else np.ones(len(T["uid"]), bool)

    def classes(T):
        if not isdata:
            return {"MC": np.ones(len(T["uid"]), bool)}
        c = C.classify(T, run); c.pop("_missing", None); return c
    gE, gP, gA = good_of(E), good_of(Pp), good_of(A)
    cE, cP, cA = classes(E), classes(Pp), classes(A)
    use = {428770: ["all", "J12only", "MBTSbased", "j20MBTS", "sptrk4", "other"], 427929: ["all427929", "MB427929", "sptrk4_427929"], 0: ["MC"]}[run]
    base = (np.abs(Pp["m"] - C.M0_D0) < MWIN) & (Pp["score"] > STHR) & (Pp["pt"] >= 3) & (Pp["pt"] < 30) & gP
    dm = Pp["dm"] * 1e3; pas_p = Pp["pas"].astype(float)
    win = (dm > I.SIG[0]) & (dm < I.SIG[1]); sbu = (dm > 148.5) & (dm < 165.0)
    truS = (np.isin(Pp["truth_ds"], [0, 1]) & np.isin(Pp["truth"], [0, 1])) if not isdata else None
    # tabela ev_ (zdarzenia z jakimkolwiek kandydatem h5) — tylko zliczenia łączne
    evg = ((z["ev_grl"] == 1) & (z["ev_actual_mu"] > 0)) if isdata else np.ones(len(z["ev_npu"]), bool)
    n_ev, k_ev = int(evg.sum()), int((z["ev_npu"][evg] == 0).sum())
    oth = binom(k_ev - int(E["pas"][gE].sum()), n_ev - int(gE.sum()))
    rs.update({"n_events_anycand": n_ev, "n_events_D0": int(gE.sum()), "n_pairs_base": int(base.sum()), "eps_anycand_notD0": oth})
    P(f"zdarzenia z jakimkolwiek kandydatem h5 (GRL, μ>0): {n_ev}; z kandydatem D0 (okno 1,75–2,00 GeV, pT ≥ 3 GeV): {gE.sum()}; pary D0←D* po cięciach bazowych: {base.sum()}")
    P(f"ε_veto zdarzeń z innym kandydatem h5 (bez D0): {fmt(oth)} %")
    if isdata:
        mu = E["actual_mu"][gE]
        P(f"μ (actual, zdarzenia z D0): średnia {mu.mean():.4f}, mediana {np.median(mu):.4f}, 1–99 %: {np.percentile(mu, 1):.4f}–{np.percentile(mu, 99):.4f}")
        P("rozkład nPileupVertices (zdarzenia z D0):", np.bincount(E["npu"][gE].astype(int))[:8].tolist())
        b0 = boot_eps(E["uid"][gE][:200000], np.ones(min(gE.sum(), 200000)), E["pas"][gE][:200000].astype(float)); a0 = ev_frac(E["pas"][gE][:200000])
        P(f"test błędu: bootstrap po zdarzeniach {100*b0[1]:.4f} % vs dwumianowy {100*a0[1]:.4f} % (200 tys. zdarzeń)")

    def sbs_weights(mask):
        w = np.zeros(len(dm)); w[mask] = I.sbs_w(Pp["dm"][mask], Pp["pt"][mask], pas_p[mask] > 0); return w

    def sig_eps(mask, w=None):
        """Sygnał D*: dane — SBS (wagi z fitu progowego pasm, grupy pass/fail × pT); MC — prawda w oknie Δm."""
        if not isdata:
            mm = mask & truS & win
            return boot_eps(Pp["uid"][mm], np.ones(mm.sum()), pas_p[mm])
        return boot_eps(Pp["uid"][mask], w[mask], pas_p[mask])

    def sb_eps(mask):
        mm = mask & sbu & ((~truS) if not isdata else True)
        _, fi = np.unique(Pp["uid"][mm], return_index=True); return ev_frac(Pp["pas"][mm][fi])

    # ---------------- 1. per klasa ----------------
    rows = []; W = {}
    for c in use:
        me = cE[c] & gE; mp = base & cP[c]
        W[c] = sbs_weights(mp) if (isdata and c != "other") else None
        sg = (np.nan, np.nan, 0.0) if c == "other" else sig_eps(mp, W[c])
        r = {"signal": sg, "sideband": (np.nan, np.nan, 0) if c == "other" else sb_eps(mp), "D0cand_events": ev_frac(E["pas"], me),
             "mu_mean": float(E["actual_mu"][me].mean()), "nSelBL_mean": float(E["nSelBL"][me].mean())}
        if isdata:
            r["pois_D0cand"] = mu_bins_fit(E["actual_mu"][me], E["pas"][me])
            f = r["pois_D0cand"]; r["C_PU"] = float(np.exp(-f["eps_v"] * E["actual_mu"][me]).mean())
            r["C_PU_err"] = float(f["eps_v_err"] * (E["actual_mu"][me] * np.exp(-f["eps_v"] * E["actual_mu"][me])).mean())
        rs[c] = r
        rows.append([c + (" (tylko tło)" if c == "other" else ""), f"{r['mu_mean']:.4f}", f"{r['nSelBL_mean']:.1f}", f"{sg[2]:.0f}" if c != "other" else "–", fmt(sg), fmt(r["sideband"]),
                     fmt(r["D0cand_events"]), r["D0cand_events"][2]])
    table("ε_veto [%] per klasa (zdarzenia; sygnał = D0←D* " + ("po odjęciu pasm Δm" if isdata else "z prawdy") + ")", rows,
          ["klasa", "⟨μ⟩", "⟨nSelBL⟩", "N_S", "sygnał D*", "pasma Δm", "zdarzenia z D0 (tło)", "N zd."])
    c0 = use[0]; me0 = cE[c0] & gE; mp0 = base & cP[c0]; w0 = W[c0]

    # ---------------- zależności (klasa łączna) ----------------
    def dep(name, key, edges, lab):
        rows = []; rr = rs.setdefault("dep_" + name, {})
        for a, b in zip(edges[:-1], edges[1:]):
            me = me0 & (E[key] >= a) & (E[key] < b); mp = mp0 & (Pp[key] >= a) & (Pp[key] < b)
            r = {"signal": sig_eps(mp, w0), "sideband": sb_eps(mp), "D0cand_events": ev_frac(E["pas"], me), "x_mean": float(E[key][me].mean()) if me.any() else np.nan,
                 "mu_mean": float(E["actual_mu"][me].mean()) if me.any() else np.nan,
                 "nS": float(w0[mp].sum()) if isdata else float((mp & truS & win).sum()), "nSB": int((mp & sbu).sum())}
            rr[f"{a}-{b}"] = r
            rows.append([f"{a}–{b}", f"{r['x_mean']:.4g}", f"{r['mu_mean']:.4f}", f"{r['nS']:.0f}", fmt(r["signal"]), fmt(r["sideband"]), fmt(r["D0cand_events"]), r["D0cand_events"][2]])
        table(f"ε_veto [%] w funkcji {lab} (klasa {c0})", rows, [lab, "⟨x⟩", "⟨μ⟩", "N_S", "sygnał D*", "pasma Δm", "zdarzenia z D0", "N zd."])

    if isdata:
        dep("mu", "actual_mu", MU_E[run], "μ")
        dep("lb", "lb", LB_E[run], "LB")
    dep("nsel", "nSelBL", NSEL_E, "nSelBL")
    # standaryzacja: ε zdarzeń z D0 przeważone do rozkładu nSelBL sygnału / pasm → czy sam kandydat zmienia weto przy stałej krotności
    dn = rs["dep_nsel"]; eD = np.array([v["D0cand_events"][0] for v in dn.values()]); sD = np.array([v["D0cand_events"][1] for v in dn.values()])
    for lab, key in (("signal", "nS"), ("sideband", "nSB")):
        wv = np.clip(np.array([v[key] for v in dn.values()], float), 0, None); wv /= wv.sum(); okb = np.isfinite(eD)
        pred = float((wv[okb] * eD[okb]).sum() / wv[okb].sum()); perr = float(np.sqrt(((wv[okb] * sD[okb]) ** 2).sum()) / wv[okb].sum())
        meas = rs[c0][lab]; rs[f"std_{lab}"] = {"measured": meas, "pred_from_D0events": (pred, perr), "diff": (meas[0] - pred, float(np.hypot(meas[1], perr)))}
        P(f"\nStandaryzacja nSelBL ({lab}, klasa {c0}): zmierzone {fmt(meas)} %, przewidziane ze zdarzeń z D0 przy tym samym nSelBL {100*pred:.2f} ± {100*perr:.2f} %, "
          f"różnica {100*(meas[0]-pred):+.2f} ± {100*np.hypot(meas[1], perr):.2f} %")
    # pT i L_xy kandydata
    rows = []; rr = rs.setdefault("dep_pt", {})
    for a, b in PT_BINS:
        mp = mp0 & (Pp["pt"] >= a) & (Pp["pt"] < b)
        ma = (A["pt"] >= a) & (A["pt"] < b) & gA & cA[c0]; _, fi = np.unique(A["uid"][ma], return_index=True)
        r = {"signal": sig_eps(mp, w0), "sideband": sb_eps(mp), "D0cand_events": ev_frac(A["pas"][ma][fi])}
        if isdata:
            r["pois_D0cand"] = mu_bins_fit(A["actual_mu"][ma][fi], A["pas"][ma][fi], 10)
        rr[f"{a}-{b}"] = r; rows.append([f"{a}–{b}", f"{r['signal'][2]:.0f}", fmt(r["signal"]), fmt(r["sideband"]), fmt(r["D0cand_events"]), r["D0cand_events"][2]])
    table(f"ε_veto [%] w funkcji pT(D0) [GeV] (klasa {c0})", rows, ["pT", "N_S", "sygnał D*", "pasma Δm", "zdarzenia z D0 o tym pT", "N zd."])
    rows = []; rr = rs.setdefault("dep_lxy", {})
    for a, b in zip(LXY_E[:-1], LXY_E[1:]):
        mp = mp0 & (Pp["Lxy"] >= a) & (Pp["Lxy"] < b)
        r = {"signal": sig_eps(mp, w0), "sideband": sb_eps(mp)}
        rr[f"{a}-{b}"] = r; rows.append([f"{a}–{b}", f"{r['signal'][2]:.0f}", fmt(r["signal"]), fmt(r["sideband"]), r["sideband"][2]])
    v = np.array([x["signal"][0] for x in rr.values()]); e = np.array([x["signal"][1] for x in rr.values()]); okb = np.isfinite(v) & (e > 0)
    wm = (v[okb] / e[okb] ** 2).sum() / (1 / e[okb] ** 2).sum(); chi = (((v[okb] - wm) / e[okb]) ** 2).sum(); rs["lxy_const_chi2"] = (float(chi), int(okb.sum() - 1))
    table(f"ε_veto [%] w funkcji L_xy(D0) [mm] (klasa {c0}); sygnał vs stała: χ²/ndf = {chi:.1f}/{okb.sum()-1}", rows, ["L_xy", "N_S", "sygnał D*", "pasma Δm", "N zd. pasm"])

    # ---------------- Poisson: (1−f)·exp(−ε_v·μ) ----------------
    if isdata:
        rows = []
        for c in use:
            r = rs[c]; rows.append([c, r["D0cand_events"][2]] + ff(r["pois_D0cand"]) + [f"{100*(1-r['C_PU']):.2f} ± {100*r['C_PU_err']:.2f}"])
        rr = rs.setdefault("pois_nsel", {})
        for a, b in zip(NSEL_E[:-1], NSEL_E[1:]):
            me = me0 & (E["nSelBL"] >= a) & (E["nSelBL"] < b)
            if me.sum() < 5000:
                continue
            f = mu_bins_fit(E["actual_mu"][me], E["pas"][me], 10); rr[f"{a}-{b}"] = f
            rows.append([f"{c0}, nSelBL {a}–{b}", int(me.sum())] + ff(f) + ["–"])
        for k, f in rs["dep_pt"].items():
            rows.append([f"{c0}, pT {k}", f["D0cand_events"][2]] + ff(f["pois_D0cand"]) + ["–"])
        table("Dopasowanie P(npu = 0 | μ) = (1 − f)·exp(−ε_v·μ), zdarzenia z kandydatem D0 (16 lub 10 binów kwantylowych μ)", rows,
              ["próbka", "N zd.", "f [%] (μ→0)", "ε_v", "ρ(f,ε_v)", "χ²/ndf", "strata pile-up 1 − ⟨exp(−ε_v μ)⟩ [%]"])
        dmu = rs["dep_mu"]; xs = np.array([v["mu_mean"] for v in dmu.values()]); ys = np.array([v["signal"][0] for v in dmu.values()]); es = np.array([v["signal"][1] for v in dmu.values()])
        okb = np.isfinite(ys) & (es > 0); f = rs[c0]["pois_D0cand"]
        pred = (1 - f["f"]) * np.exp(-f["eps_v"] * xs[okb]); chi = (((ys[okb] - pred) / es[okb]) ** 2).sum(); rs["signal_vs_model_mu"] = (float(chi), int(okb.sum()))
        P(f"\nSygnał D* (SBS) w binach μ wobec modelu dopasowanego do zdarzeń z D0 (bez swobodnych parametrów): χ²/n = {chi:.1f}/{okb.sum()}")
        # ---------------- wszystkie zdarzenia w drzewie (trg_), w tym bez kandydata ----------------
        zt = np.load(os.path.join(R, f"trg_{s}.npz"), allow_pickle=False)
        T = {k: zt[k] for k in ("event", "npu", "actual_mu", "grl", "lb", "nSelBL")}
        T.update({k: zt[k] for k in zt.files if k.startswith("tp_")})
        tg = (T["grl"] == 1) & (T["actual_mu"] > 0); tp = T["npu"] == 0
        incand = np.isin(T["event"], z["ev_event"]); inD0 = np.isin(T["event"], np.unique(z["all_event"]))
        tcls = C.classify(T, run)
        rows = []; rr = rs.setdefault("tree", {})
        for c in [u for u in use if u != "other"]:
            m = tcls[c] & tg; nc = m & ~incand
            r = {"all": ev_frac(tp, m), "no_cand": ev_frac(tp, nc), "D0": ev_frac(tp, m & inD0), "nsel_nocand": float(T["nSelBL"][nc].mean()),
                 "pois_nocand": mu_bins_fit(T["actual_mu"][nc], tp[nc]),
                 "pois_nocand_nsel10": mu_bins_fit(T["actual_mu"][nc & (T["nSelBL"] >= 10)], tp[nc & (T["nSelBL"] >= 10)])}
            rr[c] = r; a, b = r["pois_nocand"], r["pois_nocand_nsel10"]
            rows.append([c, r["all"][2], fmt(r["all"]), fmt(r["D0"]), fmt(r["no_cand"]), r["no_cand"][2], f"{r['nsel_nocand']:.1f}",
                         f"{100*a['f']:.2f} ± {100*a['f_err']:.2f} / {a['eps_v']:.3f} ± {a['eps_v_err']:.3f}", f"{100*b['f']:.2f} ± {100*b['f_err']:.2f} / {b['eps_v']:.3f} ± {b['eps_v_err']:.3f}"])
        table("Wszystkie zdarzenia w drzewie (GRL): ε_veto [%] z kandydatem i bez", rows,
              ["klasa", "N zd.", "wszystkie", "z D0", "bez żadnego kandydata h5", "N bez kand.", "⟨nSelBL⟩ bez kand.", "f [%] / ε_v (bez kand.)", "f [%] / ε_v (bez kand., nSelBL ≥ 10)"])
        rows = []; rr = rs.setdefault("tree_nsel", {}); m0 = tcls[c0] & tg
        for a, b in zip(NSEL_E[:-1], NSEL_E[1:]):
            mb = m0 & (T["nSelBL"] >= a) & (T["nSelBL"] < b); k = f"{a}-{b}"
            r = {"D0": ev_frac(tp, mb & inD0), "no_cand": ev_frac(tp, mb & ~incand)}
            rr[k] = r; rows.append([f"{a}–{b}", fmt(rs["dep_nsel"][k]["signal"]), fmt(rs["dep_nsel"][k]["sideband"]), fmt(r["D0"]), fmt(r["no_cand"]), r["no_cand"][2]])
        table(f"Zależność od kandydata przy stałej krotności (klasa {c0})", rows, ["nSelBL", "sygnał D*", "pasma Δm", "zdarzenia z D0", "bez kandydata", "N bez kand."])
        # populacja wielowierzchołkowa niezależna od μ (tło wiązki)
        t3 = tg & (T["npu"] >= 3)
        rs["npu_ge3"] = {"n": int(t3.sum()), "frac": float(t3.sum() / tg.sum()), "nSelBL_median": float(np.median(T["nSelBL"][t3])), "frac_in_D0_events": float((E["npu"][gE] >= 3).mean())}
        P(f"\nZdarzenia z npu ≥ 3 (GRL): {t3.sum()} = {100*t3.sum()/tg.sum():.3f} % drzewa, mediana nSelBL = {np.median(T['nSelBL'][t3]):.0f}; wśród zdarzeń z D0: {100*(E['npu'][gE] >= 3).mean():.4f} %")
        if run == 427929:
            hm = (T["grl"] == 0) & (T["actual_mu"] > 5)
            rs["highmu_nonGRL"] = {"n": int(hm.sum()), "npu_mean": float(T["npu"][hm].mean()), "mu_mean": float(T["actual_mu"][hm].mean())}
            P(f"Poza GRL, μ > 5 (informacyjnie): N = {hm.sum()}, ⟨npu⟩ = {T['npu'][hm].mean():.2f}, ⟨μ⟩ = {T['actual_mu'][hm].mean():.2f}, ⟨npu⟩/⟨μ⟩ = {T['npu'][hm].mean()/T['actual_mu'][hm].mean():.3f}")
        del T, zt
    P(f"[{s}: {time.time()-t0:.0f} s]")

# ---------------- dane vs MC ----------------
P("\n## Dane/MC\n")
rows = []; cmp_ = res["cmp_f_nsel"] = {}
for a, b in zip(NSEL_E[:-1], NSEL_E[1:]):
    k = f"{a}-{b}"; row = [f"{a}–{b}"]; cc = cmp_[k] = {}
    for s in ["mb_nd", "fd_nd"]:
        v = res[s]["dep_nsel"][k]["D0cand_events"]; cc[s] = (1 - v[0], v[1], v[2]); row.append(f"{100*(1-v[0]):.2f} ± {100*v[1]:.2f}")
    v = res["fd_nd"]["dep_nsel"][k]["signal"]; row.append(f"{100*(1-v[0]):.2f} ± {100*v[1]:.2f}" if np.isfinite(v[0]) else "–")
    for s in ["data_428770", "data_427929"]:
        f = res[s].get("pois_nsel", {}).get(k)
        if f:
            cc[s] = (f["f"], f["f_err"]); mc = cc["mb_nd"]; row.append(f"{100*f['f']:.2f} ± {100*f['f_err']:.2f}")
            row.append(f"{f['f']/mc[0]:.2f} ± {np.hypot(f['f_err']/mc[0], f['f']*mc[1]/mc[0]**2):.2f}" if mc[0] > 0 else "–")
        else:
            row += ["–", "–"]
    rows.append(row)
table("f = 1 − ε_veto bez pile-up [%] w binach nSelBL: MC wprost (μ = 0), dane z ekstrapolacji μ → 0 (zdarzenia z D0)", rows,
      ["nSelBL", "mb_nd zd. z D0", "fd_nd zd. z D0", "fd_nd sygnał D* (prawda)", "428770 (fit)", "428770/mb_nd", "427929 (fit)", "427929/mb_nd"])
rows = []
mcS = res["fd_nd"]["MC"]["signal"]; mcB = res["mb_nd"]["MC"]["D0cand_events"]
for s, cl in (("data_428770", ["all", "J12only", "MBTSbased"]), ("data_427929", ["all427929"])):
    for c in cl:
        r = res[s][c]; d = r["D0cand_events"]; sg = r["signal"]; f = r["pois_D0cand"]
        r["SF_events"] = (d[0] / mcB[0], float(np.hypot(d[1] / mcB[0], d[0] * mcB[1] / mcB[0] ** 2)))
        r["SF_signal"] = (sg[0] / mcS[0], float(np.hypot(sg[1] / mcS[0], sg[0] * mcS[1] / mcS[0] ** 2)))
        rows.append([s[5:], c, fmt(d), fmt(mcB), f"{r['SF_events'][0]:.4f} ± {r['SF_events'][1]:.4f}", fmt(sg), fmt(mcS), f"{r['SF_signal'][0]:.3f} ± {r['SF_signal'][1]:.3f}",
                     f"{100*f['f']:.2f} / {100*(1-mcB[0]):.2f}", f"{r['C_PU']:.4f} ± {r['C_PU_err']:.4f}"])
table("ε_veto dane vs MC i rozkład na człony", rows, ["run", "klasa", "dane: zd. z D0", "mb_nd: zd. z D0", "dane/MC (zd.)", "dane: sygnał D* (SBS)", "fd_nd: sygnał D* (prawda)", "dane/MC (sygnał)",
                                                     "f dane / f MC [%]", "C_PU = ⟨exp(−ε_v μ)⟩"])


def js(o):
    if isinstance(o, dict):
        return {str(k): js(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [js(v) for v in o]
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, np.integer):
        return int(o)
    return o


json.dump(js(res), open(os.path.join(OUT, "veto.json"), "w"), indent=1)
open(os.path.join(OUT, "tables.md"), "w").write("\n".join(md) + "\n")
print("done", time.time() - t0)
