#!/usr/bin/env python3
"""TAGXSEC (20.11.2026, projekt B) — C_D(SD)/C_D(ND), σ_tag^fid(D + proton w AFP) i R^fid = σ_tag/σ_incl, run 427929. Wykonuje ../KRYTERIA.md bez zmian.
Wejście (tylko odczyt): runs/lxplus/tx_mb_{sd,nd}.npz (tx_mc.py), TAGGED/out/{unblind,prep}.json, TAGGED/runs/lxplus/tg_cand_427929.npz (tylko krotność nSelBL
kandydatów — żadnych fitów), XSEC/out/xsec_results.json, SD_var/out/results.json. Wyjście: out/tagxsec.json, tabele/wyniki.md.
Użycie: ~/Documents/charm/.venv_h5/bin/python tx_ana.py"""
import json, os, sys
import numpy as np

DEF = sys.argv[1] if len(sys.argv) > 1 else "excl"     # "excl": mianownik C po prawdziwym trybie (decayType sygnału ∧ klasa braku 0/1, definicja BKGCAT); "xsec": po samym decayType (jak XSEC_2026-11-19)
assert DEF in ("excl", "xsec")

HERE = os.path.dirname(os.path.abspath(__file__)); TOP = os.path.dirname(HERE); W = os.path.dirname(TOP)
BINS = [(3, 5), (5, 8), (8, 30)]; BN = ["3-5", "5-8", "8-30"]
SP = {"D0": {"isp": 0, "pdg": 421, "win": (1.86484, 0.055), "rng": (1.70, 2.05), "BR": (0.03951, 0.00030)},
      "Dplus": {"isp": 1, "pdg": 411, "win": (1.86966, 0.045), "rng": (1.70, 2.20), "BR": (0.0931, 0.0014)}}
WP = 6                                                      # s ≥ 6 ⇔ WP70
L, dL = 0.4648, 0.0103                                      # nb^-1
A, dA = 0.97, 0.03; EAFP, dEAFP = 0.95, 0.02; ALIGN = 0.034
ETRIG, dETRIG = 1.000, 0.001; EXICAL, dEXICAL = 1.000, 0.005; ENPU, dENPU = 0.985, 0.007
rng = np.random.default_rng(20261120)
NBOOT = 1000


def ld(p):
    z = np.load(p, allow_pickle=False); return {k: z[k] for k in z.files}


SD, ND = ld(f"{TOP}/runs/lxplus/tx_mb_sd.npz"), ld(f"{TOP}/runs/lxplus/tx_mb_nd.npz")


def add_excl(T, path):
    """Dołącza g_excl (prawdziwy tryb K∓π±[π±]: klasa braku 0/1) do wierszy gen po (plik, wpis, pdg, pT)."""
    E = ld(path); n = len(T["g_pt"]); T["g_excl"] = np.zeros(n, bool); T["g_mc"] = np.full(n, -1, np.int8)
    ka = np.lexsort((np.round(T["g_pt"].astype(np.float64), 3), T["g_pdg"], T["g_entry"], T["g_fid"]))
    m = T["g_mode"][ka] == 1; ka = ka[m]
    kb = np.lexsort((np.round(E["pt"].astype(np.float64), 3), E["pdg"], E["entry"], E["fid"]))
    assert len(ka) == len(kb), (len(ka), len(kb))
    for c, d in (("g_fid", "fid"), ("g_entry", "entry"), ("g_pdg", "pdg")):
        assert np.array_equal(T[c][ka].astype(np.int64), E[d][kb].astype(np.int64)), c
    assert np.allclose(T["g_pt"][ka], E["pt"][kb], rtol=1e-4)
    T["g_mc"][ka] = E["mc"][kb]; T["g_excl"][ka] = E["mc"][kb] <= 1


add_excl(SD, f"{TOP}/runs/lxplus/ex_mb_sd.npz"); add_excl(ND, f"{TOP}/runs/lxplus/ex_mb_nd.npz")
UNB = json.load(open(f"{W}/TAGGED_2026-11-19/out/unblind.json"))["items"]
PREP = json.load(open(f"{W}/TAGGED_2026-11-19/out/prep.json"))["items"]
XS = json.load(open(f"{W}/XSEC_2026-11-19/out/xsec_results.json"))["cells"]
SDV = json.load(open(f"{W}/SD_var_2026-11-19/out/results.json"))
TG = ld(f"{W}/TAGGED_2026-11-19/runs/lxplus/tg_cand_427929.npz")
OUT = {"time": "2026-11-20", "kryteria": "KRYTERIA.md (20.11 00:10, przed liczeniem)"}
MD = []


def ibin(pt):
    return np.digitize(pt, [3, 5, 8, 30]) - 1


def gmask(T, sp, mode=True, prot=None):
    m = (np.abs(T["g_pdg"]) == SP[sp]["pdg"]) & (np.abs(T["g_eta"]) < 2.2) & (T["g_pt"] >= 3) & (T["g_pt"] < 30) & T["g_okf"]
    if mode: m &= (T["g_mode"] == 1) & (T["g_excl"] if DEF == "excl" else True)
    if prot is True: m &= (T["g_pA"] > 0) | (T["g_pC"] > 0)
    return m


def rmask(T, sp, prot=None, wp=True):
    m = (T["r_sp"] == SP[sp]["isp"]) & T["r_okf"]
    if wp: m &= T["r_s"] >= WP
    if prot is True: m &= (T["r_pA"] > 0) | (T["r_pC"] > 0)
    return m


def eff(k, n):
    k, n = np.asarray(k, float), np.asarray(n, float)
    e = k / np.maximum(n, 1); return e, np.sqrt(np.maximum(e * (1 - e), 0.25 / np.maximum(n, 1) ** 2 * 0 + e * (1 - e)) / np.maximum(n, 1))


def counts(T, sp, prot=None):
    r, g = rmask(T, sp, prot), gmask(T, sp, True, prot)
    k = np.bincount(ibin(T["r_pt"][r]), minlength=3)[:3]; n = np.bincount(ibin(T["g_pt"][g]), minlength=3)[:3]
    return k, n


RES = {}
for sp in SP:
    R = RES[sp] = {}
    # ---------------- C per bin, oszacowanie 1
    kN, nN = counts(ND, sp); kS, nS = counts(SD, sp); kP, nP = counts(SD, sp, True)
    cN, eN = eff(kN, nN); cS, eS = eff(kS, nS); cP, eP = eff(kP, nP)
    xs = [XS[f"{sp}|427929|{b}"] for b in BN]
    DEF0 = DEF; DEF = "xsec"; kX, nX = counts(ND, sp); kSx, nSx = counts(SD, sp); DEF = DEF0
    assert all(int(x["N_reco_truth_mb"]) == int(k) and int(x["N_gen_sig_mb"]) == int(n) for x, k, n in zip(xs, kX, nX)), "licznik/mianownik mb_nd ≠ XSEC"
    fex = nN / nX                                            # udział prawdziwego trybu w mianowniku XSEC per bin (1 dla DEF = xsec)
    R["f_excl_nd"] = fex; R["f_excl_nd_tot"] = float(nN.sum() / nX.sum()); R["f_excl_sd"] = nS / np.maximum(nSx, 1); R["f_excl_sd_tot"] = float(nS.sum() / nSx.sum())
    # widmo g_b: prawdziwe D (wszystkie mody) w mb_sd z protonem w oknie; wariant: wszystkie SD; wariant: widmo inkluzywne zmierzone
    ga = gmask(SD, sp, False, True); g_tag = np.bincount(ibin(SD["g_pt"][ga]), minlength=3)[:3].astype(float); n_gtag = g_tag.sum(); g_tag /= n_gtag
    gs = gmask(SD, sp, False); g_sd = np.bincount(ibin(SD["g_pt"][gs]), minlength=3)[:3].astype(float); n_gsd = g_sd.sum(); g_sd /= n_gsd
    sig_incl = np.array([x["sigma_nb"] for x in xs]) / 1e3 * fex   # µb; przy DEF = excl przeliczone na C z mianownikiem po prawdziwym trybie
    g_incl = sig_incl / sig_incl.sum()
    SF = np.array([x["SF"] for x in xs]); dSF = np.array([x["rel"]["SF"] for x in xs]) * SF

    def tot(c, g):
        return float((g * c).sum())

    def rho1(k, n, g):
        """łączne ρ₁ = Σ g C_SD / Σ g C_ND i błąd dwumianowy (SD ⊕ ND)."""
        c, e = eff(k, n); v = tot(c, g) / tot(cN, g)
        err = v * np.hypot(np.sqrt(((g * e) ** 2).sum()) / tot(c, g), np.sqrt(((g * eN) ** 2).sum()) / tot(cN, g))
        return v, err
    r1_all, r1_all_e = rho1(kS, nS, g_tag); r1_p, r1_p_e = rho1(kP, nP, g_tag)
    pooled = (kS.sum() / nS.sum()) / ((kN / nN * nS).sum() / nS.sum())          # ND przeważone do widma mianownika SD
    pooled_e = pooled / np.sqrt(kS.sum())
    R.update(kN=kN, nN=nN, kS=kS, nS=nS, kP=kP, nP=nP, cN=cN, eN=eN, cS=cS, eS=eS, cP=cP, eP=eP, g_tag=g_tag, n_gtag=n_gtag, g_sd=g_sd, n_gsd=n_gsd, g_incl=g_incl,
             rho1_bin=cS / cN, rho1_bin_err=(cS / cN) * np.hypot(eS / np.maximum(cS, 1e-9), eN / cN), rho1_all=(r1_all, r1_all_e), rho1_prot=(r1_p, r1_p_e),
             rho1_pooled=(pooled, pooled_e), SF=SF, dSF=dSF, sig_incl=sig_incl)
    # ---------------- oszacowanie 2: C_b(n) w mb_nd, przedziały nSelBL o ≥ 100 zrekonstruowanych (zachłannie od małych n)
    def edges_for(var, b):
        r = rmask(ND, sp) & (ibin(ND["r_pt"]) == b); x = np.sort(ND["r_" + var][r]); ed = [0]; i = 0
        while len(x) - i >= 200:
            j = i + 100; v = x[j - 1]
            while j < len(x) and x[j] == v: j += 1
            if len(x) - j < 100: break
            ed.append(int(v) + 1); i = j
        return np.array(ed + [10 ** 6])

    def cn_tab(var, b, boot=False):
        ed = edges_for(var, b)
        r = rmask(ND, sp) & (ibin(ND["r_pt"]) == b); g = gmask(ND, sp) & (ibin(ND["g_pt"]) == b)
        k = np.histogram(ND["r_" + var][r], ed)[0].astype(float); n = np.histogram(ND["g_" + var][g], ed)[0].astype(float)
        return ed, k, n

    def ceff(ed, k, n, nvals, w=None):
        c = k / n; ix = np.clip(np.searchsorted(ed, nvals, side="right") - 1, 0, len(c) - 1)
        w = np.ones(len(nvals)) if w is None else w
        return w.sum() / (w / c[ix]).sum()
    TAB = {(var, b): cn_tab(var, b) for var in ("nsel", "ntrk") for b in range(3)}
    R["cn"] = {f"{var}|{BN[b]}": {"edges": TAB[(var, b)][0][:-1].tolist(), "k": TAB[(var, b)][1].tolist(), "n": TAB[(var, b)][2].tolist()} for (var, b) in TAB}
    # dane: kandydaci T2 ∧ MBTS ∧ WP70 w zakresie mas (jak TAGGED tg_lib.load_cand / region_mask)
    a, c = TG["tag_A"], TG["tag_C"]
    T2 = (TG["npu"] == 0) & (~a | (TG["xical_A"] <= TG["xi_A"] + 0.02)) & (~c | (TG["xical_C"] <= TG["xi_C"] + 0.02))
    mb = (TG["trig"] & (1 << 1)) > 0
    lo, hi = SP[sp]["rng"]; m0, ww = SP[sp]["win"]
    base = TG["ok"] & (TG["ipt"] >= 0) & TG["grl"] & T2 & mb & (TG["sp"] == SP[sp]["isp"]) & (TG["s"] >= WP) & (TG["m"] >= lo) & (TG["m"] < hi)
    inw = np.abs(TG["m"] - m0) < ww
    key = f"{sp}|427929|MBTS|3-30|0.7|T2"
    assert int(base.sum()) == int(UNB[key]["n_tot"]) and abs(int((base & inw).sum()) - int(UNB[key]["n_win"])) <= 2, (base.sum(), (base & inw).sum())
    bsub = PREP[key]["bkg_in_window_extrap"] / PREP[key]["nsb_tag"]
    R["data_n"] = {"n_tot": int(base.sum()), "n_win": int((base & inw).sum()), "per_bin_win": np.bincount(TG["ipt"][base & inw], minlength=3).tolist(),
                   "per_bin_sb": np.bincount(TG["ipt"][base & ~inw], minlength=3).tolist(), "med_nsel_win": float(np.median(TG["nsel"][base & inw])),
                   "med_nsel_sb": float(np.median(TG["nsel"][base & ~inw])), "per_side_win": [int((base & inw & a & ~c).sum()), int((base & inw & c & ~a).sum()), int((base & inw & a & c).sum())],
                   "per_side_sb": [int((base & ~inw & a & ~c).sum()), int((base & ~inw & c & ~a).sum()), int((base & ~inw & a & c).sum())]}
    for nm, T in (("mb_nd", ND), ("mb_sd", SD)):
        r = rmask(T, sp); R["data_n"][f"med_nsel_reco_{nm}"] = float(np.median(T["r_nsel"][r])) if r.any() else None

    def rho2_sets(sets, tabs, g):
        """sets[b] = (nvals, w) per bin pT; → ρ₂ per bin i łącznie (wagi g · C_ND)."""
        rb = np.ones(3)
        for b in range(3):
            nv, w = sets[b]
            if len(nv) == 0 or w.sum() <= 0: continue
            ed, k, n = tabs[b]; rb[b] = ceff(ed, k, n, nv, w) / (k.sum() / n.sum())
        return rb, float((g * cN * rb).sum() / (g * cN).sum())

    def dsets(sel, wfun=None):
        return [(TG["nsel"][sel & (TG["ipt"] == b)].astype(float), (np.ones(int((sel & (TG["ipt"] == b)).sum())) if wfun is None else wfun(sel & (TG["ipt"] == b)))) for b in range(3)]
    tabs_nsel = [TAB[("nsel", b)] for b in range(3)]; tabs_ntrk = [TAB[("ntrk", b)] for b in range(3)]
    S_win = dsets(base & inw); S_sb = dsets(base & ~inw)
    r2w_b, r2w = rho2_sets(S_win, tabs_nsel, g_tag); r2s_b, r2s = rho2_sets(S_sb, tabs_nsel, g_tag)
    # okno − pasma: histogram w przedziałach krotności, wagi −bsub dla pasm
    r2d_b = np.ones(3); neg = False
    for b in range(3):
        ed, k, n = tabs_nsel[b]; cc = k / n
        hw = np.histogram(S_win[b][0], ed)[0].astype(float); hs = np.histogram(S_sb[b][0], ed)[0].astype(float); h = hw - bsub * hs
        if (h < 0).any(): neg = True
        if h.sum() > 0 and (h / cc).sum() > 0: r2d_b[b] = h.sum() / (h / cc).sum() / (k.sum() / n.sum())
    r2d = float((g_tag * cN * r2d_b).sum() / (g_tag * cN).sum())
    # bootstrap: dane (Poisson po kandydatach okna) × MC (k ~ dwumian)
    bs = []
    for _ in range(NBOOT):
        tb = []
        for b in range(3):
            ed, k, n = tabs_nsel[b]; kk = rng.binomial(n.astype(int), k / n).astype(float); kk[kk == 0] = 0.5; tb.append((ed, kk, n))
        ss = [(nv, rng.poisson(1.0, len(nv)).astype(float)) for nv, _ in S_win]
        bs.append(rho2_sets(ss, tb, g_tag)[1])
    r2w_stat = float(np.std(bs, ddof=1))
    # domknięcie na MC SD: rozkład krotności zrekonstruowanego sygnału mb_sd
    cl = {}
    for var, tabs in (("nsel", tabs_nsel), ("ntrk", tabs_ntrk)):
        rs = rmask(SD, sp)
        ss = [(SD["r_" + var][rs & (ibin(SD["r_pt"]) == b)].astype(float), np.ones(int((rs & (ibin(SD["r_pt"]) == b)).sum()))) for b in range(3)]
        rb, rt = rho2_sets(ss, tabs, g_tag); cl[var] = (rb, rt)
    nonclos = r1_all / cl["nsel"][1] - 1; nonclos_e = r1_all_e / cl["nsel"][1]
    s_var = max(abs(r2s - r2w), abs(r2d - r2w)) if not neg else abs(r2s - r2w)
    s_clos = max(abs(nonclos), nonclos_e) * r2w
    rho, drho = r2w, float(np.sqrt(r2w_stat ** 2 + s_var ** 2 + s_clos ** 2))
    tension = abs(r1_all - r2w) / np.hypot(r1_all_e, np.hypot(r2w_stat, s_var))
    R.update(rho2_win=(r2w_b, r2w), rho2_sb=(r2s_b, r2s), rho2_diff=(r2d_b, r2d), rho2_diff_neg=neg, rho2_stat=r2w_stat, closure=cl, nonclosure=(nonclos, nonclos_e),
             rho=rho, drho=drho, drho_parts={"stat": r2w_stat, "warianty": s_var, "niedomkniecie": s_clos}, tension_sigma=float(tension), rho_b=r2w_b)
    # ε_wierzch: udział prawdziwych D (fid., z protonem) w mb_sd ze zdarzeń o ntrk < 5
    gp = gmask(SD, sp, False, True)
    R["f_ntrk_lt5"] = float((SD["g_ntrk"][gp] < 5).mean()); R["f_ntrk_le2"] = float((SD["g_ntrk"][gp] <= 2).mean())
    R["med_ntrk_gen_sd_prot"] = float(np.median(SD["g_ntrk"][gp])); R["med_ntrk_gen_nd"] = float(np.median(ND["g_ntrk"][gmask(ND, sp)]))
    # widmo: przeliczenia dla literatury (pT ≥ 3,5 → ≥ 3; |η| < 2,5 → < 2,2), z prawdy mb_sd (wszystkie mody, wszystkie zdarzenia SD)
    e = (np.abs(SD["g_pdg"]) == SP[sp]["pdg"]) & SD["g_okf"]
    R["lit_conv"] = {"pt3_over_pt35_eta25": float((e & (SD["g_pt"] >= 3) & (np.abs(SD["g_eta"]) < 2.5)).sum() / max((e & (SD["g_pt"] >= 3.5) & (np.abs(SD["g_eta"]) < 2.5)).sum(), 1)),
                     "eta22_over_eta25_pt3": float((e & (SD["g_pt"] >= 3) & (np.abs(SD["g_eta"]) < 2.2)).sum() / max((e & (SD["g_pt"] >= 3) & (np.abs(SD["g_eta"]) < 2.5)).sum(), 1)),
                     "f_win_sd_fid": float(n_gtag / n_gsd)}

# =================================================================================================== złożenie
def assemble(sp, N, st_lo, st_hi, shape, label):
    R = RES[sp]; key = f"{sp}|427929|MBTS|3-30|0.7|T2"
    BR, dBR = SP[sp]["BR"]; Nacc, dNacc = UNB[key]["N_acc"], UNB[key]["N_acc_err"]
    g, cN, SF, rb = R["g_tag"], R["cN"], R["SF"], R["rho_b"]
    # wagi per strona (tylko D0: fity K7 istnieją)
    eA, eC = 0.968, 0.965; em = 0.5 * (eA + eC); edA, edC = 1 - 0.012, 1 - 0.003; edm = 1 - 0.0075
    AE = {"A": A * edA / edm * EAFP * eA / em, "C": A * edC / edm * EAFP * eC / em}
    if sp == "D0":
        nA, nC = UNB["D0|427929|MBTS_A|3-30|0.7|T2"]["fit"]["N"], UNB["D0|427929|MBTS_C|3-30|0.7|T2"]["fit"]["N"]
        inv = (nA / AE["A"] + nC / AE["C"]) / (nA + nC); AEeff = 1 / inv; side = {"N_A": nA, "N_C": nC, "AE_A": AE["A"], "AE_C": AE["C"]}
    else:
        AEeff = A * EAFP; side = {"uwaga": "brak fitów per strona dla D+ → średnia"}
    Cbar = float((g * cN * SF * rb).sum())
    den = Cbar * BR * L * AEeff * ETRIG * EXICAL * ENPU                       # nb^-1
    s = (N - Nacc) / den / 1e3                                                # µb
    contrib = g * cN * SF * rb / Cbar
    wincl = np.array([x["N"] for x in PREP[key]["incl"]]); wincl = wincl / wincl.sum()
    xs = [XS[f"{sp}|427929|{b}"] for b in BN]
    # warianty g_b
    Cbar_incl = float((R["g_incl"] * cN * SF * rb).sum()); Cbar_sd = float((R["g_sd"] * cN * SF * rb).sum())
    bud = {
        "kształt fitu tagowanego (TAGGED)": shape / (N - Nacc) if shape else None,
        "N_acc": dNacc / (N - Nacc),
        "obciążenie fitu (domknięcie inkluzywne, liniowo)": float((wincl * np.array([x["rel"]["fit_bias"] for x in xs])).sum()),
        "SF klasyfikatora (3–5 GeV: scenariusz)": float(np.sqrt(((contrib * R["dSF"] / SF) ** 2).sum())),
        "ρ = C_D(tag)/C_D(ND)": R["drho"] / R["rho"],
        "widmo pT (g_b: MC SD z protonem wobec inkluzywnego)": abs(Cbar_incl / Cbar - 1),
        "stat. MC C_b (mb_nd)": float(np.sqrt(((contrib * R["eN"] / cN) ** 2).sum())),
        "nieprompt (z XSEC)": float((contrib * np.array([x["rel"]["nonprompt"] for x in xs])).sum()),
        "A (geometria)": dA / A, "ε_AFP": dEAFP / EAFP, "ustawienie / optyka": ALIGN,
        "ε_trig": dETRIG, "ε_ξcal (jednostronnie +)": dEXICAL, "ε_npu": dENPU / ENPU, "ε_wierzch (udział ntrk < 5)": R["f_ntrk_lt5"],
        "BR": dBR / BR, "L (2,2 %)": dL / L}
    keys_noshape = [k for k in bud if bud[k] is not None and not k.startswith("kształt")]
    syst_noshape = float(np.sqrt(sum(bud[k] ** 2 for k in keys_noshape)))
    syst = float(np.hypot(syst_noshape, bud["kształt fitu tagowanego (TAGGED)"] or 0))
    # ---- stosunek do inkluzywnego
    s_incl = float(R["sig_incl"].sum()); Rf = s / s_incl
    # reszta SF i C w stosunku: wariacja SF_b wspólna dla licznika i mianownika
    def ratio_with(SFv, cv):
        cb = float((g * cv * SFv * rb).sum()); si = float(sum(x["N"] / (cv[i] * SFv[i]) for i, x in enumerate(xs)))
        return (1 / cb) / si
    r0 = ratio_with(SF, cN); dsf = []; dc = []
    for i in range(3):
        v = SF.copy(); v[i] += R["dSF"][i]; dsf.append(ratio_with(v, cN) / r0 - 1)
        v = cN.copy(); v[i] += R["eN"][i]; dc.append(ratio_with(SF, v) / r0 - 1)
    inc_stat = float(np.sqrt(sum((R["sig_incl"][i] * xs[i]["rel"]["stat"]) ** 2 for i in range(3))) / s_incl)
    inc_model = float(np.sqrt(sum((R["sig_incl"][i] * xs[i]["rel"]["fit_model"]) ** 2 for i in range(3))) / s_incl)
    budR = {
        "kształt fitu tagowanego (TAGGED)": bud["kształt fitu tagowanego (TAGGED)"], "N_acc": bud["N_acc"],
        "stat. N_incl": inc_stat, "model tła N_incl": inc_model,
        "SF (reszta po skróceniu)": float(np.sqrt(sum(x ** 2 for x in dsf))), "stat. MC C_b (reszta)": float(np.sqrt(sum(x ** 2 for x in dc))),
        "ρ = C_D(tag)/C_D(ND)": bud["ρ = C_D(tag)/C_D(ND)"], "widmo pT (g_b)": bud["widmo pT (g_b: MC SD z protonem wobec inkluzywnego)"],
        "A (geometria)": dA / A, "ε_AFP": dEAFP / EAFP, "ustawienie / optyka": ALIGN, "ε_trig (stosunek)": 0.005, "ε_ξcal (jednostronnie +)": dEXICAL, "ε_npu": dENPU / ENPU,
        "ε_wierzch (udział ntrk < 5)": R["f_ntrk_lt5"]}
    systR_noshape = float(np.sqrt(sum(v ** 2 for k, v in budR.items() if v is not None and not k.startswith("kształt"))))
    systR = float(np.hypot(systR_noshape, budR["kształt fitu tagowanego (TAGGED)"] or 0))
    return {"label": label, "N": N, "N_acc": Nacc, "stat_lo": st_lo, "stat_hi": st_hi, "shape": shape, "Cbar": Cbar, "Cbar_g_incl": Cbar_incl, "Cbar_g_sd": Cbar_sd,
            "contrib": contrib.tolist(), "AEeff": AEeff, "side": side, "den_nb": den, "sigma_ub": s, "sig_stat_lo": s * st_lo / (N - Nacc), "sig_stat_hi": s * st_hi / (N - Nacc),
            "syst_rel": syst, "syst_rel_noshape": syst_noshape, "budget": bud, "sigma_incl_ub": s_incl, "R": Rf, "R_stat_lo": Rf * st_lo / (N - Nacc), "R_stat_hi": Rf * st_hi / (N - Nacc),
            "R_syst_rel": systR, "R_syst_rel_noshape": systR_noshape, "budget_R": budR}


FIN = {}
for sp in SP:
    key = f"{sp}|427929|MBTS|3-30|0.7|T2"; u, p = UNB[key], PREP[key]
    nom = assemble(sp, u["fit"]["N"], u["stat_lo"], u["stat_hi"], u["syst"]["razem"], "nominał (fit TAGGED)")
    nv = (u["n_win"] - p["bkg_in_window_extrap"]) / u["sig_win_frac"]
    ev = np.sqrt(u["n_win"] + p["bkg_in_window_extrap"] ** 2 / p["nsb_tag"]) / u["sig_win_frac"]
    var = assemble(sp, nv, ev, ev, None, "wariant okno − pasma")
    FIN[sp] = {"nom": nom, "var": var}

# ---- przewidywania Pythii (poziom cząstek, 3–30 GeV z R per bin × σ_incl Pythii A3 per bin z XSEC T5: 1230 / 481 / 160 µb dla D0; D+: 641 / 254 / 85)
PY_INCL = {"D0": np.array([1230., 481., 160.]), "Dplus": np.array([641., 254., 85.])}
PRED = {}
for sp, kk in (("D0", "D0"), ("Dplus", "Dp")):
    w = PY_INCL[sp] / PY_INCL[sp].sum(); un = SDV["uncertainty"][kk]
    r_soft = float((np.array(SDV["scenarios"]["A3off"][kk]["R_pt"][1:]) * w).sum()); r_hard = float((np.array(SDV["scenarios"]["hard_t4"][kk]["R_pt"][1:]) * w).sum())
    allr = {s: float((np.array(v[kk]["R_pt"][1:]) * w).sum()) for s, v in SDV["scenarios"].items() if kk in v and v[kk].get("R_pt")}
    PRED[sp] = {"R_soft": r_soft, "R_soft_lo": r_soft * un["factor_dn"], "R_soft_hi": r_soft * un["factor_up"], "R_hard": r_hard,
                "R_hard_lo": r_hard * un["hard"]["factor_dn"], "R_hard_hi": r_hard * un["hard"]["factor_up"], "per_scen": allr, "py_incl_ub": float(PY_INCL[sp].sum())}


def verdict(lo68, hi68, lo95, hi95, plo, phi):
    if hi68 >= plo and lo68 <= phi: return "zgodny"
    if hi95 >= plo and lo95 <= phi: return "napięcie"
    return "niezgodny"


for sp in SP:
    for v in FIN[sp].values():
        Rf = v["R"]; sl = np.hypot(v["R_stat_lo"], Rf * v["R_syst_rel"]); sh = np.hypot(v["R_stat_hi"], Rf * v["R_syst_rel"])
        v["R_68"] = [Rf - sl, Rf + sh]; v["R_95"] = [Rf - 1.96 * sl, Rf + 1.96 * sh]
        P = PRED[sp]
        v["verdict_soft"] = verdict(*v["R_68"], *v["R_95"], P["R_soft_lo"], P["R_soft_hi"]); v["verdict_hard"] = verdict(*v["R_68"], *v["R_95"], P["R_hard_lo"], P["R_hard_hi"])
        v["pull_soft"] = (Rf - P["R_soft"]) / sl; v["pull_hard"] = (Rf - P["R_hard"]) / sl


def cv(o):
    if isinstance(o, dict): return {str(k): cv(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)): return [cv(x) for x in o]
    if isinstance(o, np.ndarray): return o.tolist()
    if isinstance(o, (np.floating, np.integer, np.bool_)): return o.item()
    return o


OUT.update(mc=cv(RES), wynik=cv(FIN), pred=cv(PRED), stale={"L": [L, dL], "A": [A, dA], "eAFP": [EAFP, dEAFP], "align": ALIGN, "eps_npu": [ENPU, dENPU]})
os.makedirs(f"{TOP}/out", exist_ok=True); os.makedirs(f"{TOP}/tabele", exist_ok=True)
OUT["DEF"] = DEF
json.dump(OUT, open(f"{TOP}/out/tagxsec_{DEF}.json", "w"), indent=1, ensure_ascii=False)

# =================================================================================================== tabele
f3 = lambda x: f"{x:.3f}".replace(".", ",")  # noqa: E731
f2 = lambda x: f"{x:.2f}".replace(".", ",")  # noqa: E731
f1 = lambda x: f"{x:.1f}".replace(".", ",")  # noqa: E731
pc = lambda x: f"{100 * x:.1f}".replace(".", ",")  # noqa: E731
MD.append("# TAGXSEC — tabele (generowane przez skrypty/tx_ana.py; nie edytować ręcznie)\n")
MD.append("**Definicja mianownika C: " + ("prawdziwy tryb K∓π±[π±] (decayType sygnału ∧ klasa brakującej masy 0/1, jak BKGCAT)" if DEF == "excl" else "sam decayType 10 / 4 (jak XSEC_2026-11-19)") + f"** (DEF = {DEF}).\n")
for sp in SP:
    R = RES[sp]
    MD.append(f"\n## {sp}\n\n### X1. C_D per bin pT (WP70, drzewo H): mb_nd, mb_sd (wszystkie), mb_sd z prawdziwym protonem 0,03 ≤ ξ ≤ 0,10\n")
    MD.append("| bin pT | mb_nd k/n | C_ND | mb_sd k/n | C_SD | ρ₁ = C_SD/C_ND | mb_sd+p k/n | C_SD+p | g_b (SD+p) | g_b (SD) | g_b (inkl. zmierzone) | SF |")
    MD.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for b in range(3):
        MD.append(f"| {BN[b]} | {R['kN'][b]}/{R['nN'][b]} | {f3(R['cN'][b])} ± {f3(R['eN'][b])} | {R['kS'][b]}/{R['nS'][b]} | {f3(R['cS'][b])} ± {f3(R['eS'][b])} | {f2(R['rho1_bin'][b])} ± {f2(R['rho1_bin_err'][b])} | "
                  f"{R['kP'][b]}/{R['nP'][b]} | {f3(R['cP'][b])} ± {f3(R['eP'][b])} | {f3(R['g_tag'][b])} | {f3(R['g_sd'][b])} | {f3(R['g_incl'][b])} | {f3(R['SF'][b])} ± {f3(R['dSF'][b])} |")
    MD.append(f"\nUdział prawdziwego trybu w mianowniku po decayType: mb_nd {' / '.join(f3(x) for x in R['f_excl_nd'])} (łącznie {f3(R['f_excl_nd_tot'])}); mb_sd {' / '.join(f3(x) for x in R['f_excl_sd'])} (łącznie {f3(R['f_excl_sd_tot'])}).")
    MD.append(f"\nPrawdziwe {sp} (wszystkie mody, obszar fid.) w mb_sd: {int(R['n_gsd'])}, z protonem w oknie: {int(R['n_gtag'])} (f_win = {pc(R['lit_conv']['f_win_sd_fid'])} %).")
    MD.append(f"\n### X2. ρ łącznie (3–30 GeV, wagi g_b z MC SD z protonem)\n\n| oszacowanie | ρ | uwagi |\n|---|---|---|")
    MD.append(f"| 1a: mb_sd wprost, wszystkie zdarzenia SD | {f3(R['rho1_all'][0])} ± {f3(R['rho1_all'][1])} | {int(R['kS'].sum())} zrekonstruowanych |")
    MD.append(f"| 1b: mb_sd wprost, proton w oknie | {f3(R['rho1_prot'][0])} ± {f3(R['rho1_prot'][1])} | {int(R['kP'].sum())} zrekonstruowanych |")
    MD.append(f"| 1c: mb_sd wprost, zliczenie łączne (ND przeważone do widma mianownika SD) | {f3(R['rho1_pooled'][0])} ± {f3(R['rho1_pooled'][1])} | |")
    MD.append(f"| 2: krotność nSelBL, dane — okno sygnału (NOMINAŁ) | {f3(R['rho2_win'][1])} ± {f3(R['rho2_stat'])} (bootstrap) | per bin: {' / '.join(f3(x) for x in R['rho2_win'][0])}; kandydatów w oknie per bin {R['data_n']['per_bin_win']} |")
    MD.append(f"| 2: dane — pasma boczne | {f3(R['rho2_sb'][1])} | per bin: {' / '.join(f3(x) for x in R['rho2_sb'][0])}; kandydatów {R['data_n']['per_bin_sb']} |")
    MD.append(f"| 2: dane — okno − pasma | {f3(R['rho2_diff'][1])} | per bin: {' / '.join(f3(x) for x in R['rho2_diff'][0])}; ujemne przedziały: {'tak (tylko opisowo)' if R['rho2_diff_neg'] else 'nie'} |")
    MD.append(f"| domknięcie: metoda 2 na krotności sygnału reco mb_sd (nSelBL) | {f3(R['closure']['nsel'][1])} | per bin: {' / '.join(f3(x) for x in R['closure']['nsel'][0])} |")
    MD.append(f"| domknięcie: to samo, nTracksPt500Eta25 | {f3(R['closure']['ntrk'][1])} | per bin: {' / '.join(f3(x) for x in R['closure']['ntrk'][0])} |")
    MD.append(f"| niedomknięcie ρ₁(1a)/ρ₂(SD) − 1 | {pc(R['nonclosure'][0])} % ± {pc(R['nonclosure'][1])} % | |")
    MD.append(f"| **przyjęte ρ** | **{f3(R['rho'])} ± {f3(R['drho'])}** | stat {f3(R['drho_parts']['stat'])}, warianty {f3(R['drho_parts']['warianty'])}, niedomknięcie {f3(R['drho_parts']['niedomkniecie'])}; ρ₁(1a) wobec ρ₂(dane): {f1(R['tension_sigma'])}σ |")
    d = R["data_n"]
    MD.append(f"\nMediana nSelBL: dane tagowane okno {f1(d['med_nsel_win'])}, pasma {f1(d['med_nsel_sb'])}; sygnał reco mb_nd {f1(d['med_nsel_reco_mb_nd'])}, mb_sd {f1(d['med_nsel_reco_mb_sd'])}. "
              f"Strony kandydatów (tylko A / tylko C / obie): okno {d['per_side_win']}, pasma {d['per_side_sb']}.")
    MD.append(f"Prawdziwe D w mb_sd z protonem: mediana nTracksPt500Eta25 = {f1(R['med_ntrk_gen_sd_prot'])} (mb_nd: {f1(R['med_ntrk_gen_nd'])}); udział zdarzeń z < 5 torami {pc(R['f_ntrk_lt5'])} %, ≤ 2 torami {pc(R['f_ntrk_le2'])} %.")
    MD.append(f"\n### X3. C_b(n) w mb_nd (nSelBL; przedziały o ≥ 100 zrekonstruowanych)\n")
    for b in range(3):
        t = R["cn"][f"nsel|{BN[b]}"]; c = np.array(t["k"]) / np.array(t["n"])
        MD.append(f"- {BN[b]} GeV: " + "; ".join(f"n ≥ {e}: {f3(x)}" for e, x in zip(t["edges"], c)))
    for tag in ("nom", "var"):
        v = FIN[sp][tag]
        MD.append(f"\n### X4-{tag}. {sp}: {v['label']}{' — OPISOWO (kontrolki D⁺ niespełnione)' if sp == 'Dplus' else ''}\n")
        MD.append(f"N = {f1(v['N'])} −{f1(v['stat_lo'])} +{f1(v['stat_hi'])}{'' if v['shape'] is None else ' ± ' + f1(v['shape']) + ' (kształt)'}; N_acc = {f2(v['N_acc'])}; C̄_tag = {f3(v['Cbar'])} "
                  f"(z g_b inkluzywnym {f3(v['Cbar_g_incl'])}; udziały binów w C̄: {' / '.join(pc(x) for x in v['contrib'])} %); ⟨A·ε_AFP⟩ = {f3(v['AEeff'])}; mianownik = {f3(v['den_nb'] * 1e3)} µb⁻¹.")
        MD.append(f"\n**σ_tag^fid = {f1(v['sigma_ub'])} −{f1(v['sig_stat_lo'])} +{f1(v['sig_stat_hi'])} (stat) ± {f1(v['sigma_ub'] * v['syst_rel'])} (syst) µb** "
                  f"(syst {pc(v['syst_rel'])} %; bez kształtu fitu tagowanego {pc(v['syst_rel_noshape'])} %); L 7 % zamiast 2,2 %: syst {pc(np.sqrt(v['syst_rel'] ** 2 - (dL / L) ** 2 + 0.07 ** 2))} %.")
        MD.append(f"\n**R^fid = σ_tag^fid/σ_incl(3–30) = {f2(100 * v['R'])} −{f2(100 * v['R_stat_lo'])} +{f2(100 * v['R_stat_hi'])} (stat) ± {f2(100 * v['R'] * v['R_syst_rel'])} (syst) %** "
                  f"(σ_incl = {f1(v['sigma_incl_ub'])} µb; syst {pc(v['R_syst_rel'])} %, bez kształtu {pc(v['R_syst_rel_noshape'])} %); przedział 68 % [{f2(100 * v['R_68'][0])}; {f2(100 * v['R_68'][1])}] %, 95 % [{f2(100 * v['R_95'][0])}; {f2(100 * v['R_95'][1])}] %.")
        MD.append("\n| pozycja budżetu | σ_tag^fid [%] | R^fid [%] |\n|---|---|---|")
        allk = list(dict.fromkeys(list(v["budget"].keys()) + list(v["budget_R"].keys())))
        for k in allk:
            a_, b_ = v["budget"].get(k), v["budget_R"].get(k)
            MD.append(f"| {k} | {'–' if a_ is None else pc(a_)} | {'–' if b_ is None else pc(b_)} |")
        P = PRED[sp]
        MD.append(f"\nPythia (poziom cząstek, 3–30 GeV): miękka A3 R = {f2(100 * P['R_soft'])} % ({f2(100 * P['R_soft_lo'])}–{f2(100 * P['R_soft_hi'])}) → {v['verdict_soft']} ({v['pull_soft']:+.1f}σ od centralnej); "
                  f"twarda R = {f2(100 * P['R_hard'])} % ({f2(100 * P['R_hard_lo'])}–{f2(100 * P['R_hard_hi'])}) → {v['verdict_hard']} ({v['pull_hard']:+.1f}σ). "
                  f"σ_pred = R × σ_incl zmierzone: {f1(P['R_soft'] * v['sigma_incl_ub'])} ({f1(P['R_soft_lo'] * v['sigma_incl_ub'])}–{f1(P['R_soft_hi'] * v['sigma_incl_ub'])}) µb; "
                  f"× σ_incl Pythii A3 ({f1(P['py_incl_ub'])} µb): {f1(P['R_soft'] * P['py_incl_ub'])} µb; twarda: {f1(P['R_hard'] * v['sigma_incl_ub'])} ({f1(P['R_hard_lo'] * v['sigma_incl_ub'])}–{f1(P['R_hard_hi'] * v['sigma_incl_ub'])}) µb.")
    lc = R["lit_conv"]
    MD.append(f"\nPrzeliczenia do literatury z prawdy mb_sd ({sp}): N(pT ≥ 3)/N(pT ≥ 3,5) przy |η| < 2,5 = {f3(lc['pt3_over_pt35_eta25'])}; N(|η| < 2,2)/N(|η| < 2,5) przy pT ≥ 3 = {f3(lc['eta22_over_eta25_pt3'])}.")
open(f"{TOP}/tabele/wyniki_{DEF}.md", "w").write("\n".join(MD) + "\n")
print("\n".join(MD))
