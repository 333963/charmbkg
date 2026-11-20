"""EST_spread — zestawienie estymatorów tła pod pikiem (ε_B / N_B w SR) z plików wynikowych projektu B + nowe liczby (est_new.py).
Każdy wiersz ma ścieżkę źródła. Nic nie jest zgadywane: wartości czytane z json; brak błędu w źródle → puste pole.
Wyjście: out/tabela_estymatory.csv (długa), out/rozrzut.csv (rozrzut per komórka), out/tables.md (fragmenty do raportu), out/collect.json
"""
import csv, glob, json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
W = os.path.abspath(os.path.join(HERE, ".."))
OUT = os.path.join(HERE, "out")
rows = []
BINS = ["3-5", "5-8", "8-30"]


def J(p):
    return json.load(open(os.path.join(W, p)))


def add(group, sp, pt, wp, est, cls, val, err, qty, sample, mode, src, new=0, note=""):
    rows.append(dict(group=group, species=sp, pt=pt, wp=int(round(float(wp) * 100)) if float(wp) <= 1 else int(wp), estimator=est, cls=cls,
                     value=float(val), err=("" if err is None or not np.isfinite(err) else float(err)), quantity=qty, sample=sample,
                     mode=mode, source=src, new=new, note=note))


# ---------------------------------------------------------------- G1–G3: NOWE (bt_*.npz)
N = json.load(open(os.path.join(OUT, "est_new.json")))
WPSN = N["meta"]["wps"]
SRC_N = "EST_spread/out/est_new.json (z B_trig/runs/lxplus/bt_*.npz; wagi: B_trig/BKG_427929/out/s3_w_*.npz)"
SAMP = {"d427": "run 427929 (MBTS_1 & GRL), v2.0, 30 % zdarzeń, klasyfikator H; tylko D0 inkluzywne",
        "d428": "run 428770 (82 % L1_J12), v2.0, 30 % zdarzeń, klasyfikator H; tylko D0 inkluzywne"}
CLS_SB = {"pol2": "SB: ML exp(pol2)", "pol1": "SB: ML exp(pol1)", "pol3": "SB: ML exp(pol3)", "sym": "SB: zliczenie 75–125 MeV",
          "lr": "SB: interpolacja liniowa L/P", "far": "SB: exp(pol2), tylko |dm|>=90 MeV"}
for run in ("d427", "d428"):
    g = f"G1_D0_dane_{run[1:]}" if run == "d427" else f"G2_D0_dane_{run[1:]}"
    for b in BINS:
        for e, d in N["data"][f"{run}|{b}"].items():
            for j, wp in enumerate(WPSN):
                add(g, "D0", b, wp, e, CLS_SB[e], d["val"][j], d["err"][j], "eps_B(SR)", SAMP[run], "dane", SRC_N, 1)
        for e, d in N["mc_est"][b].items():
            if e.startswith("mcV") and not e.endswith(run):
                continue
            nm = e.replace("_" + run, "")
            cl = {"mcraw_t2": "MC surowe truth>=2", "mcraw_t8": "MC surowe truth>=8", "mcV1_t2": "MC przeważone bez score (V1)",
                  "mcV2_t2": "MC przeważone ze score (V2, typu SALAD)"}[nm]
            for j, wp in enumerate(WPSN):
                add(g, "D0", b, wp, nm, cl, d["val"][j], d["err"][j], "eps_B(SR)", SAMP[run] + "; MC mb_nd v2.0 30 %", "dane",
                    SRC_N, 1, "błąd = stat MC (bez wariancji treningu wag)")
for tn in ("t2", "t8"):
    for b in BINS:
        d = N["mc_closure"][f"{tn}|{b}"]
        for e in N["meta"]["sb_est"]:
            for j, wp in enumerate(WPSN):
                add(f"G3_D0_MC_bt_{tn}", "D0", b, wp, e, CLS_SB[e], d[e]["ratio"][j], d[e]["ratio_err"][j], "eps_B est/prawda",
                    f"MC mb_nd v2.0 30 %, tło truth>={tn[1:]}, klasyfikator H", "domknięcie MC", SRC_N, 1,
                    f"prawda eps={d['true']['val'][j]:.5f}")

# ---------------------------------------------------------------- G4: D+ dane podpróbka S8
for kind in ("hgb", "sk_mlp"):
    p = f"S8_salad/work/out/data_{kind}.json"
    d = J(p)
    smp = "run 428770, h5 v1.1, podpróbka data_p0 (parzystość 1), D+ po wecie D*; klasyfikator etapu S8"
    for wp in ("0.8", "0.7", "0.5"):
        f = d[f"effB_data_fit@{wp}"]
        for b in BINS + ["all"]:
            if kind == "hgb":
                add("G4_Dp_dane_S8sub", "Dplus", b, wp, "pol2_fit", "SB: fit pol2 danych przed/po cięciu", f[b], f["err" if b == "all" else b + "_err"],
                    "eps_B(SR)", smp, "dane", p)
                for k, nm, cl in (("effB_mcraw_SR", "mcraw_t8", "MC surowe truth>=8"), ("effB_dataSB", "sb_asym", "SB: zliczenie, całe (asymetryczne)"),
                                  ("effB_dataSBsym", "sb_sym", "SB: zliczenie symetryczne 70–170 MeV"),
                                  ("effB_dataSB_left", "sb_left", "SB: tylko lewe"), ("effB_dataSB_right", "sb_right", "SB: tylko prawe")):
                    add("G4_Dp_dane_S8sub", "Dplus", b, wp, nm, cl, d[f"{k}@{wp}"][b], None, "eps_B(SR)", smp, "dane", p)
            add("G4_Dp_dane_S8sub", "Dplus", b, wp, f"salad_{kind}", f"SALAD ({kind})", d[f"effB_salad_SR@{wp}"][b], None, "eps_B(SR)", smp, "dane", p)

# ---------------------------------------------------------------- G5/G6: pełna lumi S8 (lxplus)
for sp, g in (("Dplus", "G5_Dp_dane_S8full"), ("D0", "G6_D0_dane_S8full")):
    p = f"S8_salad/runs/lxplus/salad_s8/bkg_shape_salad_{sp}.json"
    d = J(p)
    smp = f"run 428770, h5 v1.1, pełna lumi (n_data={d['n_data']}), nominalny H; SALAD {d['kind']}, MC truth>=8"
    for b in BINS + ["all"]:
        for wp in ("80", "70", "50"):
            for k, nm, cl in (("effB_salad_SR", "salad_sk_mlp", "SALAD (sk_mlp)"), ("effB_mcraw_SR", "mcraw_t8", "MC surowe truth>=8"),
                              ("effB_data_SB", "sb_asym", "SB: zliczenie, całe (asymetryczne)")):
                add(g, sp, b, wp, nm, cl, d[k][b][wp], None, "eps_B(SR)", smp, "dane", p)

# ---------------------------------------------------------------- G7: D0 dane B1
for fn, tmin in (("data_eff_D0.json", 8), ("data_eff_D0_mc2.json", 2)):
    p = f"B1_salad_v1/out_D0/{fn}"
    d = J(p)
    smp = "run 428770, h5 v1.1, podpróbka data_p1 (n=%d), D0 po wecie refleksji 25 MeV, pT 3–30; klasyfikator etapu B1" % d["n_data_p1"]
    for wp in ("0.8", "0.7", "0.5"):
        x = d[f"wp{wp}"]
        if tmin == 8:
            add("G7_D0_dane_B1", "D0", "all", wp, "pol2_fit", "SB: fit pol2 danych przed/po cięciu", x["eff_data_SR"], x["err"], "eps_B(SR)", smp, "dane", p)
            add("G7_D0_dane_B1", "D0", "all", wp, "sb_copy", "SB: zliczenie (kopia SB)", x["eff_sbcopy"], None, "eps_B(SR)", smp, "dane", p)
        add("G7_D0_dane_B1", "D0", "all", wp, f"mcraw_t{tmin}", f"MC surowe truth>={tmin}", x["eff_mcraw_SR"], None, "eps_B(SR)", smp, "dane", p)
        for k in ("hgb", "sk_mlp"):
            add("G7_D0_dane_B1", "D0", "all", wp, f"salad_{k}_t{tmin}", f"SALAD ({k}), MC truth>={tmin}", x[f"eff_salad_{k}_SR"], None, "eps_B(SR)", smp, "dane", p)

# ---------------------------------------------------------------- G8/G9: S4 (FETA) — dane i domknięcie MC
p = "S4_feta/runs/lxplus/full/results.json"
d = J(p)
smp = "run 428770, h5 v1.1, n_data=%d (bieg lxplus full), D0 po wecie refleksji, pT 3–30; XGB etapu S4, fit na połowie B" % d["P"]["n_data"]
for wp in ("80", "70", "50"):
    x = d["E_dane"]["Rwin"][f"WP{wp}"]
    for k, nm, cl in (("epsB_SR_z_fitu", "pol2_fit", "SB: fit Gauss+exp(pol2) danych"), ("epsB_SR_FETA", "feta", "generator: FETA"),
                      ("epsB_SR_SALAD_mlp", "salad_mlp", "SALAD (mlp)"), ("epsB_SR_SB_kopia", "sb_copy", "SB: zliczenie (kopia SB)")):
        add("G8_D0_dane_S4", "D0", "all", wp, nm, cl, x[k], None, "eps_B(SR)", smp, "dane", p)
    tv = d["B_szablon"]["prawda_epsB_SR"][f"WP{wp}"]
    for k, nm, cl in (("FETA", "feta", "generator: FETA"), ("MC_surowe", "mcraw_t8", "MC surowe truth>=8"), ("SALAD_mlp", "salad_mlp", "SALAD (mlp)"),
                      ("SALAD_hgb", "salad_hgb", "SALAD (hgb)"), ("SALAD_hgb_nom", "salad_hgb_nom", "SALAD (hgb bez m)"),
                      ("SB_kopia(B')", "sb_copy", "SB: zliczenie (kopia SB)"), ("flow_S2(CATHODE)", "flow", "generator: flow p(x|m)"),
                      ("idealne(A')", "ideal", "kontrola: druga kopia prawdy (skala fluktuacji)")):
        add("G9_D0_MC_S4", "D0", "all", wp, nm, cl, d["B_szablon"][k][f"WP{wp}"], None, "eps_B est/prawda",
            "MC mb_nd v1.1 (lxplus full), pseudo-dane = MC połowa B ze znanym gładkim zniekształceniem; n_pseudo=%d" % d["P"]["n_pseudo"],
            "domknięcie MC", p, 0, f"prawda eps={tv:.5f}")

# ---------------------------------------------------------------- G10/G11: S6 candrot (szablony kształtu masy, N_pred/N_true) + ε ratio
for sp in ("D0", "Dplus"):
    p = f"S6_kand/work/out/s6_{sp}.json"
    d = J(p)
    g = "G10_D0_MC_S6" if sp == "D0" else "G11_Dp_MC_S6"
    smp = f"MC mb_nd v1.1 podpróbki, połowa B, prawda SR n={d['n']['true_sr']}, pT 3–30; klasyfikator etapu S6"
    for wp in ("0.8", "0.7", "0.5"):
        add(g, sp, "all", wp, "pol2", "SB: fit pol2 (χ²) prawdy", d["pol2_pred_ratio"][wp], None, "N_B(SR) pred/prawda", smp, "domknięcie MC", p)
        for k, c in d["controls"].items():
            nm, srcc, vm = k.split("|") if "|" in k else (k, "", "")
            if srcc and srcc != "sb":
                continue                        # tylko źródło „sb” (realizowalne na danych); all/bkg pominięte
            tag = nm + ("_resample" if vm == "resample" else "")
            if "pred_ratio_wp" in c and vm == "keep":
                add(g, sp, "all", wp, tag, {"rot_pi": "rotacja toru o π", "rot_rand": "rotacja toru losowa", "flip_eta": "odbicie η toru",
                                             "swap": "zamiana K↔π (D+: pseudo-WS)", "mix": "tor z innego zdarzenia"}[nm] + " — szablon masy, norm. w SB",
                    c["pred_ratio_wp"][wp], c["pred_ratio_wp"][wp] * c["pred_err_wp"][wp], "N_B(SR) pred/prawda", smp, "domknięcie MC", p)
            add(g + "_eps", sp, "all", wp, tag, "ε_B obiektu kontrolnego w SR / prawda", c["eps_ratio"][wp], None, "eps_B est/prawda", smp, "domknięcie MC", p)

# ---------------------------------------------------------------- G12: VAL D+ MC (i dane)
p = "VAL_validate/runs/smoke_vm/validation/Dplus.json"
d = J(p)["variants"]["H"]
for scope, dd in [("all", d["sculpt_mcbkg"]["all"])] + [(b, d["sculpt_mcbkg"]["bins"][b]) for b in BINS]:
    for wp in ("80", "70", "50"):
        x = dd[wp]
        t = x["eps"]["SR"]
        smp = "MC mb_nd v1.1 (bieg smoke_vm), tło bez sygnału, n_SR=%d; nominalny H" % t["n"]
        add("G12_Dp_MC_VAL", "Dplus", scope, wp, "pol2", "SB: fit pol2", x["pol2_closure"], x["pol2_closure_err"], "eps_B est/prawda", smp, "domknięcie MC", p,
            0, f"prawda eps={t['eps']:.5f}")
        for k, nm, cl in (("SB_asym", "sb_asym", "SB: zliczenie, całe (asymetryczne)"), ("SB_sym", "sb_sym", "SB: zliczenie symetryczne"),
                          ("inner_sym", "inner", "pasmo wewnętrzne (zliczenie)")):
            e = x["eps"][k]
            r = e["eps"] / t["eps"]
            er = r * np.sqrt(1 / max(e["n_pass"], 1) + 1 / max(t["n_pass"], 1))
            add("G12_Dp_MC_VAL", "Dplus", scope, wp, nm, cl, r, er, "eps_B est/prawda", smp, "domknięcie MC", p, 0,
                "stosunek i błąd dwumianowy policzone tu z n_pass w źródle" + ("; pasmo wewn. zachodzi na SR — nie jest niezależne" if k == "inner_sym" else ""))

# ---------------------------------------------------------------- G13/G14: domknięcia SALAD (S8, B1) — D+ i D0
def closure_rows(g, sp, p, smp, pt_bins=True):
    d = J(p)
    for wp in ("0.8", "0.7", "0.5"):
        tr = d[f"eff_true_SR@{wp}"]
        rel = d.get(f"eff_true_relerr@{wp}")
        for b in (BINS + ["all"]) if pt_bins else ["all"]:
            if tr[b] <= 0:
                continue
            for k, nm, cl in (("eff_salad_SR", f"salad_{d['kind']}", f"SALAD ({d['kind']})"), ("eff_mcraw_SR", "mcraw_t8", "MC surowe"),
                              ("eff_sbcopy", "sb_copy", "SB: zliczenie (kopia SB)")):
                add(g, sp, b, wp, nm, cl, d[f"{k}@{wp}"][b] / tr[b], (rel if b == "all" else None), "eps_B est/prawda", smp, "domknięcie MC", p, 0,
                    "błąd = σ_stat prawdy (tylko 'all')")


for kind in ("hgb", "sk_mlp"):
    closure_rows("G13_Dp_MC_S8smooth", "Dplus", f"S8_salad/work/out/closure_{kind}.json",
                 "MC mb_nd v1.1 podpróbki S6, pseudo-dane = MC B gładko zniekształcone w m; klasyfikator etapu S8")
for sp, dirn in (("Dplus", "out"), ("D0", "out_D0")):
    for f in sorted(glob.glob(os.path.join(W, "B1_salad_v1", dirn, "closure_*.json"))):
        rel = os.path.relpath(f, W)
        mode = os.path.basename(f)[len("closure_"):-len(".json")]
        closure_rows(f"G14_{'Dp' if sp == 'Dplus' else 'D0'}_MC_B1|{mode}", sp, rel,
                     f"MC mb_nd v1.1 podpróbki, pseudo-dane = MC parz. 1, tryb B1 '{mode}'; klasyfikator etapu B1", pt_bins=False)

# ---------------------------------------------------------------- G15/G16: generatory S (kopuła, flow) na MC, D0
p = "S_bkggen/work/out/s2_mc.json"
d = J(p)
for k, nm, cl in (("sb_copy", "sb_both", "SB: zliczenie, oba pasma"), ("sb_near_copy", "sb_near", "SB: zliczenie, bliskie pasmo"),
                  ("copula", "copula", "generator: kopuła Gaussa"), ("flow", "flow", "generator: flow p(x|m), 40 epok CPU")):
    for wp in ("80", "70", "50"):
        add("G15_D0_MC_S2gen", "D0", "all", wp, nm, cl, d[k][f"effB_ratio_{wp}"], None, "eps_B est/prawda",
            "MC mb_nd v1.1 podpróbki S, pT 3–30, bez weta refleksji; klasyfikator etapu S", "domknięcie MC", p, 0, f"prawda eps={d[k][f'effB_{wp}_ref']:.5f}")
p = "S_bkggen/work/out/s2b_veto.json"
d = J(p)
for k, nm, cl in (("sb_veto", "sb_both_veto", "SB: oba pasma + weto refleksji"), ("sb_far_veto", "sb_far_veto", "SB: dalekie pasmo + weto refleksji"),
                  ("sb_near_veto", "sb_near_veto", "SB: bliskie pasmo + weto"), ("copula_veto", "copula_veto", "generator: kopuła + weto")):
    for j, wp in enumerate(("80", "70", "50")):
        add("G15b_D0_MC_S2veto", "D0", "all", wp, nm, cl, d[k]["effB_ratio"][j], None, "eps_B est/prawda",
            "MC mb_nd v1.1 podpróbki S, pT 3–30, prawda po wecie refleksji; klasyfikator etapu S", "domknięcie MC", p)
p = "S2_flow/runs/cloud_v1/results.json"
d = J(p)
for wp in ("0.8", "0.7", "0.5"):
    x = d["B_szablon_masy"]["wp"][wp]
    for k, nm, cl in (("eps_sr_flow", "flow", "generator: flow p(x|m) (S2)"), ("eps_sr_copula", "copula", "generator: kopuła (S2)")):
        add("G16_D0_MC_S2flow", "D0", "all", wp, nm, cl, x[k] / x["eps_sr_true"], None, "eps_B est/prawda", "MC mb_nd v1.1 podpróbki (S2 cloud_v1), pT 3–30",
            "domknięcie MC", p, 0, f"prawda eps={x['eps_sr_true']:.5f}")

# ================================================================ zapis długiej tabeli
cols = ["group", "species", "pt", "wp", "estimator", "cls", "value", "err", "quantity", "sample", "mode", "new", "source", "note"]
with open(os.path.join(OUT, "tabela_estymatory.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols)
    w.writeheader()
    for r in rows:
        w.writerow({k: (f"{r[k]:.6g}" if isinstance(r[k], float) else r[k]) for k in cols})
print("wierszy:", len(rows), "grup:", len({r["group"] for r in rows}))

# ================================================================ rozrzut
REP = np.load(os.path.join(OUT, "est_new_rep.npz"))
SBE = N["meta"]["sb_est"]
MCN = ["mcraw_t2", "mcraw_t8", "mcV1_t2_d427", "mcV2_t2_d427", "mcV1_t2_d428", "mcV2_t2_d428"]
BUD = {"Dplus": {80: 5, 70: 10, 50: 25}, "D0": {80: 5, 70: 10, 50: 15}}


def stats(v):
    v = np.asarray(v, float)
    m = v.mean()
    return dict(n=len(v), mean=m, rms=float(np.sqrt(((v - m) ** 2).mean()) / m), std1=float(v.std(ddof=1) / m) if len(v) > 1 else np.nan,
                rng=float((v.max() - v.min()) / m), half=float((v.max() - v.min()) / 2 / m), vmin=float(v.min()), vmax=float(v.max()))


def chi2_cov(v, cov):
    """χ² hipotezy „wszystkie estymatory równe” z pełną kowariancją (kontrasty względem pierwszego)."""
    n = len(v)
    D = np.zeros((n - 1, n)); D[:, 0] = -1; D[np.arange(n - 1), np.arange(1, n)] = 1
    d = D @ v
    C = D @ cov @ D.T
    return float(d @ np.linalg.pinv(C, rcond=1e-10) @ d), n - 1


def chi2_indep(v, e):
    v, e = np.asarray(v, float), np.asarray(e, float)
    w = 1 / e ** 2
    m = (w * v).sum() / w.sum()
    return float((w * (v - m) ** 2).sum()), len(v) - 1


spread = []
SETS = {"A_interp5": ["pol2", "pol1", "pol3", "sym", "lr"],
        "B_klasy4": ["pol2", "sym", "lr", "mcV2_t2"],
        "C_klasy4+MCraw": ["pol2", "sym", "lr", "mcV2_t2", "mcraw_t2"],
        "E_trio": ["pol2", "sym", "mcV2_t2"]}
for run in ("d427", "d428"):
    for b in BINS:
        Rd = REP[f"{run}|{b}"]             # [rep, est, wp]
        Rm = REP[f"mcest|{b}"]
        for j, wp in enumerate(WPSN):
            for sn, names in SETS.items():
                v, blocks = [], []
                di = [SBE.index(x) for x in names if x in SBE]
                mi = [MCN.index(x if x.startswith("mcraw") else f"{x}_{run}") for x in names if x not in SBE]
                v = [N["data"][f"{run}|{b}"][x]["val"][j] for x in names if x in SBE] + \
                    [N["mc_est"][b][MCN[i]]["val"][j] for i in mi]
                n = len(v)
                cov = np.zeros((n, n))
                cov[:len(di), :len(di)] = np.cov(Rd[:, di, j].T)
                if mi:
                    cov[len(di):, len(di):] = np.atleast_2d(np.cov(Rm[:, mi, j].T))
                c2, ndf = chi2_cov(np.array(v), cov)
                s = stats(v)
                spread.append(dict(group=f"{'G1' if run == 'd427' else 'G2'}_D0_dane_{run[1:]}", species="D0", pt=b, wp=int(float(wp) * 100), set=sn,
                                   estimators="+".join(names), mode="dane", chi2=c2, ndf=ndf, chi2_kind="pełna kowariancja (bootstrap)",
                                   stat_nom=N["data"][f"{run}|{b}"]["pol2"]["err"][j] / N["data"][f"{run}|{b}"]["pol2"]["val"][j], **s))
# domknięcie MC bt
for tn in ("t2", "t8"):
    for b in BINS:
        Rr = REP[f"mcclos_{tn}|{b}"]
        for j, wp in enumerate(WPSN):
            names = SETS["A_interp5"]
            ii = [SBE.index(x) for x in names]
            v = np.array([N["mc_closure"][f"{tn}|{b}"][x]["ratio"][j] for x in names])
            c2, ndf = chi2_cov(v, np.cov(Rr[:, ii, j].T))
            s = stats(v)
            e = np.array([N["mc_closure"][f"{tn}|{b}"][x]["ratio_err"][j] for x in names])
            ib = int(np.argmin(np.abs(v - 1)))
            spread.append(dict(group=f"G3_D0_MC_bt_{tn}", species="D0", pt=b, wp=int(float(wp) * 100), set="A_interp5", estimators="+".join(names),
                               mode="domknięcie MC", chi2=c2, ndf=ndf, chi2_kind="pełna kowariancja (bootstrap)", best=names[ib], best_dev=float(v[ib] - 1),
                               best_err=float(e[ib]), nom_dev=float(v[0] - 1), nom_err=float(e[0]), mean_dev=float(s["mean"] - 1),
                               truth_in_range=int(s["vmin"] <= 1 <= s["vmax"]), stat_nom=float(e[0]), **s))


def generic(group, sp, names, mode, label=None, pts=None, excl_note=""):
    by = {}
    for r in rows:
        if r["group"] == group and r["estimator"] in names:
            by.setdefault((r["pt"], r["wp"]), {})[r["estimator"]] = (r["value"], r["err"])
    for (pt, wp), dd in sorted(by.items(), key=lambda x: (str(x[0][0]), -x[0][1])):
        if pts and pt not in pts:
            continue
        nm = [n for n in names if n in dd]
        if len(nm) < 3:
            continue
        v = np.array([dd[n][0] for n in nm])
        s = stats(v)
        rec = dict(group=group, species=sp, pt=pt, wp=wp, set=label or "wszystkie", estimators="+".join(nm), mode=mode, **s)
        ee = [dd[n][1] for n in nm]
        if all(x != "" for x in ee):
            rec["chi2"], rec["ndf"] = chi2_indep(v, ee)
            rec["chi2_kind"] = "błędy niezależne (przybliżenie)"
        if mode == "domknięcie MC":
            ib = int(np.argmin(np.abs(v - 1)))
            rec.update(best=nm[ib], best_dev=float(v[ib] - 1), nom_dev=float(v[0] - 1), mean_dev=float(s["mean"] - 1),
                       truth_in_range=int(s["vmin"] <= 1 <= s["vmax"]))
        spread.append(rec)


generic("G4_Dp_dane_S8sub", "Dplus", ["pol2_fit", "sb_sym", "salad_hgb", "salad_sk_mlp"], "dane", "B_klasy4 (fit, SBsym, SALAD×2)")
generic("G4_Dp_dane_S8sub", "Dplus", ["pol2_fit", "sb_sym", "salad_hgb", "salad_sk_mlp", "mcraw_t8"], "dane", "C_klasy4+MCraw")
generic("G4_Dp_dane_S8sub", "Dplus", ["pol2_fit", "sb_sym", "salad_hgb", "salad_sk_mlp", "mcraw_t8", "sb_asym"], "dane", "D_wszystkie6 (z obciążonym SB asym)")
generic("G7_D0_dane_B1", "D0", ["pol2_fit", "sb_copy", "salad_hgb_t2", "salad_sk_mlp_t2"], "dane", "B_klasy4 (fit, SB, SALAD×2 truth>=2)")
generic("G7_D0_dane_B1", "D0", ["pol2_fit", "sb_copy", "salad_hgb_t8", "salad_sk_mlp_t8"], "dane", "B_klasy4 (SALAD truth>=8)")
generic("G7_D0_dane_B1", "D0", ["pol2_fit", "sb_copy", "salad_hgb_t2", "salad_sk_mlp_t2", "mcraw_t2"], "dane", "C_klasy4+MCraw")
generic("G8_D0_dane_S4", "D0", ["pol2_fit", "sb_copy", "salad_mlp", "feta"], "dane", "fit, SB, SALAD mlp, FETA")
generic("G9_D0_MC_S4", "D0", ["sb_copy", "salad_mlp", "salad_hgb", "feta", "flow"], "domknięcie MC", "SB, SALAD×2, FETA, flow")
generic("G9_D0_MC_S4", "D0", ["sb_copy", "salad_mlp", "salad_hgb"], "domknięcie MC", "SB, SALAD×2")
generic("G10_D0_MC_S6", "D0", ["pol2", "rot_rand", "flip_eta", "mix", "swap", "rot_pi"], "domknięcie MC", "pol2 + 5 szablonów candrot")
generic("G10_D0_MC_S6", "D0", ["pol2", "rot_rand", "mix"], "domknięcie MC", "ALICE-3: pol2, rotacja, mieszanie")
generic("G11_Dp_MC_S6", "Dplus", ["pol2", "rot_rand", "flip_eta", "mix", "swap", "rot_pi"], "domknięcie MC", "pol2 + 5 szablonów candrot")
generic("G11_Dp_MC_S6", "Dplus", ["pol2", "rot_rand", "mix"], "domknięcie MC", "ALICE-3: pol2, rotacja, mieszanie")
generic("G12_Dp_MC_VAL", "Dplus", ["pol2", "sb_sym", "sb_asym"], "domknięcie MC", "pol2, SBsym, SBasym")
generic("G13_Dp_MC_S8smooth", "Dplus", ["sb_copy", "salad_hgb", "salad_sk_mlp", "mcraw_t8"], "domknięcie MC", "SB, SALAD×2, MC surowe")
generic("G15_D0_MC_S2gen", "D0", ["sb_both", "sb_near", "copula", "flow"], "domknięcie MC", "SB×2, kopuła, flow")
generic("G15b_D0_MC_S2veto", "D0", ["sb_both_veto", "sb_far_veto", "sb_near_veto", "copula_veto"], "domknięcie MC", "SB×3 + weto, kopuła")

sc = ["group", "species", "pt", "wp", "set", "mode", "n", "estimators", "mean", "rms", "std1", "rng", "half", "vmin", "vmax", "chi2", "ndf", "chi2_kind",
      "stat_nom", "best", "best_dev", "best_err", "nom_dev", "nom_err", "mean_dev", "truth_in_range"]
with open(os.path.join(OUT, "rozrzut.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=sc)
    w.writeheader()
    for r in spread:
        w.writerow({k: (f"{r[k]:.5g}" if isinstance(r.get(k), float) else r.get(k, "")) for k in sc})
json.dump({"rows": rows, "spread": spread}, open(os.path.join(OUT, "collect.json"), "w"))
print("rozrzut:", len(spread))
