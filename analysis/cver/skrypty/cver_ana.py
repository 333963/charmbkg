#!/usr/bin/env python3
"""CVER — analiza lokalna: N_gen wg dwoch definicji, ulamek ekskluzywny per bin pT, C poprawne, czynniki korekcyjne.
Wejscia: in/part_grp_*.npz (cver_truth.py, 78 plikow mb_nd), XSEC out/xsec_results.json (tylko odczyt), BKGCAT in/bkgcat.npz (licznik WP70 per tryb).
Uzycie: .venv_h5/bin/python cver_ana.py  → out/cver.json, out/tabele.md"""
import glob, json, math, os, collections
import numpy as np
H = os.path.dirname(os.path.abspath(__file__)); B = os.path.dirname(H); CW = os.path.dirname(B)
BINS = [(3, 5), (5, 8), (8, 30)]; BN = ["3-5", "5-8", "8-30"]
PDG = {"D0": (0.03936, 0.00030, 0.000148), "Dplus": (0.0931, 0.0014, 0.0)}   # CF, blad, DCS (D0)
SP = {"D0": (421, 10, 101, (1.70, 2.05)), "Dplus": (411, 4, 201, (1.70, 2.20))}  # pdg, decayType, kod wlasny CF, zakres mas histogramu WPSCAN
parts = sorted(glob.glob(os.path.join(B, "in", "part_grp_*.npz")))
Z = [np.load(p) for p in parts]
off, T, R = 0, collections.defaultdict(list), collections.defaultdict(list)
for z in Z:
    for k in z.files:
        if k.startswith("r_"):
            R[k].append(z[k] + (off if k == "r_trow" else 0))
        elif k not in ("nev", "nfiles", "bad_idx", "all_pdg", "other_dau"):
            T[k].append(z[k])
    off += len(z["pdg"])
T = {k: np.concatenate(v) for k, v in T.items()}; R = {k: np.concatenate(v) for k, v in R.items()}
nev = int(sum(z["nev"][0] for z in Z)); nf = int(sum(z["nfiles"][0] for z in Z))
oth = collections.Counter()
for z in Z:
    for a, b in z["other_dau"]: oth[int(a)] += int(b)
X = json.load(open(os.path.join(CW, "XSEC_2026-11-19", "out", "xsec_results.json")))
out = {"n_events": nev, "n_files": nf, "other_daughter_pdg": dict(oth), "bad_idx": [int(sum(z["bad_idx"][i] for z in Z)) for i in (0, 1)], "sp": {}}
L = [f"# CVER — tabele (mb_nd, {nf} plikow, {nev} zdarzen)", "", f"Corki spoza {{K±, π±, π⁰}} na listach TruthDmesons (D⁰, D⁺, pT ≥ 3, |η| < 2,5): {dict(oth) or 'brak'}", ""]
apdg, aeta = np.abs(T["pdg"]), np.abs(T["eta"])
bz = np.load(os.path.join(CW, "BKGCAT_2026-11-20", "in", "bkgcat.npz"), allow_pickle=True)
bkeys, btab = bz["keys"], bz["tab"]
def binom(k, n): return math.sqrt(max(k * (n - k), 0) / n ** 3) if n else 0.0
for isp, (sp, (pdg, dty, cf, mr)) in enumerate(SP.items()):
    s = apdg == pdg
    isdt = T["dty"] == dty
    # definicja wlasna, niezalezna od decayType: lista corek = dokladnie K-pi+ (CF) [pi+]; D0 takze DCS K+pi- (kod 1010)
    topo_cf = T["code"] == cf
    topo_dcs = (T["code"] == 1010) if sp == "D0" else np.zeros(len(s), bool)
    topo = topo_cf | topo_dcs
    dm = T["dmvis"]
    ex0 = topo & (np.abs(dm) < 1.0)                 # bez fotonu (|m_D - m_wid| < 1 MeV)
    exF = topo & (np.abs(dm) >= 1.0) & (np.abs(T["mm2"]) < 100.0 ** 2) & (T["emiss"] > 0)   # brak zgodny z fotonami (|m_brak| < 100 MeV)
    mmiss = np.sign(T["mm2"]) * np.sqrt(np.abs(T["mm2"]))
    k0 = topo & ~ex0 & ~exF & (np.abs(mmiss - 497.611) < 3); eta_ = topo & ~ex0 & ~exF & (np.abs(mmiss - 547.862) < 3)
    rest = topo & ~ex0 & ~exF & ~k0 & ~eta_
    o = {"bins": {}, "check": {}}
    fid0 = s & (aeta < 2.2)
    o["check"] = {"decayType_rows": int((fid0 & isdt).sum()), "decayType_and_topo": int((fid0 & isdt & topo).sum()), "topo_not_decayType": int((fid0 & topo & ~isdt).sum()),
                  "decayType_not_topo": int((fid0 & isdt & ~topo).sum()), "nD_eq_nlist_in_decayType": float((T["nD"] == T["nlist"])[fid0 & isdt].mean()),
                  "nD_values_decayType": {int(a): int(b) for a, b in zip(*np.unique(T["nD"][fid0 & isdt], return_counts=True))},
                  "nD_values_excl": {int(a): int(b) for a, b in zip(*np.unique(T["nD"][fid0 & isdt & (ex0 | exF)], return_counts=True))},
                  "nD_values_nonexcl": {int(a): int(b) for a, b in zip(*np.unique(T["nD"][fid0 & isdt & ~(ex0 | exF)], return_counts=True))}}
    # histogram m_D - m_wid dla decayType sygnalu
    hh, he = np.histogram(dm[fid0 & isdt], bins=[-1e9, -1, -0.1, 0.1, 1, 5, 30, 100, 300, 480, 515, 530, 565, 800, 1e9])
    o["check"]["dmvis_hist_MeV"] = {f"{he[i]:g}..{he[i+1]:g}": int(hh[i]) for i in range(len(hh))}
    L += [f"## {sp}", "", f"Kontrola definicji (pT ≥ 3, |η| < 2,2): wiersze decayType {dty}: {o['check']['decayType_rows']}; z lista corek dokladnie K∓π±{'π±' if sp == 'Dplus' else ''}: {o['check']['decayType_and_topo']}; "
          f"ta lista corek poza decayType {dty}: {o['check']['topo_not_decayType']}; nDaughters: wszystkie {o['check']['nD_values_decayType']}, ekskluzywne {o['check']['nD_values_excl']}, nieekskluzywne {o['check']['nD_values_nonexcl']}",
          f"Histogram m_D − m_wid [MeV] (decayType {dty}): {o['check']['dmvis_hist_MeV']}", "",
          "| bin pT | N_gen wszystkie D | N_gen decayType (XSEC: N_gen_sig) | XSEC zapisane | ekskl. bez fotonu | ekskl. + foton | K⁰ | η | inne | f_ekskl (z FSR) [%] | f_ekskl bez FSR [%] | w tym DCS | BR_MC ekskl. [%] (z FSR) | BR_MC bez FSR [%] | PDG [%] | MC/PDG |", "|" + "---|" * 16]
    for b, bn in zip(BINS + [(3, 1e9)], BN + ["≥3"]):
        inb = fid0 & (T["pt"] >= b[0]) & (T["pt"] < b[1])
        n_all, n_dt = int(inb.sum()), int((inb & isdt).sum())
        c = {k: int((inb & isdt & m).sum()) for k, m in (("ex0", ex0), ("exF", exF), ("k0", k0), ("eta", eta_), ("rest", rest))}
        n_ex = c["ex0"] + c["exF"]; n_dcs = int((inb & isdt & topo_dcs & (ex0 | exF)).sum())
        f, fe = n_ex / n_dt, binom(n_ex, n_dt)
        f0 = c["ex0"] / n_dt
        br, bre = n_ex / n_all, math.sqrt(n_ex) / n_all
        pd_ = PDG[sp][0] + PDG[sp][2]
        xs = X["mc"][sp]["bins"].get(bn, {})
        o["bins"][bn] = {"N_gen_all": n_all, "N_gen_decayType": n_dt, "N_gen_sig_XSEC": xs.get("N_gen_sig"), "N_gen_all_XSEC": xs.get("N_gen_all"), **c, "N_gen_excl": n_ex, "N_dcs": n_dcs,
                         "f_excl": f, "f_excl_err": fe, "f_excl_noFSR": f0, "BR_mc_excl": br, "BR_mc_excl_err": bre, "BR_mc_noFSR": c["ex0"] / n_all, "BR_mc_decayType": n_dt / n_all, "BR_pdg": pd_, "ratio_pdg": br / pd_}
        L.append(f"| {bn} | {n_all} | {n_dt} | {xs.get('N_gen_sig', '—')} / {xs.get('N_gen_all', '—')} | {c['ex0']} | {c['exF']} | {c['k0']} | {c['eta']} | {c['rest']} | {100 * f:.2f} ± {100 * fe:.2f} | {100 * f0:.2f} | {n_dcs} | "
                 f"{100 * br:.3f} ± {100 * bre:.3f} | {100 * c['ex0'] / n_all:.3f} | {100 * pd_:.3f} | {br / pd_:.3f} ± {bre / pd_:.3f} |")
    # reco przypisane (recoMatchIndices), preselekcja: decayType reco = kanal, Lxy >= 0; kategoria reco i masa wg klasy prawdy
    tr = R["r_trow"]; ok = s[tr] & (aeta[tr] < 2.2) & isdt[tr] & (R["r_lxy"] >= 0)
    cls = np.where((ex0 | exF)[tr], 0, 1)
    inr = (R["r_m"] >= mr[0]) & (R["r_m"] < mr[1])
    L += ["", f"Kandydaci reco przypisani do prawdziwych D z decayType {dty} (recoMatchIndices; preselekcja ntupli, Lxy ≥ 0, bez ciecia score):", "",
          "| klasa prawdy | bin pT reco | wszystkich | w zakresie mas fitu | z tego kategoria reco 0/1 (Sig) | kategoria 6/7 (DMM) | inne kategorie | mediana masy [GeV] | kwantyle 5–95 % [GeV] |", "|" + "---|" * 9]
    o["reco"] = {}
    for ic, cn in enumerate(("ekskluzywne", "nieekskluzywne")):
        for b, bn in zip(BINS, BN):
            m = ok & (cls == ic) & (R["r_pt"] >= b[0]) & (R["r_pt"] < b[1])
            mi = m & inr
            cat = R["r_cat"][mi]
            q = np.quantile(R["r_m"][m], [0.05, 0.5, 0.95]) if m.any() else [float("nan")] * 3
            d = {"n": int(m.sum()), "n_in_range": int(mi.sum()), "sig": int(np.isin(cat, (0, 1)).sum()), "dmm": int(np.isin(cat, (6, 7)).sum()),
                 "cats": {int(a): int(b_) for a, b_ in zip(*np.unique(R["r_cat"][m], return_counts=True))}, "q": [float(x) for x in q],
                 "cat_mismatch_truthside": int((R["r_cat"][m] != R["r_cat2"][m]).sum())}
            o["reco"][f"{cn}|{bn}"] = d
            L.append(f"| {cn} | {bn} | {d['n']} | {d['n_in_range']} | {d['sig']} | {d['dmm']} | {d['n_in_range'] - d['sig'] - d['dmm']} | {q[1]:.3f} | {q[0]:.3f}–{q[2]:.3f} |")
    L += ["", "Kategorie reco (wszystkie masy) per klasa: " + "; ".join(f"{k}: {v['cats']}" for k, v in o["reco"].items()), ""]
    # licznik WP70 per tryb z BKGCAT (s >= indeks WP70), porownanie z XSEC N_reco_truth
    L += ["| bin | N_reco Sig WP70 (BKGCAT, suma) | XSEC N_reco_truth | w tym tryb ekskl. (+X0) | z brakujaca masa (K⁰, η, X) | C_XSEC [%] | C_popr = N_reco/N_gen_ekskl [%] | C_popr/C_XSEC | blad MC stosunku | C' = N_reco/N_gen_wszystkie [10⁻³] | C_popr × BR_PDG [10⁻³] | C' / (C_popr × BR_PDG) | σ_popr/σ_XSEC (def. 1) | σ_popr/σ_XSEC (def. 2: C′) |", "|" + "---|" * 14]
    for jb, bn in enumerate(BN):
        cell = X["cells"][f"{sp}|427929|{bn}"]
        nreco = cell["N_reco_truth_mb"]
        t = btab[(btab[:, 0] == isp) & (btab[:, 1] == jb)]
        thr = bz[f"thr|{sp}|{bn}"]
        # indeks WP70: najmniejsze s0 takie, ze suma Sig (bez Refl) po s >= s0 odtwarza N_reco_truth XSEC
        bymode = collections.defaultdict(lambda: np.zeros(32))
        for k in np.unique(t[:, 2]):
            nm = str(bkeys[k])
            if nm.startswith("Sig"):
                tt = t[t[:, 2] == k]; bymode[nm.split("|", 1)[1]] += np.bincount(tt[:, 3], weights=tt[:, 5], minlength=32)
        tot = sum(bymode.values()); cum = tot[::-1].cumsum()[::-1]
        s0 = int(np.argmin(np.abs(cum - nreco)))
        n_b = int(cum[s0]); n_ex_r = int(sum(v[s0:].sum() for m_, v in bymode.items() if not (m_.endswith("+K0") or m_.endswith("+eta") or m_.endswith("+X") or m_ == "nolink")))
        g = o["bins"][bn]
        C_x = cell["C"]; C_p = nreco / g["N_gen_excl"]; Cprime = nreco / g["N_gen_all"]
        ratio = C_p / C_x; rerr = ratio * binom(g["N_gen_excl"], g["N_gen_decayType"]) / g["f_excl"]
        k1 = C_x / C_p                                             # sigma_popr / sigma_XSEC
        k2 = (C_x * cell["BR"]) / Cprime
        o["bins"][bn].update({"N_reco_truth_XSEC": nreco, "N_reco_sig_BKGCAT_wp70": n_b, "s0": s0, "N_reco_excl": n_ex_r, "C_XSEC": C_x, "C_corr": C_p, "C_corr_err_rel": math.sqrt((1 - C_p) / nreco),
                              "Cprime": Cprime, "ratio_C": ratio, "ratio_C_err": rerr, "corr_sigma": k1, "corr_sigma_err": k1 * binom(g["N_gen_excl"], g["N_gen_decayType"]) / g["f_excl"], "corr_sigma_def2": k2,
                              "corr_sigma_def2_err": k2 * math.sqrt(1 / g["N_gen_excl"] - 1 / g["N_gen_all"]) if False else k2 * math.sqrt(max(1 / g["N_gen_decayType"] - 1 / g["N_gen_all"], 0))})
        L.append(f"| {bn} | {n_b} | {nreco} | {n_ex_r} | {n_b - n_ex_r} | {100 * C_x:.2f} | {100 * C_p:.2f} | {ratio:.4f} | {rerr:.4f} | {1e3 * Cprime:.3f} | {1e3 * C_p * cell['BR']:.3f} | {Cprime / (C_p * cell['BR']):.3f} | {k1:.4f} ± {o['bins'][bn]['corr_sigma_err']:.4f} | {k2:.4f} ± {o['bins'][bn]['corr_sigma_def2_err']:.4f} |")
    L.append("")
    out["sp"][sp] = o
# przekroje poprawione (427929) i stosunek D+/D0
lum = X["lumi"]["427929"][0] if "427929" in X["lumi"] else None
L += ["## Przekroje 427929 (wartosci XSEC z xsec_results.json × czynnik; nic nie przeliczano w XSEC)", "", f"L = {lum} nb⁻¹", "",
      "| gatunek | bin | σ_XSEC [µb] | czynnik (def. 1) | σ_popr [µb] | czynnik (def. 2, C′) | σ_popr′ [µb] | Pythia A3 inel. [µb] | dane/Pythia XSEC | dane/Pythia popr. |", "|" + "---|" * 10]
out["xsec"] = {}
for sp in SP:
    for bn in BN:
        c = X["cells"][f"{sp}|427929|{bn}"]; g = out["sp"][sp]["bins"][bn]
        sx = c["sigma_nb"] / 1e3 if not isinstance(c["sigma_nb"], dict) else None
        py = X["mc"][sp]["bins"][bn]["sigma_pythia_inel_ub"]
        out["xsec"][f"{sp}|{bn}"] = {"sigma_XSEC_ub": sx, "k1": g["corr_sigma"], "sigma_corr_ub": sx * g["corr_sigma"], "k2": g["corr_sigma_def2"], "sigma_corr2_ub": sx * g["corr_sigma_def2"], "pythia_ub": py}
        L.append(f"| {sp} | {bn} | {sx:.1f} | {g['corr_sigma']:.4f} ± {g['corr_sigma_err']:.4f} | {sx * g['corr_sigma']:.1f} | {g['corr_sigma_def2']:.4f} | {sx * g['corr_sigma_def2']:.1f} | {py:.1f} | {sx / py:.2f} | {sx * g['corr_sigma'] / py:.2f} |")
L += ["", "| bin | R = D⁺/D⁰ XSEC (427929, SF⁺ = SF⁰) | czynnik k⁺/k⁰ (def. 1) | R popr. | czynnik (def. 2) | Pythia A3 |", "|---|---|---|---|---|---|"]
out["ratio"] = {}
for bn in BN:
    r = X["ratio_DpD0"][f"427929|{bn}"] if f"427929|{bn}" in X["ratio_DpD0"] else None
    gp, g0 = out["sp"]["Dplus"]["bins"][bn], out["sp"]["D0"]["bins"][bn]
    k = gp["corr_sigma"] / g0["corr_sigma"]; ke = k * math.hypot(gp["corr_sigma_err"] / gp["corr_sigma"], g0["corr_sigma_err"] / g0["corr_sigma"])
    k2 = gp["corr_sigma_def2"] / g0["corr_sigma_def2"]
    out["ratio"][bn] = {"R_XSEC": r["R_sf_equal"] if r else None, "k": k, "k_err": ke, "R_corr": r["R_sf_equal"] * k if r else None, "k_def2": k2, "pythia": r["pythia_A3_mbnd"] if r else None}
    L.append(f"| {bn} | {r['R_sf_equal']:.3f} | {k:.4f} ± {ke:.4f} | {r['R_sf_equal'] * k:.3f} | {k2:.4f} | {r['pythia_A3_mbnd']:.3f} |" if r else f"| {bn} | — | {k:.4f} ± {ke:.4f} | — | {k2:.4f} | — |")
os.makedirs(os.path.join(B, "out"), exist_ok=True)
json.dump(out, open(os.path.join(B, "out", "cver.json"), "w"), indent=1, ensure_ascii=False)
open(os.path.join(B, "out", "tabele.md"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
