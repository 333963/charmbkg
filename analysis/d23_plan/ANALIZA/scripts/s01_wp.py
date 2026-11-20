"""D23/ANALIZA krok 0+1: kontrola wejść i wybór punktu pracy na MC + pasmach bocznych danych (okno sygnału ZAŚLEPIONE).
Zapis: out/wp.json.  Użycie: python s01_wp.py"""
import json
import os
import sys

import numpy as np
import pandas as pd
from iminuit import Minuit

HERE = os.path.dirname(os.path.abspath(__file__))
TOP = os.path.dirname(HERE)
CW = os.path.dirname(os.path.dirname(TOP))
IN = os.path.join(os.path.dirname(TOP), "runs", "lxplus")
SCAN = os.path.join(CW, "DLZ_dplus", "SCAN_2026-11-17")
sys.path.insert(0, os.path.join(SCAN, "K1c", "scripts"))
sys.path.insert(0, os.path.join(CW, "DLZ_dplus", "DATA_2026-11-17", "scripts"))
import fastfit as FF  # noqa: E402
import fitlib as FL  # noqa: E402

M0 = FF.M0
RUNS = {"427929": 1, "428770": 0}       # klasa histogramu
Y35 = {"427929": 4331.0, "428770": 5917.0}   # N/eps 3-5 GeV WP70 pol3 (LOWPT_calib T1)
GRID = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
BLIND = 0.060


def sb_fit(h, order=3):
    """exp(pol) do pasm bocznych (biny 5 MeV, okno ±60 MeV wyłączone) -> B w ±45 MeV, chi2/ndf pasm."""
    mask = np.abs(FF.C - M0) > BLIND
    names = ["nb"] + [f"b{i+1}" for i in range(order)]

    def nll(*p):
        lam = np.clip(p[0] * FF.bkg_bins(p[1:]), 1e-9, None)
        return float(np.sum((lam - h * np.log(lam))[mask]))
    m = Minuit(nll, h[mask].sum() / mask.mean(), -0.5, *([0.0] * (order - 1)), name=names)
    m.errordef = 0.5
    m.limits["nb"] = (0, None)
    m.migrad(ncall=20000)
    m.hesse()
    p = [m.values[n] for n in names]
    lam = p[0] * FF.bkg_bins(p[1:])
    win = np.abs(FF.C - M0) < FF.WIN
    chi2 = float((((h - lam) ** 2 / np.maximum(lam, 1))[mask]).sum())
    return float(lam[win].sum()), chi2, int(mask.sum() - len(names)), bool(m.valid)


def main():
    out = {"inputs": {}, "mc": {}, "wp": {}, "choice": {}}
    S = json.load(open(os.path.join(IN, "summary.json")))
    # ---------- MC
    mc = pd.read_parquet(os.path.join(IN, "dp_fd_nd.parquet"))
    hm = np.load(os.path.join(IN, "hist_fd_nd.npz"))
    se, me, pe = hm["s_edges"], hm["m_edges"], hm["pt_edges"]
    print("pt_edges", pe, "cls MC", list(hm["cls"]), "dp_h", hm["dp_h"].shape)
    sig = mc[(mc.is_sig == 1) & (mc.dstar_veto == 0)]
    out["inputs"]["fd_nd"] = dict(rows=len(mc), files=int(mc.file_idx.nunique()), pt=[float(mc.pt.min()), float(mc.pt.max())], n_sig=int(len(sig)),
                                  n_sig_veto=int(((mc.is_sig == 1) & (mc.dstar_veto == 1)).sum()), sum_w_sig=float(sig.weight.sum()),
                                  w_minmax=[float(sig.weight.min()), float(sig.weight.max())], hist_sig=float(hm["dp_h"][1].sum()),
                                  sumwps_noveto=float(mc.w_ps[mc.dstar_veto == 0].sum()), hist_all=float(hm["dp_h"][0].sum()),
                                  truth_counts={str(k): int(v) for k, v in sig.truth.value_counts().items()})
    print("MC:", out["inputs"]["fd_nd"])
    s, w = sig.score.to_numpy(float), sig.weight.to_numpy(float)
    W = w.sum()
    # eps przy progach 3-5
    thr35 = json.load(open(os.path.join(SCAN, "cand", "train", "thr.json")))["H_s12345"]["3-5"]
    out["mc"]["eps_at_thr35"] = {k: float(w[s >= v].sum() / W) for k, v in thr35.items()}
    ptv = sig.pt.to_numpy(float)
    out["mc"]["eps_at_thr35_pt"] = {f"{lo}-{hi}": {k: float(w[(s >= v) & (ptv >= lo) & (ptv < hi)].sum() / w[(ptv >= lo) & (ptv < hi)].sum()) for k, v in thr35.items()}
                                    for lo, hi in ((2.0, 2.5), (2.5, 3.0))}
    sa = sig.score_avg.to_numpy(float)
    out["mc"]["eps_at_thr35_scoreavg"] = {k: float(w[sa >= v].sum() / W) for k, v in thr35.items()}
    out["mc"]["score_sd_mean_sig"] = float(sig.score_sd.mean())
    print("eps_S MC 2-3 przy progach WP z 3-5:", out["mc"]["eps_at_thr35"])
    print("  per pT:", out["mc"]["eps_at_thr35_pt"])
    print("  score_avg:", out["mc"]["eps_at_thr35_scoreavg"])
    # siatka WP: kwantyle -> krawędzie score
    o = np.argsort(-s)
    c = np.cumsum(w[o]) / W
    wps = []
    for x in GRID:
        t = s[o][min(np.searchsorted(c, x), len(o) - 1)]
        ie = int(np.argmin(np.abs(se - t)))
        te = float(se[ie])
        ps = s >= te
        sh = FL.mc_shape(sig.m.to_numpy(float)[ps], w[ps])
        fw = float(np.diff(FL.sig_cdf(np.array([M0 - FF.WIN, M0 + FF.WIN]), sh, 0, 1))[0])
        wps.append(dict(wp=x, thr_q=float(t), thr=te, ie=ie, eps=float(w[ps].sum() / W), n_mc=int(ps.sum()), shape=sh, f_win=fw))
        print(f"WP{int(x*100)}: thr_q {t:.4f} -> {te:.3f}  eps {wps[-1]['eps']:.4f} n_mc {ps.sum()}  mu {sh['mu']*1e3:.2f} rms_eff {sh['rms_eff']*1e3:.2f} MeV f_win {fw:.3f}")
    sh_all = FL.mc_shape(sig.m.to_numpy(float), w)
    out["mc"]["shape_all"] = sh_all
    out["mc"]["shape35_wp70"] = json.load(open(os.path.join(SCAN, "K1c", "out", "prep.json")))["shapes"]["H|3-5|0.7"]
    out["wp"]["grid"] = wps
    # zgodność hist MC sygnału z wierszami (po wecie, niew.)
    hs = np.histogram(s, se)[0]
    out["inputs"]["fd_nd"]["sig_rows_vs_hist_score_maxdiff"] = float(np.abs(hs - hm["dp_h"][1].sum((0, 2))).max())
    print("MC sygnał: wiersze vs hist (score) max|Δ| =", out["inputs"]["fd_nd"]["sig_rows_vs_hist_score_maxdiff"])

    # ---------- dane (zaślepione)
    for run, cls in RUNS.items():
        d = pd.read_parquet(os.path.join(IN, f"dp_data_{run}.parquet"),
                            columns=["run", "event", "m", "pt", "score", "dstar_veto", "ps_tier", "w_ps", "grl", "trig_HLT_mb_mbts_L1MBTS_1", "trig_HLT_noalg_L1J12"])
        hd = np.load(os.path.join(IN, f"hist_data_{run}.npz"))
        H = hd["dp_h"]
        ss = S["samples"][f"data_{run}"]
        ok = d.dstar_veto == 0
        if cls == 1:
            selc = (d.grl == 1) & (d["trig_HLT_mb_mbts_L1MBTS_1"] == 1)
        else:
            selc = np.ones(len(d), bool)
        info = dict(rows=len(d), rows_summary=ss["dp_rows"], files=ss["n_files"], pt=[float(d.pt.min()), float(d.pt.max())], runs=[int(x) for x in d.run.unique()],
                    tiers=[int((d.ps_tier == i).sum()) for i in range(4)], w_ps_vals=sorted(float(x) for x in d.w_ps.unique()),
                    hist_all=float(H[0].sum()), hist_cls=[float(H[i].sum()) for i in range(4)], cls=[str(x) for x in hd["cls"]],
                    sumwps_all=float(d.w_ps[ok].sum()))
        # Σw_ps vs hist: pasma (tier 0, poza oknem 60 MeV, score<0.93) i score>=0.93 (pełne) — okno sygnału tylko zsumowane w jedną liczbę
        mi = np.clip(np.searchsorted(me, d.m.to_numpy(), side="right") - 1, 0, 499)
        cm = 0.5 * (me[1:] + me[:-1])
        blind = np.abs(cm - M0) < BLIND
        i93 = int(np.argmin(np.abs(se - 0.93)))
        Hc = H[cls].sum(0)   # [score, m]
        dsel = d[ok & selc]
        lo = dsel.score < se[i93]
        outw = np.abs(dsel.m - M0) >= BLIND
        info["side_lo_sumw"] = float(dsel.w_ps[lo & outw].sum())
        info["side_lo_hist"] = float(Hc[:i93][:, ~blind].sum())
        info["win_lo_sumw"] = float(dsel.w_ps[lo & ~outw].sum())
        info["win_lo_hist"] = float(Hc[:i93][:, blind].sum())
        info["hi_rows"] = int((~lo).sum())
        info["hi_hist"] = float(Hc[i93:].sum())
        info["hi_w_ps_all1"] = bool((dsel.w_ps[~lo] == 1).all())
        info["n_events_hi"] = int(dsel[~lo][["run", "event"]].drop_duplicates().shape[0])
        info["frac_cls"] = float(Hc.sum() / H[0].sum())
        info["frac_L1J12"] = float(H[2].sum() / H[0].sum())
        info["frac_mbts"] = float(H[1].sum() / H[0].sum())
        out["inputs"][run] = info
        print(f"\n{run}: wiersze {len(d)} (summary {ss['dp_rows']}), pliki {ss['n_files']}, pT {info['pt']}, tiers {info['tiers']}, w_ps {info['w_ps_vals']}")
        print(f"  hist klasy {info['hist_cls']}  udział klasy analizy {info['frac_cls']:.3f}, mbts {info['frac_mbts']:.3f}, L1J12 {info['frac_L1J12']:.3f}")
        print(f"  Σw_ps pasma score<0.93: {info['side_lo_sumw']:.0f} vs hist {info['side_lo_hist']:.0f} (stosunek {info['side_lo_sumw']/info['side_lo_hist']:.4f}; 1/√n = {1/np.sqrt(info['side_lo_sumw']/100):.4f})")
        print(f"  Σw_ps okno  score<0.93: {info['win_lo_sumw']:.0f} vs hist {info['win_lo_hist']:.0f} (stosunek {info['win_lo_sumw']/info['win_lo_hist']:.4f})")
        print(f"  score>=0.93: wiersze {info['hi_rows']} vs hist {info['hi_hist']:.0f}; w_ps==1: {info['hi_w_ps_all1']}; zdarzeń {info['n_events_hi']}")
        # pasma boczne per WP
        rows = []
        Hm5 = Hc.reshape(Hc.shape[0], 100, 5).sum(2)
        for wpd in wps:
            h = Hm5[wpd["ie"]:].sum(0).astype(float)
            hb = h.copy()
            B, chi2, ndf, valid = sb_fit(hb)
            B2 = sb_fit(hb, 2)[0]
            B4 = sb_fit(hb, 4)[0]
            Sx = 0.6 * Y35[run] * wpd["eps"] * wpd["f_win"]
            Z = Sx / np.sqrt(Sx + B)
            n_sb = float(h[np.abs(FF.C - M0) > BLIND].sum())
            rows.append(dict(wp=wpd["wp"], thr=wpd["thr"], eps=wpd["eps"], n_side=n_sb, B=B, B_pol2=B2, B_pol4=B4, chi2=chi2, ndf=ndf, valid=valid, S_exp=Sx, Z_exp=float(Z), SB_exp=float(Sx / B)))
            print(f"  WP{int(wpd['wp']*100)} thr {wpd['thr']:.3f} eps {wpd['eps']:.3f}: pasma {n_sb:.0f}  B(±45) {B:.0f} (pol2 {B2:.0f}, pol4 {B4:.0f}; chi2 {chi2:.0f}/{ndf})  S_exp {Sx:.0f}  S/B {Sx/B:.4f}  Z_exp {Z:.2f}")
        out["wp"][run] = rows
        zs = np.array([r["Z_exp"] for r in rows])
        ib = int(np.argmax(zs))
        # remis < 2 % -> luźniejszy
        for j in range(ib):
            if zs[j] > 0.98 * zs[ib]:
                ib = j
                break
        out["choice"][run] = dict(wp=rows[ib]["wp"], thr=rows[ib]["thr"], eps=rows[ib]["eps"], Z_exp=rows[ib]["Z_exp"], S_exp=rows[ib]["S_exp"], B=rows[ib]["B"],
                                  argmax_wp=rows[int(np.argmax(zs))]["wp"])
        print(f"  WYBÓR {run}: WP{int(rows[ib]['wp']*100)} (argmax WP{int(rows[int(np.argmax(zs))]['wp']*100)})")
    json.dump(out, open(os.path.join(TOP, "out", "wp.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
