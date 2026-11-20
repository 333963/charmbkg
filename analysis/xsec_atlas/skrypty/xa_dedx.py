#!/usr/bin/env python3
"""XSEC_atlas C — dE/dx: wplyw przesuniecia MC (charmpy training.yaml: dedx_mc_shift = -0,023, addytywnie do surowego dE/dx corek) na wydajnosc
drzewa H przy stalym progu WP70, per gatunek i bin pT. Tylko MC (sygnal prawdy fd_nd + mb_nd z tabel INCAL)."""
import json, os, sys
import numpy as np
H = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(H, "skrypty")); sys.path.insert(0, os.path.join(os.path.dirname(H), "INCAL_2026-11-19", "scripts"))
import xa_lib as X
import inc_lib as IL
THR = {"D0": {"3-5": 0.822, "5-8": 0.8719, "8-30": 0.8631}, "Dplus": {"3-5": 0.9415, "5-8": 0.9517, "8-30": 0.9346}}
SHIFT = -0.023
out = {}
L = ["| gatunek | bin | N sygnału MC | ε nominalna | max \\|Δscore\\| odtw. | ε(dE/dx − 0,023)/ε | ε(dE/dx + 0,023)/ε | systematyka ± [%] | ⟨dE/dx K⟩ |", "|---|---|---|---|---|---|---|---|---|"]
for sp in ("D0", "Dplus"):
    T = {}
    for s in ("fd_nd", "mb_nd"):
        t = X.load(s, f"{sp}__sig_")
        for k, v in t.items():
            T.setdefault(k, []).append(v)
    T = {k: np.concatenate(v) for k, v in T.items()}
    sc0 = IL.score(sp, T["uid"], IL.feats(T, sp))
    dmax = float(np.nanmax(np.abs(sc0 - T["score"])))
    def shifted(d):
        O = dict(T)
        for i in range(IL.ND[sp]):
            x = T[f"d{i}_dedx"].astype(float)
            O[f"d{i}_dedx"] = np.where(x > 0, x + d, x).astype(np.float32)
        return IL.score(sp, T["uid"], IL.feats(O, sp))
    scm, scp = shifted(SHIFT), shifted(-SHIFT)
    for bn, (lo, hi) in (("3-5", (3, 5)), ("5-8", (5, 8)), ("8-30", (8, 30))):
        m = IL.sigmask(T, sp) & (T["pt"] >= lo) & (T["pt"] < hi)
        thr = THR[sp][bn]
        p0, pm, pp = (sc0[m] >= thr), (scm[m] >= thr), (scp[m] >= thr)
        e0 = p0.mean(); rm, rp = pm.mean() / e0, pp.mean() / e0
        # blad stosunku: tylko zdarzenia zmieniajace decyzje
        em = np.sqrt((p0 != pm).sum()) / p0.sum(); ep = np.sqrt((p0 != pp).sum()) / p0.sum()
        sy = max(abs(rm - 1), abs(rp - 1))
        out[f"{sp}|{bn}"] = {"n": int(m.sum()), "eps": float(e0), "r_minus": float(rm), "r_plus": float(rp), "r_minus_err": float(em), "r_plus_err": float(ep), "syst": float(sy)}
        L.append(f"| {sp} | {bn} | {m.sum()} | {e0:.3f} | {dmax:.1e} | {rm:.4f} ± {em:.4f} | {rp:.4f} ± {ep:.4f} | {100 * sy:.1f} | {np.nanmean(np.where(T['d0_dedx'][m] > 0, T['d0_dedx'][m], np.nan)):.3f} |")
print("\n".join(L))
json.dump(out, open(os.path.join(H, "out/dedx.json"), "w"), indent=1)
open(os.path.join(H, "tabele/dedx.md"), "w").write("# dE/dx: przesunięcie MC ∓0,023 (charmpy `training.yaml`: `dedx_mc_shift: -0.023`) → wydajność drzewa H przy stałym progu WP70\n\n"
    "Sygnał prawdy (fd_nd + mb_nd, tabele INCAL `__sig_`), mianownik: preselekcja nominalna (L_xy, weto D*). Systematyka = max |ε(±)/ε − 1| (pełna zmiana; nasze MC nie ma przesunięcia).\n\n" + "\n".join(L) + "\n")
