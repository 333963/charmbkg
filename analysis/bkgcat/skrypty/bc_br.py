#!/usr/bin/env python3
"""BKGCAT: tablica rozpadów mb_nd z prawdy generatora (TruthDmesons, |η| < 2,5) vs PDG 2026 + wagi BR_PDG/BR_MC per tryb (do przeważenia szablonu).
PDG: przepisane z listingów otwartych 19/20.11.2026 (pdg.lbl.gov/2026/listings/rpp2026-list-{D-zero,D-plus-minus,D-star-2010-plus-minus,eta}.pdf, pdftotext);
numer Γ_i podany przy każdej pozycji. Lista córek w TruthDmesons zawiera wyłącznie π±, K±, π0 (i D0) → „tryb” = topologia widzialna + klasa brakującej masy:
'' (m_D − m_wid < 30 MeV), +X0 (m_brak < 100 MeV: foton/y, ℓν przy małej masie), +K0, +eta (±3 MeV), +X (reszta). POST HOC.
Wyjście: tabele/br.md, out/br.json (wagi)."""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bc_lib as L  # noqa: E402

PDGS = (421, 411, 413)
MCN = ["", "+X0", "+K0", "+eta", "+X"]
ALL = (0, 1, 2, 3, 4)
# (hadron, topologia, klasy sumowane, PDG [%], błąd [%], opis ze źródłem)
PDG = [
    (421, "K- pi+", (0, 1), 3.936, 0.030, "Γ47 K⁻π⁺ (PDG zawiera FSR)"),
    (421, "K+ pi-", (0, 1), 0.0148, 0.0007, "K⁺π⁻ DCS + mieszanie (BR_frag §1, PDG 2026)"),
    (421, "K- pi+ pi0", (0, 1), 14.4, 0.5, "Γ68 K⁻π⁺π⁰"),
    (421, "K- pi+ pi+ pi-", (0, 1), 8.20, 0.14, "Γ87 K⁻2π⁺π⁻"),
    (421, "K- pi+ pi0 pi0", (0, 1), 8.86, 0.23, "Γ101 K⁻π⁺2π⁰"),
    (421, "K- pi+ pi+ pi- pi0", (0, 1), 4.3, 0.4, "Γ104 K⁻2π⁺π⁻π⁰"),
    (421, "K- K+", (0, 1), 0.407, 0.006, "Γ266 K⁺K⁻"),
    (421, "pi+ pi-", (0, 1), 0.1451, 0.0024, "Γ151 π⁺π⁻"),
    (421, "pi+ pi- pi0", (0, 1), 1.49, 0.06, "Γ153 π⁺π⁻π⁰"),
    (421, "pi+ pi- pi0 pi0", (0, 1), 1.002, 0.031, "Γ212 π⁺π⁻2π⁰"),
    (421, "K- K+ pi0", (0, 1), 0.341, 0.014, "Γ289 K⁺K⁻π⁰"),
    (421, "K- K+ pi+ pi-", (0, 1), 0.278, 0.017, "Γ303 K⁺K⁻π⁺π⁻"),
    (421, "pi0 pi0", (0, 1), 0.0825, 0.0025, "Γ152 2π⁰"),
    (421, "K-", ALL, 6.955, 0.025, "Γ20 + Γ21 K⁻e⁺ν + K⁻μ⁺ν (topologia „samo K⁻”, wszystkie klasy)"),
    (421, "pi-", ALL, 3.405, 0.07, "Γ39 + Γ40 + Γ25 + Γ29: π⁻ℓ⁺ν + K̄⁰π⁻ℓ⁺ν (topologia „samo π⁻”)"),
    (421, "pi+ pi-", (2,), 5.68, 0.32, "2 × Γ53 K⁰_Sπ⁺π⁻ (założenie K⁰_L = K⁰_S)"),
    (421, "pi0", (2,), 2.215, 0.04, "Γ48 + Γ49 K⁰_Sπ⁰ + K⁰_Lπ⁰"),
    (421, "K- pi+", (2,), 0.68, 0.10, "2 × Γ268 K⁰_SK⁻π⁺ (założenie K⁰_L = K⁰_S)"),
    (421, "K- pi+", (3,), 0.740, 0.020, "Γ113 K⁻π⁺η × B(η→2γ) = 1,88 % × 39,36 %"),
    (411, "K- pi+ pi+", (0, 1), 9.31, 0.14, "Γ68 K⁻2π⁺"),
    (411, "K- pi+ pi+ pi0", (0, 1), 6.35, 0.08, "Γ93 K⁻2π⁺π⁰"),
    (411, "K- K+ pi+", (0, 1), 0.961, 0.017, "Γ177 K⁺K⁻π⁺"),
    (411, "K- K+ pi+ pi0", (0, 1), 0.662, 0.032, "Γ186 K⁺K⁻π⁺π⁰"),
    (411, "pi+ pi+ pi-", (0, 1), 0.327, 0.009, "Γ131 2π⁺π⁻"),
    (411, "pi+ pi0", (0, 1), 0.1243, 0.0033, "Γ130 π⁺π⁰"),
    (411, "pi+ pi+ pi- pi0", (0, 1), 1.165, 0.030, "Γ147 2π⁺π⁻π⁰"),
    (411, "K- pi+ pi+ pi+ pi-", (0, 1), 0.57, 0.05, "Γ117 K⁻3π⁺π⁻"),
    (411, "K- pi+", ALL, 7.62, 0.37, "Γ22 + Γ33 K⁻π⁺e⁺ν + K⁻π⁺μ⁺ν (topologia „K⁻π⁺”)"),
    (411, "pi+", (2,), 3.014, 0.06, "Γ66 + Γ67 K⁰_Sπ⁺ + K⁰_Lπ⁺"),
    (411, "pi+ pi0", (2,), 14.62, 0.40, "2 × Γ78 K⁰_Sπ⁺π⁰ (założenie K⁰_L = K⁰_S)"),
    (411, "pi+ pi+ pi-", (2,), 6.16, 0.16, "2 × Γ107 K⁰_S2π⁺π⁻ (założenie K⁰_L = K⁰_S)"),
    (411, "K+", (2,), 0.623, 0.018, "Γ170 + Γ171 K⁰_SK⁺ + K⁰_LK⁺"),
    (411, "(nic)", ALL, 17.49, 0.12, "Γ17 + Γ18 K̄⁰e⁺ν + K̄⁰μ⁺ν (dolne ograniczenie topologii „nic widzialnego”)"),
]


def mode_str(code):
    code = int(code)
    n = [(code // 10 ** i) % 10 for i in range(6)]
    s = "".join(x * k for x, k in zip(("D0 ", "K- ", "K+ ", "pi+ ", "pi- ", "pi0 "), (n[5], n[0], n[1], n[2], n[3], n[4]))).strip()
    return s if s else "(nic)"


def main():
    fn = sys.argv[sys.argv.index("--in") + 1] if "--in" in sys.argv else L.IN
    z = np.load(fn, allow_pickle=False)
    meta = json.loads(str(z["meta"]))
    tabs = {}
    for r in (0, 1, 2):
        k, n = z[f"br{r}_key"].astype(np.int64), z[f"br{r}_n"].astype(np.float64)
        tabs[r] = {"pi": k % 4, "fromB": (k // 4) % 2, "mc": (k // 8) % 8, "code": k // 64, "n": n}
    out = ["# BKGCAT — rozpady D w prawdzie generatora mb_nd vs PDG 2026 (POST HOC)", "",
           f"mb_nd: {meta['n_files']} plików, {sum(p['nev'] for p in meta['parts'])} zdarzeń, {sum(p['n_truthD'] for p in meta['parts'])} wpisów TruthDmesons; "
           f"hadrony inne niż D⁰, D⁺, D*⁺: {meta['pdg_other']} (**D_s i Λ_c nie są zapisywane w TruthDmesons — tablicy dla nich nie da się zbudować**).",
           "Region: |η| < 2,5 (klasa brakującej masy wymaga precyzji float; ułamki bez cięcia η i dla pT ≥ 3 GeV, |η| < 2,2 — kolumny kontrolne). Prompt + z b razem.",
           "Lista córek zawiera tylko π±, K±, π⁰ (rezonanse rozwinięte; K⁰, η→γγ, leptony, ν, γ niewidoczne) → tryb = topologia widzialna + klasa braku.", ""]
    J = {"modes": {}, "weights": {}}
    for ip, P in enumerate(PDGS):
        fr = {}
        for r in (0, 1, 2):
            t = tabs[r]
            m = t["pi"] == ip
            tot = t["n"][m].sum()
            d = {}
            for c, mc, n in zip(t["code"][m], t["mc"][m], t["n"][m]):
                d[(int(c), int(mc))] = d.get((int(c), int(mc)), 0.0) + n
            fr[r] = (d, tot)
        d1, tot1 = fr[1]
        # udział z b
        t = tabs[1]
        fb = t["n"][(t["pi"] == ip) & (t["fromB"] == 1)].sum() / tot1
        out += [f"## |pdg| = {P}: {tot1:.0f} wpisów w |η| < 2,5 (udział z b: {100 * fb:.1f} %)", ""]
        if P != 413:
            out += ["### Porównanie z PDG", "", "| tryb (topologia, klasy) | MC [%] | ± stat | PDG [%] | ± | MC/PDG | źródło PDG | MC bez cięcia η | MC pT ≥ 3, \\|η\\| < 2,2 |", "|---|---|---|---|---|---|---|---|---|"]
            for (pp, topo, cls, val, err, src) in PDG:
                if pp != P:
                    continue
                def f(r):
                    d, tot = fr[r]
                    n = sum(v for (c, mc), v in d.items() if mode_str(c) == topo and mc in cls)
                    return n / tot, np.sqrt(max(n, 1)) / tot
                x, ex = f(1)
                lab = topo + ("" if cls == (0, 1) else " [wszystkie klasy]" if cls == ALL else MCN[cls[0]])
                out.append(f"| {lab} | {100 * x:.3f} | {100 * ex:.3f} | {val:g} | {err:g} | {100 * x / val:.3f} | {src} | {100 * f(0)[0]:.3f} | {100 * f(2)[0]:.3f} |")
                J["modes"][f"{P}|{lab}"] = {"mc": x, "mc_err": ex, "pdg": val / 100, "pdg_err": err / 100, "ratio_mc_pdg": 100 * x / val}
                for mc in cls:
                    J["weights"][f"{P}|{topo}{MCN[mc]}"] = val / (100 * x) if x > 0 else 1.0
            out.append("")
        else:
            nd0 = sum(v for (c, mc), v in d1.items() if (c // 100000) % 10 > 0)
            out += [f"D*⁺ → D⁰π⁺ (lista zawiera D⁰): {100 * nd0 / tot1:.2f} ± {100 * np.sqrt(nd0) / tot1:.2f} %; PDG Γ1 67,7 ± 0,5 % (MC/PDG = {nd0 / tot1 / 0.677:.3f}); "
                    f"bez D⁰ (D⁺π⁰, D⁺γ): {100 * (1 - nd0 / tot1):.2f} %; PDG Γ2 + Γ3 = 32,3 %.", ""]
            J["modes"]["413|D0 pi+"] = {"mc": nd0 / tot1, "pdg": 0.677, "ratio_mc_pdg": nd0 / tot1 / 0.677}
        out += ["### Wszystkie tryby ≥ 0,2 % (topologia × klasa braku)", "", "| tryb | MC [%] | N |", "|---|---|---|"]
        for (c, mc), n in sorted(d1.items(), key=lambda kv: -kv[1]):
            if n / tot1 < 0.002:
                break
            out.append(f"| {mode_str(c)}{MCN[mc]} | {100 * n / tot1:.3f} | {n:.0f} |")
        out.append("")
    # decayType sygnału
    for dty, P, nm_, br in ((10, 421, "D⁰, decayType 10 („K π”)", 3.936), (4, 411, "D⁺, decayType 4 („K π π”)", 9.31)):
        v = z[f"dty{dty}_mc"].astype(float)
        d1, tot1 = None, None
        ip = PDGS.index(P)
        t = tabs[1]
        tot1 = t["n"][t["pi"] == ip].sum()
        f_all = v.sum() / tot1
        out += [f"## {nm_}: skład klas braku (|η| < 2,5)", "",
                f"Wpisów: {v.sum():.0f} = {100 * f_all:.3f} % wszystkich; topologie w tej klasie: {[mode_str(c) for c in z[f'dty{dty}_codes']]}.", "",
                "| klasa | N | udział | × ułamek decayType = BR_MC [%] |", "|---|---|---|---|"]
        for i, nme in enumerate(["ekskluzywny (m_D − m_wid < 30 MeV)", "+X0 (radiacyjny, m_brak < 100 MeV)", "+K0", "+eta", "+X (inne)"]):
            out.append(f"| {nme} | {v[i]:.0f} | {100 * v[i] / v.sum():.2f} % | {100 * f_all * v[i] / v.sum():.3f} |")
        fx = (v[0] + v[1]) / v.sum()
        e = np.sqrt(fx * (1 - fx) / v.sum())
        out += ["", f"**Ułamek prawdziwego trybu sygnału (ekskl. + radiacyjny) w decayType: {100 * fx:.2f} ± {100 * e:.2f} %; BR_MC(tryb sygnału) = {100 * f_all * fx:.3f} % wobec PDG {br} % "
                f"(MC/PDG = {100 * f_all * fx / br:.3f}); „cały decayType”/PDG = {100 * f_all / br:.3f}.**", ""]
        J[f"dty{dty}"] = {"frac_of_all": f_all, "f_true_mode": fx, "f_true_mode_err": e, "classes": v.tolist()}
    # skutek dla C z XSEC: mianownik = cały decayType, licznik = Sig w zakresie mas przy WP70
    B = L.BK(fn)
    out += ["## Skutek uboczny dla XSEC (C = N_reco^prawda / N_gen(decayType sygnału)) — do weryfikacji w etapie XSEC; XSEC nie zmieniany", "",
            "f_gen = ułamek prawdziwego trybu sygnału (ekskl. + radiacyjny) w mianowniku (pT ≥ 3 GeV, |η| < 2,2); f_reco = to samo w liczniku (Sig w zakresie mas, WP70).",
            "Jeżeli σ ma dotyczyć trybu z BR PDG, to σ_popr/σ_XSEC = f_gen/f_reco.", "", "| gat. | bin | f_gen | f_reco (WP70) | N_reco Sig: ekskl. + rad. / z brakiem (K0, η, X) | σ_popr/σ_XSEC |", "|---|---|---|---|---|---|"]
    t = tabs[2]
    J["xsec"] = {}
    for sp, P, topo in (("D0", 421, "K- pi+"), ("Dplus", 411, "K- pi+ pi+")):
        m = (t["pi"] == PDGS.index(P)) & np.array([mode_str(c) == topo for c in t["code"]])
        ntot = t["n"][m].sum()
        nex = t["n"][m & (t["mc"] <= 1)].sum()
        fg = nex / ntot
        for bn in L.BINS:
            c = B.cat(sp, bn, 0.7)
            a = sum(h.sum() for k, h in c.items() if L.group(k) == "Sig")
            ex = sum(h.sum() for k, h in c.items() if L.group(k) == "Sig" and k.split("|")[1] in (topo, topo + "+X0"))
            out.append(f"| {sp} | {bn} | {fg:.4f} | {ex / a:.4f} | {ex:.0f} / {a - ex:.0f} | {fg / (ex / a):.3f} |")
            J["xsec"][f"{sp}|{bn}"] = {"f_gen": fg, "f_reco": ex / a, "ratio": fg / (ex / a)}
    out.append("")
    # histogram m_brak dla K- pi+ (D0) — piki
    k, n = z["mm_key"].astype(np.int64), z["mm_n"].astype(float)
    mb, pi, code = k % 400, (k // 400) % 4, k // 1600
    for P, topo in ((421, "K- pi+"), (411, "K- pi+ pi+")):
        m = (pi == PDGS.index(P)) & np.array([mode_str(c) == topo for c in code])
        h = np.zeros(400)
        np.add.at(h, mb[m], n[m])
        top = np.argsort(-h)[:8]
        out += [f"m_brak dla topologii „{topo}” ({P}), klasy ≠ ekskluzywna: największe biny 5 MeV: " + ", ".join(f"{5 * i}–{5 * i + 5} MeV: {h[i]:.0f}" for i in sorted(top)) + f"; razem {h.sum():.0f}.", ""]
    open(os.path.join(L.BASE, "tabele", "br.md"), "w").write("\n".join(out) + "\n")
    json.dump(J, open(os.path.join(L.BASE, "out", "br.json"), "w"), indent=1)
    print("zapisano tabele/br.md, out/br.json")


if __name__ == "__main__":
    main()
