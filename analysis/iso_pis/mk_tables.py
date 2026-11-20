#!/usr/bin/env python3
"""ISO_pis — tabele markdown z iso_pis.json (sam odczyt). Użycie: mk_tables.py iso_pis.json > tables.md"""
import json, sys
r = json.load(open(sys.argv[1]))
f = lambda x, n=3: ("%.*f ± %.*f" % (n, x[0], n, x[1])) if x else "–"
g = lambda x, n=3: ("%.*f" % (n, x)) if x is not None else "–"
for s, S in r["samples"].items():
    print("\n## %s\n" % s)
    print("sygnał inkluzywny %d, próbka tagu %d, D0 z tagiem w tabeli inkluzywnej %d; własny stożek (tokeny) = gałąź w 1 MeV: %.3f; gałąź − pT(π_s) < 0: %.4f\n"
          % (S["n_sig_incl"], S["n_tag"], S["n_incl_tagged"], S["own_vs_branch_sig_incl"]["abs_lt_1MeV"], S["branch_minus_pis_in40"]["frac_negative_gt_5MeV"]))
    print("### T1. Wolny pion w stożku (próbka tagu)\n")
    print("| pT D0 | n tag | pT(π_s) q10/50/90 [GeV] | π_s pT < 0,2 (brak w Tracks) | ΔR < 0,4 | ΔR < 0,4 i pT > 0,5 | ΔR < 0,4 i pT > 1 | ΔR < 1,0 | ΔR q50/q99 | wkład do ptcone40 q10/50/90 [GeV] | udział w ptcone40 q25/50/75 | stożek = sam π_s | ⟨ptcone40⟩ tag / bez π_s / reszta b | π_s przenosi przez 2 GeV | udział b-hadronów tag / b |")
    print("|" + "---|" * 15)
    for b, o in S["bins"].items():
        p = o["pis"]
        print("| %s | %d | %s / %s / %s | %s | %s | %s | %s | %s | %s / %s | %s / %s / %s | %s / %s / %s | %s | %s / %s / %s | %s | %s / %s |" % (
            b, o["n_tag"], g(p["pt_q10_25_50_75_90"][0], 2), g(p["pt_q10_25_50_75_90"][2], 2), g(p["pt_q10_25_50_75_90"][4], 2), g(p["frac_pt_lt_0.2"]), g(p["in40"]),
            g(p["in40_pt_gt_0.5"]), g(p["in40_pt_gt_1.0"]), g(p["in100"]), g(p["dr_q50_90_99"][0], 2), g(p["dr_q50_90_99"][2], 2),
            g(p["contrib40_GeV_q10_25_50_75_90"][0], 2), g(p["contrib40_GeV_q10_25_50_75_90"][2], 2), g(p["contrib40_GeV_q10_25_50_75_90"][4], 2),
            g(p["share_of_ptcone40_q25_50_75"][0], 2), g(p["share_of_ptcone40_q25_50_75"][1], 2), g(p["share_of_ptcone40_q25_50_75"][2], 2), g(p["cone_only_pis"]),
            g(p["ptcone40_mean_tag"], 2), g(p["ptcone40_mean_tag_nopis"], 2), g(p["ptcone40_mean_b"], 2), g(p["pis_crosses_2GeV"]), g(o["fb_tag"]), g(o["fb_b"])))
    print("\n### T2. Wydajności cięć izolacyjnych\n")
    print("| cięcie | pT D0 | ε_a (tag D*) | ε_b (reszta) | ε_c (tag bez π_s) | ε_i (inkl.) | a/b | a/i | c/a | c/b | c/b prompt | c/b z b | am/bm (maska bez prawdy) | a/b po przeważeniu pT |")
    print("|" + "---|" * 14)
    for c in r["cuts"]:
        for b, o in S["bins"].items():
            q = o["cuts"][c]
            print("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |" % (c, b, f(q["a"]), f(q["b"]), f(q["c"]), f(q["i"]), f(q["a/b"]), f(q["a/i"]), f(q["c/a"]),
                  f(q["c/b"]), f(q["c/b_t0"]), f(q["c/b_t1"]), f(q["am/bm"]), f(q["a_rw/b"])))
    print("\n### T3. Efekt drugiego rzędu: SF_inkl / SF_tag, gdy aktywność spoza π_s przesunięta tak, by SF_tag = cel (π_s bez zmian)\n")
    print("| cięcie | pT D0 | model | " + " | ".join("SF_tag %s" % t for t in ("0.98", "0.95", "0.90", "0.76", "0.51")) + " |")
    print("|" + "---|" * 8)
    for c in r["cuts"]:
        for b, o in S["bins"].items():
            so = o["cuts"][c]["second_order"]
            for m in ("scale", "shift"):
                print("| %s | %s | %s | " % (c, b, {"scale": "u → k·u", "shift": "u → u + δ"}[m]) + " | ".join(
                    ("%.3f" % so[m][t]["SF_incl/SF_tag"]) if so[m][t] else "–" for t in ("0.98", "0.95", "0.90", "0.76", "0.51")) + " |")
    print("\n### T4. Wrażliwość ε inkluzywnej (ptcone40 < 2 GeV) na wydajność rekonstrukcji π_s\n")
    print("| pT D0 | ε_i | ε_i bez torów maski | Δε_i/ε_i na 10 % utraty π_s | udział tagu w inkl. | udział b1 w b |")
    print("|---|---|---|---|---|---|")
    for b, o in S["bins"].items():
        x = o["pis_eff_sens_ptcone40lt2"]
        print("| %s | %.4f | %.4f | %+.4f | %.3f | %.3f |" % (b, x["eps_i"], x["eps_i_all_mask_removed"], x["d_eps_rel_per_10pct_pis_loss"], o["tag_frac_of_incl"], o["b1_frac_of_b"]))
