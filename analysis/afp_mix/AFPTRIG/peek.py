import json, sys
run = sys.argv[1]
O = json.load(open(f"results_raw_{run}.json"))[run]
def f(e): return f"{e['n_acc']:.1f}±{e['err']:.1f}"
for st, Os in O.items():
    print("=====", run, st, "pool", Os["nominal_pool"])
    for s in "AC":
        o = Os[s]; n = o["nominal"]
        print(f" {s}: SR {n['n_sr']:.0f} CR {n['n_cr']:.0f} Nacc {f(n)} f={n['f_acc']:.4f}±{n['f_acc_err_stat']:.4f} pool {n['n_pool']} prot {n['n_protons']}")
        print("   pool strat  :", {k: f(v) for k, v in o["pool"].items()})
        print("   pool nostrat:", {k: f(v) for k, v in o["pool_nostrat"].items()})
        print("   protons     :", {k: f(v) for k, v in o["protons"].items()})
        print("   strat       :", {k: f(v) for k, v in o["strat"].items()})
        print("   margin      :", {k: f(v) for k, v in o["margin"].items()})
        for k, v in o["closure"].items():
            print("   closure", k, {p: f"obs {e['n_sr']:.0f} pred {e['n_acc']:.1f}±{e['err']:.1f} r={e['n_sr']/max(e['n_acc'],1e-9):.3f}" for p, e in v.items()})
        for k, v in o["dep"].items():
            print("   dep", k, [f"{r['lo']:.3g}-{r['hi']:.3g}: {r['f_acc']:.4f}±{r['f_acc_err_stat']:.4f}" if r['f_acc'] is not None else None for r in v])
