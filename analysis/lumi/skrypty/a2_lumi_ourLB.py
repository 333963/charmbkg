"""a2: świetlność w LB obecnych w naszych plikach (GRL analizy), z oficjalnego ilumicalc (cvmfs, OflLumi-Run3-003, livetrigger L1_EM22VHI)."""
import pandas as pd, numpy as np, glob
T = "/home/k/Documents/charm/charm_work/TAGREG_2026-11-18/runs/lxplus"
off = {f.split("tabele/")[1][:14] + ("_SIT" if "SIT__" in f else "_SIT_A" if "SIT_A" in f else "_SIT_C" if "SIT_C" in f else ""): pd.read_csv(f) for f in sorted(glob.glob("tabele/cvmfs_*__LumiMetaData.csv"))}
for k, a in off.items():
    print("==", k)
    for run in (427929, 428770):
        g = a[a.RunNbr == run]; lb = pd.read_csv(f"{T}/lb_summary_{run}.csv"); ours = lb[lb.n_grl > 0].lumiBlock
        m = g[g.LBStart.isin(ours)]
        miss = sorted(set(ours) - set(g.LBStart))
        print(f"  {run}: oficjalne LB {len(g)} ({g.LBStart.min()}-{g.LBStart.max()}), L_all = {g.IntLumi.sum()/1e3:.3f} nb-1 | nasze LB {len(ours)}, znalezione {len(m)}, brak {miss[:10]}"
              f" | L(nasze LB) = {m.IntLumi.sum()/1e3:.4f} nb-1, livefrac śr.waż. {np.average(m.Livefrac, weights=m.IntLumi/m.Livefrac):.4f}, min {m.Livefrac.min():.3f}, <mu> {np.average(m.AvergeInteractionPerXing, weights=m.IntLumi):.4f}, LArfrac min {m.LArfrac.min():.3f}, czas {m.DeltaT.sum():.0f} s")
        if run == 428770:
            for nm, sel in (("LB 4-25 + 210-324 (GRL analizy)", g[((g.LBStart >= 4) & (g.LBStart <= 25)) | ((g.LBStart >= 210) & (g.LBStart <= 324))]), ("LB 4-26 + 207-325 (GRL oficjalny)", g[((g.LBStart >= 4) & (g.LBStart <= 26)) | ((g.LBStart >= 207) & (g.LBStart <= 325))]), ("wszystkie w pliku", g)):
                print(f"      428770 {nm}: n_LB {len(sel)}, L = {sel.IntLumi.sum()/1e3:.4f} nb-1")
