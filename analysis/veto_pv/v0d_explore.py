import numpy as np
R = "/home/k/Documents/charm/charm_work/B_trig/runs/lxplus/"
for s in ["fd_nd", "mb_nd", "data_428770", "data_427929"]:
    f = np.load(R + f"bt_{s}.npz")
    k = f["ev_file_idx"].astype(np.int64) * (1 << 40) + f["ev_event"]
    print(s, "ev", len(k), "uniq (file,event)", len(np.unique(k)), "uniq event", len(np.unique(f["ev_event"])), "nfile_idx", len(np.unique(f["ev_file_idx"])))
    ks = np.sort(k)
    for pre in ["", "all_"]:
        kk = f[pre + "file_idx"].astype(np.int64) * (1 << 40) + f[pre + "event"]
        i = np.clip(np.searchsorted(ks, kk), 0, len(ks) - 1); ok = ks[i] == kk
        o = np.argsort(k); j = o[i]
        print("  ", pre or "pairs", "matched", ok.mean(), "npu agree", (f["ev_npu"][j] == f[pre + "e_npu"])[ok].mean(), "nSelBL agree", (f["ev_nSelBL"][j] == f[pre + "e_nSelBL"])[ok].mean())
