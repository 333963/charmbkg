#!/usr/bin/env python3
"""p0: lista gałęzi AFP/Truth w ROOT v2.0 (MC mb_sd i dane) + kilka zdarzeń. Tylko odczyt."""
import sys, glob, uproot, numpy as np
for s in sys.argv[1:]:
    fs = sorted(glob.glob(f"/eos/user/k/klyczak/charm_raw_v20/{s}/*dpt3to150*.root"))
    f = uproot.open(fs[0]); print("##", s, len(fs), fs[0].split("/")[-1], [k for k in f.keys()][:8])
    tn = [k for k in f.keys() if f[k].classname.startswith("TTree")][0]; t = f[tn]
    print("tree", tn, t.num_entries)
    for b in t.keys():
        bl = b.lower()
        if any(x in bl for x in ("afp", "truthproton", "xi_cal", "truthd", "mcchannel", "pileup", "nprim", "grl", "lumiblock", "actual_mu")) and "trig" not in bl and "hlt" not in bl and "l1_" not in bl:
            print("  ", b, t[b].typename)
    print("  n trig-like branches:", sum(1 for b in t.keys() if "HLT" in b or "L1_" in b))
    print("  L1 AFP:", [b for b in t.keys() if "L1_AFP" in b][:30])
