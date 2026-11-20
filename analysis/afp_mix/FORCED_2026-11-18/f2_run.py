"""FORCED f2 — bieg AFP-MIX v2 (run_data.py BEZ ZMIAN, importowany) z dodatkowym wykluczeniem zdarzeń u biorców i dawców.
Jedyna ingerencja: opakowanie run_data.trigger_classes — zdarzenia z maski wariantu tracą status „tryger nie-AFP” (wypadają z selekcji).
Użycie: f2_run.py <W1a|W1|W2>  → <wariant>/RAPORT_DANE_v2.md, result_v2.json. Tylko frakcje tagu przypadkowego i liczniki."""
import os, sys, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); AM = os.path.dirname(HERE); sys.path.insert(0, AM)
import run_data as RD
var = sys.argv[1]
ex = pd.read_parquet(f"{HERE}/out/excl_428770.parquet"); exen = ex.event_number.to_numpy()[ex[var].to_numpy()]
_orig = RD.trigger_classes
def patched(ev):
    nonafp, afp, tab, cols = _orig(ev)
    m = np.isin(ev["event_number"].to_numpy(), exen)
    print(f"[FORCED {var}] wykluczam {int(m.sum())} zdarzeń (z nich nie-AFP wg 11 bitów: {int((m & nonafp).sum())})", flush=True)
    return nonafp & ~m, afp, tab, cols
RD.trigger_classes = patched
sys.exit(RD.main(["--in", f"{AM}/runs/lxplus", "--out", f"{HERE}/{var}", "--samples", "data_428770"]))
