"""PEAK_pos: kompaktowy wyciąg D0←D* (HYP_d0/runs/lxplus_dstar) — tylko potrzebne kolumny."""
import os
import pyarrow.parquet as pq
SRC = "/home/k/Documents/charm/charm_work/HYP_d0/runs/lxplus_dstar"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
COLS = ["m", "pt", "eta", "dm", "score", "soft_pi_charge", "hypo_from_tag", "hypo_swap", "truth", "dstar_conflict", "dstar_tag",
        "pis_pt", "weight", "d_charge0", "evblk", "file"]
for s in ["data_428770", "fd_nd"]:
    t = pq.read_table(os.path.join(SRC, f"{s}.parquet"), columns=COLS)
    t.to_pandas().to_parquet(os.path.join(OUT, f"dstar_{s}.parquet"))
    print(s, t.num_rows)
