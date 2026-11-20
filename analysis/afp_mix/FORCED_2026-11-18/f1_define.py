"""FORCED f1 — definicje wykluczenia zdarzeń „w drzewie tylko przez łańcuch AFP” (428770) i liczniki. Tylko bity i liczniki.
Wyjście: out/excl_428770.parquet (event_number + maski wariantów), out/f1.log, out/t1_liczniki.csv, out/t2_ptag_nselbl.csv."""
import os, sys, json, numpy as np, pandas as pd
W = os.path.expanduser("~/Documents/charm/charm_work"); sys.path.insert(0, f"{W}/TAGREG_2026-11-18/v2")
import v2common as v2
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = f"{HERE}/out"
RUN = int(sys.argv[1]) if len(sys.argv) > 1 else 428770
ev = pd.read_parquet(f"{W}/AFP_mix/runs/lxplus/data_{RUN}_events.parquet")
n = len(ev); en = ev.event_number.to_numpy()
print("zdarzeń", n, "unikalnych event_number", len(np.unique(en)))
assert len(np.unique(en)) == n
tc = [c for c in ev.columns if c.startswith("trigPassed_")]
P = {c[11:]: ev[c].to_numpy() == 1 for c in tc}
nonafp = np.any([P[c] for c in P if "L1AFP" not in c], axis=0); afp11 = np.any([P[c] for c in P if "L1AFP" in c], axis=0)
grl = ev.passGRL_MinBiasAFP.to_numpy() == 1; npu0 = ev.nPileupVertices.to_numpy() == 0
sel = grl & nonafp & npu0                                   # próbka mieszania (dawcy)
S13 = P["HLT_noalg_L1J12"] | P["HLT_j20_L1MBTS_2"] | P["HLT_noalg_L1RD0_FILLED"]   # łańcuchy nie-AFP z S wśród 11 zapisanych
tbpafp = ev.trigTBP_L1_AFP_A_OR_C.to_numpy() == 1
tagAC = (ev.afp_tag_A.to_numpy() > 0) | (ev.afp_tag_C.to_numpy() > 0)
anyp = (ev.afp_nProtons_A.to_numpy() > 0) | (ev.afp_nProtons_C.to_numpy() > 0)
ns = ev.nSelBL.to_numpy(); hasD = ev.n_cand_sel.to_numpy() > 0 if "n_cand_sel" in ev else np.zeros(n, bool)
# pełne bity (tylko zdarzenia tagowane)
b = v2.load(RUN); H, hlt, bev = b["H"], b["hlt"], b["ev"]
ben = bev.event_number.to_numpy()
print("bits: wierszy", len(ben), "unikalnych", len(np.unique(ben)))
assert len(np.unique(ben)) == len(ben)
idx = pd.Index(en).get_indexer(ben); assert (idx >= 0).all(), "zdarzenie z bits bez wiersza w events"
# zgodność klucza: LB, nSelBL, npu muszą się zgadzać
assert (ev.lumiBlock.to_numpy()[idx] == bev.lumiBlock.to_numpy()).all() and (ns[idx] == bev.nSelBL.to_numpy()).all() and (ev.nPileupVertices.to_numpy()[idx] == bev.npu.to_numpy()).all()
for c in ("HLT_noalg_L1J12", "HLT_j20_L1MBTS_2", "HLT_mb_sptrk_pt4_L1MBTS_2", "HLT_noalg_L1AFP_A_OR_C"):
    assert (P[c][idx] == H[:, hlt.index(c)]).all(), c
print("klucz event_number: jednoznaczny w obu tabelach; LB, nSelBL, npu i 4 bity kontrolne zgodne dla wszystkich", len(ben))
inbits = np.zeros(n, bool); inbits[idx] = True
print("tag (afp_tag_A|C) w events:", int(tagAC.sum()), " w bits:", int(inbits.sum()), " tag ∧ ¬bits:", int((tagAC & ~inbits).sum()), " bits ∧ ¬tag:", int((inbits & ~tagAC).sum()))
print("sel:", int(sel.sum()), " sel ∧ bits:", int((sel & inbits).sum()), " sel ∧ ≥1 proton ∧ ¬bits:", int((sel & anyp & ~inbits).sum()))
Sf = f"{W}/TAGREG_2026-11-18/v2/out/strumien_S_{RUN}.csv"; S = pd.read_csv(Sf); Sset = set(S.lancuch[S.w_S])
inS = np.array([c in Sset for c in hlt]); isafp = np.array(["AFP" in c for c in hlt])
print("S nie-AFP:", [c for c, a, f in zip(hlt, inS, isafp) if a and not f]); print("S AFP:", [c for c, a, f in zip(hlt, inS, isafp) if a and f])
Snon = np.zeros(n, bool); Safp = np.zeros(n, bool); Hany = np.zeros(n, bool)
Snon[idx] = H[:, inS & ~isafp].any(1); Safp[idx] = H[:, inS & isafp].any(1); Hany[idx] = H.any(1)
# --- definicje
forced_tag = inbits & ~Snon & Safp                       # tagowane: pełne bity — brak łańcucha nie-AFP z S, jest łańcuch AFP z S
noS_tag = inbits & ~Snon & ~Safp                         # tagowane bez żadnego łańcucha S (nierozpoznane)
unt_cand = ~inbits & ~S13 & tbpafp                       # nietagowane: brak bitu S13, L1_AFP przed preskalą = 1 (mogły wejść łańcuchem AFP)
V = {"W1a": forced_tag,                                  # tylko tagowane (definicja przeglądu), nietagowane bez zmian
     "W1": forced_tag | unt_cand,                        # + nietagowane-kandydaci (aproksymacja TBP)
     "W2": ~S13}                                         # symetrycznie: wymóg bitu łańcucha nie-AFP z S wśród 11 (J12 ∨ j20_L1MBTS_2), tagowane i nietagowane
rows = []
def cnt(name, m):
    r = dict(maska=name, N=int(m.sum()), w_sel=int((m & sel).sum()), w_sel_tagAC=int((m & sel & inbits).sum()), w_sel_zD=int((m & sel & hasD).sum()),
             w_sel_zD_tagAC=int((m & sel & hasD & inbits).sum()), w_sel_nietag=int((m & sel & ~inbits).sum()), w_sel_nSel_lt5=int((m & sel & (ns < 5)).sum()))
    rows.append(r); return r
cnt("sel (nominał: GRL ∧ nie-AFP(11) ∧ npu0)", np.ones(n, bool))
cnt("tagowane, wymuszone (pełne bity: ¬S_nieAFP ∧ S_AFP)", forced_tag); cnt("tagowane, bez żadnego łańcucha S", noS_tag)
cnt("tagowane, ¬S13 ale z innym łańcuchem nie-AFP z S (j10f, EM, MU…)", inbits & ~S13 & Snon)
cnt("nietagowane, ¬S13", ~inbits & ~S13); cnt("nietagowane, ¬S13 ∧ TBP L1_AFP", unt_cand); cnt("nietagowane, ¬S13 ∧ TBP L1_AFP ∧ ≥1 proton AFP", unt_cand & anyp)
for k, m in V.items(): cnt(f"wykluczone {k}", m)
t = pd.DataFrame(rows); t.to_csv(f"{OUT}/t1_liczniki_{RUN}.csv", index=False)
pd.set_option("display.width", 250, "display.max_colwidth", 80); print(t.to_string())
# W2: ile tagowanych z prawdziwym łańcuchem nie-AFP z S traci (nadmiarowe wykluczenie) i czy W2 ⊇ wymuszone
print("wymuszone ⊂ W2:", bool((forced_tag & sel & ~V["W2"]).sum() == 0), " wymuszone∧sel∧S13:", int((forced_tag & sel & S13).sum()))
# --- częstość tagu vs nSelBL (tag = afp_tag_A|C, jak w przeglądzie)
edges = [0, 5, 10, 15, 20, 30, 40, 60, 10000]; R = []
for k, m in [("nominal", np.zeros(n, bool))] + list(V.items()):
    s = sel & ~m; row = {"wariant": k}
    for lo, hi in zip(edges[:-1], edges[1:]):
        a = s & (ns >= lo) & (ns < hi); row[f"N_{lo}-{hi}"] = int(a.sum()); row[f"Ntag_{lo}-{hi}"] = int((a & inbits).sum()); row[f"P_{lo}-{hi}"] = (a & inbits).sum() / max(a.sum(), 1)
    row["wzrost_lt5_vs_40-60"] = row["P_0-5"] / row["P_40-60"]; row["wzrost_5-10_vs_40-60"] = row["P_5-10"] / row["P_40-60"]
    k5, k40 = row["Ntag_0-5"], row["Ntag_40-60"]; row["wzrost_err"] = row["wzrost_lt5_vs_40-60"] * np.sqrt(1 / max(k5, 1) + 1 / max(k40, 1))
    R.append(row)
t2 = pd.DataFrame(R); t2.to_csv(f"{OUT}/t2_ptag_nselbl_{RUN}.csv", index=False)
print(t2[["wariant"] + [c for c in t2.columns if c.startswith(("Ntag_", "wzrost"))]].to_string())
print(t2[["wariant"] + [c for c in t2.columns if c.startswith("P_")]].to_string())
if RUN == 428770:
    pd.DataFrame({"event_number": en, **{k: m for k, m in V.items()}}).loc[np.any(list(V.values()), axis=0)].to_parquet(f"{OUT}/excl_{RUN}.parquet", index=False)
