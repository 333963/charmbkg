"""AFPSTREAM — wspólne: klasy trygera obejmujące strumień AFP (wyłącznie bity trigPassed), wczytanie par / all_* / ev_*.
Tylko odczyt z ../runs/lxplus/*.npz. Estymatory z ../ISO_origin/iso_common.py i ../BKG_427929/bkg_common.py (import, bez kopii)."""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
BT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(BT, "ISO_origin"))
sys.path.insert(0, os.path.join(BT, "BKG_427929"))
sys.path.insert(0, BT)
OUT = os.path.join(HERE, "out")
R = os.path.join(BT, "runs", "lxplus")
MBTS_CHAINS = ["HLT_j20_L1MBTS_2", "HLT_mb_sptrk_pt4_L1MBTS_2", "HLT_mb_sptrk_pt2_L1MBTS_2", "HLT_mb_mbts_L1MBTS_1",
               "HLT_noalg_mb_L1MBTS_1", "HLT_noalg_mb_L1MBTS_2", "HLT_mb_sptrk_L1MBTS_1"]
J12 = "HLT_noalg_L1J12"
AN = "HLT_noalg_L1AFP_A_OR_C"
AS = "HLT_mb_sptrk_pt4_L1AFP_A_OR_C"
RD0 = "HLT_noalg_L1RD0_FILLED"
ALLCH = MBTS_CHAINS + [J12, AN, AS, RD0]
# klasy rozłączne 428770 (kolejność = kolejność w tabelach)
DISJ8 = ["AFPnoalg", "AFPsptrk4", "AFPboth", "AFP&J12", "AFP&MBTS", "J12only!A", "MBTS!A", "none"]
AFP_STREAM = ["AFPnoalg", "AFPsptrk4", "AFPboth", "AFP&J12", "AFP&MBTS"]


def bits(z, prefix):
    """{chain: bool} z kolumn <prefix>tp_<chain> (== 1; −1/brak = False)."""
    return {c: (z[prefix + "tp_" + c] == 1) if (prefix + "tp_" + c) in z else None for c in ALLCH}


def classes8(tp):
    n = len(tp[J12])
    M = np.zeros(n, bool)
    for c in MBTS_CHAINS:
        M |= tp[c]
    J, an, as_ = tp[J12], tp[AN], tp[AS]
    A = an | as_
    cl = {"all": np.ones(n, bool),
          "AFPnoalg": an & ~as_ & ~J & ~M, "AFPsptrk4": as_ & ~an & ~J & ~M, "AFPboth": an & as_ & ~J & ~M,
          "AFP&J12": A & J & ~M, "AFP&MBTS": A & M, "J12only!A": J & ~M & ~A, "MBTS!A": M & ~A, "none": ~J & ~M & ~A}
    # sumy / klasy dotychczasowe / diagnostyczne (nierozłączne)
    cl["AFPonly"] = A & ~J & ~M
    cl["AFPsptrk4_any"] = as_ & ~J & ~M            # sptrk4AFP z noalg lub bez (bez J12/MBTS)
    cl["AFPstream"] = A
    cl["J12only"] = J & ~M                          # dotychczasowa (zawiera AFP&J12)
    cl["MBTSbased"] = M                             # dotychczasowa
    cl["other"] = ~J & ~M                           # dotychczasowa (= AFPonly ∪ none)
    cl["j20!J12"] = tp["HLT_j20_L1MBTS_2"] & ~J
    cl["_A"] = A; cl["_J"] = J; cl["_M"] = M; cl["_an"] = an; cl["_as"] = as_; cl["_rd0"] = tp[RD0]
    return cl


def classes9(tp, grl):
    mb = tp["HLT_mb_mbts_L1MBTS_1"] & grl
    J, j20 = tp[J12], tp["HLT_j20_L1MBTS_2"]
    A = tp[AN] | tp[AS]
    return {"MB427929": mb, "MB9&!j20&!J12": mb & ~j20 & ~J, "MB9&j20&!J12": mb & j20 & ~J,
            "MB9&AFPsptrk4": mb & tp[AS], "MB9&AFPnoalg": mb & tp[AN], "MB9&AFPsptrk4&!j20&!J12": mb & tp[AS] & ~j20 & ~J,
            "MB9&!AFP": mb & ~A, "MB9&sptrk4MBTS": mb & tp["HLT_mb_sptrk_pt4_L1MBTS_2"], "_A": A}


def atlas_label(ax, extra="", x=0.04, y=0.95, fs=11):
    t = ax.text(x, y, "ATLAS", transform=ax.transAxes, fontweight="bold", fontstyle="italic", fontsize=fs, va="top")
    ax.annotate("", xy=(1, 0), xycoords=t, xytext=(4, 0), textcoords="offset points", fontsize=fs, va="bottom")
    if extra:
        ax.text(x, y - 0.075, extra, transform=ax.transAxes, fontsize=fs - 3, va="top")


def jdump(o, path):
    def conv(x):
        if isinstance(x, (np.floating, np.integer)):
            return x.item()
        if isinstance(x, np.ndarray):
            return x.tolist()
        if isinstance(x, np.bool_):
            return bool(x)
        raise TypeError(type(x))
    json.dump(o, open(path, "w"), indent=1, default=conv, ensure_ascii=False)
