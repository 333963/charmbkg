"""PREZENTACJA_2026-11-20 — wspolny styl rysunkow (projekt B)."""
import json
import os
import sys

sys.dont_write_bytecode = True
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
CW = os.path.dirname(BASE)
FIG = os.path.join(BASE, "fig")
TAGDIR = os.path.join(BASE, "local_tagged")
J = lambda *p: json.load(open(os.path.join(CW, *p)))  # noqa: E731

# paleta kategoryczna (zweryfikowana walidatorem CVD), stala kolejnosc
C = {"blue": "#2a78d6", "orange": "#eb6834", "aqua": "#1baf7a", "yellow": "#eda100", "magenta": "#e87ba4", "green": "#008300", "violet": "#6250d6", "red": "#e34948"}
INK, INK2, GRID = "#0b0b0b", "#52514e", "#d9d8d2"
DATA = INK                                   # dane zawsze czarne
SEQ = "Blues"
LABEL = ""
BINS = ["3-5", "5-8", "8-30"]
BINLAB = {"3-5": "3–5", "5-8": "5–8", "8-30": "8–30"}
EDGES = {"3-5": (3, 5), "5-8": (5, 8), "8-30": (8, 30)}

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 11, "axes.labelsize": 12, "axes.titlesize": 12, "legend.fontsize": 9.5, "xtick.labelsize": 10.5, "ytick.labelsize": 10.5,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2, "text.color": INK, "axes.linewidth": 0.8,
    "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True, "legend.frameon": False,
    "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white", "lines.linewidth": 1.8, "errorbar.capsize": 0,
    "axes.grid": False, "grid.color": GRID, "grid.linewidth": 0.6, "pdf.fonttype": 42, "mathtext.default": "regular"})


def lumi_run(run="427929"):
    """(L [nb^-1], blad) z LUMI_2026-11-19/lumi.json dla trygera MBTS (427929: nasza kopia plikow)."""
    lj = J("LUMI_2026-11-19", "lumi.json")
    return lj


def header(fig, run_text, y=0.985, x=0.01, extra=None):
    """Etykieta na gorze rysunku: status + energia, run, L."""
    fig.text(x, y, LABEL, ha="left", va="top", fontsize=10.5, fontweight="bold", color=INK)
    t = r"$\sqrt{s}$ = 13.6 TeV, " + run_text
    if extra:
        t += ", " + extra
    fig.text(x, y - 0.03 * 6.0 / fig.get_figheight(), t, ha="left", va="top", fontsize=10, color=INK2)


def save(fig, name, tagged=False):
    d = TAGDIR if tagged else FIG
    fig.savefig(os.path.join(d, name + ".png"), dpi=200)
    fig.savefig(os.path.join(d, name + ".pdf"))
    plt.close(fig)
    print("zapisano", os.path.join(d, name + ".png"))


def L427():
    """L = 0,4648 nb^-1 z final.json (L_nb, L_rel) — run 427929, nasza kopia plikow."""
    R = J("XSEC_final_2026-11-20", "v2", "out", "final.json")
    return R["L_nb"], R["L_nb"] * R["L_rel"]


def run427_text():
    L, e = L427()
    return f"run 427929 (minimum bias), L = {L:.3f} ± {e:.3f} nb$^{{-1}}$"
