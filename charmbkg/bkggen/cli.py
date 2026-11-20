"""CLI: `charmbkg bkggen --method ...` (etap S, moduł addytywny; dispatch w charmbkg/cli.py jak `pairs`).

  charmbkg bkggen --method sideband    -c cfg [--sample data_none] [--which near|far|both]
  charmbkg bkggen --method reflection  -c cfg [--sample data_none]          # D0, szablon K↔π z danych
  charmbkg bkggen --method dstar_wrong -c cfg [--sample data_none]          # D*: zły Δm + losowy π_s (D0 + tor z innego kand.)
  charmbkg bkggen --method mix         -c cfg [--sample data_none]          # pary D0–D0: to samo zdarzenie vs mixing
  charmbkg bkggen --method generate    -c cfg [--sample data_none] [--model flow|copula] [--n N]
  charmbkg bkggen --method validate    -c cfg                               # wszystkie metody vs MC truth≥8 (incl_mc)
Wejście: próbka z YAML (`--sample`, pełne h5 przez charmbkg) albo podpróbka npz (`--npz 'sub/data_p*.npz'`, extract_sub.py).
Wyjście: <out_dir>/bkggen/<method>_<sample>.h5 (+ .json z metrykami); generator: model .pkl.
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import pickle
import sys
import time

import numpy as np

log = logging.getLogger("charmbkg")


def _load(cfg, a, species, decay_type=None, nd=None):
    from . import data as D
    if decay_type is None and a.decay_type is not None:
        decay_type = a.decay_type
        nd = nd or cfg["species"][species]["n_daughters"]
    if a.npz:
        return D.load_npz(cfg, a.npz, species, decay_type=decay_type, nd=nd, tag=a.sample)
    return D.load_h5(cfg, a.sample, species, decay_type=decay_type, nd=nd, keep_frac=a.keep_frac)


def _save(df, path, attrs):
    from ..store import Appender
    os.makedirs(os.path.dirname(path), exist_ok=True)
    keep = [c for c in df.columns if df[c].dtype != object or c in ("ctrl", "side", "sample")]
    app = Appender(path, {k: (json.dumps(v) if isinstance(v, (dict, list)) else v) for k, v in attrs.items()})
    app.add(df[keep])
    app.close()
    log.info("zapisano %s (%d wierszy)", path, len(df))


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "bkggen":
        argv = argv[1:]
    if "--method" in argv and argv[argv.index("--method") + 1] == "candrot":   # AA: S6 jako metoda bkggen
        from .candrot import main as _cr
        rest = [x for i, x in enumerate(argv) if x != "candrot" and not (x == "--method")]
        return _cr(rest)
    ap = argparse.ArgumentParser(prog="charmbkg bkggen", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--method", required=True,
                    choices=["sideband", "reflection", "dstar_wrong", "mix", "generate", "validate", "candrot"])
    ap.add_argument("-c", "--config", default=None)
    ap.add_argument("--species", default="D0", choices=["D0", "Dplus"])
    ap.add_argument("--sample", default="data_none")
    ap.add_argument("--npz", default=None, help="glob podpróbki npz zamiast h5")
    ap.add_argument("--which", default="both", choices=["near", "far", "both"])
    ap.add_argument("--model", default="flow", choices=["flow", "copula"])
    ap.add_argument("--n", type=int, default=None, help="liczba generowanych (domyślnie = oczekiwane tło w SR)")
    ap.add_argument("--n-train", type=int, default=300000)
    ap.add_argument("--epochs", type=int, default=25)
    ap.add_argument("--keep-frac", type=float, default=1.0)
    ap.add_argument("--decay-type", type=int, default=None,
                    help="inny decayType niż gatunek (np. 110 = D0 WS, 210 = rotacja w produkcji bkgctrl)")
    ap.add_argument("--veto-refl", action="store_true",
                    help="D0: odrzuć kandydatów z |m_swap−m0|<25 MeV (refleksje K↔π) — zalecane dla sidebandów/generatora")
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--h5-dir", default=None)
    a = ap.parse_args(argv)

    from .. import config as C
    from ..logs import setup_log as _setup_log
    from . import control as K
    from . import data as D
    ov = {"paths": {}}
    if a.out_dir:
        ov["paths"]["out_dir"] = a.out_dir
    if a.h5_dir:
        ov["paths"]["h5_dir"] = a.h5_dir
    cfg = C.load(a.config, ov, require_paths=True)
    od = os.path.join(cfg["paths"]["out_dir"], "bkggen")
    _setup_log(od)
    sp = cfg["species"][a.species]
    t0 = time.time()
    tag = f"{a.method}_{a.sample}_{a.species}"
    if a.method == "sideband":
        L = _load(cfg, a, a.species)
        L = L[D.in_range(L, sp)]
        if a.veto_refl and a.species == "D0":
            L = L[np.abs(D.swapped_mass(L) - sp["m0"]) > 0.025]
        out = K.sideband(L, sp, a.species, a.which)
        _save(out, os.path.join(od, f"{tag}_{a.which}.h5"), {"regions": K.REGIONS[a.species]})
    elif a.method == "reflection":
        if a.species != "D0":
            raise SystemExit("reflection: tylko D0")
        L = _load(cfg, a, "D0")
        L = L[D.in_range(L, sp)]
        R = K.reflection(L, sp)
        h, e = K.reflection_template(R, *sp["mass_range"], cfg["fit"]["nbins"])
        _save(R, os.path.join(od, f"{tag}.h5"), {"template": h.tolist(), "edges": e.tolist()})
        np.savez(os.path.join(od, f"{tag}_template.npz"), h=h, edges=e)
    elif a.method == "dstar_wrong":
        Ls = _load(cfg, a, "D0", decay_type=cfg["dstar_control"]["dstar_decay_type"], nd=3)
        W = K.dstar_wrong_dm(Ls, tuple(cfg["dstar_control"]["dm_sideband"]))
        _save(W, os.path.join(od, f"{tag}_dm.h5"), {})
        L0 = _load(cfg, a, "D0")
        Rp = K.dstar_random_pis(L0[np.abs(L0["m"] - sp["m0"]) < 0.03], sp)
        Rp = Rp[Rp["pis_pt"] > float(np.nanpercentile(Ls["d_pt2"], 1))]
        _save(Rp, os.path.join(od, f"{tag}_randpis.h5"), {})
    elif a.method == "mix":
        L = _load(cfg, a, a.species)
        L = L[D.in_range(L, sp)].reset_index(drop=True)
        I, J = K.same_event_pairs(L)
        Im, Jm = K.mixed_pairs(L, I, J, factor=2)
        _save(K.pair_frame(L, I, J, "same"), os.path.join(od, f"{tag}_same.h5"), {})
        _save(K.pair_frame(L, Im, Jm, "mixed"), os.path.join(od, f"{tag}_mixed.h5"), {})
    elif a.method == "generate":
        from . import gen as G
        from .. import features as F
        L = _load(cfg, a, a.species)
        L = L[D.in_range(L, sp)].reset_index(drop=True)
        feats = [f for f in F.usable(L.head(100000), F.feature_list(cfg, "H")) if f != "abs_impact_sig"]
        sbw = K.REGIONS[a.species]["near"][0]
        sb = np.abs(L["m"].to_numpy() - sp["m0"]) > sbw
        if a.veto_refl and a.species == "D0":
            sb &= np.abs(D.swapped_mass(L) - sp["m0"]) > 0.025
        SB = L[sb]
        if len(SB) > a.n_train:
            SB = SB.sample(a.n_train, random_state=0)
        m0 = sp["m0"]
        gm = G.BkgGenerator(feats, m0, lambda c: np.abs(np.asarray(c) - m0) > sbw, sp["mass_range"],
                            model=a.model, epochs=a.epochs).fit(SB, m_all_for_density=L["m"].to_numpy())
        w = K.REGIONS[a.species]["sr"][1]
        n = a.n or int(round(gm.mass.expected(m0 - w, m0 + w)))
        gen = gm.sample(n, (m0 - w, m0 + w))
        gen["abs_impact_sig"] = np.abs(gen["impact_sig"])
        _save(gen, os.path.join(od, f"{tag}_{a.model}.h5"), {"n_train": len(SB), "features": feats})
        with open(os.path.join(od, f"{tag}_{a.model}.pkl"), "wb") as f:
            pickle.dump(gm, f)
    elif a.method == "validate":
        raise SystemExit("validate: uruchom skrypty studiów S1/S2 (S_bkggen/work/s1_control.py, s2_gen_mc.py) — "
                         "wymagają klasyfikatora; patrz README_S.md")
    log.info("bkggen %s gotowe (%.0f s)", a.method, time.time() - t0)
    return 0


if __name__ == "__main__":
    sys.exit(main())
