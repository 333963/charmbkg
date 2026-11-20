"""charmbkg — wejście konsolowe projektu B (tło podobne do sygnału).

  charmbkg bkggen  --method sideband|reflection|dstar_wrong|mix|generate|validate|candrot ...   (= charmbkg-bkggen ...)
  charmbkg s8salad -c CFG --species Dplus ...      SALAD: reweighting MC tła → dane w SB, kształt tła pod pikiem
  charmbkg hyp-extract --sample ... --h5 ... --out ...   hipotezy D0 K↔π (arbiter PID) → parquet
  charmbkg synth --h5-dir DIR [--max-rows N]       syntetyczne h5 v1.1 (tylko testy/przykłady, NIE fizyka)
"""
import argparse
import sys


def _bkggen(argv):
    from .bkggen.cli import main
    return main(argv)


def _s8salad(argv):
    from .s8salad import main
    return main(argv)


def _hyp_extract(argv):
    from .hyp_extract import main
    return main(argv[1:])


def _synth(argv):
    ap = argparse.ArgumentParser(prog="charmbkg synth")
    ap.add_argument("synth")
    ap.add_argument("--h5-dir", required=True)
    ap.add_argument("--out-dir", default=None, help="ignorowane (zgodność z `charmtag synth`)")
    ap.add_argument("--max-rows", type=int, default=20000)
    a = ap.parse_args(argv)
    from .synth import make
    print(make(a.h5_dir, n=a.max_rows))
    return 0


_SUB = {"bkggen": _bkggen, "s8salad": _s8salad, "hyp-extract": _hyp_extract, "synth": _synth}


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] not in _SUB:
        print(__doc__)
        return 0 if argv and argv[0] in ("-h", "--help") else 2
    return _SUB[argv[0]](argv)


def bkggen_main(argv=None):
    """`charmbkg-bkggen --method ...` (= `charmbkg bkggen ...`)."""
    return _bkggen(list(sys.argv[1:] if argv is None else argv))
