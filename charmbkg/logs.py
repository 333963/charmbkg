"""Logowanie CLI (kopia `charmtag.cli._setup_log` z charmtag 2.0.3; plik logu `charmbkg.log` w out_dir)."""
import logging
import os
import sys

log = logging.getLogger("charmbkg")


def setup_log(out_dir):
    os.makedirs(out_dir, exist_ok=True)
    fmt = logging.Formatter("%(asctime)s %(levelname)s %(message)s", "%H:%M:%S")
    log.setLevel(logging.INFO)
    if not log.handlers:
        h = logging.StreamHandler(sys.stdout)
        h.setFormatter(fmt)
        log.addHandler(h)
        fh = logging.FileHandler(os.path.join(out_dir, "charmbkg.log"))
        fh.setFormatter(fmt)
        log.addHandler(fh)
