"""Ścieżki wyjść (kopia helperów `charmtag.pipeline.mdir/sdir` z charmtag 2.0.3, ten sam układ katalogów)."""
import os


def mdir(cfg, species, variant):  # S5: paths.models_dir (opcjonalnie: modele z innego biegu, np. nominalnego)
    return os.path.join(cfg["paths"].get("models_dir") or os.path.join(cfg["paths"]["out_dir"], "models"), f"{species}_{variant}")


def sdir(cfg, species, variant):  # S5: paths.scores_dir (opcjonalnie: wyniki z innego biegu)
    return os.path.join(cfg["paths"].get("scores_dir") or os.path.join(cfg["paths"]["out_dir"], "scores"), f"{species}_{variant}")
