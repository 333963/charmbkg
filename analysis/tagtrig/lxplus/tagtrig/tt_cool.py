#!/usr/bin/env python
"""TAGTRIG / węzeł lumi: klucze preskal per zakres LB z COOL i wartości preskal L1/HLT z TriggerDB — TYLKO ODCZYT, best-effort.
Uruchamiany w środowisku Athena (tt_lumi.sh). NIETESTOWANE lokalnie (API TrigConfStorage / TrigConfIO z pamięci) — każdy krok w try/except,
linie „COOL: …” trafiają do lumi_STATUS.txt; surowe struktury → --raw (json), żeby dało się je odczytać nawet przy zmianie API.
Wyjście --out (csv): run, level (L1 | HLT), key, lb_first, lb_last, name, prescale, enabled.
"""
import argparse, csv, json, sys, traceback


def note(*a):
    print("COOL:", *a, flush=True)


def ranges(d, run):
    """Z odpowiedzi TriggerCoolUtil ({run: {nazwa: [(klucz, lb_od, lb_do), …] | wartość}}) → {nazwa: [(klucz, od, do)]}."""
    out = {}
    r = d.get(run, d.get(str(run), {})) if isinstance(d, dict) else {}
    for k, v in (r.items() if isinstance(r, dict) else []):
        if isinstance(v, (list, tuple)) and v and all(isinstance(x, (list, tuple)) and len(x) >= 3 for x in v):
            out[k] = [(int(x[0]), int(x[1]), int(x[2])) for x in v]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", type=int, required=True); ap.add_argument("--out", required=True); ap.add_argument("--raw", required=True)
    ap.add_argument("--hlt", nargs="*", default=[]); ap.add_argument("--l1", nargs="*", default=[])
    a = ap.parse_args()
    raw, rows, dbalias = {}, [], "TRIGGERDB_RUN3"
    try:
        from TrigConfStorage.TriggerCoolUtil import TriggerCoolUtil
        db = TriggerCoolUtil.GetConnection("CONDBR2")
        rr = [[a.run, a.run]]
        for nm in ("getHLTConfigKeys", "getHLTPrescaleKeys", "getL1ConfigKeys", "getBunchGroupKey"):
            try:
                raw[nm] = getattr(TriggerCoolUtil, nm)(db, rr)
                note(f"{nm}: {str(raw[nm])[:300]}")
            except Exception as e:
                note(f"{nm}: BŁĄD {type(e).__name__}: {e}")
        cfg = raw.get("getHLTConfigKeys", {})
        c = cfg.get(a.run, {}) if isinstance(cfg, dict) else {}
        if isinstance(c, dict) and isinstance(c.get("DB"), str) and c["DB"].startswith("TRIGGERDB"):
            dbalias = c["DB"].split(";")[0].split()[0]
    except Exception as e:
        note(f"TriggerCoolUtil: BŁĄD {type(e).__name__}: {e}")
        traceback.print_exc()
    json.dump(raw, open(a.raw, "w"), indent=1, default=str)
    hr = ranges(raw.get("getHLTPrescaleKeys", {}), a.run)
    lr = ranges(raw.get("getL1ConfigKeys", {}), a.run)
    hk = next((v for k, v in hr.items() if "PSK" in k.upper()), [])
    lk = next((v for k, v in lr.items() if "PSK" in k.upper()), [])
    note(f"zakresy kluczy: HLTPSK {hk[:6]}{'…' if len(hk) > 6 else ''}; L1PSK {lk[:6]}{'…' if len(lk) > 6 else ''}; dbalias {dbalias}")
    cache = {}

    def l1ps(key):
        if ("L1", key) not in cache:
            from TrigConfIO.L1TriggerConfigAccess import L1PrescalesSetAccess
            cache[("L1", key)] = L1PrescalesSetAccess(dbalias=dbalias, l1pskey=key)
        return cache[("L1", key)]

    def hltps(key):
        if ("HLT", key) not in cache:
            from TrigConfIO.HLTTriggerConfigAccess import HLTPrescalesSetAccess
            cache[("HLT", key)] = HLTPrescalesSetAccess(dbalias=dbalias, hltpskey=key)
        return cache[("HLT", key)]

    def val(acc, level, name):
        try:
            return float(acc.prescale(name)), int(bool(acc.enabled(name)))
        except Exception:
            pass
        if level == "L1":       # zapas: surowy „cut” (preskala = 0xFFFFFF / (0x1000000 − cut))
            e = acc["cutValues"][name]
            return 0xFFFFFF / float(0x1000000 - int(e["cut"])), int(bool(e.get("enabled", True)))
        e = acc["prescales"][name]
        return float(e["prescale"]), int(bool(e.get("enabled", True)))

    for level, keys, names, get in (("L1", lk, a.l1, l1ps), ("HLT", hk, a.hlt, hltps)):
        bad = 0
        for key, lb0, lb1 in keys:
            try:
                acc = get(key)
            except Exception as e:
                note(f"{level} PSK {key}: BŁĄD odczytu TriggerDB {type(e).__name__}: {e}"); bad += 1
                continue
            for nme in names:
                try:
                    p, en = val(acc, level, nme)
                    rows.append([a.run, level, key, lb0, lb1, nme, p, en])
                except Exception as e:
                    rows.append([a.run, level, key, lb0, lb1, nme, "", ""]); bad += 1
                    if bad <= 5:
                        note(f"{level} PSK {key} {nme}: BŁĄD {type(e).__name__}: {e}")
        note(f"{level}: {len(keys)} zakresów kluczy, błędów {bad}")
    if rows:
        with open(a.out, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["run", "level", "key", "lb_first", "lb_last", "name", "prescale", "enabled"]); w.writerows(rows)
    note(f"wierszy preskal: {len(rows)}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as e:      # niekrytyczne
        note(f"BŁĄD ogólny {type(e).__name__}: {e}")
        traceback.print_exc()
        sys.exit(0)
