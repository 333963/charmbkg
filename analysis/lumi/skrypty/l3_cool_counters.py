"""l3: odczyt COOL (Frontier, tylko odczyt): menu L1 (CTP id), liczniki L1 per LB (TBP/TAP/TAV), czasy LB.
Użycie: python l3_cool_counters.py <run> <item1,item2,...>  -> l1counters_<run>.csv, lblb_<run>.csv"""
import sys, csv
from PyCool import cool
from CoolConvUtilities.AtlCoolLib import indirectOpen
run = int(sys.argv[1]); items = sys.argv[2].split(",")
def P(*a): print(*a, flush=True)
db = indirectOpen("COOLONL_TRIGGER/CONDBR2", True, True); P("otwarte COOLONL_TRIGGER/CONDBR2")
since, until = (run << 32), (run << 32) + 0xFFFFFFFF
# 1) menu L1 w COOL
ids = {}
try:
    f = db.getFolder("/TRIGGER/LVL1/Menu"); it = f.browseObjects(since, until, cool.ChannelSelection.all()); n = 0
    names = {}
    while it.goToNext():
        o = it.currentRef(); names[o.channelId()] = str(o.payload()["ItemName"]); n += 1
    it.close(); P("COOL /TRIGGER/LVL1/Menu: obiektów", n, "; przykłady", list(names.items())[:8])
    ids = {v: k for k, v in names.items() if v in items}
    P("znalezione w COOL:", ids)
except Exception as e: P("ERR Menu", repr(e))
# 2) menu L1 z TriggerDB (Run 3)
if len(ids) < len(items):
    try:
        from TrigConfIO.L1TriggerConfigAccess import L1MenuAccess
        smk = int(sys.argv[3]) if len(sys.argv) > 3 else 3094
        m = L1MenuAccess(dbalias="TRIGGERDB_RUN3", smkey=smk)
        allit = m.items()
        for nm in items:
            if nm in allit: ids[nm] = int(allit[nm]["ctpid"])
        P("TriggerDB SMK", smk, "ctpid:", ids)
    except Exception as e: P("ERR TriggerDB", repr(e))
# 3) czasy LB
rows = []
f = db.getFolder("/TRIGGER/LUMI/LBLB"); it = f.browseObjects(since, until, cool.ChannelSelection(0))
while it.goToNext():
    o = it.currentRef(); p = o.payload(); rows.append((o.since() & 0xFFFFFFFF, int(p["StartTime"]), int(p["EndTime"])))
it.close(); P("LBLB wierszy", len(rows))
with open(f"lblb_{run}.csv", "w") as fo:
    w = csv.writer(fo); w.writerow(["LB", "start_ns", "end_ns"]); w.writerows(rows)
# 4) liczniki
out = []
for fn in ("/TRIGGER/LUMI/LVL1COUNTERS",):
    try:
        f = db.getFolder(fn); P(fn, "payload:", [str(s.name()) for s in f.payloadSpecification()] if hasattr(f.payloadSpecification(), "__iter__") else f.payloadSpecification().size())
        for nm, cid in ids.items():
            it = f.browseObjects(since, until, cool.ChannelSelection(cid)); n = 0
            while it.goToNext():
                o = it.currentRef(); p = o.payload()
                out.append((run, nm, cid, o.since() & 0xFFFFFFFF, int(p["BeforePrescale"]), int(p["AfterPrescale"]), int(p["L1Accept"]))); n += 1
            it.close(); P(" ", nm, "ctpid", cid, "wierszy", n)
    except Exception as e: P("ERR", fn, repr(e))
with open(f"l1counters_{run}.csv", "w") as fo:
    w = csv.writer(fo); w.writerow(["run", "item", "ctpid", "LB", "TBP", "TAP", "TAV"]); w.writerows(out)
P("koniec", len(out))
