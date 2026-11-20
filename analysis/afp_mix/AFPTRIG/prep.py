import lib, sys, json
for r in lib.RUNS:
    d = lib.load(r, force=True); print(r, len(d['lb']), json.loads(str(d['info'])))
