"""DNET etap 3 — pętla treningu z DODATKOWYM członem dekorelacji od „aktywności” otoczenia. Kod wspólny (DisCo, HP, modele) z
../dnet1/dn_models.py i ../dnet2/dn2_models.py — te pliki NIE są zmieniane.

fit3 = dn_models.fit (ta sama kolejność operacji i to samo zużycie globalnego generatora torch → człon masy bez zmian; przy lam_act = 0
wynik identyczny z dn_models.fit, co sprawdza `dn3_train.py selftest`) + człon
    lam_act · L_act,   L_act = średnia po zmiennych aktywności j z dCor²(sigmoid(z), act_j)   (ważone, DisCo jak dla masy),
liczony na wierszach `amask` wsadu (np. sygnał MC albo tło MC), co najwyżej HP["dec_max"] wierszy; losowanie podpróbki z WŁASNEGO
generatora (nie rusza strumienia losowego członu masy). Epoka wybierana po pełnej stracie walidacyjnej (BCE + λ·masa + λ_act·aktywność)."""
import copy
import time

import numpy as np
import torch
from torch import nn

import dn_models as DM


def fit3(net, inputs, y, w, mass, dmask, dw, val, dec, lam, seed, log, device="cpu", act=None, amask=None, aw=None, lam_act=0.0):
    """Argumenty jak dn_models.fit; act [n, k] = zmienne aktywności, amask = wiersze członu aktywności, aw = ich wagi, lam_act = λ_act."""
    HP = DM.HP
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    dev = torch.device(device)
    net = net.to(dev)
    T = [torch.as_tensor(np.asarray(x), dtype=torch.float32) for x in inputs]
    yt, wt = torch.as_tensor(y, dtype=torch.float32), torch.as_tensor(w, dtype=torch.float32)
    mt, dt, dwt = torch.as_tensor(mass, dtype=torch.float32), torch.as_tensor(dmask, dtype=torch.bool), torch.as_tensor(dw, dtype=torch.float32)
    use_act = act is not None and lam_act > 0
    if use_act:
        at = torch.as_tensor(np.asarray(act), dtype=torch.float32).reshape(len(y), -1)
        amt, awt = torch.as_tensor(amask, dtype=torch.bool), torch.as_tensor(aw, dtype=torch.float32)
    g_tr, g_va = torch.Generator(), torch.Generator()
    g_tr.manual_seed(seed + 7919)
    tr, va = np.flatnonzero(~val), np.flatnonzero(val)
    nb = max(1, int(np.ceil(len(tr) / HP["batch"])))
    opt = torch.optim.AdamW(net.parameters(), lr=HP["lr"], weight_decay=HP["wd"])
    sch = torch.optim.lr_scheduler.OneCycleLR(opt, max_lr=HP["lr"], total_steps=HP["epochs"] * nb, pct_start=0.15)
    bce = nn.BCEWithLogitsLoss(reduction="none")
    wmean = float(wt[tr].mean())

    def loss_of(idx, gen):
        idx = torch.as_tensor(idx)
        z = net(*[t[idx].to(dev) for t in T])
        lb = (bce(z, yt[idx].to(dev)) * wt[idx].to(dev)).mean() / wmean
        ld = z.new_zeros(())
        if dec != "n" and lam > 0:                                          # człon masy — 1:1 z dn_models.fit
            k = torch.nonzero(dt[idx]).reshape(-1)
            if len(k) > HP["dec_max"]:
                k = k[torch.randperm(len(k))[:HP["dec_max"]]]
            if len(k) > 16:
                ii = idx[k]
                ld = DM.dec_loss(dec, torch.sigmoid(z[k.to(dev)]), mt[ii].to(dev), dwt[ii].to(dev))
        la = z.new_zeros(())
        if use_act:                                                         # człon aktywności (nowy)
            k = torch.nonzero(amt[idx]).reshape(-1)
            if len(k) > HP["dec_max"]:
                k = k[torch.randperm(len(k), generator=gen)[:HP["dec_max"]]]
            if len(k) > 16:
                ii = idx[k]
                s, ww = torch.sigmoid(z[k.to(dev)]), awt[ii].to(dev)
                la = torch.stack([DM.disco(s, at[ii, j].to(dev), ww) for j in range(at.shape[1])]).mean()
        return lb + lam * ld + lam_act * la, lb, ld, la
    best, best_w, best_ep, hist, t0 = float("inf"), None, -1, [], time.time()
    for ep in range(HP["epochs"]):
        net.train()
        order = rng.permutation(tr)
        for i in range(nb):
            ls = loss_of(order[i * HP["batch"]:(i + 1) * HP["batch"]], g_tr)[0]
            opt.zero_grad(set_to_none=True); ls.backward(); opt.step(); sch.step()
        net.eval()
        with torch.no_grad():
            tot = np.zeros(4); n = 0
            torch.manual_seed(seed + 1000)                                  # ta sama podpróbka dekorelacji w każdej epoce
            g_va.manual_seed(seed + 1000)
            for i in range(0, len(va), 4 * HP["batch"]):
                idx = va[i:i + 4 * HP["batch"]]
                tot += np.array([float(v) for v in loss_of(idx, g_va)]) * len(idx); n += len(idx)
            tot /= max(n, 1)
            torch.manual_seed(seed + ep + 1)
        hist.append([float(x) for x in tot])
        if tot[0] < best:
            best, best_w, best_ep = tot[0], copy.deepcopy(net.state_dict()), ep
    net.load_state_dict(best_w)
    net = net.to("cpu").eval()
    log(f"    epoka {best_ep}/{HP['epochs']}: val {best:.4f} (BCE {hist[best_ep][1]:.4f}, masa {hist[best_ep][2]:.5f}, aktywność {hist[best_ep][3]:.5f}), "
        f"{time.time() - t0:.0f}s")
    return net, {"best_epoch": best_ep, "hist": hist, "sec": round(time.time() - t0, 1), "hist_cols": ["total", "bce", "dcor2_mass", "dcor2_act"]}
