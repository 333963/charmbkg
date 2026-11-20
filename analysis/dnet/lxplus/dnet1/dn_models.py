"""DNET — modele i straty (torch): MLP tabelaryczny z przetwarzaniem wstępnym W GRAFIE (→ ONNX bez zewnętrznego skalera),
dekorelacja od masy: DisCo (ważone dCor², Kasieczka–Shih 2001.05310) i MoDe[n] (Kitouni i in. 2010.09745; tu wariant z gładkim
indykatorem), pętla treningu wspólna dla MLP (etap 1) i modeli na zbiorach torów (etap 2: dnet2/dn2_models.py)."""
import copy
import time

import numpy as np
import torch
from torch import nn

HP = {"lr": 2e-3, "wd": 1e-5, "epochs": 30, "batch": 4096, "dec_max": 1024, "val_frac": 0.1, "mode_bins": 8, "mode_thr": 24, "mode_tau": 0.03}


# ------------------------------------------------------------------------------------------------------------- dekorelacja
def disco(s, m, w=None, eps=1e-9):
    """Ważone dCor²(s, m) ∈ [0, 1] (Székely); w — wagi (normalizowane do 1)."""
    s, m = s.reshape(-1, 1).float(), m.reshape(-1, 1).float()
    n = s.shape[0]
    if n < 4:
        return s.new_zeros(())
    w = torch.full((n,), 1.0 / n, device=s.device) if w is None else (w.reshape(-1).float() / w.sum())
    a, b = torch.cdist(s, s, p=1), torch.cdist(m, m, p=1)

    def ctr(d):
        r = (d * w[None, :]).sum(1, keepdim=True)
        return d - r - r.t() + (r.reshape(-1) * w).sum()
    A, B = ctr(a), ctr(b)
    ww = w[:, None] * w[None, :]
    ab, aa, bb = (A * B * ww).sum(), (A * A * ww).sum(), (B * B * ww).sum()
    return ab / torch.sqrt(torch.clamp(aa * bb, min=eps))


def mode(s, m, w=None, order=1, nb=8, nthr=24, tau=0.03):
    """MoDe[order] z gładkim indykatorem: F_j(t) = dystrybuanta score w j-tym binie masy (równoliczne biny); kara = ważona średnia
    kwadratu reszty F_j(t) po odjęciu dopasowania wielomianem stopnia ≤ order w masie (order 0 → płasko, 1 → liniowo dozwolone)."""
    s, m = s.reshape(-1).float(), m.reshape(-1).float()
    n = s.shape[0]
    if n < 8 * nb:
        return s.new_zeros(())
    w = torch.ones(n, device=s.device) if w is None else w.reshape(-1).float()
    q = torch.quantile(m, torch.linspace(0, 1, nb + 1, device=s.device))
    j = torch.clamp(torch.bucketize(m, q[1:-1].contiguous()), 0, nb - 1)
    oh = torch.nn.functional.one_hot(j, nb).float() * w[:, None]                    # [n, nb]
    Wj = oh.sum(0) + 1e-9
    t = torch.quantile(s.detach(), torch.linspace(0.02, 0.98, nthr, device=s.device))
    ind = torch.sigmoid((t[None, :] - s[:, None]) / tau)                           # [n, nthr]
    Fj = (oh.t() @ ind) / Wj[:, None]                                              # [nb, nthr]
    mc = (oh * m[:, None]).sum(0) / Wj
    x = 2 * (mc - mc.min()) / torch.clamp(mc.max() - mc.min(), min=1e-9) - 1
    P = torch.stack([x ** k for k in range(order + 1)], 1)                         # [nb, order+1]
    PW = P * Wj[:, None]
    beta = torch.linalg.solve(P.t() @ PW + 1e-6 * torch.eye(order + 1, device=s.device), PW.t() @ Fj)
    R = Fj - P @ beta
    return ((R ** 2).mean(1) * Wj).sum() / Wj.sum()


def dec_loss(kind, s, m, w):
    if kind == "d":
        return disco(s, m, w)
    if kind == "m":
        return mode(s, m, w, 1, HP["mode_bins"], HP["mode_thr"], HP["mode_tau"])
    raise ValueError(kind)


# ------------------------------------------------------------------------------------------------------------- MLP
class Prep(nn.Module):
    """NaN/±inf → mediana; (x − mediana)/IQR; asinh. Stałe = bufory modelu (liczone z próbki treningowej)."""

    def __init__(self, med, scale):
        super().__init__()
        self.register_buffer("med", torch.as_tensor(med, dtype=torch.float32))
        self.register_buffer("scale", torch.as_tensor(scale, dtype=torch.float32))

    def forward(self, x):
        bad = torch.isnan(x) | torch.isinf(x)
        x = torch.where(bad, self.med.expand_as(x), x)
        z = (x - self.med) / self.scale
        return torch.sign(z) * torch.log(torch.abs(z) + torch.sqrt(z * z + 1.0))        # asinh z operatorów podstawowych (eksporter ONNX nie ma aten::asinh)


def prep_stats(X):
    X = np.where(np.isfinite(X), X, np.nan)
    med = np.nanmedian(X, axis=0)
    q1, q3 = np.nanpercentile(X, 25, axis=0), np.nanpercentile(X, 75, axis=0)
    sc = (q3 - q1) / 1.349
    sd = np.nanstd(X, axis=0)
    sc = np.where(sc > 0, sc, np.where(sd > 0, sd, 1.0))
    return np.where(np.isfinite(med), med, 0.0).astype(np.float32), sc.astype(np.float32)


class ResBlock(nn.Module):
    def __init__(self, h):
        super().__init__()
        self.n, self.l = nn.LayerNorm(h), nn.Linear(h, h)

    def forward(self, x):
        return x + self.l(torch.nn.functional.silu(self.n(x)))


class MLP(nn.Module):
    """arch m64: 3×64 ReLU (przepis CharmNet z gałęzi minimum-bias charmpy); r128: Linear→3 bloki rezydualne (LayerNorm, SiLU)→Linear. Wyjście: logit."""

    def __init__(self, nin, arch, med, scale):
        super().__init__()
        self.prep = Prep(med, scale)
        if arch == "m64":
            self.body = nn.Sequential(nn.Linear(nin, 64), nn.ReLU(), nn.Linear(64, 64), nn.ReLU(), nn.Linear(64, 64), nn.ReLU(), nn.Linear(64, 1))
        elif arch == "r128":
            self.body = nn.Sequential(nn.Linear(nin, 128), ResBlock(128), ResBlock(128), ResBlock(128), nn.LayerNorm(128), nn.SiLU(), nn.Linear(128, 1))
        else:
            raise ValueError(arch)

    def forward(self, x):
        return self.body(self.prep(x)).squeeze(-1)


class Sig(nn.Module):
    """Opakowanie do eksportu: wyjście = sigmoid(logit) ∈ (0, 1), jak CharmNet w charmpp."""

    def __init__(self, net):
        super().__init__()
        self.net = net

    def forward(self, *x):
        return torch.sigmoid(self.net(*x))


def nparams(net):
    return int(sum(p.numel() for p in net.parameters()))


# ------------------------------------------------------------------------------------------------------------- trening
def fit(net, inputs, y, w, mass, dmask, dw, val, dec, lam, seed, log, device="cpu"):
    """inputs: krotka tablic numpy (pierwszy wymiar = wiersze); y bool; w = wagi BCE (0 → wiersz tylko do dekorelacji);
    dmask/dw = wiersze i wagi członu dekorelacji; val = maska walidacji (wybór epoki). dec ∈ {n, d, m}, lam = λ."""
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    dev = torch.device(device)
    net = net.to(dev)
    T = [torch.as_tensor(np.asarray(x), dtype=torch.float32) for x in inputs]
    yt, wt = torch.as_tensor(y, dtype=torch.float32), torch.as_tensor(w, dtype=torch.float32)
    mt, dt, dwt = torch.as_tensor(mass, dtype=torch.float32), torch.as_tensor(dmask, dtype=torch.bool), torch.as_tensor(dw, dtype=torch.float32)
    tr, va = np.flatnonzero(~val), np.flatnonzero(val)
    nb = max(1, int(np.ceil(len(tr) / HP["batch"])))
    opt = torch.optim.AdamW(net.parameters(), lr=HP["lr"], weight_decay=HP["wd"])
    sch = torch.optim.lr_scheduler.OneCycleLR(opt, max_lr=HP["lr"], total_steps=HP["epochs"] * nb, pct_start=0.15)
    bce = nn.BCEWithLogitsLoss(reduction="none")
    wmean = float(wt[tr].mean())

    def loss_of(idx):
        idx = torch.as_tensor(idx)
        z = net(*[t[idx].to(dev) for t in T])
        lb = (bce(z, yt[idx].to(dev)) * wt[idx].to(dev)).mean() / wmean
        ld = z.new_zeros(())
        if dec != "n" and lam > 0:
            k = torch.nonzero(dt[idx]).reshape(-1)
            if len(k) > HP["dec_max"]:
                k = k[torch.randperm(len(k))[:HP["dec_max"]]]
            if len(k) > 16:
                ii = idx[k]
                ld = dec_loss(dec, torch.sigmoid(z[k.to(dev)]), mt[ii].to(dev), dwt[ii].to(dev))
        return lb + lam * ld, lb, ld
    best, best_w, best_ep, hist, t0 = float("inf"), None, -1, [], time.time()
    for ep in range(HP["epochs"]):
        net.train()
        order = rng.permutation(tr)
        for i in range(nb):
            ls, _, _ = loss_of(order[i * HP["batch"]:(i + 1) * HP["batch"]])
            opt.zero_grad(set_to_none=True); ls.backward(); opt.step(); sch.step()
        net.eval()
        with torch.no_grad():
            tot = np.zeros(3); n = 0
            torch.manual_seed(seed + 1000)                                  # ta sama podpróbka dekorelacji w każdej epoce
            for i in range(0, len(va), 4 * HP["batch"]):
                idx = va[i:i + 4 * HP["batch"]]
                tot += np.array([float(v) for v in loss_of(idx)]) * len(idx); n += len(idx)
            tot /= max(n, 1)
            torch.manual_seed(seed + ep + 1)
        hist.append([float(x) for x in tot])
        if tot[0] < best:
            best, best_w, best_ep = tot[0], copy.deepcopy(net.state_dict()), ep
    net.load_state_dict(best_w)
    net = net.to("cpu").eval()
    log(f"    epoka {best_ep}/{HP['epochs']}: val {best:.4f} (BCE {hist[best_ep][1]:.4f}, dec {hist[best_ep][2]:.5f}), {time.time() - t0:.0f}s")
    return net, {"best_epoch": best_ep, "hist": hist, "sec": round(time.time() - t0, 1)}


def predict(net, inputs, chunk=200_000, device="cpu"):
    n = len(inputs[0])
    out = np.empty(n, np.float32)
    net = net.to(device).eval()
    with torch.no_grad():
        for a in range(0, n, chunk):
            z = net(*[torch.as_tensor(np.asarray(x[a:a + chunk]), dtype=torch.float32).to(device) for x in inputs])
            out[a:a + chunk] = torch.sigmoid(z).cpu().numpy()
    net.to("cpu")
    return out


def export_onnx(net, example, path, names):
    """ONNX (opset 17), dynamiczny wymiar wsadu; kontrola onnxruntime vs torch na wsadzie INNEJ wielkości niż eksport. Zwraca dict (ok, max_abs_diff | error)."""
    m = Sig(net).eval()
    ex = tuple(torch.as_tensor(np.asarray(x), dtype=torch.float32) for x in example)
    try:
        kw = dict(input_names=list(names), output_names=["score"], dynamic_axes={n: {0: "batch"} for n in list(names) + ["score"]}, opset_version=17, verbose=False)
        try:
            torch.onnx.export(m, tuple(x[:32] for x in ex), path, dynamo=False, **kw)      # eksport na 32 wierszach, kontrola na wszystkich
        except TypeError:
            torch.onnx.export(m, tuple(x[:32] for x in ex), path, **kw)
        import onnxruntime as ort
        s = ort.InferenceSession(path, providers=["CPUExecutionProvider"])
        o = s.run(None, {n: x.numpy() for n, x in zip(names, ex)})[0].reshape(-1)
        with torch.no_grad():
            r = m(*ex).numpy().reshape(-1)
        return {"ok": True, "max_abs_diff": float(np.max(np.abs(o - r))), "n": int(len(r))}
    except Exception as e:  # noqa: BLE001
        return {"ok": False, "error": f"{type(e).__name__}: {e}"[:300]}
