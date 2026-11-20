"""DNET etap 2 — modele na zbiorze torów (córki z pełnymi parametrami + tory z otoczenia) z cechami globalnymi kandydata (H):
  ds — DeepSets (Zaheer i in. 1703.06114; w HEP: Energy Flow Networks 1810.05165, ATLAS DIPS): φ na tor → średnia + maksimum → ρ;
  tr — mały transformer z uwagą (wzorzec GN2, ATLAS 2505.19689): 2 warstwy, d = 64, 4 głowy, token globalny (cechy H) jako odczyt.
Uwaga własna (bez nn.MultiheadAttention) → pewny eksport ONNX. Wyjście: logit. Strata i pętla treningu: dnet1/dn_models.py (DisCo)."""
import math

import torch
from torch import nn

import dn_models as DM

TOK = ["log_pt", "deta", "dphi", "q_rel", "dedx", "npts", "d0", "z0", "is_K", "is_pi", "is_env"]
NCONT = 8          # pierwsze 8 cech tokenu standaryzowane (stałe z próbki treningowej), flagi bez zmian


class Block(nn.Module):
    def __init__(self, d, heads):
        super().__init__()
        self.h, self.n1, self.n2 = heads, nn.LayerNorm(d), nn.LayerNorm(d)
        self.qkv, self.proj = nn.Linear(d, 3 * d), nn.Linear(d, d)
        self.ff = nn.Sequential(nn.Linear(d, 2 * d), nn.SiLU(), nn.Linear(2 * d, d))

    def forward(self, x, m):
        B, T, d = x.shape
        q, k, v = self.qkv(self.n1(x)).reshape(B, T, 3, self.h, d // self.h).permute(2, 0, 3, 1, 4)
        att = (q @ k.transpose(-2, -1)) / math.sqrt(d // self.h) + (m[:, None, None, :] - 1.0) * 1e4
        x = x + self.proj((torch.softmax(att, -1) @ v).transpose(1, 2).reshape(B, T, d))
        return x + self.ff(self.n2(x))


class SetNet(nn.Module):
    def __init__(self, nglob, arch, tmean, tstd, gmed, gscale, d=64):
        super().__init__()
        self.arch = arch
        self.gprep = DM.Prep(gmed, gscale)
        self.register_buffer("tmean", torch.as_tensor(tmean, dtype=torch.float32))
        self.register_buffer("tstd", torch.as_tensor(tstd, dtype=torch.float32))
        self.embed = nn.Sequential(nn.Linear(len(TOK), d), nn.SiLU(), nn.Linear(d, d))
        self.gl = nn.Linear(nglob, d)
        if arch == "tr":
            self.blocks = nn.ModuleList([Block(d, 4), Block(d, 4)])
            self.norm = nn.LayerNorm(d)
            hin = d
        elif arch == "ds":
            hin = 3 * d
        else:
            raise ValueError(arch)
        self.head = nn.Sequential(nn.Linear(hin, 128), nn.SiLU(), nn.Linear(128, 64), nn.SiLU(), nn.Linear(64, 1))

    def forward(self, tok, mask, glob):
        x = torch.cat([(tok[..., :NCONT] - self.tmean) / self.tstd, tok[..., NCONT:]], -1)
        h = self.embed(x)
        g = self.gl(self.gprep(glob))
        if self.arch == "tr":
            h = torch.cat([g[:, None, :], h], 1)
            m = torch.cat([torch.ones_like(mask[:, :1]), mask], 1)
            for b in self.blocks:
                h = b(h, m)
            o = self.norm(h[:, 0])
        else:
            mm = mask[..., None]
            mean = (h * mm).sum(1) / torch.clamp(mm.sum(1), min=1.0)
            mx = (h + (mm - 1.0) * 1e4).max(1).values
            o = torch.cat([mean, mx, torch.nn.functional.silu(g)], -1)
        return self.head(o).squeeze(-1)
