"""DNET etap 5a — modele dodatkowe. Kod etapów 1–4 (../dnet1, ../dnet2, ../dnet4) NIE jest zmieniany.
  trp — transformer etapu 2 (dn2_models.SetNet „tr”) + cechy par w uwadze wg ParT (Qu, Li, Qian 2202.03772): ln Δ, ln k_T, ln z, ln m² dla par
        otoczenie–otoczenie i córka–otoczenie, osadzone małym MLP i dodane do logitów uwagi obu warstw (jedno U wspólne, jak w ParT).
        Para córka–córka jest WYZEROWANA (m² pary córek = masa kandydata → przeciek masy), token globalny nie ma cech par.
  Tagger otoczenia E nie ma własnej klasy: to SetNet „tr” z nglob = 2 (log pT, |η| kandydata) na samych tokenach otoczenia (dn5_lib.inputs)."""
import math

import torch
from torch import nn

import dn2_models as D2

NPAIR = 4
PAIR_SCALE = (2.0, 2.0, 1.0, 4.0)        # stała skala cech par (ln Δ, ln k_T, ln z, ln m²) — rząd wielkości; resztę dopasowuje pierwsza warstwa


def pair_feats(tok, mask):
    """tok [B, T, 11] SUROWE (log pT, Δη, Δφ względem D, ...; flagi is_K, is_pi, is_env), mask [B, T] → (cechy [B, T, T, 4], maska par [B, T, T]).
    Cząstki bezmasowe: m² = 2 pT_a pT_b (cosh Δη − cos Δφ); Δ = hypot(Δη, Δφ); k_T = min(pT) Δ; z = min(pT)/(pT_a + pT_b)."""
    lpt, eta, phi = tok[..., 0], tok[..., 1], tok[..., 2]
    de = eta[:, :, None] - eta[:, None, :]
    dp = torch.remainder(phi[:, :, None] - phi[:, None, :] + math.pi, 2 * math.pi) - math.pi
    lmin = torch.minimum(lpt[:, :, None], lpt[:, None, :])
    pa, pb = torch.exp(lpt)[:, :, None], torch.exp(lpt)[:, None, :]
    d2 = de * de + dp * dp
    ld = 0.5 * torch.log(d2 + 1e-8)
    lkt = lmin + ld
    lz = lmin - torch.log(pa + pb)
    lm2 = torch.log(torch.clamp(2 * pa * pb * (torch.cosh(de) - torch.cos(dp)), min=1e-8))
    f = torch.stack([ld, lkt, lz, lm2], -1) / torch.as_tensor(PAIR_SCALE, dtype=tok.dtype, device=tok.device)
    dau = tok[..., 8] + tok[..., 9]                                   # 1 dla córek (is_K lub is_pi)
    pm = mask[:, :, None] * mask[:, None, :] * (1.0 - dau[:, :, None] * dau[:, None, :])
    pm = pm * (1.0 - torch.eye(tok.shape[1], dtype=tok.dtype, device=tok.device))[None]
    return f, pm


class BlockP(D2.Block):
    def forward(self, x, m, U=None):
        B, T, d = x.shape
        q, k, v = self.qkv(self.n1(x)).reshape(B, T, 3, self.h, d // self.h).permute(2, 0, 3, 1, 4)
        att = (q @ k.transpose(-2, -1)) / math.sqrt(d // self.h) + (m[:, None, None, :] - 1.0) * 1e4
        if U is not None:
            att = att + U
        x = x + self.proj((torch.softmax(att, -1) @ v).transpose(1, 2).reshape(B, T, d))
        return x + self.ff(self.n2(x))


class SetNetP(D2.SetNet):
    """SetNet „tr” + osadzenie cech par dodawane do logitów uwagi. Sygnatura forward jak SetNet (tok, mask, glob)."""

    def __init__(self, nglob, tmean, tstd, gmed, gscale, d=64, heads=4, dp=32):
        super().__init__(nglob, "tr", tmean, tstd, gmed, gscale, d)
        self.arch = "trp"
        self.blocks = nn.ModuleList([BlockP(d, heads), BlockP(d, heads)])
        self.pair = nn.Sequential(nn.Linear(NPAIR, dp), nn.SiLU(), nn.Linear(dp, dp), nn.SiLU(), nn.Linear(dp, heads))

    def forward(self, tok, mask, glob):
        f, pm = pair_feats(tok, mask)
        U = (self.pair(f) * pm[..., None]).permute(0, 3, 1, 2)          # [B, h, T, T]
        U = torch.nn.functional.pad(U, (1, 0, 1, 0))                    # wiersz/kolumna 0 = token globalny, bez cech par
        x = torch.cat([(tok[..., :D2.NCONT] - self.tmean) / self.tstd, tok[..., D2.NCONT:]], -1)
        h = self.embed(x)
        g = self.gl(self.gprep(glob))
        h = torch.cat([g[:, None, :], h], 1)
        m = torch.cat([torch.ones_like(mask[:, :1]), mask], 1)
        for b in self.blocks:
            h = b(h, m, U)
        return self.head(self.norm(h[:, 0])).squeeze(-1)


def build(arch, nglob, tmean, tstd, gmed, gscale):
    if arch == "trp":
        return SetNetP(nglob, tmean, tstd, gmed, gscale)
    return D2.SetNet(nglob, arch, tmean, tstd, gmed, gscale)
