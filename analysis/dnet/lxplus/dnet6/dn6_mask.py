"""DNET etap 6 — maska toru zgodnego z wolnym pionem D*→D0π w tokenach otoczenia. BEZ prawdy MC: identyczna w danych i MC, dla sygnału i tła,
w treningu i w zastosowaniu. Tor otoczenia jest usuwany (token wyzerowany → maska 0 w dn2_train.tokens), gdy
    Δm = m(D0 + tor z masą pionu) − m(D0) < DM_MAX   oraz   ładunek toru = −ładunek K kandydata (D*+ → D0 π+, D0 → K− π+).
m(D0) = masa kandydata Kπ w jego przypisaniu; czteropęd D0 z (pT, η, φ, m) kandydata; tor z (pT, η_D + Δη, φ_D + Δφ) tokenu (float16).
Sprawdzenie wzoru i dobór DM_MAX: dn6_mask.py study (MC: Δm tokenu prawdziwego wolnego pionu wobec p_dm tabeli; udział usuniętych; strata torów)."""
import numpy as np

MPI = 0.13957
DM_MAX = float(__import__("os").environ.get("DN6_DM_MAX", "0.155"))     # GeV; uzasadnienie: DNET/KRYTERIA_etap6.md §1
SIGN = int(__import__("os").environ.get("DN6_SIGN", "1"))                # 1 = tylko ładunek zgodny (q_tor · q_K = −1); 0 = oba ładunki


def token_dm(E, pt, eta, phi, m):
    """Δm [n, NT] (float64; NaN dla pustych tokenów). E: (pT, Δη, Δφ, q, dE/dx, npts) [n, NT, 6]."""
    E = E.astype(np.float64)
    pt, eta, phi, m = (np.asarray(x, np.float64)[:, None] for x in (pt, eta, phi, m))
    tpt = E[:, :, 0]
    teta = eta + E[:, :, 1]
    ED = np.sqrt(m ** 2 + (pt * np.cosh(eta)) ** 2)
    Et = np.sqrt(MPI ** 2 + (tpt * np.cosh(teta)) ** 2)
    m2 = m ** 2 + MPI ** 2 + 2 * (ED * Et - pt * tpt * np.cos(E[:, :, 2]) - pt * np.sinh(eta) * tpt * np.sinh(teta))
    dm = np.sqrt(np.clip(m2, 0, None)) - m
    return np.where(tpt > 0, dm, np.nan)


def veto(E, Dd, pt, eta, phi, m, dm_max=None, sign=None):
    """→ bool [n, NT]: tokeny do usunięcia."""
    dm_max = DM_MAX if dm_max is None else dm_max
    sign = SIGN if sign is None else sign
    dm = token_dm(E, pt, eta, phi, m)
    v = np.isfinite(dm) & (dm < dm_max)
    if sign:
        qK = np.nan_to_num(Dd[:, 0, 7].astype(np.float64))
        v &= (E[:, :, 3].astype(np.float64) * qK[:, None]) < 0
    return v


def apply(E, Dd, pt, eta, phi, m, chunk=500_000):
    """→ (E z wyzerowanymi tokenami [kopia], liczba usuniętych na wiersz [int8])."""
    out = np.array(E, copy=True)
    nrm = np.zeros(len(out), np.int8)
    for a0 in range(0, len(out), chunk):
        sl = slice(a0, a0 + chunk)
        v = veto(out[sl], Dd[sl], pt[sl], eta[sl], phi[sl], m[sl])
        blk = out[sl]
        blk[v] = 0
        out[sl] = blk
        nrm[sl] = v.sum(1)
    return out, nrm


def mask_tables(S, log=None):
    """W miejscu: S[próbka][a_E / p_E] ← po masce; dopisuje a_nrm / p_nrm (liczba usuniętych) i a_E0 NIE jest trzymane (pamięć)."""
    for s, T in S.items():
        for pre in ("a_", "p_"):
            if pre + "E" not in T:
                continue
            T[pre + "E"], T[pre + "nrm"] = apply(T[pre + "E"], T[pre + "Dd"], T[pre + "pt"], T[pre + "eta"], T[pre + "phi"], T[pre + "m"])
            if log:
                E = T[pre + "E"]
                log(f"  maska D* ({s} {pre}): wierszy z ≥ 1 usuniętym torem {float((T[pre + 'nrm'] > 0).mean()):.4f}, usuniętych torów/wiersz "
                    f"{float(T[pre + 'nrm'].mean()):.4f} (Δm < {DM_MAX}, znak {SIGN})")
    return S
