"""Shared LD score regression engine for the UK Biobank scan.

Cross-trait LD score regression with a delete-one-block jackknife computed from
block sufficient statistics, and the contamination index

    excess = rg(D,C) - rg(D,R) rg(R,C)

re-formed inside every jackknife replicate (its three rg terms share SNPs and are
not independent, so their standard errors cannot simply be propagated).

Jackknife blocks are CONTIGUOUS runs of SNPs in genomic order (chromosome, then
position), as in reference LDSC. Neighbouring SNPs are in LD; random SNP blocks
would put correlated SNPs in different blocks and understate standard errors.
The initial release (v0.1) shuffled SNPs before blocking; see CHANGELOG.md.

Used by 25_ukb_ldsc_fast.py, 28_ukb_validation.py and 29_ukb_scan.py.
"""
import os
import numpy as np, pandas as pd

ROOT = os.environ.get("ROOT", "/tmp")
BL = 200


def load_ldscores(path=None):
    """LD scores sorted in genomic order. Uses BP if present (22c writes it);
    otherwise keeps the within-chromosome file order, which is positional."""
    ld = pd.read_parquet(path or f"{ROOT}/ld/ldscores.parquet").rename(columns={'SNP': 'rsid'})
    ld = ld.drop_duplicates('rsid')
    keys = ['CHR', 'BP'] if 'BP' in ld.columns else ['CHR']
    return ld.sort_values(keys, kind='mergesort').reset_index(drop=True)


def contiguous_blocks(n, nblocks=BL):
    """Block label for each of n SNPs already in genomic order."""
    lab = np.empty(n, dtype=np.int32)
    for i, p in enumerate(np.array_split(np.arange(n), nblocks)):
        lab[p] = i
    return lab


def block_sums(x, y, w, lab, nblocks=BL):
    """The five weighted sums that determine a WLS slope, per block: shape (5, B)."""
    return np.stack([np.bincount(lab, weights=v, minlength=nblocks)
                     for v in (w, w * x, w * x * x, w * y, w * x * y)])


def slope(S):
    A, B, C, Dd, E = S
    return (A * E - B * Dd) / (A * C - B * B)


def full_and_loo(S):
    """Full-data slope and the B leave-one-block-out slopes."""
    T = S.sum(axis=1)
    return slope(T), slope(T[:, None] - S)


def jk_se(loo):
    B = len(loo)
    return float(np.sqrt((B - 1) / B * ((loo - loo.mean()) ** 2).sum()))


class Store:
    def __init__(self, names, ld=None, ssdir=None):
        ld = load_ldscores() if ld is None else ld
        self.M = len(ld)
        self.rsid = ld.rsid.values
        self.l2 = ld.L2.values.astype(float)
        self.w = 1 / np.maximum(self.l2, 1)
        self.lab = contiguous_blocks(len(self.l2))
        ssdir = ssdir or f"{ROOT}/ukb/ss"
        self.z = {}; self.neff = {}
        for n in names:
            d = pd.read_parquet(f"{ssdir}/{n}.parquet", columns=['rsid', 'z', 'n_eff'])
            d = d.drop_duplicates('rsid')
            z = pd.Series(d.z.values, index=d.rsid.values).reindex(self.rsid).values
            self.z[n] = np.where(np.abs(z) < 30, z, np.nan)       # chi2 < 900
            self.neff[n] = float(d.n_eff.iloc[0])
        self._h = {}; self._r = {}

    # Regression weights follow reference LDSC: 1/max(l,1) for LD between regression
    # SNPs, times the inverse of the expected variance of the test statistic given a
    # first-pass estimate (heteroskedasticity weights). v0.2 used 1/max(l,1) only,
    # which is unbiased but inefficient; with contiguous jackknife blocks that
    # inefficiency shows up as much wider standard errors.
    def _h2_weights(self, n, ok):
        x = self.l2[ok] * self.neff[n] / self.M; y = self.z[n][ok] ** 2
        b = slope(block_sums(x, y, self.w[ok], np.zeros(ok.sum(), dtype=np.int32), 1).sum(axis=1))
        h = float(np.clip(b, 0, 1))
        return x, y, self.w[ok] / (2 * (1 + self.neff[n] * h * self.l2[ok] / self.M) ** 2), h

    def h2S(self, n):
        if n not in self._h:
            ok = np.isfinite(self.z[n])
            x, y, w, h = self._h2_weights(n, ok)
            self._h[n] = block_sums(x, y, w, self.lab[ok]); self._hfirst = getattr(self, '_hfirst', {})
            self._hfirst[n] = h
        return self._h[n]

    def rhoS(self, a, b):
        k = tuple(sorted((a, b)))
        if k not in self._r:
            za, zb = self.z[a], self.z[b]; ok = np.isfinite(za) & np.isfinite(zb)
            l = self.l2[ok]; Na, Nb = self.neff[a], self.neff[b]
            x = l * np.sqrt(Na * Nb) / self.M; y = za[ok] * zb[ok]
            r0 = slope(block_sums(x, y, self.w[ok], np.zeros(ok.sum(), dtype=np.int32), 1).sum(axis=1))
            self.h2S(a); self.h2S(b)
            ha, hb = self._hfirst[a], self._hfirst[b]
            var = (1 + Na * ha * l / self.M) * (1 + Nb * hb * l / self.M) + (np.sqrt(Na * Nb) * r0 * l / self.M) ** 2
            self._r[k] = block_sums(x, y, self.w[ok] / var, self.lab[ok])
        return self._r[k]

    def h2(self, n):
        f, l = full_and_loo(self.h2S(n)); return f, jk_se(l)

    def rho(self, a, b):
        f, l = full_and_loo(self.rhoS(a, b))
        za, zb = self.z[a], self.z[b]
        return f, jk_se(l), int((np.isfinite(za) & np.isfinite(zb)).sum())

    def _rg_full_loo(self, a, b):
        rf, rl = full_and_loo(self.rhoS(a, b))
        af, al = full_and_loo(self.h2S(a)); bf, bl = full_and_loo(self.h2S(b))
        return rf / np.sqrt(abs(af * bf)), rl / np.sqrt(np.abs(al * bl))

    def rg(self, a, b):
        g, gl = self._rg_full_loo(a, b); return g, jk_se(gl)

    def excess(self, Dp, Rp, Cp, winsor=1.25):
        dc, dcl = self._rg_full_loo(Dp, Cp)
        dr, drl = self._rg_full_loo(Dp, Rp)
        rc, rcl = self._rg_full_loo(Rp, Cp)
        cl = (lambda v: np.clip(v, -winsor, winsor)) if winsor else (lambda v: v)
        e = cl(dc) - cl(dr) * cl(rc)
        el = cl(dcl) - cl(drl) * cl(rcl)
        return e, jk_se(el), dict(rg_DC=dc, rg_DR=dr, rg_RC=rc)
