"""Exact delete-one-block jackknife for the contamination index.

Kept under its original number for the pipeline order. The implementation now
lives in ldsc3.py (imported by 28_ukb_validation.py and 29_ukb_scan.py), with
contiguous genomic jackknife blocks. Running this file prints a self-check.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ldsc3 import *  # noqa: F401,F403

if __name__ == '__main__':
    import numpy as np
    rng = np.random.default_rng(1)
    x = rng.uniform(1, 20, 5000); w = 1 / x; y = 0.3 * x + 2 + rng.normal(0, 1, 5000)
    lab = contiguous_blocks(len(x))
    f, loo = full_and_loo(block_sums(x, y, w, lab))
    W = np.diag(w); X = np.column_stack([np.ones_like(x), x])
    ref = np.linalg.solve(X.T @ W @ X, X.T @ W @ y)[1]
    print(f"five-sum slope {f:.6f} vs direct WLS {ref:.6f}; jackknife SE {jk_se(loo):.4f}")
