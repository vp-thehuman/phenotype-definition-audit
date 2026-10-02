"""Fast checks of the statistical core. No downloads.  Run: python3 tests/test_core.py"""
import os, sys, importlib.util
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "code"))
from ldsc3 import contiguous_blocks, block_sums, full_and_loo, jk_se, slope

def load_12():
    spec = importlib.util.spec_from_file_location("mvmr_src", os.path.join(HERE, "..", "code", "12_mvmr.py"))
    src = open(spec.origin).read().split("rows=[]")[0]      # helpers only, no data loading
    ns = {"__file__": spec.origin}; exec(src, ns); return ns

def test_five_sum_slope_matches_wls():
    rng = np.random.default_rng(0)
    x = rng.uniform(1, 30, 4000); w = 1 / x; y = 0.4 * x + 1 + rng.normal(0, 2, 4000)
    f, _ = full_and_loo(block_sums(x, y, w, contiguous_blocks(len(x))))
    X = np.column_stack([np.ones_like(x), x]); W = np.diag(w)
    assert abs(f - np.linalg.solve(X.T @ W @ X, X.T @ W @ y)[1]) < 1e-10

def test_loo_matches_refit():
    rng = np.random.default_rng(1)
    x = rng.uniform(1, 30, 2000); w = 1 / x; y = 0.2 * x + rng.normal(0, 1, 2000)
    lab = contiguous_blocks(len(x), 20)
    _, loo = full_and_loo(block_sums(x, y, w, lab, 20))
    k = lab != 7
    S = block_sums(x[k], y[k], w[k], np.zeros(k.sum(), dtype=np.int32), 1)
    assert abs(loo[7] - slope(S.sum(axis=1))) < 1e-10

def test_blocks_are_contiguous():
    lab = contiguous_blocks(1000, 200)
    assert np.all(np.diff(lab) >= 0) and lab.max() == 199

def test_random_blocks_understate_se_under_local_correlation():
    """Why v0.2 changed the blocks: with locally correlated data, shuffled blocks give smaller SEs."""
    rng = np.random.default_rng(2); n = 20000
    smooth = lambda v: np.convolve(v, np.ones(50) / 50, mode="valid")
    x = 10 + 40 * smooth(rng.normal(size=n + 49))       # LD scores vary smoothly along the genome
    e = smooth(rng.normal(size=n + 49))                  # and so do the test statistics
    w = np.ones(n); y = 0.1 * x + 5 * e * x
    lab_c = contiguous_blocks(n)
    lab_r = lab_c.copy(); rng.shuffle(lab_r)
    se_c = jk_se(full_and_loo(block_sums(x, y, w, lab_c))[1])
    se_r = jk_se(full_and_loo(block_sums(x, y, w, lab_r))[1])
    assert se_r < se_c

def test_conditional_f_includes_reference_error():
    ns = load_12(); rng = np.random.default_rng(3)
    bk = rng.normal(0, 0.05, 60); bj = 0.5 * bk + rng.normal(0, 0.03, 60)
    sj = np.full(60, 0.01); sk = np.full(60, 0.02)
    f_sw = ns["sw_condF"](bj, sj, bk, sk); f_zero = ns["sw_condF"](bj, sj, bk, np.zeros(60))
    assert f_sw < f_zero              # ignoring the reference SE overstates strength

if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn(); print("ok", name)
