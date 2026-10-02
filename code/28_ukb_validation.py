"""Validation: reproduce known genetic correlations before scanning anything new.

Expected (released build): self-reported asthma with hayfever about 0.38;
self-reported depression with neuroticism about 0.74; touchscreen and
self-reported hypertension close to 1.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ldsc3 import Store

tests = [('20002_1111', '20002_1387', 'asthma SR vs hayfever SR'),
         ('20002_1286', '20127_irnt', 'depression SR vs neuroticism'),
         ('20002_1111', '6152_8', 'asthma SR vs touchscreen asthma'),
         ('21001_irnt', '2178', 'BMI vs overall health rating'),
         ('20002_1065', '6150_4', 'hypertension SR vs touchscreen HBP')]
S = Store(sorted({x for a, b, _ in tests for x in (a, b)}))
for a, b, lab in tests:
    g, se = S.rg(a, b)
    print(f"rg {lab:36s} {g: .3f}  (se {se:.3f}, k={S.rho(a, b)[2]:,})", flush=True)
