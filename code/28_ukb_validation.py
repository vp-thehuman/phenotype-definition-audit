import sys; sys.path.insert(0,'/tmp/ukb')
from ldsc import load,h2,rho
import math
tests=[('20002_1111','20002_1387','asthma SR vs hayfever SR'),
       ('20002_1286','20127_irnt','depression SR vs neuroticism'),
       ('20002_1111','6152_8','asthma SR vs touchscreen asthma'),
       ('21001_irnt','2178','BMI vs overall health rating'),
       ('20002_1065','6150_4','hypertension SR vs touchscreen HBP')]
d={}; H={}
for a,b,_ in tests:
    for x in (a,b):
        if x not in d:
            d[x]=load(x); H[x]=h2(d[x])[0]
for a,b,lab in tests:
    r,se,k=rho(d[a],d[b]); s=math.sqrt(abs(H[a]*H[b]))
    print(f"rg {lab:36s} {r/s: .3f}  (se {se/s:.3f}, k={k:,})",flush=True)
