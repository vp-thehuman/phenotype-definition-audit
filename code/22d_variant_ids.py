"""Map 1000 Genomes variants (chr:pos:ref:alt, GRCh37) to rsIDs.

The 1000 Genomes phase 3 v5b VCFs carry no rsIDs (ID column is '.'), so 22b names
each variant chr:pos:ref:alt and 22c converts names to rsIDs with this table, built
from the Neale lab round-2 variants manifest (also GRCh37, same naming).
"""
import os, gzip, subprocess
import pandas as pd
ROOT = os.environ.get("ROOT", "/tmp")
URL = "https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/variants.tsv.bgz"
out = f"{ROOT}/ld/id2rsid.parquet"
if not os.path.exists(out):
    proc = subprocess.Popen(['curl', '-sL', '--retry', '5', URL], stdout=subprocess.PIPE)
    f = gzip.open(proc.stdout, 'rt')
    h = f.readline().rstrip('\n').split('\t'); iv, ir = h.index('variant'), h.index('rsid')
    rows = []
    for line in f:
        q = line.split('\t', ir + 1)
        if q[ir].startswith('rs'):
            rows.append((q[iv], q[ir]))
    proc.wait()
    d = pd.DataFrame(rows, columns=['id', 'rsid']).drop_duplicates('id')
    d.to_parquet(out, index=False)
    print(f"id2rsid: {len(d):,} variants with rsIDs")
