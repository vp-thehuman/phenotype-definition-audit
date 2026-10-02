"""Whitelist of Neale round-2 variants that are in the LD-score SNP set.

24_ukb_stream.py keeps only these rows from each streamed GWAS. Columns:
variant (Neale ID chr:pos:ref:alt), rsid, AF (alt-allele frequency, UK Biobank).
"""
import os, gzip, subprocess
import pandas as pd
ROOT = os.environ.get("ROOT", "/tmp")
URL = "https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/variants.tsv.bgz"
snps = set(pd.read_parquet(f"{ROOT}/ld/ldscores.parquet").SNP.astype(str))
os.makedirs(f"{ROOT}/ukb", exist_ok=True)
proc = subprocess.Popen(['curl', '-sL', '--retry', '5', URL], stdout=subprocess.PIPE)
f = gzip.open(proc.stdout, 'rt')
h = f.readline().rstrip('\n').split('\t'); ix = {c: i for i, c in enumerate(h)}
rows = []
for line in f:
    q = line.rstrip('\n').split('\t')
    if q[ix['rsid']] in snps:
        rows.append((q[ix['variant']], q[ix['rsid']], float(q[ix['AF']])))
proc.wait()
W = pd.DataFrame(rows, columns=['variant', 'rsid', 'AF']).drop_duplicates('rsid')
W = W[(W.AF > 0) & (W.AF < 1)]
W.to_parquet(f"{ROOT}/ukb/whitelist.parquet", index=False)
print(f"whitelist: {len(W):,} variants of {len(snps):,} LD-score SNPs")
