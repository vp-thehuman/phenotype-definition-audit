#!/bin/bash
# Compute LD scores from 1000G EUR (n=503) directly, so the pipeline has no
# dependency on the (currently unavailable) precomputed eur_w_ld_chr files.
# Common SNPs (MAF>=0.05), 1 Mb window. Density is lower than the canonical
# HapMap3 set, so scores are uniformly deflated; genetic CORRELATION is a ratio
# and is therefore unaffected. Heritability on this scale is not interpreted.
set -e
mkdir -p /tmp/ldsc/out
for c in $(seq 1 22); do
  /tmp/plink2 --bfile /tmp/ld/EUR --chr $c --maf 0.05 --thin-count 20000 --seed 42 \
    --r2-unphased --ld-window-kb 1000 --ld-window 999999 --ld-window-r2 0 \
    --out /tmp/ldsc/c$c --silent 2>/dev/null || true
  python3 - "$c" <<'PY'
import sys,pandas as pd,numpy as np
c=sys.argv[1]; N=503
d=pd.read_csv(f"/tmp/ldsc/c{c}.vcor",sep="\t")
d.columns=[x.replace('#','') for x in d.columns]
a,b,r2=d['ID_A'],d['ID_B'],d['UNPHASED_R2']
r2a=r2-(1-r2)/(N-2)                      # bias-corrected r^2
s=pd.concat([pd.Series(r2a.values,index=a.values),pd.Series(r2a.values,index=b.values)]).groupby(level=0).sum()
ids=pd.unique(pd.concat([a,b]))
out=pd.DataFrame({'SNP':ids}).assign(L2=lambda x:x.SNP.map(s).fillna(0)+1.0, CHR=int(c))
out.to_parquet(f"/tmp/ldsc/out/ld_chr{c}.parquet",index=False)
print(f"chr{c}: {len(out):,} SNPs, mean L2 {out.L2.mean():.2f}",flush=True)
PY
  rm -f /tmp/ldsc/c$c.vcor /tmp/ldsc/c$c.log
done
python3 - <<'PY'
import pandas as pd,glob
d=pd.concat([pd.read_parquet(f) for f in glob.glob("/tmp/ldsc/out/ld_chr*.parquet")],ignore_index=True)
d.to_parquet("/tmp/ldsc/ldscores.parquet",index=False)
print("TOTAL",len(d),"SNPs; mean L2",round(d.L2.mean(),2))
PY
