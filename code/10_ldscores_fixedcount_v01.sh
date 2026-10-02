#!/bin/bash
ROOT=${ROOT:-/tmp}; export ROOT
# Compute LD scores from 1000G EUR (n=503) directly, so the pipeline has no
# dependency on the (currently unavailable) precomputed eur_w_ld_chr files.
# Common SNPs (MAF>=0.05), 1 Mb window. Density is lower than the canonical
# HapMap3 set, so scores are uniformly deflated; genetic CORRELATION is a ratio
# and is therefore unaffected. Heritability on this scale is not interpreted.
# NOTE: this is the build used for the AD family (Tables 4-7). It samples a fixed
# COUNT per chromosome; the UK Biobank scan uses the fixed-FRACTION build in 22b.
set -e
mkdir -p ${ROOT}/ldsc/out
for c in $(seq 1 22); do
  ${ROOT}/plink2 --bfile ${ROOT}/ld/EUR --chr $c --maf 0.05 --thin-count 20000 --seed 42 \
    --r2-unphased --ld-window-kb 1000 --ld-window 999999 --ld-window-r2 0 \
    --out ${ROOT}/ldsc/c$c --silent 2>/dev/null || true
  python3 - "$c" <<'PY'
import sys,os,pandas as pd,numpy as np
ROOT=os.environ["ROOT"]; c=sys.argv[1]; N=503
d=pd.read_csv(f"{ROOT}/ldsc/c{c}.vcor",sep="\t")
d.columns=[x.replace('#','') for x in d.columns]
a,b,r2=d['ID_A'],d['ID_B'],d['UNPHASED_R2']
r2a=r2-(1-r2)/(N-2)                      # bias-corrected r^2
s=pd.concat([pd.Series(r2a.values,index=a.values),pd.Series(r2a.values,index=b.values)]).groupby(level=0).sum()
ids=pd.unique(pd.concat([a,b]))
bim=pd.read_csv(f"{ROOT}/ld/EUR.bim",sep=r"\s+",header=None,usecols=[1,3],names=['SNP','BP'])
out=pd.DataFrame({'SNP':ids}).assign(L2=lambda x:x.SNP.map(s).fillna(0)+1.0, CHR=int(c))
out['BP']=out.SNP.map(bim.drop_duplicates('SNP').set_index('SNP').BP)   # for contiguous jackknife blocks
out.to_parquet(f"{ROOT}/ldsc/out/ld_chr{c}.parquet",index=False)
print(f"chr{c}: {len(out):,} SNPs, mean L2 {out.L2.mean():.2f}",flush=True)
PY
  rm -f ${ROOT}/ldsc/c$c.vcor ${ROOT}/ldsc/c$c.log
done
python3 - <<'PY'
import os,pandas as pd,glob
ROOT=os.environ["ROOT"]
d=pd.concat([pd.read_parquet(f) for f in glob.glob(f"{ROOT}/ldsc/out/ld_chr*.parquet")],ignore_index=True)
d=d.sort_values(["CHR","BP"]).reset_index(drop=True)
d.to_parquet(f"{ROOT}/ldsc/ldscores.parquet",index=False)
print("TOTAL",len(d),"SNPs; mean L2",round(d.L2.mean(),2))
PY
