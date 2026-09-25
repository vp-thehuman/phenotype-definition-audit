#!/bin/bash
# Identical clumping for every exposure: p1=5e-8, r2=0.001, kb=10000, 1000G EUR (n=503)
set -e
cd /tmp/mr
mkdir -p clump
for E in E1_EAGLE_criteria E2_UKB_selfreport E3_UKB_ICD E4_allergic_broad E5_BUDU_meta; do
  python3 -c "
import pandas as pd
d=pd.read_parquet('data/$E.parquet')
d[['rsid','p']].rename(columns={'rsid':'ID','p':'P'}).to_csv('clump/$E.assoc',sep='\t',index=False)
"
  /tmp/plink2 --bfile /tmp/ld/EUR \
    --clump clump/$E.assoc --clump-p1 5e-8 --clump-r2 0.001 --clump-kb 10000 \
    --clump-id-field ID --clump-p-field P \
    --out clump/$E --silent 2>/dev/null || true
  n=$( [ -f clump/$E.clumps ] && tail -n +2 clump/$E.clumps | grep -c . || echo 0 )
  echo "$E  independent_instruments=$n"
done
