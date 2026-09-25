#!/usr/bin/env python3
"""LD-score SNP set for the depression family and its candidate contaminant."""
import pandas as pd, numpy as np, os
G="/tmp/gwas"; D="/tmp/mr/data/ldsc"; os.makedirs(D,exist_ok=True)
snps=set(pd.read_parquet("/tmp/ldsc/ldscores.parquet").SNP.astype(str))
A={'hm_rsid':'rsid','hm_effect_allele':'ea','hm_other_allele':'oa','hm_beta':'beta',
   'hm_effect_allele_frequency':'eaf','standard_error':'se','p_value':'p'}
B={'rsid':'rsid','effect_allele':'ea','other_allele':'oa','beta':'beta',
   'effect_allele_frequency':'eaf','standard_error':'se','p_value':'p'}
JOBS=[("DEP1_one.h.tsv.gz","DEP_1sym"),("DEP2_two.h.tsv.gz","DEP_2sym"),
      ("DEP3_three.h.tsv.gz","DEP_3sym"),("DEP4_fourfive.h.tsv.gz","DEP_45sym"),
      ("DEP5_cidi.h.tsv.gz","DEP_cidi"),("NEUROTICISM.h.tsv.gz","CON_neuroticism")]
for fn,name in JOBS:
    p=f"{D}/{name}.parquet"
    if os.path.exists(p): print("skip",name,flush=True); continue
    out=[]
    for ch in pd.read_csv(f"{G}/{fn}",sep="\t",compression='gzip',chunksize=1_000_000,low_memory=False):
        M=A if 'hm_rsid' in ch.columns else B
        src=[c for c in M if c in ch.columns]
        ch=ch[src].rename(columns=M); ch=ch.loc[:,~ch.columns.duplicated()]
        ch['rsid']=ch['rsid'].astype(str)
        ch=ch[ch.rsid.isin(snps)]
        if len(ch): out.append(ch)
    d=pd.concat(out,ignore_index=True)
    for c in ['beta','se','p','eaf']: d[c]=pd.to_numeric(d[c],errors='coerce')
    d['ea']=d.ea.str.upper(); d['oa']=d.oa.str.upper()
    d=d.dropna(subset=['beta','se']); d=d[d.se>0].drop_duplicates('rsid')
    d['z']=d.beta/d.se
    v=d.dropna(subset=['eaf']); v=v[(v.eaf>0.05)&(v.eaf<0.95)]
    d['n_eff']=float(np.median(1/(2*v.eaf*(1-v.eaf)*v.se**2)))
    d[['rsid','ea','oa','eaf','beta','se','z','p','n_eff']].to_parquet(p,index=False)
    print(f"{name:18s} {len(d):>7,} SNPs",flush=True)
print("done")
