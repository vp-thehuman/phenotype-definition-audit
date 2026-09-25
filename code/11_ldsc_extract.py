#!/usr/bin/env python3
"""Extract the LD-score SNP set from every GWAS for LD score regression."""
import pandas as pd, numpy as np, os
G="/tmp/gwas"; D="/tmp/mr/data/ldsc"; os.makedirs(D,exist_ok=True)
snps=set(pd.read_parquet("/tmp/ldsc/ldscores.parquet").SNP.astype(str))
print("LD-score SNPs:",len(snps),flush=True)
FG={'rsids':'rsid','ref':'oa','alt':'ea','pval':'p','beta':'beta','sebeta':'se','af_alt':'eaf'}
HM={'hm_rsid':'rsid','hm_effect_allele':'ea','hm_other_allele':'oa','hm_beta':'beta',
    'hm_effect_allele_frequency':'eaf','standard_error':'se','p_value':'p'}
JOBS=[
 ("EAGLE_AD.txt","AD_criteria",{'rsID':'rsid','reference_allele':'ea','other_allele':'oa','eaf':'eaf',
   'beta':'beta','se':'se','p.value':'p'},None,"\t",False),
 ("UKB_SELFREP_ECZ.h.tsv.gz","AD_selfreport",HM,'gzip',"\t",False),
 ("UKB_ICD_ECZDERM.h.tsv.gz","AD_icd",HM,'gzip',"\t",False),
 ("ALLERGIC_COMPOSITE.h.tsv.gz","AD_allergic_broad",HM,'gzip',"\t",False),
 ("BUDU_AD.tsv","AD_meta",{'RSID':'rsid','EA':'ea','OA':'oa','effect_allele_frequency':'eaf',
   'beta':'beta','standard_error':'se','p_value':'p'},None,"\t",True),
 ("FG_ASTHMA.gz","OUT_asthma",FG,'gzip',"\t",False),
 ("FG_ALLERG_RHINITIS.gz","OUT_rhinitis",FG,'gzip',"\t",False),
 ("FG_J10_COPD.gz","OUT_copd",FG,'gzip',"\t",False),
 ("FG_K11_IBD_STRICT.gz","OUT_ibd",FG,'gzip',"\t",False),
 ("FG_G6_MIGRAINE.gz","OUT_migraine",FG,'gzip',"\t",False),
 ("FG_DEP.gz","OUT_depression",FG,'gzip',"\t",False),
 ("FG_ANX.gz","OUT_anxiety",FG,'gzip',"\t",False),
]
for fn,name,ren,comp,sep,quo in JOBS:
    p=f"{D}/{name}.parquet"
    if os.path.exists(p): print("skip",name,flush=True); continue
    out=[]
    for ch in pd.read_csv(f"{G}/{fn}",sep=sep,compression=comp,chunksize=1_000_000,low_memory=False):
        if quo:
            ch.columns=[c.strip('"') for c in ch.columns]
            for c in ch.columns:
                if ch[c].dtype==object: ch[c]=ch[c].str.strip('"')
        src=[c for c in ren if c in ch.columns]
        ch=ch[src].rename(columns=ren); ch=ch.loc[:,~ch.columns.duplicated()]
        ch['rsid']=ch['rsid'].astype(str).str.split(',').str[0]
        ch=ch[ch.rsid.isin(snps)]
        if len(ch): out.append(ch)
    d=pd.concat(out,ignore_index=True)
    for c in ['beta','se','p','eaf']:
        if c in d.columns: d[c]=pd.to_numeric(d[c],errors='coerce')
    d['ea']=d.ea.astype(str).str.upper(); d['oa']=d.oa.astype(str).str.upper()
    d=d.dropna(subset=['beta','se']); d=d[d.se>0].drop_duplicates('rsid')
    d['z']=d.beta/d.se
    # effective N from SE
    v=d.dropna(subset=['eaf']); v=v[(v.eaf>0.05)&(v.eaf<0.95)]
    d['n_eff']=float(np.median(1/(2*v.eaf*(1-v.eaf)*v.se**2)))
    d[['rsid','ea','oa','eaf','beta','se','z','p','n_eff']].to_parquet(p,index=False)
    print(f"{name:20s} {len(d):>7,} SNPs  N_eff={d.n_eff.iloc[0]:,.0f}",flush=True)
print("done")
