#!/usr/bin/env python3
"""MVMR needs every instrument's effect in EVERY exposure, including SNPs that do
not reach significance in that exposure. Re-stream the full files for the union set."""
import os
ROOT=os.environ.get("ROOT","/tmp")
import pandas as pd, numpy as np, os
G=f"{ROOT}/gwas"; D=f"{ROOT}/mr/data"
snps=set(pd.read_csv(f"{D}/instrument_union.txt",header=None)[0].astype(str))
print("union:",len(snps),flush=True)
SPECS=[
 ("EAGLE_AD.txt","X1_EAGLE",{'rsID':'rsid','chromosome':'chr','position':'pos','reference_allele':'ea',
   'other_allele':'oa','eaf':'eaf','beta':'beta','se':'se','European_N':'n','p.value':'p'},None,'logOR',"\t"),
 ("UKB_SELFREP_ECZ.h.tsv.gz","X2_UKB_selfreport",{'hm_rsid':'rsid','hm_effect_allele':'ea','hm_other_allele':'oa',
   'hm_beta':'beta','hm_effect_allele_frequency':'eaf','standard_error':'se','p_value':'p'},'gzip','linear',"\t"),
 ("UKB_ICD_ECZDERM.h.tsv.gz","X3_UKB_ICD",{'hm_rsid':'rsid','hm_effect_allele':'ea','hm_other_allele':'oa',
   'hm_beta':'beta','hm_effect_allele_frequency':'eaf','standard_error':'se','p_value':'p'},'gzip','linear',"\t"),
 ("ALLERGIC_COMPOSITE.h.tsv.gz","X4_allergic_broad",{'hm_rsid':'rsid','hm_effect_allele':'ea','hm_other_allele':'oa',
   'hm_beta':'beta','hm_effect_allele_frequency':'eaf','standard_error':'se','p_value':'p'},'gzip','logOR',"\t"),
 ("BUDU_AD.tsv","X5_BUDU",{'RSID':'rsid','EA':'ea','OA':'oa','effect_allele_frequency':'eaf',
   'beta':'beta','standard_error':'se','p_value':'p'},None,'logOR',"\t"),
]
NFIX={"X2_UKB_selfreport":461199,"X3_UKB_ICD":484598,"X4_allergic_broad":360838,"X5_BUDU":864982}
for fn,name,ren,comp,scale,sep in SPECS:
    out=[]
    for ch in pd.read_csv(f"{G}/{fn}",sep=sep,compression=comp,chunksize=1_000_000,low_memory=False):
        if name=="X5_BUDU":
            ch.columns=[c.strip('"') for c in ch.columns]
            for c in ch.columns:
                if ch[c].dtype==object: ch[c]=ch[c].str.strip('"')
        src=[c for c in ren if c in ch.columns]
        ch=ch[src].rename(columns=ren); ch=ch.loc[:,~ch.columns.duplicated()]
        ch['rsid']=ch['rsid'].astype(str)
        ch=ch[ch.rsid.isin(snps)]
        if len(ch): out.append(ch)
    d=pd.concat(out,ignore_index=True)
    for c in ['beta','se','p','eaf']: d[c]=pd.to_numeric(d[c],errors='coerce')
    d['n']=NFIX.get(name,np.nan)
    if scale=='linear':
        full=pd.read_parquet(f"{D}/{'E2_UKB_selfreport' if name.startswith('X2') else 'E3_UKB_ICD'}.parquet")
        v=full.dropna(subset=['eaf','se','n']); v=v[(v.eaf>0.05)&(v.eaf<0.95)]
        mu=float(np.median((v.se**2)*v.n*2*v.eaf*(1-v.eaf)))
        d['beta']/=mu; d['se']/=mu
    d['ea']=d.ea.astype(str).str.upper(); d['oa']=d.oa.astype(str).str.upper()
    d=d.dropna(subset=['beta','se']).drop_duplicates('rsid')
    d[['rsid','ea','oa','eaf','beta','se','p','n']].to_parquet(f"{D}/{name}.parquet",index=False)
    print(f"{name:22s} matched={len(d):>4}",flush=True)
