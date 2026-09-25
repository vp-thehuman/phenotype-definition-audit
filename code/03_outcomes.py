#!/usr/bin/env python3
"""Extract instrument SNPs from the three outcome GWAS (streamed)."""
import pandas as pd, numpy as np, os, glob
G="/tmp/gwas"; D="/tmp/mr/data"
EXPS=["E1_EAGLE_criteria","E2_UKB_selfreport","E3_UKB_ICD","E4_allergic_broad","E5_BUDU_meta"]
snps=set()
for e in EXPS:
    c=pd.read_csv(f"/tmp/mr/clump/{e}.clumps",sep=r"\s+")
    snps |= set(c['ID'].astype(str))
print("union instruments:",len(snps),flush=True)
pd.Series(sorted(snps)).to_csv(f"{D}/instrument_union.txt",index=False,header=False)

def grab(path,name,ren,comp,scale,beta_from_or=False,sep="\t",eaf_col=None,nfix=None):
    out=[]
    kw=dict(sep=sep,compression=comp,chunksize=1_000_000)
    if sep=="\t": kw['low_memory']=False
    else: kw['engine']='python'
    for ch in pd.read_csv(path,**kw):
        ch=ch.rename(columns=ren)
        ch['rsid']=ch['rsid'].astype(str).str.split(',').str[0]
        ch=ch[ch.rsid.isin(snps)]
        if len(ch): out.append(ch)
    d=pd.concat(out,ignore_index=True) if out else pd.DataFrame()
    if beta_from_or: d['beta']=np.log(pd.to_numeric(d['OR'],errors='coerce'))
    if eaf_col and eaf_col in d.columns: d['eaf']=pd.to_numeric(d[eaf_col],errors='coerce')
    if 'eaf' not in d.columns: d['eaf']=np.nan
    if nfix: d['n']=nfix
    if 'n' not in d.columns: d['n']=np.nan
    for c in ['beta','se','p','eaf','n']:
        if c in d.columns: d[c]=pd.to_numeric(d[c],errors='coerce')
    d['ea']=d.ea.astype(str).str.upper(); d['oa']=d.oa.astype(str).str.upper()
    d['scale']=scale
    d=d.dropna(subset=['beta','se']).drop_duplicates('rsid')
    d[['rsid','ea','oa','eaf','beta','se','p','n','scale']].to_parquet(f"{D}/{name}.parquet",index=False)
    print(f"{name:20s} matched={len(d):>5,}",flush=True)

grab(f"{G}/FG_DEP.gz","O1_FG_depression",
     {'rsids':'rsid','ref':'oa','alt':'ea','pval':'p','sebeta':'se','af_alt':'eaf'},'gzip','logOR',nfix=None)
grab(f"{G}/FG_ANX.gz","O2_FG_anxiety",
     {'rsids':'rsid','ref':'oa','alt':'ea','pval':'p','sebeta':'se','af_alt':'eaf'},'gzip','logOR',nfix=None)
grab(f"{G}/PGC_MDD.gz","O3_PGC_MDD",
     {'SNP':'rsid','A1':'ea','A2':'oa','SE':'se','P':'p','Neff':'n','FRQ_U_113154':'eaf'},'gzip','logOR',
     beta_from_or=True,sep=r"\s+")
