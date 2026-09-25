#!/usr/bin/env python3
"""Extract instruments from the extra FinnGen outcomes spanning a range of
genetic relatedness to atopy."""
import pandas as pd, numpy as np, os
G="/tmp/gwas"; D="/tmp/mr/data"
snps=set(pd.read_csv(f"{D}/instrument_union.txt",header=None)[0].astype(str))
FG={'rsids':'rsid','ref':'oa','alt':'ea','pval':'p','beta':'beta','sebeta':'se','af_alt':'eaf'}
for fn,name in [("FG_ALLERG_RHINITIS.gz","O5_FG_rhinitis"),("FG_J10_COPD.gz","O6_FG_copd"),
                ("FG_K11_IBD_STRICT.gz","O7_FG_ibd"),("FG_G6_MIGRAINE.gz","O8_FG_migraine")]:
    p=f"{D}/{name}.parquet"
    if os.path.exists(p): print("skip",name,flush=True); continue
    out=[]
    for ch in pd.read_csv(f"{G}/{fn}",sep="\t",compression='gzip',chunksize=1_000_000,low_memory=False):
        ch=ch[[c for c in FG if c in ch.columns]].rename(columns=FG)
        ch['rsid']=ch['rsid'].astype(str).str.split(',').str[0]
        ch=ch[ch.rsid.isin(snps)]
        if len(ch): out.append(ch)
    d=pd.concat(out,ignore_index=True)
    for c in ['beta','se','p','eaf']: d[c]=pd.to_numeric(d[c],errors='coerce')
    d['ea']=d.ea.str.upper(); d['oa']=d.oa.str.upper(); d['n']=np.nan; d['scale']='logOR'
    d=d.dropna(subset=['beta','se']).drop_duplicates('rsid')
    d[['rsid','ea','oa','eaf','beta','se','p','n','scale']].to_parquet(p,index=False)
    print(f"{name:18s} matched={len(d)}",flush=True)
