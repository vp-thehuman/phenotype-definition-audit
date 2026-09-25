#!/usr/bin/env python3
"""Pass 1: stream each GWAS, keep only SNPs with p < 1e-5 (instrument candidates).
Standard schema: rsid, chr, pos, ea, oa, eaf, beta, se, p, n, scale
Match key throughout is rsID (avoids b37/b38 liftover: FinnGen R11 is GRCh38).
"""
import pandas as pd, numpy as np, os
G="/tmp/gwas"; OUT="/tmp/mr/data"; os.makedirs(OUT,exist_ok=True)
PTHRESH=1e-5
COLS=['rsid','chr','pos','ea','oa','eaf','beta','se','p','n','scale']

def stream(path, name, ren, scale, sep="\t", n=None, comp=None,
           beta_from_or=False, strip_quotes=False, eaf_col=None):
    keep=[]; tot=0
    for ch in pd.read_csv(path, sep=sep, compression=comp, chunksize=1_000_000,
                          low_memory=False, engine='c' if sep!="\s+" else 'python'):
        if strip_quotes:
            ch.columns=[c.strip('"') for c in ch.columns]
            for c in ch.columns:
                if ch[c].dtype==object: ch[c]=ch[c].str.strip('"')
        src=[c for c in ren if c in ch.columns]
        extra=[c for c in ['OR',eaf_col] if c and c in ch.columns]
        ch=ch[src+[c for c in extra if c not in src]].rename(columns=ren)
        ch=ch.loc[:,~ch.columns.duplicated()]
        if beta_from_or: ch['beta']=np.log(pd.to_numeric(ch['OR'],errors='coerce'))
        if eaf_col and eaf_col in ch.columns: ch['eaf']=pd.to_numeric(ch[eaf_col],errors='coerce')
        for c in ['pos','eaf','beta','se','p']:
            if c in ch.columns: ch[c]=pd.to_numeric(ch[c],errors='coerce')
        if 'eaf' not in ch.columns: ch['eaf']=np.nan
        if 'n' not in ch.columns: ch['n']=n
        else: ch['n']=pd.to_numeric(ch['n'],errors='coerce').fillna(n if n else np.nan)
        ch['scale']=scale
        tot+=len(ch)
        keep.append(ch.loc[ch.p<PTHRESH, COLS])
    d=pd.concat(keep,ignore_index=True)
    d=d.dropna(subset=['rsid','ea','oa','beta','se','p'])
    d=d[(d.se>0)&np.isfinite(d.beta)]
    d['rsid']=d.rsid.astype(str).str.split(',').str[0]
    d=d[d.rsid.str.startswith('rs')]
    d['ea']=d.ea.astype(str).str.upper(); d['oa']=d.oa.astype(str).str.upper()
    d['chr']=d['chr'].astype(str)
    d=d.drop_duplicates('rsid')
    d.to_parquet(f"{OUT}/{name}.parquet",index=False)
    print(f"{name:22s} scanned={tot:>11,}  p<1e-5={len(d):>7,}  gws(5e-8)={(d.p<5e-8).sum():>6,}  min_p={d.p.min():.2e}",flush=True)

# EXPOSURES
stream(f"{G}/EAGLE_AD.txt","E1_EAGLE_criteria",
       {'rsID':'rsid','chromosome':'chr','position':'pos','reference_allele':'ea',
        'other_allele':'oa','eaf':'eaf','beta':'beta','se':'se','European_N':'n','p.value':'p'},'logOR')

stream(f"{G}/UKB_SELFREP_ECZ.h.tsv.gz","E2_UKB_selfreport",
       {'hm_rsid':'rsid','hm_chrom':'chr','hm_pos':'pos','hm_effect_allele':'ea',
        'hm_other_allele':'oa','hm_beta':'beta','hm_effect_allele_frequency':'eaf',
        'standard_error':'se','p_value':'p'},'linear',comp='gzip',n=461199)

stream(f"{G}/UKB_ICD_ECZDERM.h.tsv.gz","E3_UKB_ICD",
       {'hm_rsid':'rsid','hm_chrom':'chr','hm_pos':'pos','hm_effect_allele':'ea',
        'hm_other_allele':'oa','hm_beta':'beta','hm_effect_allele_frequency':'eaf',
        'standard_error':'se','p_value':'p'},'linear',comp='gzip',n=484598)

stream(f"{G}/ALLERGIC_COMPOSITE.h.tsv.gz","E4_allergic_broad",
       {'hm_rsid':'rsid','hm_chrom':'chr','hm_pos':'pos','hm_effect_allele':'ea',
        'hm_other_allele':'oa','hm_beta':'beta','hm_effect_allele_frequency':'eaf',
        'standard_error':'se','p_value':'p'},'logOR',comp='gzip',n=360838)

stream(f"{G}/BUDU_AD.tsv","E5_BUDU_meta",
       {'RSID':'rsid','chromosome':'chr','base_pair_location':'pos','EA':'ea','OA':'oa',
        'effect_allele_frequency':'eaf','beta':'beta','standard_error':'se','p_value':'p'},'logOR',
       strip_quotes=True,n=864982)
print("exposures done",flush=True)
