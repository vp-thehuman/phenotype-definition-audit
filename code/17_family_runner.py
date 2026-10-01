#!/usr/bin/env python3
"""Generalised runner: apply the definition audit to any trait family.

A family is a set of GWAS of the same trait that differ in case definition,
one of which is nominated as the reference (strictest / best ascertained),
plus a candidate contaminating condition.
"""
import os
ROOT=os.environ.get("ROOT","/tmp")
import pandas as pd, numpy as np, os, subprocess, sys, json
G=f"{ROOT}/gwas"; D=f"{ROOT}/mr/data"; C=f"{ROOT}/mr/clump"
HM_A={'hm_rsid':'rsid','hm_chrom':'chr','hm_pos':'pos','hm_effect_allele':'ea','hm_other_allele':'oa',
      'hm_beta':'beta','hm_effect_allele_frequency':'eaf','standard_error':'se','p_value':'p'}
HM_B={'rsid':'rsid','chromosome':'chr','base_pair_location':'pos','effect_allele':'ea','other_allele':'oa',
      'beta':'beta','effect_allele_frequency':'eaf','standard_error':'se','p_value':'p'}
def pick_map(cols):
    return HM_A if 'hm_rsid' in cols else HM_B
FAMILY={
 "DEP1_one":      ("DEP1_one.h.tsv.gz",      "1 endorsed symptom measure",  57321+232552),
 "DEP2_two":      ("DEP2_two.h.tsv.gz",      "2 endorsed symptom measures", 21468+232552),
 "DEP3_three":    ("DEP3_three.h.tsv.gz",    "3 endorsed symptom measures",  9738+232552),
 "DEP4_fourfive": ("DEP4_fourfive.h.tsv.gz", "4-5 endorsed symptom measures",4887+232552),
 "DEP5_cidi":     ("DEP5_cidi.h.tsv.gz",     "CIDI lifetime depression",    28982+232552),
 "NEUROTICISM":   ("NEUROTICISM.h.tsv.gz",   "Neuroticism (contaminant)",   393411),
}
def stream(fn,name,n,pthresh=1e-5):
    p=f"{D}/{name}.parquet"
    if os.path.exists(p): print("skip",name,flush=True); return
    keep=[]; tot=0
    for ch in pd.read_csv(f"{G}/{fn}",sep="\t",compression='gzip',chunksize=1_000_000,low_memory=False):
        M=pick_map(ch.columns)
        src=[c for c in M if c in ch.columns]
        ch=ch[src].rename(columns=M); ch=ch.loc[:,~ch.columns.duplicated()]
        for c in ['pos','eaf','beta','se','p']:
            if c in ch.columns: ch[c]=pd.to_numeric(ch[c],errors='coerce')
        ch['n']=n; ch['scale']='linear'    # all BOLT-LMM linear on 0/1 (or continuous for neuroticism)
        tot+=len(ch); keep.append(ch.loc[ch.p<pthresh])
    d=pd.concat(keep,ignore_index=True).dropna(subset=['rsid','ea','oa','beta','se','p'])
    d=d[d.se>0]; d['rsid']=d.rsid.astype(str)
    d=d[d.rsid.str.startswith('rs')]
    d['ea']=d.ea.str.upper(); d['oa']=d.oa.str.upper(); d['chr']=d['chr'].astype(str)
    d=d.drop_duplicates('rsid')
    d[['rsid','chr','pos','ea','oa','eaf','beta','se','p','n','scale']].to_parquet(p,index=False)
    print(f"{name:16s} scanned={tot:>10,}  p<1e-5={len(d):>6,}  gws={(d.p<5e-8).sum():>5,}",flush=True)

for k,(fn,lab,n) in FAMILY.items(): stream(fn,k,n)

for k in FAMILY:
    if os.path.exists(f"{C}/{k}.clumps"): print("clumped",k,flush=True); continue
    d=pd.read_parquet(f"{D}/{k}.parquet")
    d[['rsid','p']].rename(columns={'rsid':'ID','p':'P'}).to_csv(f"{C}/{k}.assoc",sep="\t",index=False)
    subprocess.run([f"{ROOT}/plink2","--bfile",f"{ROOT}/ld/EUR","--clump",f"{C}/{k}.assoc",
        "--clump-p1","5e-8","--clump-r2","0.001","--clump-kb","10000",
        "--clump-id-field","ID","--clump-p-field","P","--out",f"{C}/{k}","--silent"],
        capture_output=True)
    n=sum(1 for _ in open(f"{C}/{k}.clumps"))-1 if os.path.exists(f"{C}/{k}.clumps") else 0
    print(f"{k:16s} instruments={n}",flush=True)

u=set()
for k in FAMILY:
    if os.path.exists(f"{C}/{k}.clumps"):
        u |= set(pd.read_csv(f"{C}/{k}.clumps",sep=r"\s+")['ID'].astype(str))
pd.Series(sorted(u)).to_csv(f"{D}/instrument_union_DEP.txt",index=False,header=False)
print("DEP union instruments:",len(u),flush=True)
