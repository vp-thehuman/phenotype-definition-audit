#!/usr/bin/env python3
"""Contamination index for the depression family, and a cross-family summary.

The depression ladder (Glanville 2021) varies case definition within ONE cohort,
one array, one control set, so definition is not confounded with cohort, ancestry
or genotyping. It IS confounded with case count, which is why instrument counts
collapse for the strict definitions and MR is not attempted here. The index does
not need instruments.
"""
import os
ROOT=os.environ.get("ROOT","/tmp")
import pandas as pd, numpy as np, itertools, json, os
D=f"{ROOT}/mr/data/ldsc"; O=f"{ROOT}/mr/out"
ld=pd.read_parquet(f"{ROOT}/ldsc/ldscores.parquet").rename(columns={'SNP':'rsid'})
LDCOLS=['rsid','L2']+[c for c in ('CHR','BP') if c in ld.columns]; SORTK=[c for c in ('CHR','BP') if c in ld.columns] or ['rsid']; M=len(ld)
comp={'A':'T','T':'A','C':'G','G':'C'}
NAMES=["DEP_1sym","DEP_2sym","DEP_3sym","DEP_45sym","DEP_cidi","CON_neuroticism"]
LAB={"DEP_1sym":"1 endorsed measure","DEP_2sym":"2 endorsed measures","DEP_3sym":"3 endorsed measures",
     "DEP_45sym":"4-5 endorsed measures","DEP_cidi":"CIDI lifetime depression (reference)",
     "CON_neuroticism":"Neuroticism"}
def load(n):
    d=pd.read_parquet(f"{D}/{n}.parquet").merge(ld[LDCOLS],on='rsid').sort_values(SORTK)
    return d[d.z.abs()<30]
def wls(x,y,w):
    X=np.column_stack([np.ones(len(x)),x]); XtW=X.T*w
    return np.linalg.solve(XtW@X,XtW@y)
def jk(f,n,blocks=200,seed=0):
    idx=np.arange(n); parts=np.array_split(idx,blocks)   # contiguous genomic blocks: data are sorted by CHR, BP
    full=f(idx)
    e=np.array([f(np.concatenate([p for k,p in enumerate(parts) if k!=i])) for i in range(blocks)])
    return full,np.sqrt((blocks-1)/blocks*np.sum((e-e.mean())**2))
def h2(d):
    n=float(d.n_eff.iloc[0]); l=d.L2.values; chi=d.z.values**2; w=1/np.maximum(l,1)
    return jk(lambda ix: wls(l[ix]*n/M,chi[ix],w[ix])[1], len(d))
def rg(a,b):
    m=a.merge(b,on='rsid',suffixes=('_1','_2'))
    m=m.sort_values([k+'_1' for k in SORTK if k!='rsid'] or ['rsid']).reset_index(drop=True)
    s=np.where(m.ea_1==m.ea_2,1,np.where(m.ea_1==m.ea_2.map(comp),1,
      np.where(m.ea_1==m.oa_2,-1,np.where(m.ea_1==m.oa_2.map(comp),-1,0))))
    m=m[s!=0].copy(); zz=m.z_1.values*(m.z_2.values*s[s!=0])
    n1=float(m.n_eff_1.iloc[0]); n2=float(m.n_eff_2.iloc[0]); l=m.L2_1.values; w=1/np.maximum(l,1)
    return jk(lambda ix: wls(l[ix]*np.sqrt(n1*n2)/M,zz[ix],w[ix])[1], len(m))
dat={n:load(n) for n in NAMES}
H={n:h2(dat[n])[0] for n in NAMES}
def RG(a,b):
    if a==b: return 1.0
    r,_=rg(dat[a],dat[b]); return float(np.clip(r/np.sqrt(max(H[a],1e-9)*max(H[b],1e-9)),-1,1))
REF="DEP_cidi"; CON="CON_neuroticism"
rows=[]
for n in ["DEP_1sym","DEP_2sym","DEP_3sym","DEP_45sym","DEP_cidi"]:
    obs=RG(n,CON); ref=RG(n,REF); exp=ref*RG(REF,CON)
    rows.append(dict(family="Depression (Glanville 2021, UK Biobank)",definition=n,label=LAB[n],
                     rg_vs_reference=ref,rg_contaminant=obs,expected_if_pure=exp,excess_rg=obs-exp))
    print(f"{LAB[n]:38s} rg_ref={ref: .3f}  rg_neuroticism={obs: .3f}  expected={exp: .3f}  EXCESS={obs-exp: .3f}",flush=True)
t=pd.DataFrame(rows); t.to_csv(f"{O}/Table7_depression_family_index.csv",index=False)
print(f"\nrg(CIDI depression, neuroticism) = {RG(REF,CON):.3f}")
