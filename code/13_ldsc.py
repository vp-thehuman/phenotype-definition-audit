#!/usr/bin/env python3
"""Cross-trait LD score regression, implemented directly.

Univariate:  E[chi2_j]      = N h^2 l_j / M + 1 + a
Cross-trait: E[z1j z2j]     = sqrt(N1 N2) rho_g l_j / M + rho_e Ns / sqrt(N1 N2)
rg = rho_g / sqrt(h1^2 h2^2)

The sample-size term enters the numerator and both denominators and cancels in the
ratio, so rg is robust to the effective-N proxy used here. Heritability on this
scale is NOT interpreted, because the LD scores are computed on a thinned SNP set
and are therefore uniformly deflated.
"""
import os
ROOT=os.environ.get("ROOT","/tmp")
import pandas as pd, numpy as np, itertools, os, json
D=f"{ROOT}/mr/data/ldsc"
ld=pd.read_parquet(f"{ROOT}/ldsc/ldscores.parquet").rename(columns={'SNP':'rsid'})
M=len(ld)
comp={'A':'T','T':'A','C':'G','G':'C'}

def load(name):
    d=pd.read_parquet(f"{D}/{name}.parquet")
    d=d.merge(ld[['rsid','L2']],on='rsid')
    d=d[(d.z.abs()<30)]                      # standard chi2 < 900 filter
    return d

def wls(x,y,w):
    W=np.asarray(w); X=np.column_stack([np.ones(len(x)),x])
    XtW=X.T*W; b=np.linalg.solve(XtW@X,XtW@y)
    r=y-X@b; s2=(W*r**2).sum()/(len(x)-2)
    cov=np.linalg.inv(XtW@X)*s2*len(x)/max((W.sum()),1e-12)
    return b,np.sqrt(np.diag(np.abs(cov)))

def jackknife(f,n,blocks=200,seed=0):
    idx=np.arange(n); np.random.default_rng(seed).shuffle(idx)
    parts=np.array_split(idx,blocks); full=f(idx)
    ests=np.array([f(np.concatenate([p for k,p in enumerate(parts) if k!=i])) for i in range(blocks)])
    se=np.sqrt((blocks-1)/blocks*np.sum((ests-ests.mean())**2))
    return full,se

def h2(d):
    n=float(d.n_eff.iloc[0]); l=d.L2.values; chi=d.z.values**2
    w=1/np.maximum(l,1)
    def est(ix):
        b,_=wls(l[ix]*n/M,chi[ix],w[ix]); return b[1]
    return jackknife(est,len(d))

def rg_pair(a,b):
    m=a.merge(b,on='rsid',suffixes=('_1','_2'))
    sign=np.where(m.ea_1==m.ea_2,1,np.where(m.ea_1==m.ea_2.map(comp),1,
          np.where(m.ea_1==m.oa_2,-1,np.where(m.ea_1==m.oa_2.map(comp),-1,0))))
    m=m[sign!=0].copy(); m['z2a']=m.z_2.values*sign[sign!=0]
    n1=float(m.n_eff_1.iloc[0]); n2=float(m.n_eff_2.iloc[0])
    l=m.L2_1.values; zz=m.z_1.values*m.z2a.values; w=1/np.maximum(l,1)
    def est(ix):
        bb,_=wls(l[ix]*np.sqrt(n1*n2)/M,zz[ix],w[ix]); return bb[1]
    rho,rho_se=jackknife(est,len(m))
    return rho,rho_se,len(m)

names=[f[:-8] for f in sorted(os.listdir(D)) if f.endswith(".parquet")]
print("traits:",names,flush=True)
dat={n:load(n) for n in names}
H={}
for n in names:
    v,se=h2(dat[n]); H[n]=v
    print(f"h2* {n:20s} {v: .4f} (se {se:.4f})  [thinned-LD scale, not interpretable as h2]",flush=True)

rows=[]
for a,b in itertools.combinations(names,2):
    rho,rse,k=rg_pair(dat[a],dat[b])
    den=np.sqrt(max(H[a],1e-9)*max(H[b],1e-9))
    rg=rho/den; rgse=rse/den
    rows.append(dict(trait1=a,trait2=b,n_snp=k,rg=float(np.clip(rg,-1.5,1.5)),rg_se=float(rgse),
                     z=float(rg/rgse) if rgse>0 else np.nan))
    print(f"rg {a:20s} {b:20s} {np.clip(rg,-1.5,1.5): .3f} ({rgse:.3f})",flush=True)
pd.DataFrame(rows).to_csv(f"{ROOT}/mr/out/Table4_genetic_correlations.csv",index=False)
json.dump({k:float(v) for k,v in H.items()},open(f"{ROOT}/mr/out/h2_thinned_scale.json","w"),indent=2)
