#!/usr/bin/env python3
"""Sanity check that 'effect allele' is assigned correctly in every file.
All five exposures measure (variants of) the same trait. On shared SNPs, after
allele-matching, effect sizes must correlate POSITIVELY. A negative slope means
a column-mapping/sign error somewhere.
"""
import os
ROOT=os.environ.get("ROOT","/tmp")
import pandas as pd, numpy as np, itertools
D=f"{ROOT}/mr/data"
EXPS=["E1_EAGLE_criteria","E2_UKB_selfreport","E3_UKB_ICD","E4_allergic_broad","E5_BUDU_meta"]
d={e:pd.read_parquet(f"{D}/{e}.parquet") for e in EXPS}
# rescale linear -> logOR so slopes are interpretable
for e in EXPS:
    x=d[e]
    if x['scale'].iloc[0]=='linear':
        v=x.dropna(subset=['eaf','se','n']); v=v[(v.eaf>0.05)&(v.eaf<0.95)]
        mu=float(np.median((v.se**2)*v.n*2*v.eaf*(1-v.eaf)))
        d[e]=x.assign(beta=x.beta/mu, se=x.se/mu)
comp={'A':'T','T':'A','C':'G','G':'C'}
print(f"{'pair':<44} {'nSNP':>5} {'slope':>8} {'r':>7}  verdict")
bad=0
for a,b in itertools.combinations(EXPS,2):
    m=d[a].merge(d[b],on='rsid',suffixes=('_a','_b'))
    sign=np.where(m.ea_a==m.ea_b,1,np.where(m.ea_a==m.ea_b.map(comp),1,
          np.where(m.ea_a==m.oa_b,-1,np.where(m.ea_a==m.oa_b.map(comp),-1,0))))
    m=m[sign!=0].copy(); m['bb']=m.beta_b*sign[sign!=0]
    m=m[(m.p_a<5e-8)|(m.p_b<5e-8)]
    if len(m)<10: print(f"{a} vs {b}: too few shared"); continue
    sl=np.polyfit(m.beta_a,m.bb,1)[0]; r=np.corrcoef(m.beta_a,m.bb)[0,1]
    ok="OK" if r>0.3 else ("SIGN ERROR?" if r<-0.3 else "WEAK")
    if r<0.3: bad+=1
    print(f"{a[:20]:<21}vs {b[:20]:<21}{len(m):>5} {sl:>8.3f} {r:>7.3f}  {ok}")
print("\nsuspect pairs:",bad)
