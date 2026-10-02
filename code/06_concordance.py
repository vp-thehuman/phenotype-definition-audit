#!/usr/bin/env python3
"""Exposure-side construct test.

If two definitions index the same genetic liability but with different measurement
efficiency, instrument effects satisfy  beta_B = lambda * beta_A  with lambda a scalar
attenuation factor. Deming regression (errors in both variables) estimates lambda.
Key point: an MR ratio estimate beta_y/beta_x is invariant to lambda; an observational
association is not. That asymmetry is the mechanism this study is really about.
"""
import os
ROOT=os.environ.get("ROOT","/tmp")
import pandas as pd, numpy as np, itertools, json
D=f"{ROOT}/mr/data"
EXPS=["E1_EAGLE_criteria","E2_UKB_selfreport","E3_UKB_ICD","E4_allergic_broad","E5_BUDU_meta"]
LAB={"E1_EAGLE_criteria":"Clinician criteria (EAGLE)","E2_UKB_selfreport":"UKB hayfever/rhinitis/eczema",
     "E3_UKB_ICD":"Hospital ICD (UKB)","E4_allergic_broad":"Broad allergic composite",
     "E5_BUDU_meta":"Pooled meta-analysis"}
comp={'A':'T','T':'A','C':'G','G':'C'}
d={}
scales={}
for e in EXPS:
    x=pd.read_parquet(f"{D}/{e}.parquet")
    if x['scale'].iloc[0]=='linear':
        v=x.dropna(subset=['eaf','se','n']); v=v[(v.eaf>0.05)&(v.eaf<0.95)]
        mu=float(np.median((v.se**2)*v.n*2*v.eaf*(1-v.eaf))); scales[e]=mu
        x=x.assign(beta=x.beta/mu, se=x.se/mu)
    d[e]=x

def deming(x,y,sx,sy):
    """Deming regression through the origin, delta = var ratio."""
    delta=np.mean(sy**2)/np.mean(sx**2)
    sxx=np.sum(x**2); syy=np.sum(y**2); sxy=np.sum(x*y)
    num=syy-delta*sxx+np.sqrt((syy-delta*sxx)**2+4*delta*sxy**2)
    return num/(2*sxy)

REF="E1_EAGLE_criteria"
rows=[]
for e in EXPS:
    if e==REF: 
        rows.append(dict(definition=e,label=LAB[e],n_shared=np.nan,lambda_vs_EAGLE=1.0,pearson_r=1.0)); continue
    m=d[REF].merge(d[e],on='rsid',suffixes=('_a','_b'))
    sign=np.where(m.ea_a==m.ea_b,1,np.where(m.ea_a==m.ea_b.map(comp),1,
          np.where(m.ea_a==m.oa_b,-1,np.where(m.ea_a==m.oa_b.map(comp),-1,0))))
    m=m[sign!=0].copy(); m['bb']=m.beta_b*sign[sign!=0]
    # instruments defined by the REFERENCE only -> avoids selecting on the comparator
    m=m[m.p_a<5e-8]
    lam=deming(m.beta_a.values,m.bb.values,m.se_a.values,m.se_b.values)
    r=float(np.corrcoef(m.beta_a,m.bb)[0,1])
    rows.append(dict(definition=e,label=LAB[e],n_shared=len(m),lambda_vs_EAGLE=float(lam),pearson_r=r))
    print(f"{LAB[e]:32s} nSNP={len(m):>4}  lambda={lam:6.3f}  r={r:.3f}")
t=pd.DataFrame(rows)
t.to_csv(f"{ROOT}/mr/out/Table2_definition_concordance.csv",index=False)
l=t.lambda_vs_EAGLE.dropna()
print(f"\nlambda spans {l.min():.3f} to {l.max():.3f}  =>  {l.max()/l.min():.1f}-fold difference in "
      f"how strongly each definition captures liability")
print("mean |r| across definitions:",round(t.pearson_r.mean(),3))
json.dump(scales,open(f"{ROOT}/mr/out/scale_factors.json","w"),indent=2)

# Scale audit: mu(1-mu) is the variance of the 0/1 trait, so it reveals the case
# fraction of each linear-model file. Compare with what the file is supposed to be.
# (v0.2) For E2 this gives about 22%; UK Biobank self-reported eczema/dermatitis
# (20002_1452) is 2.6%, while 'hayfever, allergic rhinitis or eczema' (6152_9) is 23%.
EXPECTED={"E2_UKB_selfreport":"UKB 6152_9 hayfever/rhinitis/eczema is 23% (self-reported eczema 20002_1452 is 2.6%)",
          "E3_UKB_ICD":"12,176 / 484,598 = 2.5%"}
for e,v in scales.items():
    mu=(1-np.sqrt(max(1-4*v,0)))/2
    print(f"scale audit {LAB[e]:28s} mu(1-mu)={v:.4f} -> implied case fraction {mu:.1%}"
          f"   expected: {EXPECTED.get(e,'?')}")
    # v0.3: confirmed. rg(E2, UKB 6152_9 hayfever/rhinitis/eczema) = 0.98 (SE 0.04); rg(E2, 20002_1452
    # self-reported eczema) = 0.47 (SE 0.08). E2 is the touchscreen composite, labelled accordingly.
