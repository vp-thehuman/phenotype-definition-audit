#!/usr/bin/env python3
"""Two things:
1. A contamination index. Raw rg(definition, related condition) conflates genuine
   shared biology with phenotype contamination: clinician-criteria AD is itself
   correlated with rhinitis (atopic march). The excess over what a noisy measure of
   clinician-criteria AD would predict isolates contamination:
       excess = rg(D, C) - rg(D, AD_criteria) * rg(AD_criteria, C)
2. An outcome ladder: does definitional inflation scale with how genetically
   related the outcome is to the contaminating condition?
"""
import os
ROOT=os.environ.get("ROOT","/tmp")
import pandas as pd, numpy as np
from scipy import stats
O=f"{ROOT}/mr/out"; D=f"{ROOT}/mr/data"
t=pd.read_csv(f"{O}/Table4_genetic_correlations.csv")
def rg(a,b):
    if a==b: return 1.0
    r=t[((t.trait1==a)&(t.trait2==b))|((t.trait1==b)&(t.trait2==a))]
    return float(np.clip(r.rg.iloc[0],-1,1)) if len(r) else np.nan
KEY={'E1_EAGLE_criteria':'AD_criteria','E3_UKB_ICD':'AD_icd','E5_BUDU_meta':'AD_meta',
     'E2_UKB_selfreport':'AD_selfreport','E4_allergic_broad':'AD_allergic_broad'}
LAB={'E1_EAGLE_criteria':'Clinician criteria','E3_UKB_ICD':'Hospital ICD',
     'E5_BUDU_meta':'Pooled meta-analysis','E2_UKB_selfreport':'Self-report',
     'E4_allergic_broad':'Allergic composite'}
C='OUT_rhinitis'; REF='AD_criteria'
idx=[]
for e,tr in KEY.items():
    obs=rg(tr,C); via=rg(tr,REF)*rg(REF,C)
    idx.append(dict(exposure=e,label=LAB[e],trait=tr,rg_contaminant=obs,
                    rg_vs_criteria=rg(tr,REF),expected_if_pure_AD=via,excess_rg=obs-via))
ix=pd.DataFrame(idx)
mv=pd.read_csv(f"{O}/Table3_mvmr.csv"); mv=mv[mv.outcome=="Asthma"]
mvmap=dict(zip(mv.definition,mv.mvmr_OR_definition))
NAME={'Self-report':'Self-report (UK Biobank)','Hospital ICD':'Hospital ICD (UK Biobank)',
      'Allergic composite':'Broad allergic composite','Pooled meta-analysis':'Pooled meta-analysis'}
ix['mvmr_conditional_OR']=ix.label.map(lambda l: mvmap.get(NAME.get(l,''),np.nan))
ix.to_csv(f"{O}/Table5_contamination_index.csv",index=False)
print(ix[['label','rg_contaminant','expected_if_pure_AD','excess_rg','mvmr_conditional_OR']].round(3).to_string(index=False))
s=ix.dropna(subset=['mvmr_conditional_OR'])
print(f"\nexcess rg vs MVMR conditional OR: Pearson r={np.corrcoef(s.excess_rg,np.log(s.mvmr_conditional_OR))[0,1]:.3f}, "
      f"Spearman rho={stats.spearmanr(s.excess_rg,s.mvmr_conditional_OR).statistic:.2f} (k={len(s)})")

# ---- outcome ladder ----
comp={'A':'T','T':'A','C':'G','G':'C'}
OUTS={"O4_FG_asthma":("Asthma","OUT_asthma"),"O5_FG_rhinitis":("Allergic rhinitis","OUT_rhinitis"),
      "O6_FG_copd":("COPD","OUT_copd"),"O7_FG_ibd":("IBD","OUT_ibd"),
      "O8_FG_migraine":("Migraine","OUT_migraine"),"O1_FG_depression":("Depression","OUT_depression"),
      "O2_FG_anxiety":("Anxiety","OUT_anxiety")}
def harm(e,o):
    m=e.merge(o,on='rsid',suffixes=('_e','_o'))
    sign=np.where(m.ea_e==m.ea_o,1,np.where(m.ea_e==m.ea_o.map(comp),1,
          np.where(m.ea_e==m.oa_o,-1,np.where(m.ea_e==m.oa_o.map(comp),-1,0))))
    m=m[sign!=0].copy(); m['by']=m.beta_o.values*sign[sign!=0]
    m=m[(m.eaf_e.notna())&(m.eaf_e.between(0.01,0.99))]
    pal=(m.ea_e+m.oa_e).isin(['AT','TA','CG','GC'])&m.eaf_e.between(0.42,0.58)
    return m[~pal]
scales={'E2_UKB_selfreport':None,'E3_UKB_ICD':None}
exp={}
for e in KEY:
    x=pd.read_parquet(f"{D}/{e}.parquet")
    cl=pd.read_csv(f"{ROOT}/mr/clump/{e}.clumps",sep=r"\s+")
    x=x[x.rsid.isin(cl['ID'].astype(str))].copy()
    if x['scale'].iloc[0]=='linear':
        v=pd.read_parquet(f"{D}/{e}.parquet").dropna(subset=['eaf','se','n'])
        v=v[(v.eaf>0.05)&(v.eaf<0.95)]
        mu=float(np.median((v.se**2)*v.n*2*v.eaf*(1-v.eaf))); x['beta']/=mu; x['se']/=mu
    exp[e]=x
rows=[]
for ok,(olab,ldn) in OUTS.items():
    o=pd.read_parquet(f"{D}/{ok}.parquet")
    for e in KEY:
        m=harm(exp[e],o)
        if len(m)<5: continue
        bx=m.beta_e.values; by=m.by.values; sy=m.se_o.values
        s=np.where(bx<0,-1,1); bx,by=bx*s,by*s
        w=1/sy**2; b=np.sum(w*bx*by)/np.sum(w*bx**2); se=np.sqrt(1/np.sum(w*bx**2))
        q=np.sum(w*(by-b*bx)**2); df=len(bx)-1
        if df>0 and q>df: se*=np.sqrt(q/df)
        rows.append(dict(outcome=olab,outcome_trait=ldn,exposure=e,label=LAB[e],n_snp=len(m),
                         b=b,se=se,OR=float(np.exp(b))))
L=pd.DataFrame(rows)
sum_=[]
for olab,g in L.groupby('outcome'):
    w=1/g.se**2; mu=np.sum(w*g.b)/np.sum(w)
    Q=float(np.sum(w*(g.b-mu)**2)); df=len(g)-1; I2=max(0,(Q-df)/Q)*100
    ldn=g.outcome_trait.iloc[0]
    sum_.append(dict(outcome=olab,rg_with_rhinitis=rg(ldn,'OUT_rhinitis'),
        pooled_OR=float(np.exp(mu)),I2=I2,Q_p=float(stats.chi2.sf(Q,df)),
        OR_min=float(np.exp(g.b.min())),OR_max=float(np.exp(g.b.max())),
        spread=float(np.exp(g.b.max()-g.b.min()))))
S=pd.DataFrame(sum_).sort_values('rg_with_rhinitis',ascending=False)
L.to_csv(f"{O}/Table6_outcome_ladder_full.csv",index=False)
S.to_csv(f"{O}/Table6_outcome_ladder.csv",index=False)
print("\n=== OUTCOME LADDER: does definitional spread scale with relatedness to the contaminant? ===")
print(S.round(3).to_string(index=False))
print(f"\nrg(outcome, rhinitis) vs definitional spread: Pearson r={np.corrcoef(S.rg_with_rhinitis,np.log(S.spread))[0,1]:.3f} (k={len(S)})")
