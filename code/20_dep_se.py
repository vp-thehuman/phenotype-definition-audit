#!/usr/bin/env python3
"""Standard errors for the depression-family correlations and contamination index.

excess_se is now the exact delete-one-block jackknife of the whole index (all six
regressions re-fitted on the SNPs common to D, R and C, contiguous blocks), as in
the UK Biobank scan. v0.1 reported only excess_se_approx, which treats the three
rg terms as independent; both are written."""
import os
ROOT=os.environ.get("ROOT","/tmp")
import pandas as pd, numpy as np
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'19_family_index.py')).read().split('dat={n:load(n)')[0])
dat={n:load(n) for n in NAMES}
H={}; HS={}
for n in NAMES:
    v,s=h2(dat[n]); H[n]=v; HS[n]=s
def RGse(a,b):
    r,rs=rg(dat[a],dat[b]); den=np.sqrt(max(H[a],1e-9)*max(H[b],1e-9))
    return float(np.clip(r/den,-1,1)), float(rs/den)
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ldsc3 import contiguous_blocks, block_sums, full_and_loo, jk_se

def excess_exact(Dn, Rn, Cn, clip=1.0):
    """Exact jackknife SE of rg(D,C) - rg(D,R) rg(R,C) on the common SNP set."""
    keep=['rsid','ea','oa','z','n_eff','L2']+[c for c in ('CHR','BP') if c in dat[Dn].columns]
    m=dat[Dn][keep]
    for n,suf in ((Rn,'_R'),(Cn,'_C')):
        o=dat[n][['rsid','ea','oa','z','n_eff']].rename(columns=lambda c: c if c=='rsid' else c+suf)
        m=m.merge(o,on='rsid')
        sg=np.where(m.ea==m['ea'+suf],1,np.where(m.ea==m['ea'+suf].map(comp),1,
           np.where(m.ea==m['oa'+suf],-1,np.where(m.ea==m['oa'+suf].map(comp),-1,0))))
        m=m[sg!=0].copy(); m['z'+suf]=m['z'+suf].values*sg[sg!=0]
    m=m.sort_values([c for c in ('CHR','BP') if c in m.columns] or ['rsid'])
    l=m.L2.values; w=1/np.maximum(l,1); lab=contiguous_blocks(len(m))
    Z={'D':m.z.values,'R':m.z_R.values,'C':m.z_C.values}
    N={'D':float(m.n_eff.iloc[0]),'R':float(m.n_eff_R.iloc[0]),'C':float(m.n_eff_C.iloc[0])}
    h={k:full_and_loo(block_sums(l*N[k]/M,Z[k]**2,w,lab)) for k in Z}
    def g(a,b):
        rf,rl=full_and_loo(block_sums(l*np.sqrt(N[a]*N[b])/M,Z[a]*Z[b],w,lab))
        return (np.clip(rf/np.sqrt(abs(h[a][0]*h[b][0])),-clip,clip),
                np.clip(rl/np.sqrt(np.abs(h[a][1]*h[b][1])),-clip,clip))
    dc,dr,rcx=g('D','C'),g('D','R'),g('R','C')
    return float(dc[0]-dr[0]*rcx[0]), jk_se(dc[1]-dr[1]*rcx[1])

REF="DEP_cidi"; CON="CON_neuroticism"
rc,rcse=RGse(REF,CON)
print(f"rg(CIDI, neuroticism) = {rc:.3f} ({rcse:.3f})\n")
print(f"{'definition':40s}{'rg vs CIDI':>18s}{'rg vs neuroticism':>22s}{'excess':>9s}{'se':>8s}{'se v0.1':>9s}")
rows=[]
for n in ["DEP_1sym","DEP_2sym","DEP_3sym","DEP_45sym"]:
    ref,refse=RGse(n,REF); obs,obsse=RGse(n,CON)
    exc=obs-ref*rc
    # v0.1 approximation: treats the three rg terms as independent
    se_approx=float(np.sqrt(obsse**2 + (rc*refse)**2 + (ref*rcse)**2))
    exc_common,se=excess_exact(n,REF,CON)
    rows.append(dict(definition=n,label=LAB[n],rg_vs_reference=ref,rg_ref_se=refse,
                     rg_contaminant=obs,rg_con_se=obsse,excess_rg=exc,excess_se=se,
                     excess_se_approx=se_approx,excess_rg_common_snps=exc_common,
                     n_cases={"DEP_1sym":57321,"DEP_2sym":21468,"DEP_3sym":9738,"DEP_45sym":4887}[n]))
    print(f"{LAB[n]:40s}{ref:>11.3f} ({refse:.3f}){obs:>14.3f} ({obsse:.3f}){exc:>9.3f}{se:>8.3f}{se_approx:>9.3f}")
t=pd.DataFrame(rows)
t['rg_cidi_neuroticism']=rc; t['rg_cidi_neuroticism_se']=rcse
t.to_csv(f"{ROOT}/mr/out/Table7_depression_family_index.csv",index=False)
print("\nall four touchscreen definitions vs the CIDI reference:")
print(f"  mean excess = {t.excess_rg.mean():.3f}; range {t.excess_rg.min():.3f} to {t.excess_rg.max():.3f}")
print(f"  ordering by strictness? Spearman(cases, excess) = "
      f"{t[['n_cases','excess_rg']].corr(method='spearman').iloc[0,1]:.2f}")
