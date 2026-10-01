#!/usr/bin/env python3
"""Standard errors for the depression-family correlations, to see which
differences are real and which are noise."""
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
REF="DEP_cidi"; CON="CON_neuroticism"
rc,rcse=RGse(REF,CON)
print(f"rg(CIDI, neuroticism) = {rc:.3f} ({rcse:.3f})\n")
print(f"{'definition':40s}{'rg vs CIDI':>18s}{'rg vs neuroticism':>22s}{'excess':>9s}{'~se':>8s}")
rows=[]
for n in ["DEP_1sym","DEP_2sym","DEP_3sym","DEP_45sym"]:
    ref,refse=RGse(n,REF); obs,obsse=RGse(n,CON)
    exc=obs-ref*rc
    # delta-method-ish: dominated by the two rg terms
    se=float(np.sqrt(obsse**2 + (rc*refse)**2 + (ref*rcse)**2))
    rows.append(dict(definition=n,label=LAB[n],rg_vs_reference=ref,rg_ref_se=refse,
                     rg_contaminant=obs,rg_con_se=obsse,excess_rg=exc,excess_se=se,
                     n_cases={"DEP_1sym":57321,"DEP_2sym":21468,"DEP_3sym":9738,"DEP_45sym":4887}[n]))
    print(f"{LAB[n]:40s}{ref:>11.3f} ({refse:.3f}){obs:>14.3f} ({obsse:.3f}){exc:>9.3f}{se:>8.3f}")
t=pd.DataFrame(rows)
t['rg_cidi_neuroticism']=rc; t['rg_cidi_neuroticism_se']=rcse
t.to_csv(f"{ROOT}/mr/out/Table7_depression_family_index.csv",index=False)
print("\nall four touchscreen definitions vs the CIDI reference:")
print(f"  mean excess = {t.excess_rg.mean():.3f}; range {t.excess_rg.min():.3f} to {t.excess_rg.max():.3f}")
print(f"  ordering by strictness? Spearman(cases, excess) = "
      f"{t[['n_cases','excess_rg']].corr(method='spearman').iloc[0,1]:.2f}")
