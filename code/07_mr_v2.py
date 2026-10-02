#!/usr/bin/env python3
"""MR v2. Adds: MAF filter, Steiger filtering, an asthma positive control,
explicit power / minimum detectable effect, and a correlation-aware treatment of
the between-definition comparison.
"""
import os
ROOT=os.environ.get("ROOT","/tmp")
import pandas as pd, numpy as np, json, os
from scipy import stats
D=f"{ROOT}/mr/data"; OUT=f"{ROOT}/mr/out"
MAF_MIN=0.01
EXPS={"E1_EAGLE_criteria":"Clinician criteria (EAGLE)","E2_UKB_selfreport":"Self-report (UK Biobank)",
      "E3_UKB_ICD":"Hospital ICD record (UK Biobank)","E4_allergic_broad":"Broad allergic composite",
      "E5_BUDU_meta":"Pooled meta-analysis (2023)"}
OUTC={"O1_FG_depression":"Depression (FinnGen R11)","O2_FG_anxiety":"Anxiety (FinnGen R11)",
      "O3_PGC_MDD":"MDD (PGC 2018, ex-23andMe)","O4_FG_asthma":"Asthma (FinnGen R11) [positive control]"}
comp={'A':'T','T':'A','C':'G','G':'C'}
COMPSET=set('AT TA CG GC'.split())

def eff_n(d):
    """Effective sample size from logistic SEs: N_eff = 1/(2f(1-f)SE^2)."""
    v=d.dropna(subset=['eaf','se']); v=v[(v.eaf>0.05)&(v.eaf<0.95)&(v.se>0)]
    return float(np.median(1/(2*v.eaf*(1-v.eaf)*v.se**2)))

def r2_from_beta(beta,eaf,n):
    """Variance in liability explained by a SNP (log-odds scale approximation)."""
    num=2*eaf*(1-eaf)*beta**2
    return num/(num+2*eaf*(1-eaf)*(1/np.sqrt(np.maximum(n,1)))**2*n)

def harmonise(e,o):
    m=e.merge(o,on='rsid',suffixes=('_e','_o')); keep=[]
    for _,r in m.iterrows():
        ea_e,oa_e,ea_o,oa_o=r.ea_e,r.oa_e,r.ea_o,r.oa_o; b_o,eaf_o=r.beta_o,r.eaf_o
        if {ea_e,oa_e}=={ea_o,oa_o}:
            if ea_e!=ea_o: b_o=-b_o; eaf_o=1-eaf_o if pd.notna(eaf_o) else eaf_o
        else:
            try: ea2,oa2=comp[ea_o],comp[oa_o]
            except KeyError: continue
            if {ea_e,oa_e}=={ea2,oa2}:
                if ea_e!=ea2: b_o=-b_o; eaf_o=1-eaf_o if pd.notna(eaf_o) else eaf_o
            else: continue
        f=r.eaf_e
        if pd.isna(f) or min(f,1-f)<MAF_MIN: continue          # MAF filter
        if (ea_e+oa_e) in COMPSET and 0.42<f<0.58: continue     # ambiguous palindrome
        keep.append(dict(rsid=r.rsid,bx=r.beta_e,sx=r.se_e,by=b_o,sy=r.se_o,eaf=f))
    return pd.DataFrame(keep)

def ivw(bx,sx,by,sy):
    w=1/sy**2; b=np.sum(w*bx*by)/np.sum(w*bx**2); se=np.sqrt(1/np.sum(w*bx**2))
    q=np.sum(w*(by-b*bx)**2); df=len(bx)-1
    if df>0 and q>df: se*=np.sqrt(q/df)
    return b,se,q,df

def wmedian(bx,by,sx,sy,nb=1000,seed=1):
    ratio=by/bx; w=(sy**2/bx**2)**-1
    def wm(v,w):
        o=np.argsort(v); v,w=v[o],w[o]; cw=np.cumsum(w)-0.5*w; cw/=w.sum()
        k=np.searchsorted(cw,0.5)
        return v[0] if k==0 else v[k-1]+(v[k]-v[k-1])*(0.5-cw[k-1])/(cw[k]-cw[k-1])
    rng=np.random.default_rng(seed)
    return wm(ratio,w),float(np.std([wm(rng.normal(by,sy)/rng.normal(bx,sx),w) for _ in range(nb)]))

def egger(bx,by,sy):
    s=np.sign(bx); bx2,by2=bx*s,by*s
    W=np.diag(1/sy**2); X=np.column_stack([np.ones(bx2.size),bx2])
    XtWX=X.T@W@X; beta=np.linalg.solve(XtWX,X.T@W@by2)
    res=by2-X@beta; df=bx2.size-2
    sig2=(res@W@res)/df if df>0 else np.nan
    cov=np.linalg.inv(XtWX)*max(sig2,1.0)
    return beta[1],np.sqrt(cov[1,1]),beta[0],np.sqrt(cov[0,0]),df

# scale conversion for linear-scale exposures
exp_d={}; scales={}
for ek in EXPS:
    x=pd.read_parquet(f"{D}/{ek}.parquet")
    cl=pd.read_csv(f"{ROOT}/mr/clump/{ek}.clumps",sep=r"\s+")
    x=x[x.rsid.isin(cl['ID'].astype(str))].copy()
    if x['scale'].iloc[0]=='linear':
        v=pd.read_parquet(f"{D}/{ek}.parquet").dropna(subset=['eaf','se','n'])
        v=v[(v.eaf>0.05)&(v.eaf<0.95)]
        mu=float(np.median((v.se**2)*v.n*2*v.eaf*(1-v.eaf))); scales[ek]=mu
        x['beta']/=mu; x['se']/=mu
    exp_d[ek]=x

out_d={ok:pd.read_parquet(f"{D}/{ok}.parquet") for ok in OUTC}
outN={ok:eff_n(o) for ok,o in out_d.items()}
print("effective outcome N (from SEs):",{k:int(v) for k,v in outN.items()},"\n")

rows=[]
for ek,elab in EXPS.items():
    nx=eff_n(exp_d[ek])
    for ok,olab in OUTC.items():
        h=harmonise(exp_d[ek],out_d[ok])
        if len(h)<3: continue
        h=h[h.bx!=0].copy()
        # Steiger: keep SNPs explaining more variance in exposure than outcome
        r2x=2*h.eaf*(1-h.eaf)*h.bx**2; r2y=2*h.eaf*(1-h.eaf)*h.by**2
        n_fail=int((r2y>=r2x).sum()); h=h[r2x>r2y]
        bx,sx,by,sy=h.bx.values,h.sx.values,h.by.values,h.sy.values
        s=np.where(bx<0,-1,1); bx,by=bx*s,by*s
        b,se,q,df=ivw(bx,sx,by,sy); wm,wmse=wmedian(bx,by,sx,sy)
        eb,ebse,ei,eise,edf=egger(bx,by,sy)
        rows.append(dict(exposure=ek,exposure_label=elab,outcome=ok,outcome_label=olab,
            n_snp=len(h),n_steiger_dropped=n_fail,F_mean=float(np.mean((bx/sx)**2)),
            ivw_b=b,ivw_se=se,ivw_p=2*stats.norm.sf(abs(b/se)),
            ivw_or=np.exp(b),ivw_lo=np.exp(b-1.96*se),ivw_hi=np.exp(b+1.96*se),
            mdOR_80pct=float(np.exp(2.80*se)),
            Q=q,Q_df=df,Q_p=stats.chi2.sf(q,df) if df>0 else np.nan,
            wm_or=float(np.exp(wm)),wm_p=2*stats.norm.sf(abs(wm/wmse)) if wmse>0 else np.nan,
            egger_or=float(np.exp(eb)),egger_int_p=2*stats.t.sf(abs(ei/eise),edf) if edf>0 and eise>0 else np.nan))
        print(f"{elab[:31]:33s} -> {olab[:34]:36s} n={len(h):>4}(-{n_fail}) OR={np.exp(b):.3f} "
              f"({np.exp(b-1.96*se):.3f}-{np.exp(b+1.96*se):.3f}) p={2*stats.norm.sf(abs(b/se)):.2g}",flush=True)

res=pd.DataFrame(rows); res.to_csv(f"{OUT}/mr_results_v2.csv",index=False)

print("\n=== BETWEEN-DEFINITION HETEROGENEITY ===")
het=[]
for ok,olab in OUTC.items():
    s=res[res.outcome==ok]
    if len(s)<2: continue
    w=1/s.ivw_se**2; mu=np.sum(w*s.ivw_b)/np.sum(w)
    Q=float(np.sum(w*(s.ivw_b-mu)**2)); df=len(s)-1
    I2=max(0,(Q-df)/Q)*100 if Q>0 else 0
    # independent contrast: EAGLE (no UKB/FinnGen/23andMe) vs UKB self-report
    a=s[s.exposure=="E1_EAGLE_criteria"]; b_=s[s.exposure=="E2_UKB_selfreport"]
    if len(a) and len(b_):
        diff=float(a.ivw_b.iloc[0]-b_.ivw_b.iloc[0])
        sed=float(np.sqrt(a.ivw_se.iloc[0]**2+b_.ivw_se.iloc[0]**2))
        pdiff=2*stats.norm.sf(abs(diff/sed)); mdr=float(np.exp(2.80*sed))
    else: diff=sed=pdiff=mdr=np.nan
    het.append(dict(outcome=ok,outcome_label=olab,k=len(s),pooled_OR=float(np.exp(mu)),
        Q=Q,df=df,Q_p=stats.chi2.sf(Q,df),I2=I2,
        OR_min=float(np.exp(s.ivw_b.min())),OR_max=float(np.exp(s.ivw_b.max())),
        max_min_ratio=float(np.exp(s.ivw_b.max()-s.ivw_b.min())),
        indep_contrast_OR_ratio=float(np.exp(diff)),indep_contrast_p=pdiff,
        min_detectable_OR_ratio_80pct=mdr))
    print(f"{olab[:36]:38s} k={len(s)} pooledOR={np.exp(mu):.3f} I2={I2:>3.0f}% Qp={stats.chi2.sf(Q,df):.2f} "
          f"range={np.exp(s.ivw_b.min()):.2f}-{np.exp(s.ivw_b.max()):.2f} | EAGLE/UKB ratio={np.exp(diff):.3f} "
          f"(p={pdiff:.2f}, detectable ratio >{mdr:.2f})")
pd.DataFrame(het).to_csv(f"{OUT}/heterogeneity_v2.csv",index=False)
json.dump(scales,open(f"{OUT}/scale_factors.json","w"),indent=2)
