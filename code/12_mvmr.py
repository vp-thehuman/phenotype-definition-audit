#!/usr/bin/env python3
"""Multivariable MR: does a loose AD definition carry an effect on asthma that is
NOT explained by clinician-criteria AD liability?

If a loose definition is merely a noisier measure of the same AD liability, then
conditioning on clinician-criteria AD should absorb its effect entirely. If it
additionally imports non-AD atopic signal, it retains an independent effect.
That is the difference between attenuation (cancels in an MR ratio) and
non-specificity (does not).
"""
import os
ROOT=os.environ.get("ROOT","/tmp")
import pandas as pd, numpy as np, subprocess, shutil, tempfile
from scipy import stats
# v0.2: the union of the two instrument lists is re-clumped (r2 < 0.001, 10 Mb), because
# SNPs from two separately clumped lists can be in LD with each other. Set RECLUMP=0 to
# reproduce v0.1. Conditional F is now the Sanderson-Windmeijer statistic, which includes
# the reference exposure's standard errors; the v0.1 value is kept as conditional_F_v01.
RECLUMP=os.environ.get("RECLUMP","1")=="1"
PLINK2=os.environ.get("PLINK2") or shutil.which("plink2") or f"{ROOT}/plink2"
LDREF=os.environ.get("LDREF",f"{ROOT}/ld/EUR")

def reclump(ids, *tables):
    """Clump the union of instrument lists on the smallest p across exposures."""
    p=pd.concat([t[t.rsid.isin(ids)][['rsid','p']] for t in tables]).groupby('rsid').p.min()
    tmp=tempfile.mkdtemp()
    p.rename('P').rename_axis('ID').reset_index().to_csv(f"{tmp}/u.assoc",sep="\t",index=False)
    subprocess.run([PLINK2,"--bfile",LDREF,"--clump",f"{tmp}/u.assoc","--clump-p1","1",
                    "--clump-r2","0.001","--clump-kb","10000","--clump-id-field","ID",
                    "--clump-p-field","P","--out",f"{tmp}/u","--silent"],check=False)
    try: kept=set(pd.read_csv(f"{tmp}/u.clumps",sep=r"\s+")['ID'].astype(str))
    except FileNotFoundError: kept=set(ids)
    shutil.rmtree(tmp,ignore_errors=True)
    return kept

def sw_condF(bj,sj,bk,sk):
    """Sanderson-Windmeijer conditional F for exposure j given exposure k."""
    d=np.sum(bj*bk)/np.sum(bk**2)
    Q=np.sum((bj-d*bk)**2/(sj**2+d**2*sk**2))
    return float(Q/max(len(bj)-1,1))
D=f"{ROOT}/mr/data"
REF="X1_EAGLE"; REFLAB="Clinician criteria (EAGLE)"
TEST={"X2_UKB_selfreport":"UKB hayfever/rhinitis/eczema","X3_UKB_ICD":"Hospital ICD (UK Biobank)",
      "X4_allergic_broad":"Broad allergic composite","X5_BUDU":"Pooled meta-analysis"}
CLUMP={"X1_EAGLE":"E1_EAGLE_criteria","X2_UKB_selfreport":"E2_UKB_selfreport",
       "X3_UKB_ICD":"E3_UKB_ICD","X4_allergic_broad":"E4_allergic_broad","X5_BUDU":"E5_BUDU_meta"}
comp={'A':'T','T':'A','C':'G','G':'C'}
OUTS={"O4_FG_asthma":"Asthma","O1_FG_depression":"Depression"}

def align(base, other, bname, sname):
    """other: rsid, ea, oa, beta, se -> adds bname/sname aligned to base's effect allele."""
    o=other.rename(columns={'ea':'_ea','oa':'_oa','beta':'_b','se':'_s'})
    m=base.merge(o,on='rsid')
    sign=np.where(m.ea==m._ea,1,np.where(m.ea==m._ea.map(comp),1,
          np.where(m.ea==m._oa,-1,np.where(m.ea==m._oa.map(comp),-1,0))))
    m=m[sign!=0].copy()
    m[bname]=m['_b'].values*sign[sign!=0]; m[sname]=m['_s'].values
    return m.drop(columns=['_ea','_oa','_b','_s'])

rows=[]
for ok,olab in OUTS.items():
    o=pd.read_parquet(f"{D}/{ok}.parquet")
    ref=pd.read_parquet(f"{D}/{REF}.parquet")
    for tk,tlab in TEST.items():
        t=pd.read_parquet(f"{D}/{tk}.parquet")
        ivs=set(pd.read_csv(f"{ROOT}/mr/clump/{CLUMP[REF]}.clumps",sep=r"\s+")['ID'].astype(str)) | \
            set(pd.read_csv(f"{ROOT}/mr/clump/{CLUMP[tk]}.clumps",sep=r"\s+")['ID'].astype(str))
        n_union=len(ivs)
        if RECLUMP: ivs=reclump(ivs,ref,t)
        base=ref[ref.rsid.isin(ivs)][['rsid','ea','oa','eaf','beta','se']].rename(columns={'beta':'bref','se':'sref'})
        m=align(base,t[['rsid','ea','oa','beta','se']],'btest','setest')
        m=align(m,o[['rsid','ea','oa','beta','se']],'by','sout')
        m=m.dropna(subset=['bref','btest','by','sout']).query("sout>0")
        if len(m)<10: continue
        X=np.column_stack([m.bref.values,m.btest.values]); y=m.by.values; W=np.diag(1/m.sout.values**2)
        XtWX=X.T@W@X
        try: b=np.linalg.solve(XtWX,X.T@W@y)
        except np.linalg.LinAlgError: continue
        res=y-X@b; df=len(m)-2
        phi=max((res@W@res)/df,1.0)          # overdispersion
        cov=np.linalg.inv(XtWX)*phi
        se=np.sqrt(np.diag(cov))
        # conditional F for the test exposure (Sanderson-Windmeijer style)
        r=m.btest.values-m.bref.values*(np.sum(m.bref*m.btest)/np.sum(m.bref**2))
        condF_v01=float(np.sum((r/m.setest.values)**2)/max(len(m)-1,1))
        condF=sw_condF(m.btest.values,m.setest.values,m.bref.values,m.sref.values)
        condF_ref=sw_condF(m.bref.values,m.sref.values,m.btest.values,m.setest.values)
        # univariable comparison
        wu=1/m.sout.values**2
        uni=np.sum(wu*m.btest*m.by)/np.sum(wu*m.btest**2)
        rows.append(dict(outcome=olab,definition=tlab,n_snp=len(m),
            uni_OR=float(np.exp(uni)),
            mvmr_OR_definition=float(np.exp(b[1])),mvmr_lo=float(np.exp(b[1]-1.96*se[1])),
            mvmr_hi=float(np.exp(b[1]+1.96*se[1])),mvmr_p=float(2*stats.norm.sf(abs(b[1]/se[1]))),
            mvmr_OR_criteria=float(np.exp(b[0])),criteria_p=float(2*stats.norm.sf(abs(b[0]/se[0]))),
            conditional_F=condF,conditional_F_reference=condF_ref,conditional_F_v01=condF_v01,
            n_snp_union_before_reclump=n_union,reclumped=RECLUMP))
        print(f"{olab:11s} {tlab:28s} n={len(m):>4} condF={condF:6.1f} | univariable OR={np.exp(uni):.3f} "
              f"-> conditional on criteria: OR={np.exp(b[1]):.3f} "
              f"({np.exp(b[1]-1.96*se[1]):.3f}-{np.exp(b[1]+1.96*se[1]):.3f}) p={2*stats.norm.sf(abs(b[1]/se[1])):.2g}",flush=True)
t=pd.DataFrame(rows); t.to_csv(f"{ROOT}/mr/out/Table3_mvmr.csv",index=False)
