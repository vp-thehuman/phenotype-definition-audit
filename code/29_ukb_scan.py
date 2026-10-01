import os
ROOT=os.environ.get("ROOT","/tmp")
import sys, json, itertools, numpy as np, pandas as pd
sys.path.insert(0,f'{ROOT}/ukb')
import ldsc3
cfg=json.load(open(f'{ROOT}/config.json'))
names=[n for n in cfg['all']]
S=ldsc3.Store(names)
H={n:S.h2(n) for n in names}
HZ_MAIN=6.0
RANK={'hospital ICD':0,'curated endpoint':1,'doctor-diagnosed Q':2,'self-report':3,'touchscreen composite':4}
CONT=cfg['contaminants']
NCTRL=361194.0

def prevalence(a): return (a['cases'] or 0)/NCTRL

def run(hz=HZ_MAIN, winsor=1.25, ref_mode='clinical', min_prev=0.0):
    rows=[]
    conts=[c for c in CONT if H[c][0]/H[c][1]>=hz]
    for fam,arms in cfg['families'].items():
        a=[x for x in arms if H[x['pheno']][0]/H[x['pheno']][1]>=hz and prevalence(x)>=min_prev]
        if len(a)<2: continue
        a=sorted(a,key=lambda x:(RANK[x['deftype']],-(x['cases'] or 0)))
        if ref_mode=='clinical': ref=a[0]
        elif ref_mode=='curated':
            cur=[x for x in a if x['deftype']=='curated endpoint']
            if not cur: continue
            ref=cur[0]
        elif ref_mode=='largest': ref=max(a,key=lambda x:(x['cases'] or 0))
        for Dd in a:
            if Dd['pheno']==ref['pheno']: continue
            rg_dr,se_dr=S.rg(Dd['pheno'],ref['pheno'])
            if se_dr>0.25: continue
            best=None; per={}
            for c in conts:
                if c in (Dd['pheno'],ref['pheno']): continue
                e,se,comp=S.excess(Dd['pheno'],ref['pheno'],c,winsor=winsor)
                per[c]=(e,se)
                if best is None or e>best[1]: best=(c,e,se)
            if best is None: continue
            dup = rg_dr>=0.99 and max(abs(v[0]) for v in per.values())<0.02
            rows.append(dict(family=fam,def_pheno=Dd['pheno'],def_type=Dd['deftype'],def_cases=Dd['cases'],
                def_desc=Dd['desc'][:70],def_prevalence=round(prevalence(Dd),4),
                ref_pheno=ref['pheno'],ref_type=ref['deftype'],ref_cases=ref['cases'],
                rg_def_ref=rg_dr,rg_def_ref_se=se_dr,
                worst_contaminant=CONT[best[0]],worst_cont_pheno=best[0],
                excess_max=best[1],excess_max_se=best[2],duplicate_definition=dup))
    return pd.DataFrame(rows)

def summarise(T,tag):
    A=T[~T.duplicate_definition]
    return dict(analysis=tag,families=int(T.family.nunique()),comparisons=len(T),
                negative_controls=int(T.duplicate_definition.sum()),informative=len(A),
                n_excess_gt_0_1=int((A.excess_max>0.1).sum()),
                n_excess_gt_0_2=int((A.excess_max>0.2).sum()),
                n_excess_gt_AD=int((A.excess_max>0.376).sum()),
                n_excess_z_gt_2=int(((A.excess_max/A.excess_max_se)>2).sum()),
                median_excess_selfreport=round(A[A.def_type=='self-report'].excess_max.median(),3),
                median_excess_curated=round(A[A.def_type=='curated endpoint'].excess_max.median(),3),
                median_rg_def_ref=round(A.rg_def_ref.median(),3),
                n_rg_below_0_8=int((A.rg_def_ref<0.8).sum()),
                min_rg_def_ref=round(A.rg_def_ref.min(),3),
                max_control_excess=round(T[T.duplicate_definition].excess_max.max(),4) if T.duplicate_definition.any() else None)

if __name__=='__main__':
    main=run()
    main.sort_values('excess_max',ascending=False).to_csv(f'{ROOT}/ukb/Table8_ukb_definition_scan.csv',index=False)
    sens=[summarise(main,'main: h2 z>=6, clinical reference, winsor 1.25')]
    for tag,kw in [('h2 z>=4',dict(hz=4.0)),('h2 z>=8',dict(hz=8.0)),
                   ('no winsorising',dict(winsor=0)),
                   ('reference = curated endpoint',dict(ref_mode='curated')),
                   ('reference = largest arm',dict(ref_mode='largest')),
                   ('case prevalence >= 1%',dict(min_prev=0.01))]:
        try: sens.append(summarise(run(**kw),tag))
        except Exception as e: sens.append(dict(analysis=tag,families=f"ERROR {e}"))
    Ssum=pd.DataFrame(sens)
    Ssum.to_csv(f'{ROOT}/ukb/Table10_sensitivity.csv',index=False)
    pd.set_option('display.width',260)
    print(Ssum.to_string(index=False))
