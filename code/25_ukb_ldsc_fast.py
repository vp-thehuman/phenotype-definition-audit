"""Genetic correlations for every pair the scan needs (Table9).

Same WLS and 200-block jackknife as 13_ldsc.py, computed from block sufficient
statistics. Engine and contiguous jackknife blocks are in ldsc3.py."""
import os, sys, json, itertools
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ldsc3 import Store, ROOT
D=f"{ROOT}/ukb/ss"
if __name__=='__main__':
    cfg=json.load(open(f'{ROOT}/config.json'))
    names=[n for n in cfg['all'] if os.path.exists(f"{D}/{n}.parquet")]
    miss=[n for n in cfg['all'] if n not in names]
    print("loaded",len(names),"missing",miss,flush=True)
    S=Store(names)
    H={}
    for n in names:
        v,se=S.h2(n); H[n]={'h2':v,'se':se,'z':v/se if se>0 else np.nan}
    json.dump(H,open(f'{ROOT}/ukb/h2.json','w'),indent=1)
    print("h2 done",flush=True)
    need=set(); conts=[c for c in cfg['contaminants'] if c in names]
    for f,arms in cfg['families'].items():
        ps=[a['pheno'] for a in arms if a['pheno'] in names]
        for x,y in itertools.combinations(ps,2): need.add(tuple(sorted((x,y))))
        for x in ps:
            for c in conts:
                if c!=x: need.add(tuple(sorted((x,c))))
    for x,y in itertools.combinations(conts,2): need.add(tuple(sorted((x,y))))
    need=sorted(need); print("pairs:",len(need),flush=True)
    res=[]
    for i,(x,y) in enumerate(need):
        r,se,k=S.rho(x,y); res.append({'t1':x,'t2':y,'rho':r,'rho_se':se,'nsnp':k})
        if i%200==0: print(f"  {i}/{len(need)}",flush=True)
    R=pd.DataFrame(res)
    R['h1']=R.t1.map(lambda t:H[t]['h2']); R['h2_']=R.t2.map(lambda t:H[t]['h2'])
    den=np.sqrt(np.abs(R.h1*R.h2_))
    R['rg']=R.rho/den*np.sign(R.h1*R.h2_); R['rg_se']=R.rho_se/den
    R.to_csv(f'{ROOT}/ukb/rg_pairs.csv',index=False)
    print("RGDONE",len(R))
