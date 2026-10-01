"""Fast equivalent of code/13_ldsc.py: same WLS + 200-block jackknife,
computed from block sufficient statistics instead of refitting."""
import os
ROOT=os.environ.get("ROOT","/tmp")
import pandas as pd, numpy as np, json, os, itertools
LD=pd.read_parquet(f'{ROOT}/ld/ldscores.parquet').rename(columns={'SNP':'rsid'})
M=len(LD)
D=f'{ROOT}/ukb/ss'
BL=200
def blocks(n,seed=0):
    idx=np.arange(n); np.random.default_rng(seed).shuffle(idx)
    lab=np.empty(n,dtype=np.int32)
    for i,p in enumerate(np.array_split(idx,BL)): lab[p]=i
    return lab
def jk_wls(x,y,w,lab):
    sw=np.bincount(lab,weights=w,minlength=BL)
    swx=np.bincount(lab,weights=w*x,minlength=BL)
    swxx=np.bincount(lab,weights=w*x*x,minlength=BL)
    swy=np.bincount(lab,weights=w*y,minlength=BL)
    swxy=np.bincount(lab,weights=w*x*y,minlength=BL)
    T=[a.sum() for a in (sw,swx,swxx,swy,swxy)]
    def slope(A,B,C,Dd,E):
        det=A*C-B*B
        return (A*E-B*Dd)/det
    full=slope(*T)
    ests=np.array([slope(T[0]-sw[i],T[1]-swx[i],T[2]-swxx[i],T[3]-swy[i],T[4]-swxy[i]) for i in range(BL)])
    se=np.sqrt((BL-1)/BL*np.sum((ests-ests.mean())**2))
    return full,se
class Store:
    def __init__(self,names):
        ldm=LD.drop_duplicates('rsid').set_index('rsid').L2
        self.rsid=ldm.index.values
        self.l2=ldm.values.astype(float)
        self.z={}; self.neff={}
        for n in names:
            d=pd.read_parquet(f"{D}/{n}.parquet",columns=['rsid','z','n_eff'])
            z=pd.Series(d.z.values,index=d.rsid.values).reindex(self.rsid).values
            self.z[n]=np.where(np.abs(z)<30,z,np.nan)
            self.neff[n]=float(d.n_eff.iloc[0])
        self.w=1/np.maximum(self.l2,1)
        self.lab=blocks(len(self.l2))
    def h2(self,n):
        z=self.z[n]; ok=np.isfinite(z)
        x=self.l2[ok]*self.neff[n]/M; y=z[ok]**2
        return jk_wls(x,y,self.w[ok],self.lab[ok])
    def rho(self,a,b):
        za,zb=self.z[a],self.z[b]; ok=np.isfinite(za)&np.isfinite(zb)
        x=self.l2[ok]*np.sqrt(self.neff[a]*self.neff[b])/M; y=za[ok]*zb[ok]
        r,se=jk_wls(x,y,self.w[ok],self.lab[ok]); return r,se,int(ok.sum())
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
