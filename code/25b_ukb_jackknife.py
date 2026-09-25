"""Exact delete-one-block jackknife for the contamination index.

The index  excess = rg(D,C) - rg(D,R) rg(R,C)  is a function of three cross-trait
regressions and three heritabilities, all estimated on the same SNPs and the same
200 jackknife blocks.  Propagating three independent standard errors ignores the
covariance between them; here the whole statistic is re-formed inside each
delete-one-block replicate, which is exact.
"""
import pandas as pd, numpy as np, json, os
LD=pd.read_parquet('/tmp/ld/ldscores.parquet').rename(columns={'SNP':'rsid'})
M=len(LD); D='/tmp/ukb/ss'; BL=200

def _blocks(n,seed=0):
    idx=np.arange(n); np.random.default_rng(seed).shuffle(idx)
    lab=np.empty(n,dtype=np.int32)
    for i,p in enumerate(np.array_split(idx,BL)): lab[p]=i
    return lab

def _stats(x,y,w,lab):
    return np.stack([np.bincount(lab,weights=v,minlength=BL) for v in
                     (w, w*x, w*x*x, w*y, w*x*y)])          # (5,BL)

def _slope(S):
    A,B,C,Dd,E=S
    return (A*E-B*Dd)/(A*C-B*B)

def _loo(S):
    T=S.sum(axis=1,keepdims=True)
    return _slope(T-S), _slope(T)[0] if False else _slope(T.ravel()[:,None]).ravel()[0]

class Store:
    def __init__(self,names):
        ldm=LD.drop_duplicates('rsid').set_index('rsid').L2
        self.rsid=ldm.index.values; self.l2=ldm.values.astype(float)
        self.w=1/np.maximum(self.l2,1); self.lab=_blocks(len(self.l2))
        self.z={}; self.neff={}
        for n in names:
            d=pd.read_parquet(f"{D}/{n}.parquet",columns=['rsid','z','n_eff'])
            z=pd.Series(d.z.values,index=d.rsid.values).reindex(self.rsid).values
            self.z[n]=np.where(np.abs(z)<30,z,np.nan); self.neff[n]=float(d.n_eff.iloc[0])
        self._h={}; self._r={}
    def h2S(self,n):
        if n not in self._h:
            z=self.z[n]; ok=np.isfinite(z)
            self._h[n]=_stats(self.l2[ok]*self.neff[n]/M, z[ok]**2, self.w[ok], self.lab[ok])
        return self._h[n]
    def rhoS(self,a,b):
        k=tuple(sorted((a,b)))
        if k not in self._r:
            za,zb=self.z[a],self.z[b]; ok=np.isfinite(za)&np.isfinite(zb)
            self._r[k]=_stats(self.l2[ok]*np.sqrt(self.neff[a]*self.neff[b])/M,
                              za[ok]*zb[ok], self.w[ok], self.lab[ok])
        return self._r[k]
    @staticmethod
    def _full_loo(S):
        T=S.sum(axis=1)
        return _slope(T), _slope(T[:,None]-S)
    def h2(self,n):
        f,l=self._full_loo(self.h2S(n)); return f,np.sqrt((BL-1)/BL*((l-l.mean())**2).sum())
    def rg(self,a,b):
        rf,rl=self._full_loo(self.rhoS(a,b))
        af,al=self._full_loo(self.h2S(a)); bf,bl=self._full_loo(self.h2S(b))
        g=rf/np.sqrt(abs(af*bf)); gl=rl/np.sqrt(np.abs(al*bl))
        return g,np.sqrt((BL-1)/BL*((gl-gl.mean())**2).sum())
    def excess(self,Dp,Rp,Cp,winsor=1.25):
        def gl(a,b):
            rf,rl=self._full_loo(self.rhoS(a,b))
            af,al=self._full_loo(self.h2S(a)); bf,bl=self._full_loo(self.h2S(b))
            return rf/np.sqrt(abs(af*bf)), rl/np.sqrt(np.abs(al*bl))
        dc,dcl=gl(Dp,Cp); dr,drl=gl(Dp,Rp); rc,rcl=gl(Rp,Cp)
        if winsor:
            cl=lambda v: np.clip(v,-winsor,winsor)
            e=cl(dc)-cl(dr)*cl(rc); el=cl(dcl)-cl(drl)*cl(rcl)
        else:
            e=dc-dr*rc; el=dcl-drl*rcl
        se=np.sqrt((BL-1)/BL*((el-el.mean())**2).sum())
        return e,se,dict(rg_DC=dc,rg_DR=dr,rg_RC=rc)
