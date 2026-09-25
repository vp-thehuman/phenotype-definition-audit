import sys,pandas as pd,numpy as np
c=sys.argv[1]; N=503
d=pd.read_csv(f"/tmp/ld/c{c}.vcor",sep="\t")
d.columns=[x.replace('#','') for x in d.columns]
a,b,r2=d['ID_A'],d['ID_B'],d['UNPHASED_R2']
r2a=r2-(1-r2)/(N-2)
s=pd.concat([pd.Series(r2a.values,index=a.values),pd.Series(r2a.values,index=b.values)]).groupby(level=0).sum()
bim=pd.read_csv(f"/tmp/ld/c{c}.bim",sep="\t",header=None,names=['CHR','SNP','CM','BP','A1','A2'])
out=pd.DataFrame({'SNP':bim.SNP}).assign(L2=lambda x:x.SNP.map(s).fillna(0)+1.0, CHR=int(c))
out.to_parquet(f"/tmp/ld/out/ld_chr{c}.parquet",index=False)
print(f"chr{c}: {len(out):,} SNPs, mean L2 {out.L2.mean():.2f}, pairs {len(d):,}",flush=True)
