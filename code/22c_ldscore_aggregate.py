"""Per-chromosome LD scores from PLINK 2 r^2 output (called by 22b).
L2 = 1 + sum of bias-corrected r^2 within 1 Mb over ALL thinned SNPs; writes rsID (SNP), L2, CHR, BP."""
import os
ROOT=os.environ.get("ROOT","/tmp")
import sys,pandas as pd,numpy as np
c=sys.argv[1]; N=503
d=pd.read_csv(f"{ROOT}/ld/c{c}.vcor",sep="\t")
d.columns=[x.replace('#','') for x in d.columns]
a,b,r2=d['ID_A'],d['ID_B'],d['UNPHASED_R2']
r2a=r2-(1-r2)/(N-2)
s=pd.concat([pd.Series(r2a.values,index=a.values),pd.Series(r2a.values,index=b.values)]).groupby(level=0).sum()
bim=pd.read_csv(f"{ROOT}/ld/c{c}.bim",sep="\t",header=None,names=['CHR','SNP','CM','BP','A1','A2'])
out=pd.DataFrame({'ID':bim.SNP,'BP':bim.BP}).assign(L2=lambda x:x.ID.map(s).fillna(0)+1.0, CHR=int(c))
ids=pd.read_parquet(f"{ROOT}/ld/id2rsid.parquet").set_index('id').rsid   # 22d
out['SNP']=out.ID.map(ids)
out=out.dropna(subset=['SNP']).drop_duplicates('SNP')[['SNP','L2','CHR','BP']]   # LD scores use all SNPs; keep named ones
out.to_parquet(f"{ROOT}/ld/out2/ld_chr{c}.parquet",index=False)
print(f"chr{c}: {len(out):,} SNPs, mean L2 {out.L2.mean():.2f}, pairs {len(d):,}",flush=True)
