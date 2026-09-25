import sys,os,gzip,json,subprocess,pandas as pd,numpy as np
shard=int(sys.argv[1]); nsh=int(sys.argv[2])
W=pd.read_parquet('/tmp/ukb/whitelist.parquet')
wl={v:(r,a) for v,r,a in zip(W.variant,W.rsid,W.AF)}
B="https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs"
cfg=json.load(open('/tmp/config.json'))
todo=[p for i,p in enumerate(cfg['all']) if i%nsh==shard]
os.makedirs('/tmp/ukb/ss',exist_ok=True)

def fetch(p):
    url=f"{B}/{p}.gwas.imputed_v3.both_sexes.tsv.bgz"
    proc=subprocess.Popen(['curl','-sL','--retry','5','--retry-all-errors','-m','3600',url],stdout=subprocess.PIPE)
    try:
        f=gzip.open(proc.stdout,'rt',errors='replace')
        h=f.readline().rstrip('\n').split('\t'); ix={c:i for i,c in enumerate(h)}
        iv,ib,ise=ix['variant'],ix['beta'],ix['se']
        ilc=ix.get('low_confidence_variant')
        rec=[]
        for line in f:
            q=line.rstrip('\n').split('\t')
            m=wl.get(q[iv])
            if m is None: continue
            if ilc is not None and q[ilc]=='true': continue
            try: b=float(q[ib]); sd=float(q[ise])
            except: continue
            if not (sd>0): continue
            rec.append((m[0],m[1],b,sd))
    finally:
        try: proc.stdout.close()
        except: pass
        proc.wait()
    if len(rec)<50000: raise RuntimeError(f"too few rows {len(rec)}")
    return rec

for p in todo:
    out=f'/tmp/ukb/ss/{p}.parquet'
    if os.path.exists(out): continue
    rec=None
    for att in range(4):
        try:
            rec=fetch(p); break
        except Exception as e:
            print(f"  retry {p} #{att}: {e}",flush=True)
    if rec is None:
        print("FAILED",p,flush=True); continue
    d=pd.DataFrame(rec,columns=['rsid','AF','beta','se'])
    d['z']=d.beta/d.se
    d['n_eff']=float(np.nanmedian(1.0/(2*d.AF*(1-d.AF)*d.se**2)))
    d=d[np.isfinite(d.z)]
    d[['rsid','AF','beta','se','z','n_eff']].to_parquet(out+'.tmp',index=False)
    os.replace(out+'.tmp',out)
    print(f"{p:26s} snps={len(d):>7,} neff={d.n_eff.iloc[0]:>12,.0f}",flush=True)
print("SHARD",shard,"DONE",flush=True)
