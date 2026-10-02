#!/bin/bash
# LD scores from 1000 Genomes EUR: MAF >= 0.05, a uniform 6.5% sample of SNPs per
# chromosome (--thin, a FRACTION, not --thin-count), 1 Mb window, bias-corrected r^2.
ROOT=${ROOT:-/tmp}; export ROOT
HERE=$(cd "$(dirname "$0")" && pwd)
PLINK2=${PLINK2:-$(command -v plink2 || echo ${ROOT}/bin/plink2)}
python3 "$HERE/22d_variant_ids.py"          # chr:pos:ref:alt -> rsID (VCFs carry no rsIDs)
cd ${ROOT}/ld
mkdir -p out2
proc(){ c=$1
  $PLINK2 --vcf chr$c.vcf.gz --keep eur.id --snps-only --max-alleles 2 --maf 0.05 \
    --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 100 --rm-dup exclude-all --memory 2500 --make-bed --out f$c --silent || return 1
  rm -f chr$c.vcf.gz
  $PLINK2 --bfile f$c --thin 0.065 --seed 42 --memory 2500 --make-bed --out c$c --silent || return 1
  rm -f f$c.bed f$c.bim f$c.fam f$c.log
  $PLINK2 --bfile c$c --r2-unphased --ld-window-kb 1000 --ld-window 999999 \
    --ld-window-r2 0 --memory 2500 --out c$c --silent || return 1
  python3 "$HERE/22c_ldscore_aggregate.py" $c || return 1
  rm -f c$c.vcor c$c.log c$c.bed c$c.bim c$c.fam
}
for c in $(seq 1 22); do
  [ -f out2/ld_chr$c.parquet ] && continue
  [ -f chr$c.vcf.gz ] || { echo "chr$c.vcf.gz missing: re-run 22a_1000g_fetch.sh" >&2; exit 1; }
  proc $c || { echo "LD scores failed on chr$c" >&2; exit 1; }
done
python3 - <<'PY'
import os, glob, pandas as pd
ROOT=os.environ["ROOT"]
d=pd.concat([pd.read_parquet(f) for f in glob.glob(f"{ROOT}/ld/out2/ld_chr*.parquet")],ignore_index=True)
d=d.sort_values(['CHR','BP']).reset_index(drop=True)
d.to_parquet(f"{ROOT}/ld/ldscores.parquet",index=False)
print("LD scores:",len(d),"SNPs; mean L2 by chromosome",round(d.groupby('CHR').L2.mean().min(),1),"to",round(d.groupby('CHR').L2.mean().max(),1))
PY
touch ${ROOT}/ld/ld2.flag
