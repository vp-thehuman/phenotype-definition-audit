#!/bin/bash
# Genome-wide 1000 Genomes phase 3 EUR reference in PLINK format ($ROOT/ld/EUR), named by
# rsID. Used for clumping (02, 12) and for the AD-family LD scores (10). Needs
# 22d_variant_ids.py (rsID map) and the EUR sample list written by 22a.
ROOT=${ROOT:-/tmp}; export ROOT
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
PLINK2=${PLINK2:-$(command -v plink2 || echo ${ROOT}/bin/plink2)}
B=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502
mkdir -p ${ROOT}/ld/ref && cd ${ROOT}/ld/ref
[ -s ../eur.id ] || curl -sL --retry 5 "$B/integrated_call_samples_v3.20130502.ALL.panel" | awk -F'\t' 'NR>1 && $3=="EUR"{print $1}' > ../eur.id
python3 "$HERE/22d_variant_ids.py"
python3 - <<'PY'
import os, pandas as pd
R=os.environ["ROOT"]; d=pd.read_parquet(f"{R}/ld/id2rsid.parquet")
d[['id','rsid']].to_csv(f"{R}/ld/ref/id2rsid.tsv",sep="\t",index=False,header=False)
PY
one(){ c=$1
  [ -f e$c.bed ] && return 0
  for try in 1 2 3; do
    curl -sL --retry 5 -m 7200 -o chr$c.vcf.gz "$B/ALL.chr${c}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" \
      && gzip -t chr$c.vcf.gz 2>/dev/null && break
    rm -f chr$c.vcf.gz
  done
  [ -s chr$c.vcf.gz ] || { echo "download failed: chr$c" >&2; return 1; }
  $PLINK2 --vcf chr$c.vcf.gz --keep ../eur.id --snps-only --max-alleles 2 --maf 0.01 \
    --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 100 --rm-dup exclude-all \
    --memory 2000 --make-bed --out t$c --silent
  $PLINK2 --bfile t$c --update-name id2rsid.tsv --extract id2rsid.tsv --rm-dup exclude-all \
    --memory 2000 --make-bed --out e$c --silent
  rm -f chr$c.vcf.gz t$c.*
}
export -f one; export PLINK2 B
seq 1 22 | xargs -P 3 -I{} bash -c 'one {}' || true
n=$(ls e*.bed 2>/dev/null | wc -l)
[ "$n" -eq 22 ] || { echo "only $n of 22 chromosomes built; re-run this script to resume" >&2; exit 1; }
ls e*.bed | sed 's/.bed//' | sort -V > merge.txt
$PLINK2 --pmerge-list merge.txt bfile --memory 3000 --make-bed --out ../EUR --silent
rm -f e*.bed e*.bim e*.fam
echo "EUR reference: $(wc -l < ../EUR.bim) SNPs, $(wc -l < ../EUR.fam) samples"
