#!/bin/bash
# 1000 Genomes phase 3 VCFs (chr 1-22) and the list of the 503 EUR samples.
ROOT=${ROOT:-/tmp}
set -e
B=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502
mkdir -p ${ROOT}/ld && cd ${ROOT}/ld
if [ ! -s eur.id ]; then
  curl -sL --retry 5 "$B/integrated_call_samples_v3.20130502.ALL.panel" \
    | awk -F'\t' 'NR>1 && $3=="EUR"{print $1}' > eur.id
fi
echo "EUR samples: $(wc -l < eur.id)"   # expect 503
dl(){ c=$1
  [ -f chr$c.vcf.gz ] && return 0
  [ -f out2/ld_chr$c.parquet ] && return 0
  curl -sL -C - --retry 5 -m 7200 -o .chr$c.tmp "$B/ALL.chr${c}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" && mv .chr$c.tmp chr$c.vcf.gz
}
export -f dl; export B
seq 1 22 | xargs -n1 -P5 -I{} bash -c 'dl {}'
touch ${ROOT}/ld/fetch2.flag
