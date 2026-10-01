ROOT=${ROOT:-/tmp}
B=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502
cd ${ROOT}/ld
dl(){ c=$1
  [ -f chr$c.vcf.gz ] && return 0
  [ -f out2/ld_chr$c.parquet ] && return 0
  curl -sL -C - -m 7200 -o .chr$c.tmp "$B/ALL.chr${c}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" && mv .chr$c.tmp chr$c.vcf.gz
}
export -f dl; export B
seq 1 22 | xargs -n1 -P5 -I{} bash -c 'dl {}'
touch ${ROOT}/ld/fetch2.flag
