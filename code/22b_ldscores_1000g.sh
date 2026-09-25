cd /tmp/ld
mkdir -p out2
proc(){ c=$1
  /tmp/bin/plink2 --vcf chr$c.vcf.gz --keep eur.id --snps-only --max-alleles 2 --maf 0.05 \
    --rm-dup exclude-all --memory 2500 --make-bed --out f$c --silent || return 1
  rm -f chr$c.vcf.gz
  /tmp/bin/plink2 --bfile f$c --thin 0.065 --seed 42 --memory 2500 --make-bed --out c$c --silent || return 1
  rm -f f$c.bed f$c.bim f$c.fam f$c.log
  /tmp/bin/plink2 --bfile c$c --r2-unphased --ld-window-kb 1000 --ld-window 999999 \
    --ld-window-r2 0 --memory 2500 --out c$c --silent || return 1
  python3 /tmp/ld/agg2.py $c || return 1
  rm -f c$c.vcor c$c.log c$c.bed c$c.bim c$c.fam
}
while true; do
  n=$(ls out2/ld_chr*.parquet 2>/dev/null | wc -l); [ "$n" -ge 22 ] && break
  did=0
  for c in $(seq 1 22); do
    [ -f out2/ld_chr$c.parquet ] && continue
    [ -f chr$c.vcf.gz ] || continue
    proc $c && did=1
    break
  done
  [ $did -eq 0 ] && sleep 15
done
touch /tmp/ld/ld2.flag
