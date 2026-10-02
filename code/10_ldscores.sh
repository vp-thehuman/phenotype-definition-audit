#!/bin/bash
# LD scores for the AD and depression families. v0.3 uses the same fixed-FRACTION build as the
# UK Biobank scan (22a, 22b): 1000 Genomes EUR, MAF >= 0.05, 6.5% of SNPs per chromosome,
# 1 Mb window, bias-corrected r^2. The v0.1 fixed-COUNT build is kept as
# 10_ldscores_fixedcount_v01.sh for reproducing the original tables.
ROOT=${ROOT:-/tmp}; export ROOT
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
[ -s ${ROOT}/ld/ldscores.parquet ] || { bash "$HERE/22a_1000g_fetch.sh"; bash "$HERE/22b_ldscores_1000g.sh"; }
mkdir -p ${ROOT}/ldsc
cp ${ROOT}/ld/ldscores.parquet ${ROOT}/ldsc/ldscores.parquet
echo "LD scores for the AD/depression families: $(python3 -c "import pandas as pd;print(len(pd.read_parquet('${ROOT}/ldsc/ldscores.parquet')))") SNPs"
