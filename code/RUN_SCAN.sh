#!/bin/bash
ROOT=${ROOT:-/tmp}
# Reproduce the UK Biobank definition scan end to end, from nothing.
#
#   bash RUN_SCAN.sh
#
# Runtime on 2 cores / 8 GB RAM: about 2 hours, almost all of it download.
#   1000 Genomes VCFs      ~15 GB transfer,  ~15 min on a fast link
#   LD scores (22 chr)     ~45 min CPU
#   100 UK Biobank GWAS    ~60 GB streamed, ~15 min, nothing stored
#   rg + index + figures   ~5 min
# Peak disk about 25 GB; nothing is kept except the SNP-level parquet files.
#
# Requires: plink2 on PATH or at ${ROOT}/bin/plink2, python3 with
# pandas numpy pyarrow matplotlib, curl. No R, no LDSC install.
set -e
ROOT=${ROOT:-/tmp}
export ROOT
cd "$(dirname "$0")"
mkdir -p ${ROOT}/ld ${ROOT}/ukb/ss
command -v plink2 >/dev/null || [ -x ${ROOT}/bin/plink2 ] || { echo "plink2 not found on PATH or at ${ROOT}/bin/plink2" >&2; exit 1; }
python3 -c "import pandas, numpy, pyarrow, matplotlib" || { echo "pip install -r ../requirements.txt" >&2; exit 1; }

echo "[1/9] 1000 Genomes phase 3 EUR panel"
bash 22a_1000g_fetch.sh

echo "[2/9] LD scores, uniform 6.5% sample of MAF>=0.05 SNPs per chromosome"
# Uniform SAMPLING FRACTION, not a fixed count per chromosome: a fixed count
# makes SNP density, and therefore the LD score scale, depend on chromosome
# length. See 22b for the flag (--thin 0.065, not --thin-count).
bash 22b_ldscores_1000g.sh

echo "[3/9] trait families, contaminant panel and a-priori negative controls"
python3 23_ukb_families.py

echo "[4/9] whitelist of Neale variants in the LD-score SNP set"
python3 23b_ukb_whitelist.py

echo "[5/9] stream 100 GWAS, keep only LD-score SNPs (5 parallel shards)"
for s in 0 1 2 3 4; do python3 24_ukb_stream.py $s 5 & done; wait

echo "[6/9] validation: reproduce known genetic correlations"
python3 28_ukb_validation.py

echo "[7/9] all genetic correlations behind the scan (Table9)"
python3 25_ukb_ldsc_fast.py

echo "[8/9] scan with exact delete-one-block jackknife (contiguous blocks), plus sensitivity analyses"
python3 29_ukb_scan.py

echo "[9/9] figure"
python3 27_ukb_figure.py

echo
echo "outputs: Table8_ukb_definition_scan.csv, Table10_sensitivity.csv,"
echo "         Fig6_ukb_definition_scan.png/.pdf"
