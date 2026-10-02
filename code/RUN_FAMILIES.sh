#!/bin/bash
# Reproduce the AD and depression families (Tables 1-7, Figs 1-5) from public data.
#   ROOT=/path/to/workdir bash code/RUN_FAMILIES.sh
# Needs plink2, Python (requirements.txt) and about 30 GB of disk. Run RUN_SCAN.sh first,
# or let step 2 build the LD scores. Mostly download time: allow 3-4 hours.
set -e
ROOT=${ROOT:-/tmp}; export ROOT
cd "$(dirname "$0")"
mkdir -p ${ROOT}/mr/data ${ROOT}/mr/clump ${ROOT}/mr/out
echo "[1/6] public GWAS";            bash 00b_fetch_inputs.sh
echo "[2/6] EUR reference, LD scores"; bash 00a_eur_reference.sh; bash 10_ldscores.sh
echo "[3/6] AD exposures and outcomes"
python3 01_harmonise.py; bash 02_clump.sh; python3 03_outcomes.py; python3 00_validate_orientation.py
echo "[4/6] concordance, MR, MVMR"
python3 06_concordance.py; python3 07_mr_v2.py; python3 09_exposure_matrix.py; python3 12_mvmr.py
echo "[5/6] genetic correlations, index, outcome ladder, depression family"
python3 11_ldsc_extract.py; python3 13_ldsc.py; python3 14_outcome_ladder.py; python3 15_index_and_ladder.py
python3 17_family_runner.py; python3 18_dep_ldsc_extract.py; python3 19_family_index.py; python3 20_dep_se.py
echo "[6/6] figures"
python3 08_figures.py; python3 16_final_figures.py; python3 21_crossfamily_figure.py
echo "outputs in ${ROOT}/mr/out"
