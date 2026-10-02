#!/bin/bash
# Download every public GWAS used by the AD and depression families into $ROOT/gwas,
# under the file names the scripts expect. About 12 GB.
ROOT=${ROOT:-/tmp}
set -e
mkdir -p ${ROOT}/gwas && cd ${ROOT}/gwas
B=https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics
F=https://storage.googleapis.com/finngen-public-data-r11/summary_stats
while read name url; do
  [ -s "$name" ] && continue
  curl -sL --retry 5 -o "$name.tmp" "$url" && mv "$name.tmp" "$name" && echo "got $name"
done <<LIST
EAGLE_AD.txt $B/GCST003001-GCST004000/GCST003184/EAGLE_AD_no23andme_results_29072015.txt
UKB_SELFREP_ECZ.h.tsv.gz $B/GCST90029001-GCST90030000/GCST90029017/harmonised/29892013-GCST90029017-HP_0000964.h.tsv.gz
UKB_ICD_ECZDERM.h.tsv.gz $B/GCST90038001-GCST90039000/GCST90038680/harmonised/33959723-GCST90038680-HP_0000964.h.tsv.gz
ALLERGIC_COMPOSITE.h.tsv.gz $B/GCST005001-GCST006000/GCST005038/harmonised/29083406-GCST005038-EFO_0003785.h.tsv.gz
BUDU_AD.tsv $B/GCST90244001-GCST90245000/GCST90244787/GCST90244787_buildGRCh37.tsv
PGC_MDD.gz $B/GCST005001-GCST006000/GCST005839/MDD2018_ex23andMe.gz
FG_DEP.gz $F/finngen_R11_F5_DEPRESSIO.gz
FG_ANX.gz $F/finngen_R11_F5_ALLANXIOUS.gz
FG_ASTHMA.gz $F/finngen_R11_J10_ASTHMA_EXMORE.gz
FG_ALLERG_RHINITIS.gz $F/finngen_R11_ALLERG_RHINITIS.gz
FG_J10_COPD.gz $F/finngen_R11_J10_COPD.gz
FG_K11_IBD_STRICT.gz $F/finngen_R11_K11_IBD_STRICT.gz
FG_G6_MIGRAINE.gz $F/finngen_R11_G6_MIGRAINE.gz
DEP1_one.h.tsv.gz $B/GCST90014001-GCST90015000/GCST90014430/harmonised/GCST90014430.h.tsv.gz
DEP2_two.h.tsv.gz $B/GCST90014001-GCST90015000/GCST90014434/harmonised/GCST90014434.h.tsv.gz
DEP3_three.h.tsv.gz $B/GCST90014001-GCST90015000/GCST90014432/harmonised/GCST90014432.h.tsv.gz
DEP4_fourfive.h.tsv.gz $B/GCST90014001-GCST90015000/GCST90014426/harmonised/GCST90014426.h.tsv.gz
DEP5_cidi.h.tsv.gz $B/GCST90014001-GCST90015000/GCST90014428/harmonised/GCST90014428.h.tsv.gz
NEUROTICISM.h.tsv.gz $B/GCST90029001-GCST90030000/GCST90029028/harmonised/29892013-GCST90029028-EFO_0004257.h.tsv.gz
LIST
echo "inputs in ${ROOT}/gwas: $(ls | wc -l) files"
