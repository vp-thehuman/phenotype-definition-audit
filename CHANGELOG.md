# Changelog

## v0.3 (October 2026): UK Biobank scan re-run

### Fixes found while re-running
- **rsIDs for LD scores.** The 1000 Genomes phase 3 v5b VCFs have '.' in the ID column, so
  every variant shared one name and LD scores were meaningless. 22b now names variants
  chr:pos:ref:alt and 22c maps them to rsIDs with the new `22d_variant_ids.py` (Neale
  variants manifest, GRCh37).
- **LD score regression weights.** `ldsc3.py` now uses reference-LDSC heteroskedasticity
  weights (inverse expected variance from a first-pass estimate) in addition to
  1/max(LD score, 1). With contiguous jackknife blocks the unweighted estimator had
  standard errors several times larger than necessary.
- Fig6 title now reports the computed family count and control maximum.

### Results (Table8, Table8b, Table9, Table10, Fig6, ukb_h2.json regenerated)
- LD build: 383,878 SNPs; mean LD score 7.0–15.6 by chromosome. 100 GWAS streamed.
- Validation: self-reported asthma–hayfever rg 0.52 (v0.1: 0.38); self-reported
  depression–neuroticism 0.71 (v0.1: 0.74).
- Heritability standard errors are roughly 3–10 times larger than v0.1 (e.g. BMI z 19.8
  against 90.8). v0.1's shuffled blocks understated them.
- Pre-specified filter h2 z ≥ 6: 7 families (v0.1: 29), 12 comparisons, 2 negative
  controls (max excess 0.004), 10 informative; 3 with excess > 0.1, 1 > 0.2, 4 beyond
  2 SE (v0.1: 19, 11 and 34 of 35). Median excess: self-report 0.26 (2 comparisons),
  curated 0.05 (6). All three comparisons above 0.1 have a reporting or hospitalisation
  trait as worst contaminant.
- Sensitivity h2 z ≥ 4 (Table8b): 19 families, 26 informative; 11 > 0.1, 3 > 0.2,
  10 beyond 2 SE; median self-report 0.14 vs curated 0.07; healthcare-contact or
  reporting trait worst in 8 of 11.
- Tables 1–7 and Figs 1–5 (AD and depression families) were not re-run: their input
  GWAS are not distributed with the repository.


## v0.2 (branch `review-fixes`, October 2026)

Code changes only. **`results/` still holds v0.1 output** and must be regenerated
before any number in it is quoted against this code.

### Statistical corrections
- **Jackknife blocks are contiguous.** `ldsc3.py`, `13_ldsc.py` and `19_family_index.py`
  now split SNPs, sorted by chromosome and position, into 200 adjacent blocks, as in
  reference LDSC. v0.1 shuffled SNPs before blocking. With LD between neighbouring SNPs,
  shuffled blocks are not independent of each other and the jackknife understates
  standard errors (`tests/test_core.py` demonstrates the direction). Affects every rg
  SE, the index SEs, and the "34 of 35 comparisons beyond 2 SE" summary. Point estimates
  are unchanged.
- **Negative controls are defined before any result is computed.** v0.1 treated a pair
  as a definitionally identical control only if rg ≥ 0.99 *and* every excess was already
  below 0.02, so the reported pass (maximum 0.0013) was guaranteed by construction and
  could not bound the upward bias of taking a maximum over eleven contaminants. v0.2
  fixes five pairs in `23_ukb_families.py` (curated endpoint = same ICD-10 code) and
  reports their excess whatever it is. Under v0.1 touchscreen versus self-reported
  hypertension was classed as a control by the data rule; it is now an informative
  comparison. The old flag is kept as `duplicate_by_data_rule`.
- **Exact SEs in the depression family.** `20_dep_se.py` jackknifes the whole index
  (`excess_se`); the v0.1 independence approximation is kept as `excess_se_approx`.
- **Multivariable MR.** `12_mvmr.py` re-clumps the union of the reference and test
  instrument lists (r² < 0.001, 10 Mb; `RECLUMP=0` reproduces v0.1) and reports the
  Sanderson–Windmeijer conditional F, which includes the reference exposure's SEs. The
  v0.1 statistic, which ignored them, is kept as `conditional_F_v01`. The pooled
  meta-analysis arm's conditional F (5.6 in v0.1) should be re-quoted from the new run.
- **Scale audit.** `06_concordance.py` prints the case fraction implied by each
  linear-model file's scale factor. For E2 ("self-report", GCST90029017) it is about 22%,
  against 2.6% for UK Biobank self-reported eczema/dermatitis (20002_1452) and 23% for
  *hayfever, allergic rhinitis or eczema* (6152_9). The phenotype of E2 needs confirming.

### Reproducibility
- All scripts read `ROOT` from the environment (default `/tmp`).
- Added `ldsc3.py`, imported by 29 but previously missing; 25, 25b and 28 use it
  (28 previously imported a missing `ldsc` module).
- Added `23b_ukb_whitelist.py` (builds the whitelist 24 reads); 23 fetches the Neale
  phenotype manifest; 22a writes the EUR sample list; 22b calls 22c from the repository,
  fails loudly on a missing chromosome and builds `ldscores.parquet`; 22c writes positions.
- `RUN_SCAN.sh` checks dependencies and runs nine steps including the whitelist and Table9.
- Added `tests/test_core.py`.

### Numbers to correct in the abstract, memo and pre-registration
| Item | Was | Correct |
|---|---|---|
| Negative-control maximum excess | "at most 0.001" (abstract, memo); "0.013" (pre-registration D6) | 0.0013 in v0.1, and not a valid control (see above) |
| Screen versus MVMR | "rank-consistent" | Pearson r = 0.78 with log OR (k = 4); Spearman 0.6; the meta-analysis arm has excess −0.07 but MVMR OR 1.28 |
| λ range | 3.5-fold | 0.963 / 0.284 = 3.4-fold |
| Mean instrument-effect correlation | 0.96 | 0.96 includes the reference's own 1.00; the four comparisons average 0.95 |
| Asthma heterogeneity | 96%, 1.70-fold | Table6 (outcome ladder) 96.1%, 1.70; main MR (heterogeneity_v2) 93.9%, 1.45: name which |
| LD build quoted in TOOL_README | 427,480 SNPs; 0.43, 0.73 | released build 383,835 SNPs; 0.38, 0.74 |
| UK Biobank families in memo file list | 27 | 29 |

## v0.1 (initial release)
