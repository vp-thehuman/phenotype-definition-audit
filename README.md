# Phenotype-definition audit for GWAS summary statistics

A screen for phenotype-definition contamination in GWAS used as Mendelian randomisation (MR) exposures. It runs on public summary statistics and needs no instruments.

Vishnupriya Kannan and Marie Loh, Lee Kong Chian School of Medicine, Nanyang Technological University, Singapore.

> Status: research code for work under preparation (Genome Informatics 2026 abstract). Version 0.4 (2 October 2026): every table and figure in `results/` was regenerated from public data with the corrected code (contiguous jackknife blocks, reference-LDSC weights, a-priori negative controls, re-clumped MVMR with Sanderson–Windmeijer conditional F). See [CHANGELOG.md](CHANGELOG.md).

## The problem

GWAS that share one trait label, and are reused as MR exposures, often differ in case definition: clinician criteria, self-report, hospital coding, symptom questionnaires, composite phenotypes. Loosening a definition does two different things:

- **Attenuation.** A looser definition measures the same liability more noisily, so all instrument effects shrink by a scalar. In an MR ratio (outcome effect / exposure effect) the scalar cancels. Harmless for MR.
- **Non-specificity.** A looser definition also imports other conditions. Instruments then tag something broader than the trait, and estimates inflate for any outcome genetically correlated with what was imported. This does not cancel, and standard MR sensitivity analyses (F, MR-Egger, weighted median) do not detect it, because the instruments are valid for the phenotype as defined.

## The audit

Four steps, in increasing cost.

1. **Concordance (lambda).** Deming regression of instrument effects in definition D against a reference definition R. Low lambda flags a candidate but is not itself evidence of bias.
2. **Contamination index (cheap, genome-wide, no instruments).**

   `excess = rg(D, C) − rg(D, R) × rg(R, C)`

   where C is a nominated contaminant. Subtracting the reference term separates contamination from genuine shared biology.
3. **Multivariable MR (confirmatory).** Regress the outcome jointly on D and R. If D is merely noisier, conditioning absorbs it; if D retains an independent effect, it carries non-reference signal.
4. **Outcome triage.** Bias reaches only outcomes genetically related to the contaminant; `rg(outcome, C)` says which outcomes are at risk.

## Worked example: atopic dermatitis

Five European AD definitions, seven outcomes, all public data (sources in `code/01_harmonise.py`, `03`, `14`).

| Definition | lambda | excess rg | MVMR conditional OR (asthma) | conditional F | Verdict |
|---|---|---|---|---|---|
| Clinician criteria (EAGLE) | 1.00 | 0.00 | reference | — | reference |
| Hospital ICD (UK Biobank) | 0.96 | +0.08 | 1.00 (0.86–1.16), p = 0.99 | 4.4 | fully absorbed |
| Pooled meta-analysis (2023) | 0.69 | +0.04 | 1.32 (1.00–1.73), p = 0.046 | 1.0 | not identifiable (weak) |
| Hay fever, rhinitis or eczema (UK Biobank)¹ | 0.29 | +0.49 | 1.60 (1.46–1.75), p = 1e-23 | 13.1 | contaminated |
| Broad allergic composite | 0.28 | +0.49 | 2.06 (1.79–2.37), p = 2e-23 | 5.3 | contaminated |

The screen separates the two definitions MVMR retains (excess 0.49 and 0.49) from the hospital-ICD definition it absorbs (0.08). It does not flag the pooled meta-analysis (0.04), whose conditional estimate is not identifiable (conditional F 1.0). Pearson r = 0.84 between excess and log MVMR OR (Spearman 0.6, k = 4), descriptive. Outcome triage: genetic correlation of each outcome with allergic rhinitis predicts its between-definition spread (r = 0.90, k = 7).

¹ GCST90029017 (Loh 2018) is catalogued as "eczema" and was labelled self-reported eczema in v0.1. It is the UK Biobank touchscreen item *hayfever, allergic rhinitis or eczema*: its implied case fraction is 22% (the item's prevalence is 23%; self-reported eczema/dermatitis, 20002_1452, is 2.6%), and its genetic correlation is 0.98 (SE 0.04) with that item and 0.47 (SE 0.08) with self-reported eczema. It contains rhinitis cases by construction.

## Worked example: depression

The Glanville 2021 UK Biobank ladder; reference CIDI lifetime depression, contaminant neuroticism. Instrument counts collapse for the strict arms (1–8 clumped instruments), so MR is not possible, but the index still runs. rg(CIDI, neuroticism) = 0.46 (SE 0.04).

| Definition | rg vs reference | rg vs contaminant | excess (SE) |
|---|---|---|---|
| 1 endorsed measure | 0.66 | 0.72 | 0.42 (0.05) |
| 2 endorsed measures | 0.70 | 0.68 | 0.36 (0.05) |
| 3 endorsed measures | 0.55 | 0.61 | 0.35 (0.05) |
| 4–5 endorsed measures | 0.49 | 0.64 | 0.42 (0.06) |

Every touchscreen definition is contaminated relative to the interview reference, and the excess does not fall as the threshold rises: contamination belongs to the measurement instrument, not to the breadth of the definition. (v0.1 reported the strictest arm as the most contaminated, 0.59; with corrected standard errors and weights the four arms are indistinguishable.)

## Systematic scan: UK Biobank

34 candidate trait families (89 arms, ≥ 2,000 cases each) from the Neale lab round-2 GWAS and 11 candidate contaminants, with definition varying within one cohort, array and control set. Negative controls (a curated endpoint that is the same ICD-10 code as the hospital arm) are fixed in advance in `23_ukb_families.py` and reported whatever their excess.

Results (v0.3). With the pre-specified heritability filter (SNP-heritability z ≥ 6), 7 families give 12 comparisons: 2 negative controls (excess 0.004 and 0.002) and 10 informative comparisons, of which 3 carry excess above 0.1 and 4 are more than two standard errors from zero. The three above 0.1 are self-reported osteoarthritis (0.38, SE 0.07), self-reported myocardial infarction (0.14, SE 0.07) and a curated osteoarthritis endpoint (0.13, SE 0.05), and in all three the worst contaminant is illness-reporting or hospitalisation propensity. Median excess is 0.26 for self-report against 0.05 for curated endpoints. Relaxing the filter to z ≥ 4 (Table8b) gives 19 families and 26 informative comparisons: 11 above 0.1, 3 above 0.2, 10 beyond two standard errors, median excess 0.14 for self-report against 0.07 for curated endpoints, and a healthcare-contact or reporting trait as the worst contaminant in 8 of the 11. Outputs: `results/Table8_ukb_definition_scan.csv`, `results/Table8b_ukb_scan_h2z4.csv`, `results/Table10_sensitivity.csv`, `results/Fig6_ukb_definition_scan.png`.

Why v0.1 reported 29 families: its jackknife used shuffled SNP blocks, which understated standard errors roughly ten- to twentyfold, so nearly every arm passed the heritability filter and 34 of 35 comparisons appeared to exceed two standard errors.

## Running it

Dependencies: Python 3 (`pip install -r requirements.txt`), PLINK 2, curl. No R, no LDSC install.

```
pip install -r requirements.txt
python3 tests/test_core.py          # fast checks of the statistical core, no downloads
ROOT=/path/to/workdir bash code/RUN_SCAN.sh       # UK Biobank scan (Table8-10, Fig6)
ROOT=/path/to/workdir bash code/RUN_FAMILIES.sh   # AD and depression families (Tables 1-7, Figs 1-5)
```

All scripts read their working directory from `ROOT` (default `/tmp`). `plink2` must be on `PATH`, at `$ROOT/bin/plink2`, or given as `PLINK2=`.

About two hours on two cores / 8 GB RAM, mostly download (1000 Genomes ~15 GB; 100 UK Biobank GWAS ~60 GB streamed, not stored). Peak disk about 25 GB.

| Scripts | Purpose |
|---|---|
| `00`–`03` | harmonise GWAS, orientation check, PLINK 2 clumping, outcome extraction |
| `06`–`09` | concordance (lambda), univariable MR, figures, exposure matrix for MVMR |
| `10`–`13` | LD scores from 1000G EUR, LD-score SNP extraction, MVMR, cross-trait LD score regression |
| `14`–`16` | outcome ladder, contamination index, final figures |
| `17`–`21` | applying the audit to a new trait family (depression), standard errors, cross-family figure |
| `22a`–`22c` | 1000 Genomes panel, EUR sample list and LD scores (uniform sampling fraction) |
| `23`, `23b` | UK Biobank families, contaminant panel, a-priori negative controls; variant whitelist |
| `24`–`29` | streaming, all genetic correlations (Table9), validation, scan with sensitivity analyses, figure |
| `ldsc3.py` | shared LD score regression and jackknife engine used by 25, 25b, 28 and 29 |
| `tests/` | unit checks of the five-sum WLS slope, leave-one-block-out, block construction, conditional F |

## Design choices

- **Self-contained LD scores.** No dependency on the precomputed `eur_w_ld_chr` release; LD scores are computed from 1000 Genomes EUR (1 Mb window, bias-corrected r²), so the audit runs on any reference panel you can download. The phase 3 v5b VCFs carry no rsIDs, so variants are named chr:pos:ref:alt and mapped to rsIDs with the Neale variants manifest (`22d`). The v0.3 build has 383,878 SNPs, mean LD score 7.0 to 15.6 by chromosome. The thinned SNP set deflates LD scores, so heritability on this scale is not interpreted; genetic correlation is a ratio and is unaffected. The implementation reproduces published correlations (asthma–allergic rhinitis 0.834; depression–anxiety 0.882).
- **Uniform SNP density.** The released build samples a fixed fraction (`--thin 0.065`) of MAF ≥ 0.05 SNPs per chromosome rather than a fixed count, which would make LD-score magnitude depend on chromosome length. Genetic correlations moved by < 0.05 between the two builds.
- **Effective sample size from standard errors**, not reported counts, so files with incomplete metadata still run.
- **Exact standard errors.** The index is re-formed inside each of 200 delete-one-block jackknife replicates (`ldsc3.py`, and `20_dep_se.py` for the depression family), since its three rg terms share SNPs and are not independent. Blocks are contiguous runs of SNPs in genomic order, as in reference LDSC; v0.1 used shuffled SNP sets, which understates standard errors when neighbouring SNPs are in LD.
- **Negative controls fixed in advance.** Pairs where a curated endpoint is the same ICD-10 code as the hospital arm. v0.1 labelled a pair a control only if every excess was already below 0.02, so a control could not fail; that flag is still written as `duplicate_by_data_rule` for comparison.
- **Multivariable MR.** The union of the two instrument lists is re-clumped before fitting, and conditional F is the Sanderson–Windmeijer statistic, which includes the reference exposure's standard errors (`12_mvmr.py`; v0.1 values are kept in `conditional_F_v01`).
- **Winsorising** rg at 1.25 (GWAS ATLAS convention); results unchanged without it (Table10).
- **Rare binary traits.** Neale round-2 fits linear models to 0/1 phenotypes. Arms require ≥ 2,000 cases, and Table10 restricts to prevalence ≥ 1%. The ordering by definition type survives; the single most extreme result (noninfectious colitis vs hospital-coded ulcerative colitis) does not and should be quoted with that caveat.
- **The reference is the best available arm, not a gold standard.** Contamination shared by a definition and its reference is invisible to the index, biasing it towards the null. Table10 re-runs with two alternative reference rules; the ranking by definition type is preserved.

## Related work

- Cai et al. 2020, *Nature Genetics*, doi:10.1038/s41588-020-0594-5 — minimal phenotyping yields depression GWAS signals of low specificity.
- de la Fuente, Londono-Correa and Tucker-Drob 2025, *Bioinformatics* — residual genetic correlation within genomic SEM; for a single reference indicator it reduces to the index used here. This repository provides a closed form needing no model fitting or instruments, validation against multivariable MR, and a systematic scan.

## Known gaps

- **One reference panel.** LD scores and clumping use 1000 Genomes EUR (503 people); all GWAS here are European-ancestry.
- **AD and depression families (scripts 00–21)** read GWAS files from `$ROOT/gwas` under the names used in the scripts (GWAS Catalog GCST003184, GCST90029017, GCST90038680, GCST005038, GCST90244787, GCST005839, GCST90014426/28/30/32/34, GCST90029028; FinnGen R11 F5_DEPRESSIO, F5_ALLANXIOUS, J10_ASTHMA_EXMORE, ALLERG_RHINITIS, J10_COPD, K11_IBD_STRICT, G6_MIGRAINE). `00a_eur_reference.sh` builds the 1000 Genomes EUR PLINK reference used for clumping. The LD scores are the same fixed-fraction build as the scan; the v0.1 fixed-count build is kept as `10_ldscores_fixedcount_v01.sh`.
- **LD score regression weights.** All LD score regression (`ldsc3.py`, `13_ldsc.py`, `19_family_index.py`, `20_dep_se.py`) uses reference-LDSC weights: 1/max(LD score, 1) times the inverse expected variance of each statistic from a first-pass estimate. There is no two-step estimator or intercept constraint. v0.3 validation on the UK Biobank build: self-reported asthma with hayfever 0.52, self-reported depression with neuroticism 0.71.
- **Most UK Biobank arms are weakly heritable on this scale.** With honest standard errors only 7 of 34 families pass z ≥ 6. A denser LD-score SNP set (the thinned set is 6.5% of common SNPs) would increase power.

## Citation

Kannan V, Loh M. Phenotype-definition contamination is detectable from summary statistics alone: a genome-wide screen validated against multivariable Mendelian randomisation. Abstract submitted to Genome Informatics 2026, Wellcome Genome Campus.

## Licence

MIT
