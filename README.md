# Phenotype-definition audit for GWAS summary statistics

A screen for phenotype-definition contamination in GWAS used as Mendelian randomisation (MR) exposures. It runs on public summary statistics and needs no instruments.

Vishnupriya Kannan and Marie Loh, Lee Kong Chian School of Medicine, Nanyang Technological University, Singapore.

> Status: research code for work under preparation (Genome Informatics 2026 abstract). Version 0.2 corrects the jackknife blocks, the negative-control rule and the multivariable MR, and makes the UK Biobank scan run from a clean clone. **The tables and figures in `results/` were produced by v0.1 and have not yet been regenerated**; see [CHANGELOG.md](CHANGELOG.md) and [Known gaps](#known-gaps).

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

Five European AD definitions, seven outcomes, all public data.

| Definition | lambda | excess rg | MVMR conditional OR (asthma) | Verdict |
|---|---|---|---|---|
| Clinician criteria (EAGLE) | 1.00 | 0.00 | reference | clean |
| Hospital ICD (UK Biobank) | 0.96 | −0.04 | 1.03 (0.88–1.21), p = 0.69 | clean, fully absorbed |
| Pooled meta-analysis (2023) | 0.69 | −0.07 | 1.28 (1.02–1.61), p = 0.032 | partly contaminated |
| Self-report (UK Biobank)¹ | 0.29 | +0.38 | 1.57 (1.44–1.71), p = 1e-24 | contaminated |
| Broad allergic composite | 0.28 | +0.28 | 2.03 (1.79–2.31), p = 3e-28 | heavily contaminated |

¹ Under verification. The linear-to-log-odds scale factor for this file (GCST90029017, Loh 2018) implies a case fraction of about 22%. UK Biobank self-reported eczema/dermatitis (20002_1452) is 2.6%; the touchscreen item *hayfever, allergic rhinitis or eczema* (6152_9) is 23%. If the file is the combined allergy/eczema phenotype, this arm contains rhinitis cases by construction and should be relabelled. `06_concordance.py` now prints this audit for every linear-model file.

## Worked example: depression

The Glanville 2021 UK Biobank ladder; reference CIDI lifetime depression, contaminant neuroticism. Instrument counts collapse for the strict arms, so MR is not possible, but the index still runs.

| Definition | rg vs reference | rg vs contaminant | excess |
|---|---|---|---|
| 1 endorsed measure | 0.606 | 0.656 | 0.359 |
| 2 endorsed measures | 0.577 | 0.714 | 0.431 |
| 3 endorsed measures | 0.744 | 0.614 | 0.250 |
| 4–5 endorsed measures | 0.521 | 0.849 | 0.594 |

Contamination can belong to the measurement instrument rather than to the breadth of the definition: raising the threshold on a neuroticism-loaded questionnaire concentrates it instead of removing it.

## Systematic scan: UK Biobank

29 trait families from the Neale lab round-2 GWAS, 35 definition comparisons, 11 candidate contaminants, definition varying within one cohort and array. Negative controls (a curated endpoint that is the same ICD-10 code as the hospital arm) are fixed in advance in `23_ukb_families.py` and reported whatever their excess. Outputs: `results/Table8_ukb_definition_scan.csv`, `results/Table10_sensitivity.csv`, `results/Fig6_ukb_definition_scan.png`.

## Running it

Dependencies: Python 3 (`pip install -r requirements.txt`), PLINK 2, curl. No R, no LDSC install.

```
pip install -r requirements.txt
python3 tests/test_core.py          # fast checks of the statistical core, no downloads
ROOT=/path/to/workdir bash code/RUN_SCAN.sh
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

- **Self-contained LD scores.** No dependency on the precomputed `eur_w_ld_chr` release; LD scores are computed from 1000 Genomes EUR (1 Mb window, bias-corrected r²), so the audit runs on any reference panel you can download. The thinned SNP set deflates LD scores, so heritability on this scale is not interpreted; genetic correlation is a ratio and is unaffected. The implementation reproduces published correlations (asthma–allergic rhinitis 0.834; depression–anxiety 0.882).
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

- **Results not regenerated.** Everything in `results/` comes from v0.1. Re-running with v0.2 will change standard errors (contiguous blocks), the negative-control set (five a-priori pairs; touchscreen versus self-reported hypertension becomes an informative comparison), Table7 standard errors, and the MVMR instrument sets and conditional F. Point estimates of rg and of the index are unchanged by the block fix.
- **The 'self-report' AD arm** needs its phenotype definition confirmed (footnote ¹ above).
- **AD and depression families (scripts 00–21) need inputs you supply:** the GWAS files in `$ROOT/gwas` (sources in the pre-registration) and a 1000 Genomes EUR PLINK set at `$ROOT/ld/EUR` for clumping and for `10_ldscores.sh`. `10_ldscores.sh` is the fixed-count LD build used for those families; the UK Biobank scan uses the fixed-fraction build from `22b`.
- **LD score regression is simplified** relative to the reference software: weights are 1/max(LD score, 1) without the heteroskedasticity term, and there is no two-step estimator or intercept constraint. Genetic correlations reproduce published values (above).

## Citation

Kannan V, Loh M. Phenotype-definition contamination is detectable from summary statistics alone: a genome-wide screen validated against multivariable Mendelian randomisation. Abstract submitted to Genome Informatics 2026, Wellcome Genome Campus.

## Licence

MIT
