# Clinical Data Classification

**Predicting Atrial Fibrillation Following Myocardial Infarction: A Comparison of Naive Bayes and k-NN**

Academic project for Statistical Methods in Data Mining, Instituto Superior Técnico, 2026/2027.

**Group 1:** Maria Oliveira Maia, Jéssica Maria Moreira Lima, Ebba Sofia Löwenberg, Rupashree Dey and Bernardo Delgado de Brito Malheiro.

## Objective

Compare Naive Bayes and k-nearest neighbours for predicting atrial fibrillation (`FIBR_PREDS`) using admission-related variables in a myocardial infarction dataset with 1,700 patients and 124 original columns.

The project explores missing-data handling, feature filtering, class imbalance, threshold selection and the gap between validation and test performance. This is an academic experiment, not a clinically validated prediction system.

## Methodology

- Stratified 80/20 split with seed 42: 1,360 training and 340 held-out test patients (34 positives).
- Five-fold cross-validation within training data for model and preprocessing selection.
- Median/mode imputation, categorical handling and feature filtering fitted on training-fold patients only.
- Naive Bayes with missingness indicators; k-NN with Gower-distance contributions and training-derived numeric scales.
- Comparison of baseline and validation-selected preprocessing across 8,192 configurations, with 12 candidate k values.
- Classification at 0.5 and at a frozen validation-selected Youden threshold.

Outlier diagnostics were exploratory; observations were not removed or winsorised on that basis. The fitting row-missingness filter is a separate preprocessing step; validation and test patients remain in evaluation.

## Recorded results

| Model | Validation ROC-AUC | Test ROC-AUC |
| --- | ---: | ---: |
| Naive Bayes — baseline | 0.6953 | 0.7598 |
| Naive Bayes — optimized | 0.7133 | 0.7141 |
| k-NN — baseline, k=121 | 0.6291 | 0.7071 |
| k-NN — optimized, k=25 | 0.6884 | 0.6762 |

At its validation-selected threshold, baseline Naive Bayes identified 26 of 34 positive test cases: **76.5% sensitivity and 61.8% specificity**. This operating point also had low precision (18.2%), illustrating the false-positive trade-off.

**Validation improvements did not translate into higher test ROC-AUC.** Searching many configurations can make selection-data performance optimistic. Both baseline and optimized results are presented, rather than treating the optimized models as universally better. No settings should be retuned using these observed test results.

The tables in `results/` are copied from the supplied completed project. The full [project report](reports/project-report.pdf) provides the methodology and discussion.

## Portfolio package

- `R/01_data.R`: original import, variable selection and stratified split.
- `R/03_functions.R`: original preprocessing, metrics and modelling helpers.
- `R/07_search_setup.R`: original configuration grid and prediction functions.
- `run_fixed_models.R`: new portfolio entry point that fits the four already selected configurations without the archived search checkpoints.
- `results/configurations.csv`: existing fixed thresholds and selected k values.
- `results/validation_comparison.csv`: recorded validation metrics.
- `results/test_metrics.csv`: recorded test metrics.
- `data/mi_variables.csv`: variable names and descriptions used by the original import.

The three original R source files are unchanged. Raw patient records, fitted model objects containing training/reference patients, patient-level predictions and search checkpoints are excluded. This is a focused extract, not the complete historical search archive.

## Setup and run

The supplied original model-build record used R 4.3.3. Use a compatible R installation and set the working directory to this repository root.

Install the two core packages:

```r
source("install_dependencies.R", encoding = "UTF-8")
```

Obtain an authorised copy of the original course dataset and place it at `data/MI.data`. It must retain the expected original row order and 1,700 × 124 structure. The variable dictionary alone is not enough to execute the models.

Run in a fresh R session:

```r
source("run_fixed_models.R", encoding = "UTF-8")
```

This imports the dataset, recreates the training/test split, refits all four fixed configurations with their stored validation-selected thresholds and writes aggregate metrics to `outputs/report/refitted_test_metrics.csv`. It does not reconstruct the original cross-validation/search, retrain frozen submitted artifacts, or change the archived tables in `results/`.

`R/07_search_setup.R` also writes the planned grid; sourcing it does not launch the 8,192-configuration search. Some original helpers support optional EDA/outlier diagnostics requiring extra packages, but these are not invoked by the portfolio entry point.

## Verification status

Package file paths, source dependencies and recorded metric tables were checked. **R is unavailable in the packaging environment, so the new portfolio entry point has not been executed or syntax-checked with R.** The reported metrics are archived project results, not newly reproduced results. Dependency versions are not locked.

## Contact

[Maria Oliveira Maia — LinkedIn](https://www.linkedin.com/in/mariaoliveiramaia/)
