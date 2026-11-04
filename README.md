# GOSDT Mortality Classification in Severe Malaria

## Project overview

This repository contains R and Python code developed to investigate mortality classification in severe falciparum malaria using Generalized and Scalable Optimal Sparse Decision Trees (GOSDT).

The analyses explore the use of clinical and laboratory measurements to construct interpretable decision trees for classifying mortality risk. They include data preparation, multiple imputation, selection of candidate thresholds, decision tree optimisation, and evaluation of classification performance using held-out validation and repeated validation.

Analyses were conducted using both the original dataset and datasets generated through multiple imputation. The project also includes exploratory investigations of relationships between clinical measurements and procedures for selecting thresholds used in decision tree construction.

This repository documents exploratory methodological work. The resulting decision trees have not been independently validated for clinical use and should not be interpreted as established clinical prediction tools.

## Relationship to the severe malaria staging project

The analyses in this repository informed the development of the GOSDT methodology subsequently used to investigate the prognostic contribution of peripheral parasite staging in severe malaria.

The subsequent project extended the decision tree analyses to include measures of parasite developmental stage and compared classification performance with and without staging information.

The code for that project is available in the [Severe Malaria Staging repository](https://github.com/Infectious-Diseases-Data-Observatory/severe-malaria-staging-code).

The present repository is retained as a record of the earlier methodological development and exploratory analyses.

## Repository structure

The repository is organised into three principal directories:

- `data/` contains Python scripts for exploratory data analysis and examination of relationships between clinical measurements.
- `R/` contains scripts for data preparation, multiple imputation, threshold selection, decision tree analysis, and visualisation.
- `Python/` contains scripts for fitting GOSDT classifiers and evaluating classification performance.

### Data exploration (`data/`)

| File | Description |
|---|---|
| `Heatmap_correlation.py` | Generates correlation heatmaps to examine relationships between clinical and laboratory variables. |
| `Heatmap_correlation_check.py` | Provides an additional implementation for checking the correlation heatmap analysis. |
| `Linear_regression.py` | Examines relationships between selected clinical and laboratory measurements using linear regression and graphical summaries. |
| `Linear_regression_check.py` | Provides additional checks of the linear regression analysis and associated visualisations. |

### R analyses (`R/`)

| File | Description |
|---|---|
| `Impute_combined_dataset.R` | Performs multiple imputation of missing values in the combined clinical dataset. |
| `Combined_thresholds_binned.R` | Applies candidate thresholds to clinical variables in the original combined dataset to produce binary predictors for decision tree analysis. |
| `Imputed_thresholds_binned.R` | Applies candidate thresholds to clinical variables across the imputed datasets, retaining patient and imputation identifiers. |
| `Combined_analysis_script.R` | Analyses decision tree results obtained from the original combined dataset, including classification performance and patient-level predictions. |
| `Imputed_analysis_script.R` | Analyses decision tree results across imputed datasets, including classification performance and comparisons between trees. |
| `Combined_plotting_script.R` | Produces visual representations of decision trees fitted using the original combined dataset. |
| `Imputed_plotting_script.R` | Produces visual representations of decision trees fitted using imputed datasets. |
| `Thresholds_for_mortality.R` | Investigates relationships between clinical measurements and mortality to inform candidate threshold selection. |
| `Mortality_for_thresholds.R` | Examines mortality in relation to selected clinical thresholds and produces graphical summaries. |

### Python analyses (`Python/`)

| File | Description |
|---|---|
| `Combined_gosdt.py` | Fits GOSDT classifiers using thresholded predictors from the original combined dataset. |
| `Combined_gosdt_validation.py` | Evaluates GOSDT classification performance using a held-out training and test partition of the original combined dataset. |
| `Combined_gosdt_repeated_validation.py` | Evaluates classification performance across repeated training and test partitions of the original combined dataset. |
| `Imputed_gosdt.py` | Fits GOSDT classifiers using thresholded predictors from multiply imputed datasets. |
| `Imputed_gosdt_validation.py` | Evaluates classification performance using training and test partitions while preserving patient grouping across imputations. |
| `Imputed_gosdt_repeated_validation.py` | Evaluates classification performance across repeated training and test partitions of the imputed datasets, preserving patient grouping across imputations. |

The repository also contains `gosdt-severe-malaria.Rproj`, an RStudio project file, and `.gitignore`, which specifies files excluded from version control.

## Software requirements

The analyses use R and Python. The versions below were recorded from the software environment inspected in November 2026.

### R

**R version:** 4.6.1 (2026-06-24)

**Platform:** macOS, Apple Silicon (aarch64)

| Package | Version |
|---|---|
| `data.table` | 1.18.4 |
| `mice` | 3.19.0 |
| `dplyr` | 1.2.1 |
| `tidyr` | 1.3.2 |
| `ggplot2` | 4.0.3 |
| `mgcv` | 1.9.4 |
| `DiagrammeR` | 1.0.12 |
| `DiagrammeRsvg` | 0.1 |
| `rsvg` | 2.7.0 |

The packages listed above are imported or explicitly checked by the R scripts.

### Python

**Python version:** 3.13.7

| Package | Version |
|---|---|
| `gosdt` | 1.0.4 |
| `pandas` | 2.2.3 |
| `numpy` | 2.2.1 |
| `scikit-learn` | 1.6.1 |
| `matplotlib` | 3.10.0 |
| `scipy` | 1.15.1 |

SciPy is included as a dependency of the scientific computing and machine learning packages.

The Python scripts also use standard library modules, including `os`, `re`, `warnings`, and `pathlib`.

These versions describe the inspected software environment and are provided for reference. They do not establish that every analysis was originally executed using exactly these versions.

## Data availability

The underlying clinical datasets are not included in this repository because they contain restricted patient information.

Datasets generated through multiple imputation, intermediate data files, fitted decision tree objects, patient-level predictions, and generated outputs have also been excluded where they may contain or reveal restricted information.

The repository provides the analysis code but does not contain all files required for independent reproduction of the results. Access to the underlying data is subject to the relevant data sharing arrangements.

## Contact

For questions about the code or analyses in this repository, please contact:

**Mary Scott**

Email: [marypscott7@gmail.com](mailto:marypscott7@gmail.com)
