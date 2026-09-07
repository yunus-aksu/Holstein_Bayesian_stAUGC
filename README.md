# Bayesian Growth Modeling and stAUGC in Holstein Heifers

This repository contains the R scripts used in the paper: 
"A Bayesian Approach to Modeling Right-Censored Growth Data in Holstein Heifers: Introducing the Synchronized Area Under the Growth Curve (stAUGC)".

## Repository Structure
- `001_Classical_Growth_Models.R`: Fits 7 classical non-linear growth models individually.
- `002_Bayesian_Model_Selection.R`: Performs Bayesian non-linear mixed-effects modeling using `brms` (MCMC with seed 12345, 2000 iterations).
- `003_Bayesian_Logistic_PTI_stAUGC.R`: Calculates Predicted Time to Insemination (PTI) and synchronized area under the growth curve (stAUGC).
- `004_Growth_Curves_Plots.R`: Generates high-quality publication-ready plots (Spaghetti plot, comparisons, etc.).

## Reproducibility
To ensure exact reproducibility of MCMC chains and results, please run the scripts in sequential order using the provided seed values.
