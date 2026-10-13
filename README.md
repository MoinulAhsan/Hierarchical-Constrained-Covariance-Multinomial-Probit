# Hierarchical Constrained Covariance Models for Multinomial Probit Regression

R code and supporting material for the manuscript:

**A novel hierarchical constrained covariance model in multinomial probit regression for single- and multi-response datasets**

**Authors:** Md. Moinul Ahsan (Virginia Commonwealth University) and Avishek Chakraborty (University of Arkansas).

## Overview

This project develops a Bayesian covariance parameterization for multinomial probit (MNP) regression and extends it to multivariate nominal outcomes (MVMNP). The covariance matrix is decomposed as

$$
\Sigma = S R S,
$$

where `S` is a diagonal matrix of latent-variable standard deviations and `R` is a correlation matrix. Writing `phi_j = log(S_jj)`, the single-response identification constraint is

$$
\sum_{j=1}^{k}\phi_j=0,
$$

which implies that the product of the latent-variable variances equals one. For multiple responses, the zero-sum constraint is imposed separately within each response group. The model uses a constrained multivariate-normal prior for the log-standard deviations and a joint-uniform prior over valid correlation matrices.

The repository includes Markov chain Monte Carlo (MCMC) implementations, simulation studies, and single- and bivariate-response real-data analyses. A trace-constrained MNP implementation (`BN12_MCMC.R`) is included for comparison.

## Repository contents

| File or directory | Purpose |
|------------------------------------|------------------------------------|
| `Additional_functions.R` | Supporting functions used by analyses |
| `univariate_MCMC.R` | Single-response MNP MCMC |
| `Bivariate_joint_S.R` | Bivariate MNP with joint updating of scale parameters |
| `Bivariate_conditional_S.R` | Bivariate MNP with conditional updating of scale parameters |
| `BN12_MCMC.R` | Trace-constrained comparison implementation (BN12) |
| `Category_specific_X_MCMC.R` | MNP implementation with category-specific covariate information |
| `Covariate_processing.R` | Covariate-processing utilities |
| `Ordinal_probit_model.R` | Ordinal probit-related code |
| `SSVS_for_Joint_S.R` | Stochastic search variable selection for the joint-S setting |
| `Simulation_Study.Rmd` | Simulation experiments and model comparisons |
| `Real_data_Univariate.Rmd` | Single-response real-data analysis |
| `Real_data_Bivariate.Rmd` | Bivariate-response real-data analysis |
| `Real_data.csv` | Data file used by real-data scripts; check its sharing permissions before publishing |
| `Figure/` | Figures associated with the analyses |

## Methodology

The computational implementation considers two updates for the log-standard-deviation vector:

- **Joint S:** joint Metropolis–Hastings updating under the sum-to-zero constraint.
- **Conditional S:** component-wise Metropolis–Hastings proposals that preserve the same constraint.

Correlation parameters are updated while maintaining positive definiteness, with grid-based or shifted-beta proposals considered in the manuscript. The analyses report posterior summaries, predictive comparisons, and diagnostics such as HPD intervals, WAIC, BPIC, Geweke statistics, and effective sample size (ESS).

## Getting started

1.  Clone this repository:

    ``` bash
    git clone https://github.com/MoinulAhsan/Hierarchical-Constrained-Covariance-Multinomial-Probit.git
    ```

2.  Open the cloned folder in R or RStudio. Install any packages required by the scripts and R Markdown documents.

3.  Inspect and update local filesystem paths in the analysis documents as needed.

4.  Open and run `Simulation_Study.Rmd`, `Real_data_Univariate.Rmd`, or `Real_data_Bivariate.Rmd`, sourcing the required `.R` functions as indicated in those documents.

**Reproducibility note:** Package versions, computational time, and any required external inputs should be checked against the manuscript and scripts. These steps describe the intended workflow; no automated end-to-end installation or execution test is claimed here.

## Citation

Ahsan, M. M., & Chakraborty, A. *A novel hierarchical constrained covariance model in multinomial probit regression for single- and multi-response datasets.* Manuscript.

## Contact

**Md. Moinul Ahsan**\
Department of Biostatistics, Virginia Commonwealth University\
Email: [ahsanm8\@vcu.edu](mailto:ahsanm8@vcu.edu){.email}

