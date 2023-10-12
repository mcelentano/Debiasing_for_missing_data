<h1 align="center" style="margin-bottom:0px; border-bottom:0px; padding-bottom:0px">Challenges of the inconsistency regime</h1>
<h3 align="center" style="margin-bottom:0px; border-bottom:0px; padding-bottom:0px">Novel debiasing methods for missing data models</h3>
<p align="center" style="margin-bottom:0px; border-bottom:0px; padding-bottom:0px">Michael Celentano, Martin J. Wainwright</p>

<p align="center">
    <a style="text-decoration:none !important;" href="https://arxiv.org/abs/2309.01362" alt="arXiv"><img src="https://img.shields.io/badge/paper-arXiv-red" /></a>
</p>


This repository contains `R` code to reproduce the figures in "Challenges of the inconsistency regime: novel debiasing methods for missing data models."
Standard semi-parametric approaches for estimating population means when data is missing at random, or relatedly, average treatment effects, include outcome imputation, augmented inverse probability weighting (AIPW), and inverse probability weighting (IPW)/Horwitz-Thompson. These methods are inconsistent when both the outcome model and missingness/propensity model cannot be estimated consistently. We develop methods which achieve consistency for the population mean in a setting where neither the outcome nor missingness/propensity model can be estimated consistently.

#### Guide to paper's figures

* **Figures 1-3:** These plot data simulated by `simulations/standard_estimators.R`. The figures are produced by `plots_paper/make_plots_standard_estimators.R`.
* **Figure 4:** This plots data simulated by `simulations/novel_estimators.R`. The figure is produced by `plots_paper/make_plots_novel_estimators.R`.
* **Figure 5 & 6:** these plot data simulated by `simulations/lambda_dependence.R`. The figures are produced by `plots_paper/make_plots_lambda_depdences.R`.

The files in the `utils` directory define functions used by the simulations.

Figures are output to a directory `fig_paper` and simulation data is written to and read from a directory `data`.
