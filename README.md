# psg_modeling

This repository contains code that generates the results from the paper: G. Palmer, D.B. Dunson (2025+) 'Quantifying sleep apnea heterogeneity using hierarchical Bayesian modeling', as well as a Stan model implementing our approach that can be used for other data sets.

If you have any questions, find bugs, etc. please reach out to glenn.palmer@duke.edu.

# Code for results in paper

The R markdown notebooks and R files in the main directory and listed below can be used to generate the results and figures from the paper.

* `joint_model_factor.Rmd`: Fits the model to data from the APPLES study (available by request at https://sleepdata.org/datasets/apples), saves the output, summarizes fit metrics (Rhat, ESS), estimates the 2nd stage clustering, and makes some basic summary figures.
* `paper_figures.Rmd`: Using output from the above, generates the remaining figures and table results in the paper.
* `posterior_predictive_checks.Rmd`: Using the fitted model, performs the posterior predictive checks summarized in the paper and supplement.
* `clustering_supplement.Rmd`: Using the fitted model, performs and summarizes the additional clustering methodology shown in the supplement.
* `simulation_scenario[1-3].R`: Scripts that generate synthetic polysomnogram summaries, fit the model to them, and save results.

