# psg_modeling

This repository contains code that generates the results from the paper: G. Palmer, D.B. Dunson (2025+) 'Quantifying sleep apnea heterogeneity using hierarchical Bayesian modeling' (https://arxiv.org/abs/2311.16470), as well as a Stan model implementing our approach that can be used for other data sets.

If you have any questions, find bugs, etc. please reach out to glenn.palmer@duke.edu.

# Code for results in paper

The R markdown notebooks and R files in the main directory and listed below can be used to generate the results and figures from the paper.

* `joint_model_factor.Rmd`, `paper_figures.Rmd`, `clustering_supplement.Rmd`: Fits the model to data from the APPLES study (available by request at https://sleepdata.org/datasets/apples) and creates figures in the paper and supplementary information. These three notebooks should be run in the order listed.
* `simulation_scenario[1-3].R`: Scripts that generate synthetic polysomnogram summaries, fit the model to them, and save results.

