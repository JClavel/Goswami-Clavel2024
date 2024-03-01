# Morphological Evolution in a time of Phenomics
Approaches to study phenotypic/phenomic evolution used for the figures in Goswami & Clavel 2024

> "dtt_plots.R" file shows how to generate DTT with expectations from model fit (by penalized likelihood) on the high-dimensional datasets

> "stacked_pca.R" file is constructing the stacked PCA for the mammals skull dataset from Goswami et al. 2022 with expectations from both a Brownian motion and Pagel's lambda with variable rate model.

> "simulations_OU_clim.R" performs simulations under the OU climatic model (where the optimum follow a climatic curve through time) using the "fit_t_env_ou" function implemented in RPANDA (previous versions of the model are available on my gitHub). The model is compared to Brownian motion, Early Burst, OU, and BM+trend models

> "fun.R" contains wrapper functions used in the "dtt_plots.R" and the "stacked_pca.R" codes


Goswami, A., Clavel., J., 2024. Morphological Evolution in a time of Phenomics.
