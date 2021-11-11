# BTIME
Bayesian analysis of multiplex immunofluorescence data with 8 different models with various levels of zero-inflation or over-dispersion modeled:

- `fit1`: Binomial (B)
- `fit2`: Poisson (P)
- `fit3`: Negative Binomial (NB)
- `fit4`: Zero Inflated Binomial (ZIB)
- `fit5`: Zero Inflated Poisson (ZIP)
- `fit6`: Zero Inflated Negative Binomial (ZINB)
- `fit7`: Beta-Binomial (BB)
- `fit8`: Zero Inflated Beta-Binomial (ZIBB which is still pending for now)

## Instructions:

First you need to install the following libraries:

- `rstan`: For Hamiltonian MC (will install `Rcpp` and related libraries)
- `brms`: To build Bayesian GLMMs using `rstan` as a back end (will install `rstan` but I strongly recommend you install it in advanced since it is not a trivial installation)
- `loo`: To perform post-processing for model comparison using LOO-CV.
- `ggplot2`: For plotting
- `gridExtra`: For making grid plots
- `dplyr`: For data manipulation.

The script `Auxproc_miF_data_bayesian_models.R` needs to be executed first for **every** marker from the following list:

- foxp3_opal_540_positive_cells
- cd3_opal_650_positive_cells
- cd8_opal_570_positive_cells
- cd11b_opal_620_positive_cells
- cd15_opal_520_positive_cells
- cd3plus_foxp3plus_cells
- cd3plus_cd8plus_cells
- cd11bplus_cd15plus_cells

Note that to execute the code above you need to choose which of the markers above (foxp3, cd3, ..., etc.) is of your interest.  If you want the subsequent scripts to compare the models for ALL markers then you need to go ahead and find/replace the markers on the script with the exact same name you want from the list above.  Once ran for every marker you should have ended with a stored number of files like these:

- `BRMS_models_foxp3_opal_540_positive_cells.RData`
- `BRMS_models_cd3_opal_650_positive_cells.RData`
- `BRMS_models_cd8_opal_570_positive_cells.RData`
- `BRMS_models_cd11b_opal_620_positive_cells.RData`
- `BRMS_models_cd15_opal_520_positive_cells.RData`
- `BRMS_models_cd3plus_foxp3plus_cells.RData`
- `BRMS_models_cd3plus_cd8plus_cells.RData`
- `BRMS_models_cd11bplus_cd15plus_cells.RData`

When this is done the data files containing the Bayesian GLMMs listed above will be generated and stored in your local folder. After that, execute the script named `LooCVscript.R` which (using the `loo` package) will generate comparison between models for each marker and store the results in a file called `LOO_CV_Results.RData`.
Finally run the script named `Plot_ELP_diff.R` to produce a nice graphical output/plot of all the models for all the markers called `LOO_plot_BGLMM.pdf`.
