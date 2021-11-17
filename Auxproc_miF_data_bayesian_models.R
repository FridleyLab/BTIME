# BEFORE RUNNING THIS CODE READ THIS FIRST!!!!
# This is the first script that needs to be run.
# This produces Bayesian GLMMs for the markers listed below.
# The code is fully annotated and can be run all at once, HOWEVER...
# Please make sure you change the paths for the files that are read from this file
# and the files to which we are storing results into.  I STRONGLY suggest you read the 
# annotations in this code before running it so that you understand what you are doing.

# Read the data in:
# dfmif = readRDS("Peres IF_AACES_NCOCS full data 12092020.rds")
dfmif = readRDS("Sample_mIF_data.rds")
# # Subset the data by intratumoral_i_vs_peripheral_p_ == intratumoral only:
# dfmif = subset(dfmif, intratumoral_i_vs_peripheral_p_ == "Intratumoral")
# # process stage variable (this is X = group in my model)
# # remove missing values
# rmlist = which(is.na(dfmif$stage))
# dfmif = dfmif[-rmlist,]
# # now do this as per Brooke/Chris: 1 = Localized, 2 = Regional, and 3 = Distant. Group 1&2 vs 3:
# dfmif$stage_group = as.factor(ifelse(dfmif$stage == "Distant", "Stage 3", "Stage 1 or 2"))
# # My variables (total cells = total tumor + total stroma) go from variables 43:51
# dfvars = as.data.frame(names(dfmif))
# # Need to find/replace the markers below:
# dfvars[43:51,]
# [1] "total_cells"                   "foxp3_opal_540_positive_cells" "cd3_opal_650_positive_cells"  
# [4] "cd8_opal_570_positive_cells"   "cd11b_opal_620_positive_cells" "cd15_opal_520_positive_cells" 
# [7] "cd3plus_foxp3plus_cells"       "cd3plus_cd8plus_cells"         "cd11bplus_cd15plus_cells" 

# To execute the code below you need to choose which of the markers above (foxp3, cd3, ..., etc.)
# is of your interest. Then go ahead and find/replace the markers on the code below with the exact same 
# name you want from the list above.
library(brms)
# Models for cd11bplus_cd15plus_cells:
# fit1 = Binomial
start_time = Sys.time()
fit1 = brm(cd11bplus_cd15plus_cells | trials(total_cells) ~ stage_group + (1|suid),
           data = dfmif, 
           family = binomial(), cores = 8,
           iter = 14000, 
           warmup = 4000)
end_time = Sys.time()
my_time = (end_time - start_time)
print(my_time)
summary(fit1)
time_fit1 = my_time
# fit2 = Poisson
start_time = Sys.time()
fit2 = brm(cd11bplus_cd15plus_cells ~ stage_group + (1|suid),
           data = dfmif, 
           family = poisson(), cores = 8,
           iter = 14000, 
           warmup = 4000)
end_time = Sys.time()
my_time = (end_time - start_time)
print(my_time)
summary(fit2)
time_fit2 = my_time
# fit3 = Negative Binomial
start_time = Sys.time()
fit3 = brm(cd11bplus_cd15plus_cells ~ stage_group + (1|suid),
           data = dfmif, 
           family = negbinomial(), cores = 8,
           iter = 14000, 
           warmup = 4000)
end_time = Sys.time()
my_time = (end_time - start_time)
print(my_time)
summary(fit3)
time_fit3 = my_time
# fit4 = zero_inflated_binomial
start_time = Sys.time()
fit4 = brm(cd11bplus_cd15plus_cells | trials(total_cells) ~ stage_group + (1|suid),
           data = dfmif, 
           family = zero_inflated_binomial, cores = 8,
           iter = 14000, 
           warmup = 4000)
end_time = Sys.time()
my_time = (end_time - start_time)
print(my_time)
summary(fit4)
time_fit4 = my_time
# fit5 = zero_inflated_poisson
start_time = Sys.time()
fit5 = brm(cd11bplus_cd15plus_cells ~ stage_group + (1|suid),
           data = dfmif, 
           family = zero_inflated_poisson(), cores = 8,
           iter = 14000, 
           warmup = 4000)
end_time = Sys.time()
my_time = (end_time - start_time)
print(my_time)
summary(fit5)
time_fit5 = my_time
# fit6 = zero_inflated_negbinomial
start_time = Sys.time()
fit6 = brm(cd11bplus_cd15plus_cells ~ stage_group + (1|suid),
           data = dfmif, 
           family = zero_inflated_negbinomial(), cores = 8,
           iter = 14000, 
           warmup = 4000)
end_time = Sys.time()
my_time = (end_time - start_time)
print(my_time)
summary(fit6)
time_fit6 = my_time
# fit7 = beta_binomial2
# To execute this model we need to source the below file first since the beta binomial
# is not part of the BRMS suite of families but rather a custom family coded in the
# BetaBinomial2Definitions.R script file.  Please CHANGE THE PATH to wherever you stored such file.
source("BetaBinomial2Definitions.R")
start_time = Sys.time()
fit7 = brm(cd11bplus_cd15plus_cells | trials(total_cells) ~ stage_group + (1|suid),
           data = dfmif, 
           family = beta_binomial2, cores = 8,
           stanvars = stanvars,
           iter = 14000, 
           warmup = 4000)
end_time = Sys.time()
my_time = (end_time - start_time)
print(my_time)
summary(fit7)
time_fit7 = my_time
# Show total running time:
sum(time_fit1,time_fit2,time_fit3,time_fit4,time_fit5,time_fit6,time_fit7)
# Clean up variables before saving the models and necessary extra info:
remove(list=c("dfmif", "start_time", "end_time", "rmlist", "my_time", "dfvars"))
# Save the image corresponding to the models ran above using the markers names above.
# (saving locally since I like this style of work while I test things.  The complete results will be stored
# in the folder mentioned at the end of this script)
save.image("BRMS_models_cd11bplus_cd15plus_cells.RData")

# To see estimation and diagnostics results interactively do:
# shinystan::launch_shinystan(fit)

# Our results folder is stored in:
# smb://hlm/data/lab/Lab_Fridley/bayesTIME/results/