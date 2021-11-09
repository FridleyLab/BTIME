# BEFORE RUNNING THIS CODE READ THIS FIRST!!!!
# The code is fully annotated and can be run all at once, HOWEVER...
# Please make sure you change the paths for the files that are read from this file
# and the files to which we are storing results into.  I STRONGLY suggest you read the 
# annotations in this code before running it so that you understand what you are doing.

# library("bayesplot")
library("loo")
# [1] "total_cells"                   "foxp3_opal_540_positive_cells" "cd3_opal_650_positive_cells"  
# [4] "cd8_opal_570_positive_cells"   "cd11b_opal_620_positive_cells" "cd15_opal_520_positive_cells" 
# [7] "cd3plus_foxp3plus_cells"       "cd3plus_cd8plus_cells"         "cd11bplus_cd15plus_cells" 

# To execute the code below you need to choose which of the markers above (foxp3, cd3, ..., etc.)
# is of your interest. Then go ahead and find/replace the markers on the code below with the exact same 
# name you want from the list above.
load("~/Documents/MyProjects/BayesianWorld/TestStan_files/BRMS_models_cd11bplus_cd15plus_cells.RData")
source("~/Documents/MyProjects/BayesianWorld/TestStan_files/BetaBinomial2Definitions.R")

start_time = Sys.time()
loo1 = loo(fit1, save_psis = FALSE, cores = 12) #, moment_match = TRUE) crashes de R session.  No solution for this!
loo2 = loo(fit2, save_psis = FALSE, cores = 12) #, moment_match = TRUE)
loo3 = loo(fit3, save_psis = FALSE, cores = 12) #, moment_match = TRUE)
loo4 = loo(fit4, save_psis = FALSE, cores = 12) #, moment_match = TRUE)
loo5 = loo(fit5, save_psis = FALSE, cores = 12) #, moment_match = TRUE)
loo6 = loo(fit6, save_psis = FALSE, cores = 12) #, moment_match = TRUE)
expose_functions(fit7, vectorize = TRUE)
loo7 = loo(fit7, save_psis = FALSE, cores = 12) #, moment_match = TRUE)
end_time = Sys.time()
my_time = (end_time - start_time)

comp = loo_compare(loo1, loo2, loo3, loo4, loo5, loo6, loo7)
print(comp, digits = 2)
end_time = Sys.time()
my_time = (end_time - start_time)
print(my_time)
dfcomp_cd11b_cd15 = as.data.frame(comp)