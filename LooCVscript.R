# BEFORE RUNNING THIS CODE READ THIS FIRST!!!!
# This is the second script that needs to be run.
# This produces LOO-CV results for the all Bayesian GLMMs for each of the markers listed below.

# The code is fully annotated and can be run all at once, HOWEVER...
# Please make sure you change the paths for the files that are read from this file
# and the files to which we are storing results into.  I STRONGLY suggest you read the 
# annotations in this code before running it so that you understand what you are doing.

library("loo")
# [1] "total_cells"                   "foxp3_opal_540_positive_cells" "cd3_opal_650_positive_cells"  
# [4] "cd8_opal_570_positive_cells"   "cd11b_opal_620_positive_cells" "cd15_opal_520_positive_cells" 
# [7] "cd3plus_foxp3plus_cells"       "cd3plus_cd8plus_cells"         "cd11bplus_cd15plus_cells" 

# This code needs to be run AFTER all the models have been generated through
# the script named Auxproc_miF_data_bayesian_models.R
source("BetaBinomial2Definitions.R")
load("BRMS_models_foxp3_opal_540_cells.RData")
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
dfcomp_foxp3 = as.data.frame(comp)

load("BRMS_models_cd3_opal_650_cells.RData")
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
dfcomp_cd3 = as.data.frame(comp)


load("BRMS_models_cd8_opal_570_cells.RData")
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
dfcomp_cd8 = as.data.frame(comp)

load("BRMS_models_cd11b_opal_620_cells.RData")
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
dfcomp_cd11b = as.data.frame(comp)

load("BRMS_models_cd15_opal_520_cells.RData")
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
dfcomp_cd15 = as.data.frame(comp)

load("BRMS_models_cd3plus_foxp3plus_cells.RData")
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
dfcomp_cd3_foxp3 = as.data.frame(comp)

load("BRMS_models_cd3plus_cd8plus_cells.RData")
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
dfcomp_cd3_cd8 = as.data.frame(comp)

load("BRMS_models_cd11bplus_cd15plus_cells.RData")
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


save(dfcomp_foxp3,
     dfcomp_cd3,
     dfcomp_cd8,
     dfcomp_cd11b,
     dfcomp_cd15,
     dfcomp_cd3_cd8,
     dfcomp_cd3_foxp3,
     dfcomp_cd11b_cd15,
     file = "LOO_CV_Results.RData")