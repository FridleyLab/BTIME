# Read the data in:
dfmif = readRDS("~/Documents/MyProjects/BayesianWorld/TestStan_files/Peres IF_AACES_NCOCS full data 12092020.rds")
# Subset the data by intratumoral_i_vs_peripheral_p_ == intratumoral only:
dfmif = subset(dfmif, intratumoral_i_vs_peripheral_p_ == "Intratumoral")
# process stage variable (this is X = group in my model)
# remove missing values
rmlist = which(is.na(dfmif$stage))
dfmif = dfmif[-rmlist,]
# now do this as per Brooke/Chris: 1 = Localized, 2 = Regional, and 3 = Distant. Group 1&2 vs 3:
dfmif$stage_group = as.factor(ifelse(dfmif$stage == "Distant", "Stage 3", "Stage 1 or 2"))
# My variables (total cells = total tumor + total stroma) go from variables 43:51
dfvars = as.data.frame(names(dfmif))
# Need to find/replace the markers below:
dfvars[43:51,]
# [1] "total_cells"                   "foxp3_opal_540_positive_cells" "cd3_opal_650_positive_cells"  
# [4] "cd8_opal_570_positive_cells"   "cd11b_opal_620_positive_cells" "cd15_opal_520_positive_cells" 
# [7] "cd3plus_foxp3plus_cells"       "cd3plus_cd8plus_cells"         "cd11bplus_cd15plus_cells"   
library(brms)

# fit1 = Binomial
# with 2 it took 22 minutes, with 3 it took 30 minutes, Time difference of 1.314031 hours for all of them
start_time = Sys.time()
fit1 = brm(
  mvbind(foxp3_opal_540_positive_cells, 
         cd3_opal_650_positive_cells, 
         cd8_opal_570_positive_cells,
         cd11b_opal_620_positive_cells,
         cd15_opal_520_positive_cells,
         cd3plus_foxp3plus_cells,
         cd3plus_cd8plus_cells,
         cd11bplus_cd15plus_cells
         ) | trials(total_cells) ~ stage_group + (1|p|suid),
  data = dfmif, 
  family = binomial(), cores = 8,
  iter = 14000, 
  warmup = 4000)
end_time = Sys.time()
my_time = (end_time - start_time)
print(my_time)
summary(fit1)
time_fit1 = my_time

# fit6 = zero_inflated_negbinomial
# with 2 it took 6 minutes, with 3 it took 10 minutes, Time difference of 35.85713 mins for all of them
start_time = Sys.time()
fit6 = brm(
  mvbind(foxp3_opal_540_positive_cells, 
         cd3_opal_650_positive_cells, 
         cd8_opal_570_positive_cells,
         cd11b_opal_620_positive_cells,
         cd15_opal_520_positive_cells,
         cd3plus_foxp3plus_cells,
         cd3plus_cd8plus_cells,
         cd11bplus_cd15plus_cells) ~ stage_group + (1|p|suid),
           data = dfmif, 
           family = zero_inflated_negbinomial(), cores = 8,
           iter = 14000, 
           warmup = 4000)
end_time = Sys.time()
my_time = (end_time - start_time)
print(my_time)
summary(fit6)
time_fit6 = my_time

# fit3 = Negative Binomial
start_time = Sys.time()
fit3 = brm(
  mvbind(foxp3_opal_540_positive_cells, 
         cd3_opal_650_positive_cells, 
         cd8_opal_570_positive_cells,
         cd11b_opal_620_positive_cells,
         cd15_opal_520_positive_cells,
         cd3plus_foxp3plus_cells,
         cd3plus_cd8plus_cells,
         cd11bplus_cd15plus_cells) ~ stage_group + (1|suid),
  data = dfmif, 
  family = negbinomial(), cores = 8,
  iter = 14000, 
  warmup = 4000)
end_time = Sys.time()
my_time = (end_time - start_time)
print(my_time)
summary(fit3)
time_fit3 = my_time

# fit2 = Poisson
start_time = Sys.time()
fit2 = brm(
  mvbind(foxp3_opal_540_positive_cells, 
         cd3_opal_650_positive_cells, 
         cd8_opal_570_positive_cells,
         cd11b_opal_620_positive_cells,
         cd15_opal_520_positive_cells,
         cd3plus_foxp3plus_cells,
         cd3plus_cd8plus_cells,
         cd11bplus_cd15plus_cells) ~ stage_group + (1|suid),
           data = dfmif, 
           family = poisson(), cores = 8,
           iter = 14000, 
           warmup = 4000)
end_time = Sys.time()
my_time = (end_time - start_time)
print(my_time)
summary(fit2)
time_fit2 = my_time

# fit4 = zero_inflated_binomial
start_time = Sys.time()
fit4 = brm(
  mvbind(foxp3_opal_540_positive_cells, 
         cd3_opal_650_positive_cells, 
         cd8_opal_570_positive_cells,
         cd11b_opal_620_positive_cells,
         cd15_opal_520_positive_cells,
         cd3plus_foxp3plus_cells,
         cd3plus_cd8plus_cells,
         cd11bplus_cd15plus_cells
  ) | trials(total_cells) ~ stage_group + (1|p|suid),
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
fit5 = brm(
  mvbind(foxp3_opal_540_positive_cells, 
         cd3_opal_650_positive_cells, 
         cd8_opal_570_positive_cells,
         cd11b_opal_620_positive_cells,
         cd15_opal_520_positive_cells,
         cd3plus_foxp3plus_cells,
         cd3plus_cd8plus_cells,
         cd11bplus_cd15plus_cells) ~ stage_group + (1|suid),
           data = dfmif, 
           family = zero_inflated_poisson(), cores = 8,
           iter = 14000, 
           warmup = 4000)
end_time = Sys.time()
my_time = (end_time - start_time)
print(my_time)
summary(fit5)
time_fit5 = my_time

# fit7 = beta_binomial2
start_time = Sys.time()
fit7 = brm(mvbind(foxp3_opal_540_positive_cells, 
                  cd3_opal_650_positive_cells, 
                  cd8_opal_570_positive_cells,
                  cd11b_opal_620_positive_cells,
                  cd15_opal_520_positive_cells,
                  cd3plus_foxp3plus_cells,
                  cd3plus_cd8plus_cells,
                  cd11bplus_cd15plus_cells) | trials(total_cells) ~ stage_group + (1|suid),
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
dfcomp = as.data.frame(comp)

# Show total running time:
sum(time_fit1,time_fit2,time_fit3,time_fit4,time_fit5,time_fit6,time_fit7)
# Clean up variables before saving the models and necessary extra info:
remove(list=c("dfmif", "start_time", "end_time", "rmlist", "my_time", "dfvars"))
# Save the image corresponding to the models ran above using the markers names above:
save.image("~/Documents/MyProjects/BayesianWorld/TestStan_files/BRMS_Multi_Models.RData")


library(ggplot2)
library(gridExtra)
library(dplyr)
dfcomp$my_models = rownames(dfcomp) 
dfcomp = dfcomp %>% 
  mutate(the_models = case_when(
    my_models == 'fit1' ~ 'B',
    my_models == 'fit2' ~ 'P',
    my_models == 'fit3' ~ 'NB',
    my_models == 'fit4' ~ 'ZIB',
    my_models == 'fit5' ~ 'ZIP',
    my_models == 'fit6' ~ 'ZINB',
    TRUE ~ 'BB'))
p_multi = ggplot(dfcomp) +
  geom_bar(aes(x=reorder(the_models, elpd_diff), y=elpd_diff), 
           stat="identity", 
           fill="skyblue", 
           alpha=0.5) +
  geom_errorbar(aes(x=reorder(the_models, elpd_diff), 
                    ymin=elpd_diff-se_diff, ymax=elpd_diff+se_diff), 
                width=0.4, 
                colour="orange", 
                alpha=0.9, 
                size=1.3) +
  xlab("Models") +
  ggtitle("Multi Response Models") +
  coord_flip()
p_multi