# BEFORE RUNNING THIS CODE READ THIS FIRST!!!!
# This is the third script that needs to be run.
# This script produces a nice LOO-CV plot from the results stored in the
# file called "LOO_CV_Results.RData"

# The code is fully annotated and can be run all at once, HOWEVER...
# Please make sure you change the paths for the files that are read from this file
# and the files to which we are storing results into.  I STRONGLY suggest you read the 
# annotations in this code before running it so that you understand what you are doing.

# [1] "total_cells"                   "foxp3_opal_540_positive_cells" "cd3_opal_650_positive_cells"  
# [4] "cd8_opal_570_positive_cells"   "cd11b_opal_620_positive_cells" "cd15_opal_520_positive_cells" 
# [7] "cd3plus_foxp3plus_cells"       "cd3plus_cd8plus_cells"         "cd11bplus_cd15plus_cells" 

# This code needs to be run AFTER all the LOO results for all the models have been generated through
# the script named LooCVscript.R

library(ggplot2)
library(gridExtra)
library(dplyr)
load("LOO_CV_Results.RData")
dfcomp_foxp3$my_models = rownames(dfcomp_foxp3) 
dfcomp_foxp3 = dfcomp_foxp3 %>% 
  mutate(the_models = case_when(
    my_models == 'fit1' ~ 'B',
    my_models == 'fit2' ~ 'P',
    my_models == 'fit3' ~ 'NB',
    my_models == 'fit4' ~ 'ZIB',
    my_models == 'fit5' ~ 'ZIP',
    my_models == 'fit6' ~ 'ZINB',
    TRUE ~ 'BB'))
dfcomp_cd3$my_models = rownames(dfcomp_cd3) 
dfcomp_cd3 = dfcomp_cd3 %>% 
  mutate(the_models = case_when(
    my_models == 'fit1' ~ 'B',
    my_models == 'fit2' ~ 'P',
    my_models == 'fit3' ~ 'NB',
    my_models == 'fit4' ~ 'ZIB',
    my_models == 'fit5' ~ 'ZIP',
    my_models == 'fit6' ~ 'ZINB',
    TRUE ~ 'BB'))
dfcomp_cd8$my_models = rownames(dfcomp_cd8) 
dfcomp_cd8 = dfcomp_cd8 %>% 
  mutate(the_models = case_when(
    my_models == 'fit1' ~ 'B',
    my_models == 'fit2' ~ 'P',
    my_models == 'fit3' ~ 'NB',
    my_models == 'fit4' ~ 'ZIB',
    my_models == 'fit5' ~ 'ZIP',
    my_models == 'fit6' ~ 'ZINB',
    TRUE ~ 'BB'))
dfcomp_cd15$my_models = rownames(dfcomp_cd15) 
dfcomp_cd15 = dfcomp_cd15 %>% 
  mutate(the_models = case_when(
    my_models == 'fit1' ~ 'B',
    my_models == 'fit2' ~ 'P',
    my_models == 'fit3' ~ 'NB',
    my_models == 'fit4' ~ 'ZIB',
    my_models == 'fit5' ~ 'ZIP',
    my_models == 'fit6' ~ 'ZINB',
    TRUE ~ 'BB'))
dfcomp_cd11b$my_models = rownames(dfcomp_cd11b) 
dfcomp_cd11b = dfcomp_cd11b %>% 
  mutate(the_models = case_when(
    my_models == 'fit1' ~ 'B',
    my_models == 'fit2' ~ 'P',
    my_models == 'fit3' ~ 'NB',
    my_models == 'fit4' ~ 'ZIB',
    my_models == 'fit5' ~ 'ZIP',
    my_models == 'fit6' ~ 'ZINB',
    TRUE ~ 'BB'))
dfcomp_cd3_cd8$my_models = rownames(dfcomp_cd3_cd8) 
dfcomp_cd3_cd8 = dfcomp_cd3_cd8 %>% 
  mutate(the_models = case_when(
    my_models == 'fit1' ~ 'B',
    my_models == 'fit2' ~ 'P',
    my_models == 'fit3' ~ 'NB',
    my_models == 'fit4' ~ 'ZIB',
    my_models == 'fit5' ~ 'ZIP',
    my_models == 'fit6' ~ 'ZINB',
    TRUE ~ 'BB'))
dfcomp_cd3_foxp3$my_models = rownames(dfcomp_cd3_foxp3) 
dfcomp_cd3_foxp3 = dfcomp_cd3_foxp3 %>% 
  mutate(the_models = case_when(
    my_models == 'fit1' ~ 'B',
    my_models == 'fit2' ~ 'P',
    my_models == 'fit3' ~ 'NB',
    my_models == 'fit4' ~ 'ZIB',
    my_models == 'fit5' ~ 'ZIP',
    my_models == 'fit6' ~ 'ZINB',
    TRUE ~ 'BB'))
dfcomp_cd11b_cd15$my_models = rownames(dfcomp_cd11b_cd15) 
dfcomp_cd11b_cd15 = dfcomp_cd11b_cd15 %>% 
  mutate(the_models = case_when(
    my_models == 'fit1' ~ 'B',
    my_models == 'fit2' ~ 'P',
    my_models == 'fit3' ~ 'NB',
    my_models == 'fit4' ~ 'ZIB',
    my_models == 'fit5' ~ 'ZIP',
    my_models == 'fit6' ~ 'ZINB',
    TRUE ~ 'BB'))

# FoxP3
p_foxp3 = ggplot(dfcomp_foxp3) +
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
  ggtitle("FoxP3") +
  coord_flip()

# CD3
p_cd3 = ggplot(dfcomp_cd3) +
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
  ggtitle("CD3") +
  coord_flip()

# CD8
p_cd8 = ggplot(dfcomp_cd8) +
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
  ggtitle("CD8") +
  coord_flip()

# CD11b
p_cd11b = ggplot(dfcomp_cd11b) +
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
  ggtitle("CD11b") +
  coord_flip()

# CD15
p_cd15 = ggplot(dfcomp_cd15) +
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
  ggtitle("CD15") +
  coord_flip()

# CD3 + Foxp3
p_cd3_foxp3 = ggplot(dfcomp_cd3_foxp3) +
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
  ggtitle("CD3 + Foxp3") +
  coord_flip()

# CD3 + CD8
p_cd3_cd8 = ggplot(dfcomp_cd3_cd8) +
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
  ggtitle("CD3 + CD8") +
  coord_flip()

# CD11b + CD15
p_cd11b_cd15 = ggplot(dfcomp_cd11b_cd15) +
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
  ggtitle("CD11b + CD15") +
  coord_flip()



LOO_plot_BGLMM = grid.arrange(
  grobs = list(p_foxp3, p_cd3, p_cd8, 
               p_cd11b, p_cd15, p_cd3_foxp3,
               p_cd3_cd8, p_cd11b_cd15),
  widths = c(1, 1, 1),
  layout_matrix = rbind(c(1, 2, 3),
                        c(4, 5, 6),
                        c(7, 8, NA))
)

ggsave("LOO_plot_BGLMM.pdf")