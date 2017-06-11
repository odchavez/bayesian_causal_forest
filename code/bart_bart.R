rm(list = ls())
options(java.parameters = "-Xmx256g")
library(bartMachine)
library(foreach)
library(doParallel)
set.seed(123)
setwd("/work/03330/ochavez/wrangler")
source("caus_inf/code/causal_dart_functions.R")


##############################################
#                   CODE                     #
##############################################

#using simulation.2 data from Athey CAT paper
std_per = 0.2
set_val = 1
dat = read.csv("caus_inf/data/mini_test_0.2.csv")# read.csv(paste0("caus_inf/data/caus_sim_1_N=10000_P=1000_noise_percent=",std_per,"_set_",set_val,".csv"))#read.csv("caus_inf/data/mini_test_0.2.csv")#read.csv(paste0("caus_inf/data/caus_sim_1_N=10000_P=1000_noise_percent=",std_per,"_set_",set_val,".csv"))
cnames = c("X","Tau","Y0_given_X","Y1_given_X","Y","T")
train_TorC = list()
train_TorC$X = dat[,-which(names(dat) %in% cnames)]
train_TorC$Y = factor(dat$T)
NUM_AFTER_BURN_IN = 2000
NUM_BEFORE_BURN_IN = 500
TREE_NUM = 50
p_val_num = 50
show_num  = 100
## Fit bart
bart_TorC <- bartMachine(X = as.data.frame(train_TorC$X), 
                         y = train_TorC$Y, 
                         num_trees = TREE_NUM,
                         num_burn_in = NUM_BEFORE_BURN_IN,
                         num_iterations_after_burn_in = NUM_AFTER_BURN_IN,
                         mem_cache_for_speed = FALSE)
ppps = bart_machine_get_posterior(bart_TorC, train_TorC$X)

propensity_matrix = ppps$y_hat_posterior_samples[,sample(1:ncol(ppps$y_hat_posterior_samples), p_val_num, replace = TRUE)]
transformed_dat   = get_stacked_transformed_response(dat, propensity_matrix, dat$T, dat$Y)
pvars = setdiff(names(transformed_dat[[1]]), c("X","Tau","Y0_given_X","Y1_given_X","Y","T","y_transform"))
rep_num = length(transformed_dat)
cate_bart = matrix(NA, nrow = nrow(dat), ncol = rep_num)
posterior_list = list()

cores=detectCores()
cl <- makeCluster(cores[1])
registerDoParallel(cl)
cate_bart_out = foreach(i = 1:rep_num, .packages = "bartMachine") %dopar%{
#for(i in 1:rep_num){
  bart_transformed      <- bartMachine(X = transformed_dat[[i]][,which(names(transformed_dat[[i]]) %in% pvars )], 
                                       y = transformed_dat[[i]]$y_transform, 
                                       num_trees = TREE_NUM,
                                       num_burn_in = NUM_BEFORE_BURN_IN,
                                       num_iterations_after_burn_in = NUM_AFTER_BURN_IN,
                                       mem_cache_for_speed = FALSE)
  posterior_list[[i]]= bart_machine_get_posterior(bart_transformed, transformed_dat[[i]][,which(names(transformed_dat[[i]]) %in% pvars )])
  #cate_bart[,i] = posterior_list[[i]]$y_hat
  posterior_list[[i]]$y_hat
  #list()
#}
}
stopCluster(cl)
cate_bart = matrix(unlist(cate_bart_out), ncol = rep_num)
cov_ratio = est_coverage(t(cate_bart[1:show_num,]), dat[1:show_num, "Tau"])

pdf(paste0("caus_inf/plots/bart_bart_",std_per,".pdf")) 
boxplot(t(cate_bart[1:show_num,]), 
        ylim = c(-2,2), 
        main = "bart only",
        sub = paste("coverage ratio = ", cov_ratio)) 
points(1:show_num ,dat[1:show_num, "Tau"], col = 'red', pch = 4)
dev.off()






