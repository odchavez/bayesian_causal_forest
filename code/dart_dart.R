rm(list = ls())
options(java.parameters = "-Xmx256g")
library(dartMachine)
library(foreach)
library(doParallel)
set.seed(123)
setwd("/home/ochavez/bayesian_causal_forest")
source("code/causal_dart_functions.R")

##############################################
#                   CODE                     #
##############################################

#using simulation.2 data from Athey CAT paper
std_per = 0.2
set_val = 1
dat = read.csv(paste0("caus_inf/data/caus_sim_1_N=10000_P=1000_noise_percent=",std_per,"_set_",set_val,".csv"))
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

dart_TorC_2  <- bartMachine(X = as.data.frame(train_TorC$X), y = train_TorC$Y,
                            do_ard = TRUE,
                            num_trees = 200,
                            num_burn_in = 4000,
                            num_iterations_after_burn_in = NUM_AFTER_BURN_IN,
                            do_prior = TRUE,
                            mem_cache_for_speed = TRUE)
ppps_2 = bart_machine_get_posterior(dart_TorC_2, train_TorC$X)
propensity_matrix_2   <- ppps_2$y_hat_posterior_samples[,sample(1:ncol(ppps_2$y_hat_posterior_samples), p_val_num, replace = TRUE)]
transformed_dat_2 = get_stacked_transformed_response(dat, propensity_matrix_2, dat$T, dat$Y)
rep_num = length(transformed_dat_2)
cores=detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)
cate_dart_out = foreach(i = 1:rep_num, .packages = "dartMachine") %dopar%{
  
  dart_transformed_2d        = bartMachine(X = transformed_dat_2[[i]][,which(names(transformed_dat_2[[i]]) %in% paste0("V.",1:20) )], 
                                           y = transformed_dat_2[[i]]$y_transform, 
                                           mem_cache_for_speed = TRUE,
                                           do_ard = TRUE,
                                           num_trees = 200,
                                           num_burn_in = 4000, 
                                           num_iterations_after_burn_in = NUM_AFTER_BURN_IN,
                                           do_prior = TRUE)
  
  bart_machine_get_posterior(dart_transformed_2d, 
                                transformed_dat_2[[i]]
                                    [,which(names(transformed_dat_2[[i]]) 
                                        %in% paste0("V.",1:20) 
                                            )
                                    ]
                            )
}

stopCluster(cl)
cate_dart = matrix(unlist(cate_dart_out), ncol = rep_num)
cov_ratio = est_coverage(t(cate_dart[1:show_num,]), dat[1:show_num, "Tau"])

pdf(paste0("images/dart_dart_",std_per,".pdf"))  
boxplot(t(cate_dart[1:show_num,]), 
        ylim = c(-1,1), 
        main = "dart only",
        sub = paste("coverage ratio = ", cov_ratio)) 
points(1:show_num ,dat[1:show_num, "Tau"], col = 'red', pch = 4)
dev.off()







