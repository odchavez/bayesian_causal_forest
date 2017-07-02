rm(list = ls())
options(java.parameters = "-Xmx256g")
library(dartMachine)
library(foreach)
library(doParallel)
set.seed(123)
setwd("/home/ochavez/bayesian_causal_forest")
source("code/causal_dart_functions.R")

################
#temp functions#
################

causal_dart = function(dat, train_TorC, NUM_AFTER_BURN_IN, NUM_BEFORE_BURN_IN, TREE_NUM){
  
  dart_TorC_2  <- bartMachine(X = as.data.frame(train_TorC$X), y = train_TorC$Y,
                              do_ard = TRUE,
                              num_trees = TREE_NUM,
                              num_burn_in = NUM_BEFORE_BURN_IN,
                              num_iterations_after_burn_in = NUM_AFTER_BURN_IN,
                              do_prior = TRUE,
                              mem_cache_for_speed = FALSE)
  
  ppps_2 = bart_machine_get_posterior(dart_TorC_2, train_TorC$X)
  propensity_matrix_2   <- ppps_2$y_hat_posterior_samples[,sample(1:ncol(ppps_2$y_hat_posterior_samples), p_val_num, replace = TRUE)]
  transformed_dat_2 = get_stacked_transformed_response(dat, propensity_matrix_2, dat$T, dat$Y)
  rep_num = length(transformed_dat_2)
  cores=detectCores()
  cl <- makeCluster(cores[1]-1)
  registerDoParallel(cl)
  cate_dart_out = foreach(i = 1:rep_num, .packages = "dartMachine") %dopar%{
    #cate_dart_out = list()
    #for(i in 1:rep_num){
    
    dart_transformed_2d        = bartMachine(X = transformed_dat_2[[i]][,which(names(transformed_dat_2[[i]]) %in% paste0("V.",1:994) )], 
                                             y = transformed_dat_2[[i]]$y_transform, 
                                             mem_cache_for_speed = FALSE,
                                             do_ard = TRUE,
                                             num_trees = TREE_NUM,
                                             num_burn_in = NUM_BEFORE_BURN_IN, 
                                             num_iterations_after_burn_in = NUM_AFTER_BURN_IN,
                                             do_prior = TRUE)
    
    bart_machine_get_posterior(dart_transformed_2d, 
                               transformed_dat_2[[i]]
                               [,which(names(transformed_dat_2[[i]]) 
                                       %in% paste0("V.",1:994) 
                                       )
                                ]
                               )$y_hat_posterior_samples
  }
  
  stopCluster(cl)
  cate_dart = organize_causal_posterior(cate_dart_out)#matrix(unlist(cate_dart_out), ncol = rep_num)
  return(cate_dart)
}

select_percent_T = function(percent, nrow){
  
  pn = floor(percent*nrow)
  
  T_index = sample(1:nrow, pn, replace = FALSE)
  
  output = rep(0, nrow)
  output[T_index] = 1
  output = factor(output)
  return(output)
  
}

organize_causal_posterior = function(cate_dart_out){
    S = length(cate_dart_out)
    for(j in 1:S){
        if(j == 1){output = cate_dart_out[[j]]}
        else{ output = cbind(output,cate_dart_out[[j]])}
    }
    output = t(output)
    return(output)  
}
##############################################
#                   CODE                     #
##############################################

#using simulation.2 data from Athey CAT paper
cov_perc = c(0.05,0.95)
pred_num = 100
std_per = 0.64
set_val = 2
NUM_AFTER_BURN_IN = 2000
NUM_BEFORE_BURN_IN = 500
TREE_NUM = 50
p_val_num = 100
show_num  = 100
N = 10000
## Fit causal dart
squash_val = c(0.01, seq(0.05, 0.95, by = 0.05), 0.99)
coverage = rep(NA, length(squash_val))
interval_size = rep(NA, length(squash_val))
for(i in 1:length(squash_val)){
  dat = read.csv(paste0("data/caus_sim_1_N=10000_P=1000_noise_percent=",std_per,"_squash_val_",squash_val[i],".csv"))[1:N,1:(6+pred_num)]
  cnames = c("X","Tau","Y0_given_X","Y1_given_X","Y","T")
  train_TorC = list()
  train_TorC$X = dat[,-which(names(dat) %in% cnames)]
  train_TorC$Y = factor(dat$T)
  #temp = train_TorC
  #temp$Y = select_percent_T(squash_val[i], length(temp$Y))
  cate_dart        = causal_dart(dat, train_TorC, NUM_AFTER_BURN_IN, NUM_BEFORE_BURN_IN, TREE_NUM)
  coverage[i]      = est_coverage(cate_dart, dat[, "Tau"], cov_perc)
  interval_size[i] = est_coverage_size_median(cate_dart, cov_perc)

  ord_index = order(dat$Tau, decreasing = TRUE)
  ord_Tua = dat$Tau[ord_index]
  even_space_100 = seq(1, N, by = N/100)
  pdf(paste0("images/dart_dart_ind_effects_squash_val=",squash_val[i],"_perc_T=",sum(dat$T)/nrow(dat),"_N=",N,"_P=",pred_num,".pdf"))  
  boxplot(cate_dart[sample(nrow(cate_dart), 2000, replace = FALSE),ord_index[even_space_100] ], 
          ylim = c(-1,1), 
          main = "dart only")
          #sub = paste("coverage ratio = ", cov_ratio)) 
  points(1:length(even_space_100) ,ord_Tua[even_space_100], col = 'red', pch = 4)
  abline(h = 0, col = 'red')
  dev.off()    
}


pdf(paste0("images/dart_dart_coverage_vs_squash_val","_N=",N,"_P=",pred_num,".pdf"))  
par(mfrow = c(1,2))
plot(squash_val, coverage)
plot(squash_val, interval_size)
dev.off()
#pdf(paste0("images/dart_dart_ind_effects.pdf"))  
#boxplot(t(cate_dart[1:show_num,]), 
#        ylim = c(-1,1), 
#        main = "dart only")
#        #sub = paste("coverage ratio = ", cov_ratio)) 
#points(1:show_num ,dat[1:show_num, "Tau"], col = 'red', pch = 4)
#dev.off()







