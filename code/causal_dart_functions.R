###################################
#           Functions             #
###################################

transform_y = function(y, p, TorC){
  #accepts a obsrved response y, propensity 
  #score p and treatment assignment indicator TorC
  output = y*(TorC - p)/(p*(1-p))
  return(output)
}

get_matrix_transformed_response = function(propensity_matrix, TorC,response){
  
  M = ncol(propensity_matrix)
  N = nrow(propensity_matrix)
  propensity_matrix[which(propensity_matrix <= 0)] = 0.000001
  propensity_matrix[which(propensity_matrix >= 1)] = 0.999999
  transformed_matrix = matrix(NA, nrow = N, ncol = M)
  for(i in 1:M){
    p_col = propensity_matrix[,i]
    transformed_matrix[,i] = transform_y(response, p_col, TorC)
  }
  return(transformed_matrix)
}

get_stacked_transformed_response = function(dat, 
                                            propensity_matrix, 
                                            TorC, response){
  
  M = ncol(propensity_matrix)
  N = nrow(propensity_matrix)
  dat$y_transform = NA
  output = list()
  for(i in 1:M){
    print(i)
    p_col = propensity_matrix[,i]
    dat$y_transform = transform_y(response, p_col, TorC)
    output[[i]] = dat
  }
  return(output)
}

convert_to_sd_units = function(fit_matrix){
  #takes matrix of CATE's and converts to units of standard deviation
  return
}

get_train_test_cv = function(dat, fold){
  N = nrow(dat)
  start_index = round(seq(1, N, by = N/fold))
  stop_index  = rev(round(seq(N, 1, by = -N/fold)))
  output = list()
  for(i in 1:fold){
    test_index = c(start_index[i]:stop_index[i])
    output[[i]] = list()
    output[[i]]$test_index = test_index
    output[[i]]$train      = dat[-test_index,]
    output[[i]]$test       = dat[ test_index,]
  }
  return(output)
}

est_coverage = function(mat, Tau){
  MIN = apply(mat, 2, min)
  MAX = apply(mat, 2, max)
  
  index = which( Tau > MIN & Tau < MAX )
  count = length(index)
  covered_ratio = count/length(Tau)
  return(covered_ratio)
}