make_data = function(N, P, noise_std, seed){
	set.seed(seed)
	X = matrix(NA, nrow = N, ncol = P)
	colnames(X) = paste0("V.",1:P)
	T = rep(NA, N)
	#fill covariates
	for(cnum in 1:P){
		X[,cnum] = runif(N, -2,2)
		#print(paste("cnum = ", cnum))
	}
	#T|x 
	term = 0.4*X[,1]*X[,2]  + 
		0.4*X[,3]*X[,4]  + 
		0.1*X[,1]*X[,6]  + 
		0.1*X[,1]*X[,7]  + 
		0.4*X[,8]*X[,9]  + 
		0.5*X[,5]*X[,10] + 
		0.5*X[,5]*X[,5]  - 1
	prob_T = exp(term)/( 1 + exp(term))
	#hist(prob_T)
	for(i in 1:N){
		T[i] = rbinom(1,1,prob_T[i])
	}
	y1_term = 0.5*X[,3]*X[,4]
	Y1_given_X = exp(y1_term)/(1+exp(y1_term))*
				dnorm(3+0.5*X[,2]*X[,5]+0.5*X[,1]^2, 0.5) +
				1/(1 + exp(y1_term)) * 
				dnorm(-0.5+0.5*X[,2]^2 - 0.5*X[,1]*X[,3], 0.8)
	Y0_given_X = exp(-abs(X[,5]))*dnorm((0.2*rowSums(X[,1:5]))^4, 1) + 
				(1 - exp(-abs(X[,5])))*dnorm(2+0.2*rowSums(X[,1:5]^2), 1)
	s = sd(c(Y1_given_X, Y0_given_X))
	noise = rnorm(N, mean = 0, sd = s*noise_std)			
	Y = T*Y1_given_X + (1-T)*Y0_given_X + noise
	Tau = Y1_given_X - Y0_given_X
	#hist(Tau,50)
	output = cbind(Tau, Y0_given_X,Y1_given_X, Y, T, X)
	return(output)
}

make_data_w_imbalance = function(N, P, noise_std, squash_val,  seed){
  set.seed(seed)
  X = matrix(NA, nrow = N, ncol = P)
  colnames(X) = paste0("V.",1:P)
  T = rep(NA, N)
  #fill covariates
  for(cnum in 1:P){
    X[,cnum] = runif(N, -2,2)
    #print(paste("cnum = ", cnum))
  }
  #T|x 
  term = 0.4*X[,1]*X[,2]  + 
    0.4*X[,3]*X[,4]  + 
    0.1*X[,1]*X[,6]  + 
    0.1*X[,1]*X[,7]  + 
    0.4*X[,8]*X[,9]  + 
    0.5*X[,5]*X[,10] + 
    0.5*X[,5]*X[,5]  - 1
  prob_T = squash_val * exp(term)/( 1 + exp(term))
  #hist(prob_T)
  for(i in 1:N){
    T[i] = rbinom(1,1,prob_T[i])
  }
  y1_term = 0.5*X[,3]*X[,4]
  Y1_given_X = exp(y1_term)/(1+exp(y1_term))*
    dnorm(3+0.5*X[,2]*X[,5]+0.5*X[,1]^2, 0.5) +
    1/(1 + exp(y1_term)) * 
    dnorm(-0.5+0.5*X[,2]^2 - 0.5*X[,1]*X[,3], 0.8)
  Y0_given_X = exp(-abs(X[,5]))*dnorm((0.2*rowSums(X[,1:5]))^4, 1) + 
    (1 - exp(-abs(X[,5])))*dnorm(2+0.2*rowSums(X[,1:5]^2), 1)
  s = sd(c(Y1_given_X, Y0_given_X))
  noise = rnorm(N, mean = 0, sd = s*noise_std)			
  Y = T*Y1_given_X + (1-T)*Y0_given_X + noise
  Tau = Y1_given_X - Y0_given_X
  #hist(Tau,50)
  output = cbind(Tau, Y0_given_X,Y1_given_X, Y, T, X)
  return(output)
}

N = 10^4
P = 1000
per_std = 0.64
squash_val = c(0.01,seq(0.05,0.95, by = 0.05), 0.99)
#regular data
#for(d in 1:10){
#	seed = d
#	dat = make_data(N, P, per_std,seed)
#	write.csv(dat, paste0("caus_sim_1_N=",N,"_P=",P,"_noise_percent=", per_std, "_set_", d,".csv"))
#}
#imbalanced data
for(d in 1:length(squash_val)){
  seed = d
  dat = make_data_w_imbalance(N, P, per_std, squash_val[d], seed)
  write.csv(dat, paste0("data/caus_sim_1_N=",N,"_P=",P,"_noise_percent=", per_std, "_squash_val_", squash_val[d],".csv"))
}
