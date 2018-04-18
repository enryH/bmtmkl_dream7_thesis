DZA:
setwd("C:/Users/Webel/Desktop/Versuche/statisticsmaster/cancerprototyp/")
# Laptop:
setwd("C:/Users/Henry/Documents/1_Statistik/statisticsmaster/cancerprototyp/")
# Linux:
setwd("/home/enryh/Desktop/statisticsmaster/cancerprototyp/")
source("code_R/bayesian_multitask_multiple_kernel_learning_train_uncommented.R")
source("code_R/bayesian_multitask_multiple_kernel_learning_test.R")
# load(".RData")
#set the number of tasks (e.g., the number of compounds in Nature Biotechnology paper)
T <- 30  ## in Costello et. al. only 28 drugs were mentioned
#set the number of kernels (e.g., the number of views in Nature Biotechnology paper)
P <- 22
## load data and compute:
source("code_R/data_master.R")
source("code_R/bmtmkl_setparfct.R")
parameters <- setgammadisthyperparameter(1,1,1,1,1,1,1,1,1,1,200,1606)
###################################################################################################
values_to_test = list( c(10^(-10),10^(-10)) , c(1, 1 ) ,c(10^(-10),10^(10)))
sample <- c(1,8,23)
cells  <- c(1:11)
P=8
#
for (iter in c(10,50,200)){
for (gam in values_to_test) {
  parameters <- setgammadisthyperparameter(gam[1],gam[2],gam[1],gam[2],gam[1],gam[2],gam[1],gam[2],gam[1],gam[2],iter)
  state <- bayesian_multitask_multiple_kernel_learning_train(Ktrain, ytrain, parameters)
  cat("Prior for Gammas: ", gam ," ; ELBO after",parameters$iteration, " :",state$bounds[parameters$iteration],
      "\n#############################################################\n")
  #
  prediction <- bayesian_multitask_multiple_kernel_learning_test(Ktest, state)
  outofsample <- NULL
  colnames<- vector("list", 2*length(sample))
  i = 1
  for (t in sample) {
    outofsample <- cbind(outofsample, drugResponsetest[cells,t],(prediction$y[[t]]$mu[cells] * drugscale[t] + druglevel[t]) )
    colnames[(2*i-1):(2*i)] <- c(sprintf("D%dt" ,t) , sprintf("D%dp" ,t))
    i <- i+1
    # dimnames(outofsample)[2] <- c(sprintf("drug_%d_p" ,t), sprintf("drug_%d_t" ,t))
  }
  mse_ofs = 0  
  for (t in 1:T){
    mse_ofs <- mse_ofs + ((1/sum(!is.na(drugResponsetest[,t])))*sum(((drugResponsetest[,t] - (prediction$y[[t]]$mu * drugscale[t] + druglevel[t]))^2) , na.rm=TRUE) )
  }
  dimnames(outofsample)[[2]] <- as.list(colnames)
  print(outofsample)
  cat("\n#############################################################\n",
      "Mean-Squared-Error: ", mse_ofs , "\n\n")
  #print(state$parameters)
}
}
###############################################################################
source("code_R/data_master_yunscaled.R")
#
for (iter in c(50,100,200)){
for (gam in values_to_test) {
  parameters <- setgammadisthyperparameter(gam[1],gam[2],gam[1],gam[2],gam[1],gam[2],gam[1],gam[2],gam[1],gam[2],iter)
  state <- bayesian_multitask_multiple_kernel_learning_train(Ktrain, ytrain, parameters)
  cat("Prior for Gammas: ", gam ," ; ELBO after", parameters$iteration, " :",state$bounds[parameters$iteration],
      "\n#############################################################\n")
  #
  prediction <- bayesian_multitask_multiple_kernel_learning_test(Ktest, state)
  outofsample <- NULL
  colnames<- vector("list", 2*length(sample))
  i = 1
  for (t in sample) {
    outofsample <- cbind(outofsample, drugResponsetest[cells,t] , (prediction$y[[t]]$mu[cells]))
    colnames[(2*i-1):(2*i)] <- c(sprintf("D%dt" ,t) , sprintf("D%dp" ,t))
    i <- i+1
    # dimnames(outofsample)[2] <- c(sprintf("drug_%d_p" ,t), sprintf("drug_%d_t" ,t))
  }
  mse_ofs = 0  
  for (t in 1:T){
    # squared error (sq err): all observations in test sample of 18 cell lines
    mse_ofs <- mse_ofs + ((1/sum(!is.na(drugResponsetest[,t])))*sum(((drugResponsetest[,t] - prediction$y[[t]]$mu)^2) , na.rm=TRUE))
  }
  dimnames(outofsample)[[2]] <- as.list(colnames)
  print(outofsample)
  cat("\n#############################################################\n",
      "Mean-Squared-Error: ", mse_ofs , "\n\n")
  #print(state$parameters)
}
}
