#DZA:
setwd("C:/Users/Webel/Desktop/Versuche/statisticsmaster/cancerprototyp/")
# Laptop:
setwd("C:/Users/Henry/Documents/1_Statistik/statisticsmaster/cancerprototyp/")
# Fernrechner:
setwd("C:/Users/Webel/Documents/cancerprototyp/")
# Linux:
setwd("/home/enryh/Desktop/statisticsmaster/cancerprototyp/")
source("code_R/bayesian_multitask_multiple_kernel_learning_train_uncommented.R")
file.edit(("code_R/bayesian_multitask_multiple_kernel_learning_train_uncommented.R"))
source("code_R/bayesian_multitask_multiple_kernel_learning_test.R")
# load(".RData")
#set the number of tasks (e.g., the number of compounds in Nature Biotechnology paper)
T <- 30  ## in Costello et. al. only 28 drugs were mentioned
#set the number of kernels (e.g., the number of views in Nature Biotechnology paper)
P <- 22
set.seed(1606)
## load data and compute:
source("code_R/data_master_yunscaled.R")
file.edit("code_R/data_master_yunscaled.R")
source("code_R/bmtmkl_setparfct.R")
###################################################################################################
# reference value only using the mean:
MSE_insample_mean = 0
for (t in 1:T){
  # cat(colnames(y_train)[t], "\t has ", sum(!is.na(y_train[,t])), "\t observed values. \n")
  MSE_insample_mean <- ( MSE_insample_mean + 
                               1/ sum(!is.na(ic50_train[[t]])) * 
                               sum(((ic50_train[[t]] - mean(ic50_train[[t]]) )^2)))
  # cat("MSE insample based on mean prediction after drug ",t, ": ", MSE_insample_mean, "\n")
}
print(MSE_insample_mean)  # 9.105074 
# compute TSS (Total Sum of Squares -> prop. to variance)
for (t in 1:T){
  # cat(colnames(y_train)[t], "\t has ", sum(!is.na(y_train[,t])), "\t observed values. \n")
  TSS <- ( TSS + sum(((ic50_train[[t]] - mean(ic50_train[[t]]) )^2)))
  cat("TTS after drug ",t, ": ", TTS, "\n")
}
print(TSS)  # 22783.83  # to compare with RSS (Residual Sum of Squares)
TSS_scaled <- sum(unlist(ytrain)^2) # equals T
# length(unlist(ytrain)^2) -28 = 757
#
MSE_outofsample_mean = 0
for (t in c(1:T)){
  # cat(colnames(y_train)[t], "\t has ", sum(!is.na(y_train[,t])), "\t observed values. \n")
  MSE_outofsample_mean <- ( MSE_outofsample_mean + 
                           1/ sum(!is.na(drugResponsetest[,t])) * 
                           sum(((drugResponsetest[,t] - mean(ic50_train[[t]]) )^2), na.rm=TRUE))
  cat("MSE out-of-sample based on mean prediction after drug ",t, ": ", MSE_outofsample_mean, "\n")
}
# contribution of outlier (maximal value) of drug 24:
cat("contribution of outlier on drug 24:",  (1/ sum(!is.na(drugResponsetest[,10]))) * 16.776512979^2  )
# roughly 60% contribution to MSE out of sample by one single observation.
# excluding drug24: 11.56214
print(MSE_outofsample_mean)  # 29.16276
###### run model function  ########################################################################
# utils_grid_search.R
getMSE_unscaled <- function(parameters){
  # start.time <- Sys.time()
  state      <- bayesian_multitask_multiple_kernel_learning_train(Ktrain, ytrain, parameters)
  prediction <- bayesian_multitask_multiple_kernel_learning_test(Ktest, state)
  trainingfit <- bayesian_multitask_multiple_kernel_learning_test(Ktrain, state)
  ### squared error (in-sample-fit)
  MSE_insample = 0
  for (t in 1:T){
    # cat(colnames(y_train)[t], "\t has ", sum(!is.na(y_train[,t])), "\t observed values. \n")
    MSE_insample <- ( MSE_insample + 
                        1/ sum(!is.na(y_train[,t])) * 
                        sum(((ytrain[[t]] - trainingfit$y[[t]]$mu)^2)))
  } 
  ### squared error (sq err): all observations in test sample of 18 cell lines
  MSE_outofsample = 0
  for (t in 1:T){
    # cat(colnames(drugResponsetest)[t], "\t has ", sum(!is.na(drugResponsetest[,t])), "\t observed values. \n")
    MSE_outofsample <- MSE_outofsample + 
      1/ sum(!is.na(drugResponsetest[,t])) * 
      sum(((drugResponsetest[,t] - prediction$y[[t]]$mu)^2) , na.rm=TRUE)
  }
  # add sq_err to parametes
state$parameters$msetrainset          <- MSE_insample  
state$parameters$msetestset           <- MSE_outofsample    
# add final ELBO:
state$parameters$ELBOFinal   <- state$bounds[length(state$bounds)]
return(state$parameters)  # gives back a list
}
###### setup ######################################################################################
# hyperparameter grid
alpha <- c(10^(-10),    1, 10)
beta  <- c(10^(-10), 0.01, 1)
# 3^10 combinations
# https://stat.ethz.ch/R-manual/R-devel/library/base/html/expand.grid.html
grid <- expand.grid(alpha_lambda  = alpha, beta_lambda  = beta,
                    alpha_upsilon = alpha, beta_upsilon = beta, 
                    alpha_gamma   = alpha, beta_gamma   = beta,
                    alpha_omega   = alpha, beta_omega   = beta, 
                    alpha_epsilon = alpha, beta_epsilon = beta  
        )
source("code_R/bmtmkl_setparfct.R") # load function in script
i=51668
setparameters_from_grid <- function(i){
parameters <- setgammadisthyperparameter(
                           alpha_lambda  = grid[i,1], beta_lambda  = grid[i,2],
                           alpha_upsilon = grid[i,3], beta_upsilon = grid[i,4], 
                           alpha_gamma   = grid[i,5], beta_gamma   = grid[i,6],
                           alpha_omega   = grid[i,7], beta_omega   = grid[i,8], 
                           alpha_epsilon = grid[i,9], beta_epsilon = grid[i,10] ,
                           iterations= 200, seed=1606)
        return(parameters)
}
# for (i in 1:dim(grid)[1]) {if (sum(grid[i,1:10]== c(1,1,1,1,1,1,1,1,1,1))== 10) {cat("row of ones: ", i)} } 
#row of ones:  51668
parameters <- setparameters_from_grid(i)
result_i <- (getMSE_unscaled(parameters))
names(result_i)
####################################################################################
if(!require("doParallel")) install.packages("doParallel")
library("doParallel")
# 59049 / 7.889/ 60 /12 = 10.4
core_nr=12
cat("approx. run-time in h: ", dim(grid)[1]/core_nr/60/10.4)

start.time <- Sys.time()

  # parallelization
  cl <- makeCluster(core_nr)
  registerDoParallel(cl)

  result <- foreach (i = (1:dim(grid)[1]), .combine = cbind) %dopar%{
    # train and predict kernel svm
    temp_par <- setgammadisthyperparameter(alpha_lambda  = grid[i,1], beta_lambda  = grid[i,2],
                                           alpha_upsilon = grid[i,3], beta_upsilon = grid[i,4],
                                           alpha_gamma   = grid[i,5], beta_gamma   = grid[i,6],
                                           alpha_omega   = grid[i,7], beta_omega   = grid[i,8],
                                           alpha_epsilon = grid[i,9], beta_epsilon = grid[i,10] ,
                                           iterations    = 200      , seed         = 1606)
    output <- getMSE_unscaled(temp_par)
  }
  stopCluster(cl)
  # save result object to file
  file_name <- paste("results", "_",format(Sys.time(), "%g%m%d_%H-%M"),".RData", sep = "")
  rownames(result)<-NULL  #Strip out the rownames
  result <- t(result)
  save.image(file = file_name) # save whole working space in order to be able to redo everything
  print(file_name)
end.time  <- Sys.time()
print(end.time-start.time)
# ###################################################################################################
# list to dataframe:
names <- c("lambda_a","lambda_b", "upsilon_a",  "upsilon_b", "gamma_a",    "gamma_b",
           "omega_a" ,"omega_b" , "epsilon_a",  "epsilon_b",            
           "iter", "process", "seed",
           "MSE_ins", "MSE_ofs", "ELBO")
dimension_temp <- attr(result, "dim")
mat = matrix(unlist(result), nrow=dimension_temp[1], ncol=dimension_temp[2])
results_grid_ns_drugs_200iter <- as.data.frame(mat)
colnames(results_grid_ns_drugs_200iter) <- names
save(results_grid_ns_drugs_200iter, parameters, grid, file = paste("results", "_",format(Sys.time(), "%g%m%d_%H-%M"),"_df",".RData", sep = ""))
save(results_grid_ns_drugs_200iter, parameters, grid, file = paste("results_ns_drugs_200iter_df",".RData", sep = ""))
# ###################################################################################################
