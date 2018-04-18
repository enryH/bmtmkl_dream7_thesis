#DZA:
setwd("C:/Users/Webel/Desktop/Versuche/statisticsmaster/cancerprototyp/")
# Laptop:
setwd("C:/Users/Henry/Documents/1_Statistik/statisticsmaster/cancerprototyp/")
# Linux:
setwd("/home/enryh/Desktop/statisticsmaster/cancerprototyp/")
source("code_R/bayesian_multitask_multiple_kernel_learning_train_uncommented.R")
source("code_R/bayesian_multitask_multiple_kernel_learning_test.R")
# load(".RData")
#set the number of tasks (e.g., the number of compounds in Nature Biotechnology paper)
T <- 28  ## in Costello et. al. only 28 drugs were mentioned
#set the number of kernels (e.g., the number of views in Nature Biotechnology paper)
P <- 8
## load data and compute:
source("Code_R/data_master.R")
###### run model function  ########################################################################
runmodel <- function(parameters){
  # start.time <- Sys.time()
  state       <- bmtmkl(Ktrain, ytrain, parameters)
  prediction  <- bmtmkl_getpred(Ktest , state)
  trainingfit <- bmtmkl_getpred(Ktrain, state)
  ### squared error (in-sample-fit)
  MSE_insample = 0
  for (t in 1:T){
    # cat(colnames(y_train)[t], "\t has ", sum(!is.na(y_train[,t])), "\t observed values. \n")
    MSE_insample <- ( MSE_insample + 
      1/ sum(!is.na(y_train[,t])) * 
      sum(((ytrain[[t]] - (trainingfit$y[[t]]$mu * drugscale[t] + druglevel[t]))^2)))
  }  
### squared error (sq err): all observations in test sample of 18 cell lines
  MSE_outofsample = 0
  for (t in 1:T){
    # cat(colnames(drugResponsetest)[t], "\t has ", sum(!is.na(drugResponsetest[,t])), "\t observed values. \n")
    MSE_outofsample <- MSE_outofsample + 
      1/ sum(!is.na(drugResponsetest[,t])) * 
      sum(((drugResponsetest[,t] - (prediction$y[[t]]$mu * drugscale[t] + druglevel[t]))^2) , na.rm=TRUE)
  }
state$parameters$msetestset  <- MSE_outofsample  ## add sq_err to parametes
state$parameters$msetrainset <- MSE_insample
state$parameters$ELBOFinal   <- state$bounds[length(state$bounds)]
return(state$parameters)  # gives back a list
}
###### setup ######################################################################################
#fold_nr <- 1 #set to 1,...,5
core_nr <- 4 #set to nr of cores of parallelization

# hyperparameter grid
alpha <- c( 10^(-10),1, 10)
beta  <- c(10^(-10), 0.01, 1)
# https://stat.ethz.ch/R-manual/R-devel/library/base/html/expand.grid.html
grid <- expand.grid(alpha_lambda  = alpha, beta_lambda  = beta,
                    alpha_upsilon = alpha, beta_upsilon = beta, 
                    alpha_gamma   = alpha, beta_gamma   = beta,
                    alpha_omega   = alpha, beta_omega   = beta, 
                    alpha_epsilon = alpha, beta_epsilon = beta  
        )
source("code_R/bmtmkl_setparfct.R") # load function in script
i=10
parameters <- setgammadisthyperparameter(alpha_lambda  = grid[i,1], beta_lambda  = grid[i,2],
                           alpha_upsilon = grid[i,3], beta_upsilon = grid[i,4], 
                           alpha_gamma   = grid[i,5], beta_gamma   = grid[i,6],
                           alpha_omega   = grid[i,7], beta_omega   = grid[i,8], 
                           alpha_epsilon = grid[i,9], beta_epsilon = grid[i,10] ,
                           iterations= 200, seed=1606)
runmodel(parameters)
#########################################################################################
if(!require("doParallel")) install.packages("doParallel") 
library("doParallel")
#load("estimationset.RData")
#load("folds.RData")

#testset <-estimationset[folds[[fold_nr]],]
#trainset <- estimationset[-folds[[fold_nr]],]

# create downsampled balaned trainset indices
#bagging_index <- bagging.index(trainset, 5)


core_nr=2
cat("approx. run-time in h: ", dim(grid)[1]/core_nr/(3600/11))

start.time <- Sys.time()
  # loop over bagging sets
  # for (n in 1:5) {
  #   print(n)
  
  #trainset_bag <- trainset[bagging_index[[n]],]
    
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
    output <- runmodel(temp_par)        
  }
  stopCluster(cl)
  # save result object to file
  file_name <- paste("results", "_",format(Sys.time(), "%g%m%d_%H-%M"),".RData", sep = "")
  rownames(result)<-NULL  #Strip out the rownames
  result <- t(result)
  save(result, file = file_name)
  print(file_name)
# }
end.time  <- Sys.time()
print(end.time-start.time)
###################################################################################################
#list to dataframe:
names <- c("lambda_a","lambda_b", "upsilon_a",  "upsilon_b", "gamma_a",    "gamma_b",
  "omega_a" ,"omega_b" , "epsilon_a",  "epsilon_b",            
  "iter", "process", "seed",
  "MSE_ins", "MSE_ofs", "MSE_ins_resc", "ELBO")
dimension_temp <- attr(result, "dim")
mat = matrix(unlist(result), nrow=dimension_temp[1], ncol=dimension_temp[2])
df <- as.data.frame(mat)
save(df, file = paste("results", "_",format(Sys.time(), "%g%m%d_%H-%M"),"_df",".RData", sep = ""))
###################################################################################################
load()
###################################################################################################