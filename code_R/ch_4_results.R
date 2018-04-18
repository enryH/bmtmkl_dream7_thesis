source("code_R/bmtmkl_train_uncommented.R")
# file.edit("code_R/bmtmkl_train_uncommented.R")
source("code_R/bmtmkl_getpred.R")
# load(".RData")
## load data and compute:
source("code_R/data_master.R")         # load all necessary data and calculate kernels
# file.edit("code_R/data_master.R")    # inspect and understand!
source("code_R/bmtmkl_utils.R")
# file.edit("code_R/bmtmkl_setparfct.R")
if(!require("xtable")) install.packages("xtable")  # load package for latex tables
###### setup ######################################################################################
# ### build true IC50 values matrix including NAs.
drugReponseALLTRUE <-rbind (traindata$DrugResponse[1:35,] , 
                            as.matrix(read.table(filepath_testdataoutcome, header=T, row.names=1)) 
)
# check ordering and order as in traindata:
cellline_names_all == rownames(drugReponseALLTRUE)
drugReponseALLTRUE <- preprocessData(drugReponseALLTRUE, cellline_names_all)
cellline_names_all == rownames(drugReponseALLTRUE)
rownames(traindata$DrugResponse)== rownames(drugReponseALLTRUE)
# View(drugReponseALLTRUE)
#
drugResponsetest=preprocessData(drugResponsetest, cellline_names_test)
# View(drugResponsetest)
#
# summary-statistics using train+ test data:
source("code_R/ch_appendix_summary_drugs.R")
# file.edit('code_R/app_summary_drugs.R')
#
#drugnames used in Costello et.al.:
drugnames <- (as.matrix(read.csv("submission data/drugnames.csv")))
drugnames <- drugnames[,order(colnames(drugnames))]
drugnames  # drug 12, 26, 27 have been excluded
xtable((as.matrix(drugnames)))
###################################################################################################
# reference value only using the mean:

MSE_insample_mean = 0
for (t in 1:T){
  # cat(colnames(y_train)[t], "\t has ", sum(!is.na(y_train[,t])), "\t observed values. \n")
  MSE_insample_mean <- ( MSE_insample_mean + 
                           1/ sum(!is.na(ic50_train[[t]])) * 
                           sum(((ic50_train[[t]] - mean(ic50_train[[t]]) )^2)))
  flog.debug(cat("MSE insample based on mean prediction after drug ",t, ": ", MSE_insample_mean, "\n"))
}
flog.info(cat("MSE insample:", MSE_insample_mean))  # 9.119296 
# compute TSS (Total Sum of Squares -> prop. to variance)
TSS <- 0
for (t in 1:T){
  # cat(colnames(y_train)[t], "\t has ", sum(!is.na(y_train[,t])), "\t observed values. \n")
  TSS <- ( TSS + sum(((ic50_train[[t]] - mean(ic50_train[[t]]) )^2)))
  flog.debug(cat("TTS after drug ",t, ": ", TSS, "\n"))
}
print(TSS)  # 238.5528  # to compare with RSS (Residual Sum of Squares)
TSS_scaled <- sum(unlist(ytrain)^2) # equals sum(N_t) - T 
flog.info(cat("scaled drug -log IC_50 data has as a total variance: " ,TSS_scaled, "which is N - T\n"))
# length(unlist(ytrain)^2) -30 = 822
#
MSE_insample_mean_scaled = 0
for (t in 1:T) { 
  print(1/ length(ytrain[[t]]) * sum((ytrain[[t]])^2))
  MSE_insample_mean_scaled <- MSE_insample_mean_scaled + ((1 /  length(ytrain[[t]])) * sum((ytrain[[t]])^2))
}
flog.info(cat("MSE for scaled data has a reference value when using the mean as predictor: ", MSE_insample_mean_scaled ,"\n")) # 28.86858
#
MSE_outofsample_mean = 0
for (t in 1:T){
  # cat(colnames(y_train)[t], "\t has ", sum(!is.na(y_train[,t])), "\t observed values. \n")
  MSE_outofsample_mean <- (MSE_outofsample_mean + 
                              1/ sum(!is.na(drugResponsetest[,t])) * 
                              sum(((drugResponsetest[,t] - mean(ic50_train[[t]]))^2), na.rm=TRUE))
  cat("MSE out-of-sample based on mean prediction after drug ",t, ": ", MSE_outofsample_mean, "\n")
}
flog.info(cat("MSE out-of-sample", MSE_outofsample_mean))  # 29.16276
flog.info("60% of out-of-sample MSE is due to one outlier of drug 10")
# contribution of outlier (maximal value) of drug 24:
flog.info(cat("contribution of outlier on drug 24:",  
              (1/ sum(!is.na(drugResponsetest[,10]))) * 16.776512979^2 , "to MSE ofs\n" ))
flog.info("roughly 60% contribution to MSE out of sample by one single observation.")
flog.info("excluding drug24: 11.56214")
###################################################################################################
###### run model function  ########################################################################
getMSE <- function(parameters){
  # start.time <- Sys.time()
  state      <- bmtmkl_train(Ktrain, ytrain, parameters)
  prediction <- bmtmkl_test(Ktest, state)
  trainingfit<- bmtmkl_test(Ktrain, state)
  ### squared error (in-sample-fit)
  MSE_insample_rescaled = 0
  for (t in 1:T){
    # cat(colnames(y_train)[t], "\t has ", sum(!is.na(y_train[,t])), "\t observed values. \n")
    MSE_insample_rescaled <- ( MSE_insample_rescaled + 
                                 1/ sum(!is.na(y_train[,t])) * 
                               sum((((ytrain[[t]]* drugscale[t] + druglevel[t]) - 
                                       (trainingfit$y[[t]]$mu * drugscale[t] + druglevel[t]))^2))
                              )
  } 
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
    flog.debug(cat(colnames(drugResponsetest)[t], "\t has ", sum(!is.na(drugResponsetest[,t])),
                   "\t observed values. \n"))
    MSE_outofsample <- MSE_outofsample + 
      1/ sum(!is.na(drugResponsetest[,t])) * 
      sum(((drugResponsetest[,t] - (prediction$y[[t]]$mu * drugscale[t] + druglevel[t]))^2) , na.rm=TRUE)
  }
  MSE_outofsample_scaled = 0
  for (t in 1:T){
    flog.debug(cat(colnames(drugResponsetest)[t], "\t has ", sum(!is.na(drugResponsetest[,t])), 
    "\t observed values. \n"))
    MSE_outofsample_scaled <- MSE_outofsample_scaled + 
      1/ sum(!is.na(drugResponsetest[,t])) * 
      sum(((scale(drugResponsetest[,t]) - (prediction$y[[t]]$mu))^2) , na.rm=TRUE)
  } 
  # add sq_err to parametes
  state$parameters$msetrainset          <- MSE_insample  
  state$parameters$msetestset           <- MSE_outofsample    
  state$parameters$msetrainset_rescaled <- MSE_insample_rescaled
  state$parameters$msetestset_scaled    <- MSE_outofsample_scaled 
  # add final ELBO:
  state$parameters$ELBOFinal   <- state$bounds[length(state$bounds)]
  return(state$parameters)  # gives back a list
}
###################################################################################################
# hyperparameter grid used in grid:
# 
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
source("code_R/bmtmkl_utils.R") # load function in script
setparameters_from_grid <- function(i, iterations=200 , seed= 1606){
  parameters <- setgammadisthyperparameter(
    alpha_lambda  = grid[i,1], beta_lambda  = grid[i,2],
    alpha_upsilon = grid[i,3], beta_upsilon = grid[i,4], 
    alpha_gamma   = grid[i,5], beta_gamma   = grid[i,6],
    alpha_omega   = grid[i,7], beta_omega   = grid[i,8], 
    alpha_epsilon = grid[i,9], beta_epsilon = grid[i,10] ,
    iterations= iterations, seed=seed)
  return(parameters)
}
# for (i in 1:dim(grid)[1]) {if (sum(grid[i,1:10]== c(1,1,1,1,1,1,1,1,1,1))== 10) {cat("row of ones: ", i)} } 
i=51668 #row of ones:  51668
parameters <- setparameters_from_grid(i)
result_i <- (getMSE(parameters))
result_i2 <- getMSE(setparameters_from_grid(i, 40))
cbind(result_i2,result_i)
###################################################################################################
# all views
load("data/results_180113_06-50_df.RData") # default using 6 views?
results_default <- df
colnames(results_default) <- c("lambda_a","lambda_b", "upsilon_a",  "upsilon_b", "gamma_a",    "gamma_b",
                               "omega_a" ,"omega_b" , "epsilon_a",  "epsilon_b",            
                               "iter", "process", "seed",
                               "MSE_ins", "MSE_ofs", "MSE_ins_resc", "ELBO")
unique(results_default[,1])
unique(results_default[,2])
View(results_default)

load("data/results_ns_drugs_200iter_df.RData")
View(results_grid_ns_drugs_200iter)
#ELBO: maxELBO: i=  58770 ; minELBO: i=73
grid[58770,]
results_hyper_ns <- matrix(0,2,11)
results_hyper_ns[1,] <- as.matrix(results_grid_ns_drugs_200iter[58770, c(1:10,16)])
results_hyper_ns[2,] <- as.matrix(results_grid_ns_drugs_200iter[73,c(1:10,16)])
colnames(results_hyper_ns) <- colnames(results_grid_ns_drugs_200iter)[c(1:10,16)]
rownames(results_hyper_ns) <- c("maxELBO", "minELBO")
xtable(results_hyper_ns, digits=2)
###################################################################################################
###################################################################################################
# parameters, criteria considered
# predictions y (scaled and rescaled)
# mean values of hidden random variables for three models
# compare predictions with less views -> different settings of views.
# start with 6 given views and then add most predictive view?
###################################################################################################
select <- c(45278,45359,45440)  
table_min_MSE_ofs <- matrix(0, length(result_i), length(select))
i=1
for (k in select) {
  cat("grid entry: ", k , "\n")
  table_min_MSE_ofs[,i]<-  unlist(getMSE(setparameters_from_grid(k)))
  i = i + 1
}
colnames(table_min_MSE_ofs) <- select
rownames(table_min_MSE_ofs) <- names(result_i)
rm("i", "k", "select")
#
table_min_MSE_ofs
currenttable <- xtable(table_min_MSE_ofs)
print(currenttable, type="latex", file=("paper/Tables/tab_min_test_error_comp"))
#
select <- c(59041,59042,59043) #in-sample
table_min_mse_in <- matrix(0, length(result_i), length(select))
i=1
for (k in select) {
  cat("grid entry: ", k , "\n")
  table_min_mse_in[,i]<-  unlist(getMSE(setparameters_from_grid(k)))
  i = i + 1
}
colnames(table_min_mse_in) <- select
rownames(table_min_mse_in) <- names(result_i)
rm("i", "k", "select")
#
table_min_mse_in[-c(12,13,17),]
currenttable <- xtable(table_min_mse_in[-c(12,13,17),])
print(currenttable, type="latex", file=("paper/Tables/tab_min_test_error_comp"))

xtable(cbind(table_min_mse_in[-c(12,13,17),],table_min_MSE_ofs[-c(12,13,17),] ))
####################################################################################################
result_random  <- unlist(getMSE(setparameters_from_grid(51688)))
result_allones  <- unlist(getMSE(setparameters_from_grid(51668)))
result_default  <- unlist(getMSE(setparameters_from_grid(1)))
result_maxELBO <- unlist(getMSE(setparameters_from_grid(52452)))
result_minELBO <- unlist(getMSE(setparameters_from_grid(73)))
# result_minMSE_1   <- unlist(getMSE(setparameters_from_grid(59041)))
# result_minMSE_2   <- unlist(getMSE(setparameters_from_grid(59042)))
# result_minMSE_3   <- unlist(getMSE(setparameters_from_grid(59042)))
# resut_minMSE      <- unlist(getMSE(setparameters_from_grid(58983)))

# table of results:
results_hyper      <- matrix(0,5,11)
colnames(results_hyper) <- c("alpha_lambda","beta_lambda", "alpha_upsilon" ,  "beta_upsilon", "alpha_gamma",    "beta_gamma",
                             "alpha_omega" ,"beta_omega" , "alpha_epsilon" ,  "beta_epsilon", "ELBO")
# "MSE out of sample","MSE in sample" ,"MSE in sample standardized",
rownames(results_hyper) <- c("default","all ones", "max", "min", "random")
#
results_hyper[1,]  <- result_default[c(1:10, 18)]
results_hyper[2,]  <- result_allones[c(1:10, 18)]
results_hyper[3,]  <- result_maxELBO[c(1:10, 18)]
results_hyper[4,]  <- result_minELBO[c(1:10, 18)]
results_hyper[5,]  <- result_random[c(1:10,18)]
results_hyper
xtable(results_hyper, digits=c(1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0), caption="results")
##########################################################################
#14,15,16
compare_defaults <- matrix(0,4,6)
rownames(compare_defaults) <- names(result_default)[c(14,15,16,18)]
colnames(compare_defaults) <- c("mean", "default","all ones", "max", "min", "random")
compare_defaults[,1]<-c(MSE_insample_mean_scaled, MSE_outofsample_mean, MSE_insample_mean, NA)
compare_defaults[,2]<-result_default[c(14:16,18)]
compare_defaults[,3]<-result_allones[c(14:16,18)]
compare_defaults[,4]<-result_maxELBO[c(14:16,18)]
compare_defaults[,5]<-result_minELBO[c(14:16,18)]
compare_defaults[,6]<-result_random[c(14:16,18)]
compare_defaults <- compare_defaults[c(1,3,2,4),]
compare_defaults
xtable(compare_defaults, digits=2, caption="")
##########################################################################

####
# compare predicitons for drug 10
allones      <- bmtmkl_train(Ktrain, ytrain, setparameters_from_grid(51668))
default      <- bmtmkl_train(Ktrain, ytrain, setparameters_from_grid(1))
maxELBO      <- bmtmkl_train(Ktrain, ytrain, setparameters_from_grid(52452))
minELBO      <- bmtmkl_train(Ktrain, ytrain, setparameters_from_grid(73))
random       <- bmtmkl_train(Ktrain, ytrain, setparameters_from_grid(51688))
#
pred_allones <- bmtmkl_test(Ktest, allones)
pred_default <- bmtmkl_test(Ktest, default)
pred_maxELBO <- bmtmkl_test(Ktest, maxELBO)
pred_minELBO <- bmtmkl_test(Ktest, minELBO)
pred_random  <- bmtmkl_test(Ktest, random)
####
fit_allones <- bmtmkl_test(Ktrain, allones)
fit_default <- bmtmkl_test(Ktrain, default)
fit_maxELBO <- bmtmkl_test(Ktrain, maxELBO)
fit_minELBO <- bmtmkl_test(Ktrain, minELBO)
fit_random  <- bmtmkl_test(Ktrain, random)
####################################################################################
# ranking # rank largerst value with 1!
myrank <- function(x) {return(rank(-x, ties.method = "first"))}
####################################################################################
set.seed(1606)  # in order to get reproducable results
fit_sample <- sort(sample(length(ytrain[[10]]), 5) )
#
fit_y_10 <- cbind(ytrain[[10]], fit_default$y[[10]]$mu, fit_allones$y[[10]]$mu, fit_maxELBO$y[[10]]$mu, fit_minELBO$y[[10]]$mu, fit_random$y[[10]]$mu)
colnames(fit_y_10) <- c("scaled IC", "default", "ones", "max", "min", "random")
rownames(fit_y_10) <- names((y_train[,10]))[!is.na(y_train[,10])]
fit_y_10[fit_sample,]
apply(fit_y_10[fit_sample,], 2, myrank)
fit_y_10_tab1 <- cbind(fit_y_10[fit_sample,], 
      apply(fit_y_10[fit_sample,], 2, myrank)   )
fit_y_10_tab1
xtable(fit_y_10_tab1, digits=c(1,2,-2,2,-2,-2,2,0,0,0,0,0,0))
# rescale:
fit_y_10_rescaled <- fit_y_10*drugscale[10]+druglevel[10]
fit_y_10_rescaled[fit_sample,]
cbind(fit_y_10_rescaled[fit_sample,],  apply(fit_y_10_rescaled[fit_sample,], 2, myrank)) 
xtable(
  cbind(fit_y_10_rescaled[fit_sample,],  apply(fit_y_10_rescaled[fit_sample,], 2, myrank)) ,
  digits=c(1,3,3,3,3,3,3,0,0,0,0,0,0)
)
### rescale fitted values  #####
rescale_fit <- function(fit, drugno) { fit$y[[drugno]]$mu * drugscale[drugno] + druglevel[drugno]}
rescale_variance <- function(fit, drugno) {fit$y[[drugno]]$sigma * (drugscale[drugno])^2}
# compare fit to variances:
fits_and_var <- cbind(fit_allones$y[[10]]$mu , fit_allones$y[[10]]$sigma, sqrt(fit_allones$y[[10]]$sigma),
                      rescale_fit(fit_allones, 10), rescale_variance(fit_allones, 10) , sqrt(rescale_variance(fit_allones, 10)))
colnames(fits_and_var) <- c("fit", "Var","SD", "fit res.", "var res.", "SD")
rownames(fits_and_var) <- names(ytrain[[10]])
fits_and_var[fit_sample,]
xtable(fits_and_var[fit_sample,], digits=3)
####
# build out-of-sample table:
set.seed(1607) # in order to get reproducable results
pred_sample <- sample(1:18, 5)
#scale?
# results_y_ofs <- cbind(scale(drugResponsetest[,10]), pred_default$y[[10]]$mu, pred_allones$y[[10]]$mu, pred_maxELBO$y[[10]]$mu, pred_minELBO$y[[10]]$mu)
results_y_ofs <- cbind((drugResponsetest[,10]-druglevel[10])/drugscale[10], pred_default$y[[10]]$mu, pred_allones$y[[10]]$mu, pred_maxELBO$y[[10]]$mu, pred_minELBO$y[[10]]$mu, pred_random$y[[10]]$mu)
results_y_ofs_rescaled <- results_y_ofs *drugscale[10]+druglevel[10]
#results_y_ofs_rescaled <- rbind(c(NA,results_hyper[,11]),   results_y_ofs_rescaled) # add ELBO
#results_y_ofs <- rbind(c(NA,results_hyper[,11]),   results_y_ofs)          # add ELBO
colnames(results_y_ofs)   <- c("true IC", "default", "ones", "max", "min", "random")
colnames(results_y_ofs_rescaled) <- c("true IC", "default", "ones", "max", "min", "random")
rownames(results_y_ofs)          <- c(rownames(drugResponsetest))
rownames(results_y_ofs_rescaled) <- c(rownames(drugResponsetest))
results_y_ofs
results_y_ofs_rescaled
results_y_ofs_sample <- results_y_ofs[pred_sample,]
results_y_ofs_rescaled_sample <- results_y_ofs_rescaled[pred_sample,]
### latex tables:
xtable(cbind(results_y_ofs_sample, apply(results_y_ofs_sample, 2, myrank))
             , digits=c(1,3,-2,2,-2,-2,2,0,0,0,0,0,0), 
       caption="prediction given by model trained on standardized data")
xtable(cbind(results_y_ofs_rescaled_sample, apply(results_y_ofs_rescaled_sample, 2, myrank))
             , digits=c(1,3,3,3,3,3,3,0,0,0,0,0,0), 
       caption="rescaled predictions of test data given by model trained on standardized data")
##########################################################################
# same for non-standardized out-comes... missing
###### run model function  ##############
getMSE_unscaled <- function(parameters){
  # start.time <- Sys.time()
  state      <- bmtmkl_train(Ktrain, ic50_train, parameters)
  prediction <- bmtmkl_test(Ktest, state)
  trainingfit<- bmtmkl_test(Ktrain, state)
  ### squared error (in-sample-fit)
  MSE_insample = 0
  for (t in 1:T){
    # cat(colnames(y_train)[t], "\t has ", sum(!is.na(y_train[,t])), "\t observed values. \n")
    MSE_insample <- ( MSE_insample + 
                        1/ sum(!is.na(ic50_train[[t]])) * 
                        sum(((ic50_train[[t]] - trainingfit$y[[t]]$mu)^2)))
  } 
  ### squared error (sq err): all observations in test sample of 18 cell lines
  MSE_outofsample = 0
  for (t in 1:T){
    # cat(colnames(drugResponsetest)[t], "\t has ", sum(!is.na(drugResponsetest[,t])), "\t observed values. \n")
    MSE_outofsample <- MSE_outofsample + 
      1/ sum(!is.na(drugResponsetest[,t])) * 
      sum(((drugResponsetest[,t] - prediction$y[[t]]$mu )^2) , na.rm=TRUE)
  }
  # add sq_err to parametes
  state$parameters$msetrainset          <- MSE_insample  
  state$parameters$msetestset           <- MSE_outofsample    
  # add final ELBO:
  state$parameters$ELBOFinal   <- state$bounds[length(state$bounds)]
  return(state$parameters)  # gives back a list
}
##########################################################################
# comparison table:
result_allones_ns  <- unlist(getMSE_unscaled(setparameters_from_grid(51668)))
result_default_ns  <- unlist(getMSE_unscaled(setparameters_from_grid(1)))
result_maxELBO_ns  <- unlist(getMSE_unscaled(setparameters_from_grid(52452)))
result_maxELBOstar_ns <- unlist(getMSE_unscaled(setparameters_from_grid(58770)))
result_minELBO_ns  <- unlist(getMSE_unscaled(setparameters_from_grid(73)))
result_random_ns  <- unlist(getMSE_unscaled(setparameters_from_grid(51688)))
# result_minMSE_ns   <- unlist(getMSE_unscaled(setparameters_from_grid(58393)))
# resut_minMSE_ns   <- unlist(getMSE_unscaled(setparameters_from_grid(58983)))
compare_defaults <- matrix(0,3,6)
rownames(compare_defaults) <- names(result_default_ns)[c(14,15,16)]
colnames(compare_defaults) <- c("mean", "default","all ones", "max*", "min","random")
compare_defaults[,1]<-c(MSE_insample_mean, MSE_outofsample_mean,  NA)
compare_defaults[,2]<-result_default_ns[c(14:16)]
compare_defaults[,3]<-result_allones_ns[c(14:16)]
compare_defaults[,4]<-result_maxELBOstar_ns[c(14:16)]
compare_defaults[,5]<-result_minELBO_ns[c(14:16)]
compare_defaults[,6]<-result_random_ns[c(14:16)]
compare_defaults
xtable(compare_defaults, digits=2, caption="Configurations of Table 9 trained with \\emph{non-standardized} drug data")
#
#
# compare predicitons for drug 10 (non-scaled)
default_ns      <- bmtmkl_train(Ktrain, ic50_train, setparameters_from_grid(1))
allones_ns      <- bmtmkl_train(Ktrain, ic50_train, setparameters_from_grid(51668))
maxELBO_ns      <- bmtmkl_train(Ktrain, ic50_train, setparameters_from_grid(52452))
minELBO_ns      <- bmtmkl_train(Ktrain, ic50_train, setparameters_from_grid(73))
random_ns       <- bmtmkl_train(Ktrain, ic50_train, setparameters_from_grid(51688))
# predictions:
pred_default_ns <- bmtmkl_test(Ktest, default_ns)
pred_allones_ns <- bmtmkl_test(Ktest, allones_ns)
pred_maxELBO_ns <- bmtmkl_test(Ktest, maxELBO_ns)
pred_minELBO_ns <- bmtmkl_test(Ktest, minELBO_ns)
pred_random_ns <- bmtmkl_test(Ktest, random_ns)

# fitted values training data:
# predictions:
fit_default_ns <- bmtmkl_test(Ktrain, default_ns)
fit_allones_ns <- bmtmkl_test(Ktrain, allones_ns)
fit_maxELBO_ns <- bmtmkl_test(Ktrain, maxELBO_ns)
fit_minELBO_ns <- bmtmkl_test(Ktrain, minELBO_ns)
fit_random_ns <- bmtmkl_test(Ktrain, random_ns)
###
####
set.seed(1606)  # in order to get reproducable results
fit_sample <- sort(sample(length(ytrain[[10]]), 5) )
#
fit_y_10_ns <- cbind((ic50_train[[10]]), fit_default_ns$y[[10]]$mu, fit_allones_ns$y[[10]]$mu, fit_maxELBO_ns$y[[10]]$mu, fit_minELBO_ns$y[[10]]$mu, fit_random_ns$y[[10]]$mu)
colnames(fit_y_10_ns) <- c(" IC", "default", "ones", "max", "min", "random")
# rownames(fit_y_10_ns) <- names((ic50_train))
fit_y_10_ns[fit_sample,]
apply(fit_y_10_ns[fit_sample,], 2, myrank)
fit_y_10_tab1 <- cbind(fit_y_10_ns[fit_sample,], 
                       apply(fit_y_10_ns[fit_sample,], 2, myrank))
#
xtable(fit_y_10_tab1, digits=c(1,2,2,2,2,2,2,0,0,0,0,0,0), 
       caption="fitted values for 5 training cell lines given by model trained on \\emph{non-standardized data} and corresponding ranking")
#####
# build table:
results_y_ofs_ns <- cbind((drugResponsetest[,10]), pred_default_ns$y[[10]]$mu, pred_allones_ns$y[[10]]$mu, pred_maxELBO_ns$y[[10]]$mu, pred_minELBO_ns$y[[10]]$mu, pred_random_ns$y[[10]]$mu)
# tmp_elbo_ns <- c(default_ns$bounds[200], allones_ns$bounds[200], maxELBO_ns$bounds[200], minELBO_ns$bounds[200], random_ns$bounds[200])
# results_y_ofs_ns <- rbind( c(NA, tmp_elbo_ns),   results_y_ofs_ns[pred_sample,])
# add ELBO and select sample
colnames(results_y_ofs_ns)   <- c("IC50", "default", "ones", "max", "min", "random")
# rownames(results_y_ofs_ns)          <- c("ELBO", rownames(drugResponsetest))
rownames(results_y_ofs_ns)          <- rownames(drugResponsetest)
results_y_ofs_ns[pred_sample,]
#### #### ####
xtable(cbind(results_y_ofs_ns[pred_sample,], apply(results_y_ofs_ns[pred_sample,],2,myrank)) 
       , digits=c(1,3,3,3,3,3,3,0,0,0,0,0,0), 
       caption="prediction given by model trained on standardized data")
##########################################################################
# combine results for out of samples predictions
# based on model trained with scaled and as is outcomes
# outcome data has to be made comparable. 
results_y_ofs_comb <- cbind(results_y_ofs_rescaled, results_y_ofs_ns[,c(-1)])
xtable(results_y_ofs_comb, digits=3,
       caption="Out of sample predictions based on models trained on scaled and initial IC50 outcome values")
##########################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
set.seed(1606)  # in order to get reproducable results
fit_sample <- sort(sample(length(ytrain[[10]]), 5) )
#######################################################################################################################
# compare mean values of hidden random variables lambda and a :
evaluate <- function(state,t){
  name_of_model <- (substitute(state))
  mat     <- matrix(t(state$lambda[[t]]$alpha * state$lambda[[t]]$beta))
  mat     <- cbind(mat, state$a[[t]]$mu)
  rownames(mat) <- names(ytrain[[t]])
  colnames(mat) <- c(paste0("lambda_", name_of_model), paste0("a_", name_of_model))
  return(mat)
}
drugno= 10
tab_lambda_a <- cbind(evaluate(default, drugno)
                      ,evaluate(allones, drugno)
                      ,evaluate(maxELBO, drugno)
                      ,evaluate(minELBO, drugno)
                      ,evaluate(random, drugno))  # why is E[lambda] constant for minELBO?
tab_lambda_a[fit_sample,]
xtable(tab_lambda_a[fit_sample,] , digits=c(1,-1,-1,2,3,3,-1,-1,3,-1,1))
###########################################################
#todo: how to loop over model names and change persistently something?
# # just do it in algorithm
# for (model in list(default, allones, maxELBO, minELBO)) {
#   print(substitute(model))
#   # colnames(model$G[[10]]$mu) <- colnames(Ktrain[[10]])
#   # rownames(model$G[[10]]$mu) <- c(paste0("mu_g_{10,",c(1:22),"}"))
#   
# }
# per Hand:
colnames(random$G[[10]]$mu) <- colnames(Ktrain[[10]])
rownames(random$G[[10]]$mu) <- c(paste0("mu_g_{10,",c(1:22),"}"))
colnames(default$G[[10]]$mu) <- colnames(Ktrain[[10]])
rownames(default$G[[10]]$mu) <- c(paste0("mu_g_{10,",c(1:22),"}"))
colnames(allones$G[[10]]$mu) <- colnames(Ktrain[[10]])
rownames(allones$G[[10]]$mu) <- c(paste0("mu_g_{10,",c(1:22),"}"))
colnames(maxELBO$G[[10]]$mu) <- colnames(Ktrain[[10]])
rownames(maxELBO$G[[10]]$mu) <- c(paste0("mu_g_{10,",c(1:22),"}"))
colnames(minELBO$G[[10]]$mu) <- colnames(Ktrain[[10]])
rownames(minELBO$G[[10]]$mu) <- c(paste0("mu_g_{10,",c(1:22),"}"))
#
expectation_upsilon <- function(model, drugno){print(model$upsilon$alpha[drugno]*model$upsilon$beta[drugno])}
#
xtable(t(default$G[[10]]$mu[1:6,])[fit_sample,] , digits=3, caption="Intermediate Outputs of default model")
expectation_upsilon(default, 10)
xtable(t(allones$G[[10]]$mu[1:6,])[fit_sample,] , digits=3, caption="Intermediate outputs of ones model")
expectation_upsilon(allones, 10)
xtable(t(maxELBO$G[[10]]$mu[1:6,])[fit_sample,] , digits=3, caption="Intermediate outputs of max ELBO model")
xtable(t(minELBO$G[[10]]$mu[1:6,])[fit_sample,] , digits=3, caption="Intermediate outputs of min ELBO model")
xtable(t(random$G[[10]]$mu[1:6,])[fit_sample,] , digits=1, caption="Intermediate Outputs of random model")
expectation_upsilon(random, 10)
t(default$G[[10]]$mu[1:6,])
library("pastecs")
xtable(stat.desc(t(default$G[[10]]$mu[1:6,]), basic=TRUE)[c(4:6,8,9, 13),],
       digits=4, caption="intermediate outputs summary statistics - default model" )  # default
xtable(rbind(default$be$mu[31:36], (default$omega$alpha * default$omega$beta)[1:6]),digits=-2)

xtable(stat.desc(t(allones$G[[10]]$mu[1:6,]), basic=TRUE)[c(4:6,8,9,13),],
       digits=3, caption="intermediate outputs summary statistics - ones model" )  # allones
xtable(rbind(allones$be$mu[31:36] , (allones$omega$alpha * allones$omega$beta)[1:6]),digits=4)

xtable(stat.desc(t(maxELBO$G[[10]]$mu[1:6,]), basic=TRUE)[c(4:6,8,9,13),],
       digits=3, caption="intermediate outputs summary statistics - maxELBO" )  # max
xtable(stat.desc(t(minELBO$G[[10]]$mu[1:6,]), basic=TRUE)[c(4:6,8,9,13),],
       digits=3, caption="intermediate outputs summary statistics - minELBO" ) # min
xtable(stat.desc(t(random$G[[10]]$mu[1:6,]), basic=TRUE)[c(4:6,8,9,13),],
       digits=1, caption="intermediate outputs summary statistics - random model" ) # random
xtable(rbind(random$be$mu[31:36], (random$omega$alpha * random$omega$beta)[1:6]),digits=-2)
###########################################################
table_coef <- function(state, i){
  table <- c(state$be$mu[i],     # mean of bias of drug 10
           state$be$mu[(T+1):(T+P)],    # mean of kernel coefficients for 22 views/input kernels
           state$gamma$alpha[i] * state$gamma$beta[i] , # precision of bias
           state$omega$alpha * state$omega$beta , # precisions for kernel coefficients
           state$epsilon$alpha[i] * state$epsilon$beta[i])  # precision for outcome
  #
  names(table)<- c(paste0("bias_drug_", i), paste0("e_", c(1:22)),"E_gamma", paste0("E_omega", c(1:22)), paste0("E_epsilon_", i))
  #
  return(table)
}
table           <- matrix(0,47,5)
colnames(table) <- c("default", "ones", "minELBO", "maxELBO", "random")
i=10
rownames(table) <- c(paste0("bias_drug_", i), paste0("e_", c(1:22)),"E_gamma", paste0("E_omega", c(1:22)), "E_epsilon")
table[,1] <- table_coef(default, 10)
table[,2] <- table_coef(allones, 10)
table[,3] <- table_coef(minELBO, 10)
table[,4] <- table_coef(maxELBO, 10)
table[,5] <- table_coef(random, 10) 
View(table)
xtable(t(table[c(1:7),]),digits=-2)
xtable(table[c(1:7,46,47),], digits=c(1,-2,4,-2,-2))
# why all kernel coeff. have the same precision? bug?
# biases?
table_biases <- function(state){
  table        <- c(state$be$mu[1:T],
                       ( state$gamma$alpha * state$gamma$beta)
                        )
  names(table) <- c(paste0("b_{", c(1:T)), paste0("E_gamma", c(1:T)))
  return(table)
  }
# how to not repeat building the table? 
table           <- matrix(0,2*T,5)
colnames(table) <- c("default", "ones", "minELBO", "maxELBO", "random")
rownames(table) <- c(paste0("\\ExpectationVarDist{b_{", c(1:T), "}}"), paste0("E_gamma", c(1:T)))
table[,1] <- table_biases(default)
table[,2] <- table_biases(allones)
table[,3] <- table_biases(minELBO)
table[,4] <- table_biases(maxELBO)
table[,5] <- table_biases(random) 
View(table)
xtable(table,digits=c(1,-3,3,-3,-3,3), caption=c("Biases and their precisions for models in grid"))
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
source("code_R/data_master.R")
source("code_R/bmtmkl_train_mse.R")
source("code_R/bmtmkl_train_mse_yunscaled.R")
#
random      <- bmtmkl_train_mse(Ktrain, ytrain, setparameters_from_grid(51688))
allones     <- bmtmkl_train_mse(Ktrain, ytrain, setparameters_from_grid(51668))
default     <- bmtmkl_train_mse(Ktrain, ytrain, setparameters_from_grid(1))
#
last_decreasing_iteration <- function(MSE_vector, Iter){
        change_sign          <- sign(MSE_vector[12:Iter] - MSE_vector[11:(Iter-1)])
        last_decreasing_iter <- head(which(change_sign==1), n=1) + 10  # augment by eleven since first ten discarded.
        return(last_decreasing_iter)
}
#test last_decreasing_iteration():
last_decreasing_iteration(c(20:1,2:4), length(c(20:1,2:4))) # result should be twenty
#
find_min_mses <- function(model){
      Iter <- model$parameters$iteration
      last_decreasing_iter_MSE_insample <- last_decreasing_iteration(model$MSE_insample, Iter)
      last_decreasing_iter_MSE_insample_rescaled <- last_decreasing_iteration(model$MSE_insample_rescaled, Iter)
      last_decreasing_iter_MSE_outofsample  <- last_decreasing_iteration(model$MSE_outofsample, Iter)
      cat(" MSE insample \t\t decreases until iteration: ", last_decreasing_iter_MSE_insample , "\n",
          "MSE insample rescaled \t decreases until iteration: ", last_decreasing_iter_MSE_insample_rescaled , "\n",
          "MSE out-of-sample\t decreases until iteration: ", last_decreasing_iter_MSE_outofsample , "\n")
}
find_min_mses(random)
find_min_mses(allones)
find_min_mses(default)
###
random_ns      <- bmtmkl_train_mse_yunscaled(Ktrain, ic50_train, setparameters_from_grid(51688, 400))
allones_ns     <- bmtmkl_train_mse_yunscaled(Ktrain, ic50_train, setparameters_from_grid(51668, 400))
default_ns     <- bmtmkl_train_mse_yunscaled(Ktrain, ic50_train, setparameters_from_grid(1, 500))
#here
find_min_mses_yunscaled <- function(model){
  Iter <- model$parameters$iteration
  last_decreasing_iter_MSE_insample <- last_decreasing_iteration(model$MSE_insample, Iter)
  #last_decreasing_iter_MSE_insample_rescaled <- last_decreasing_iteration(model$MSE_insample_rescaled, Iter)
  last_decreasing_iter_MSE_outofsample  <- last_decreasing_iteration(model$MSE_outofsample, Iter)
  cat(" MSE insample \t\t decreases until iteration: ", last_decreasing_iter_MSE_insample , "\n",
      #"MSE insample rescaled \t decreases until iteration: ", last_decreasing_iter_MSE_insample_rescaled , "\n",
      "MSE out-of-sample\t decreases until iteration: ", last_decreasing_iter_MSE_outofsample , "\n")
}
find_min_mses_yunscaled(random_ns)
find_min_mses_yunscaled(allones_ns)
find_min_mses_yunscaled(default_ns)
#Graphs illustrating trad-off between ELBO and MSEs 
### default:
s=2
x <- s: length(default$bounds)
y1 <- default$MSE_insample[s:length(default$bounds)]
y2 <- default$MSE_insample_rescaled[s:length(default$bounds)]
y3 <- default$MSE_outofsample[s: length(default$bounds)]
g4 <- default$bounds[s: length(default$bounds)]
par(mar=c(5,4,1,5)+.1) # set margin
plot(x,y1 ,type="p",col="red", xlab="Iteration", ylab="MSE", ylim=c(0,30))
par(new=TRUE)
lines(x, y2 ,type="p",col="blue")
par(new=TRUE)
lines(x, y3 ,type="p",col="grey")
par(new=TRUE)
plot(x, g4,type="p",col="green",xaxt="n",yaxt="n",xlab="",ylab="")
axis(4, at=c(0,30))
mtext("ELBO",side=4,line=3)
legend("right",col=c("red","blue", "grey", "green"),lty=1,legend=c("MSE in sample", "MSE in sample rescaled", "MSE out of sample","ELBO"))
### ones model:
s=2
x <- s: length(allones$bounds)
y1 <- allones$MSE_insample[s:length(allones$bounds)]
y2 <- allones$MSE_insample_rescaled[s:length(allones$bounds)]
y3 <- allones$MSE_outofsample[s: length(allones$bounds)]
g4 <- allones$bounds[s: length(allones$bounds)]
par(mar=c(5,4,1,5)+.1) # set margin
axis(2, at=c(0,30))
plot(x, y2 ,type="p",col="blue", xlab="Iteration", ylab="MSE" , ylim=c(0,30))
# par(new=TRUE)
lines(x,y1 ,type="p",col="red")
par(new=TRUE)
lines(x, y3 ,type="p",col="grey")
par(new=TRUE)
plot(x, g4,type="p",col="green",xaxt="n",yaxt="n",xlab="",ylab="")
axis(4)
mtext("ELBO",side=4,line=3)
legend("right",col=c("red","blue", "grey", "green"),lty=1,legend=c("MSE in sample", "MSE in sample rescaled", "MSE out of sample","ELBO"))
### 
### random model:
s=5
x <- s: length(random$bounds)
y1 <- random$MSE_insample[s:length(random$bounds)]
y2 <- random$MSE_insample_rescaled[s:length(random$bounds)]
y3 <- random$MSE_outofsample[s: length(random$bounds)]
g4 <- random$bounds[s: length(random$bounds)]
par(mar=c(5,4,1,5)+.1) # set margin
plot(x,y1 ,type="p",col="red", xlab="Iteration", ylab="MSE")
par(new=TRUE)
lines(x, y2 ,type="p",col="blue")
par(new=TRUE)
lines(x, y3 ,type="p",col="grey")
par(new=TRUE)
plot(x, g4,type="p",col="green",xaxt="n",yaxt="n",xlab="",ylab="")
axis(4)
mtext("ELBO",side=4,line=3)
legend("right",col=c("red","blue", "grey", "green"),lty=1,legend=c("MSE in sample", "MSE in sample rescaled", "MSE out of sample","ELBO"))
#######################################################################################################################
# ELBO-Cutoff of relative change of 10^(-4) would not have been reached for model of ones:
allones     <- bmtmkl_train_mse(Ktrain, ytrain, setparameters_from_grid(51668, 500))
I <- length(allones$bounds)
stop_criterium <- 1- (allones$bounds[2:I] / allones$bounds[1:(I-1)])
cat("number of times stop-criterium has been reached: ", sum(stop_criterium < 10^(-4)),"\n")
which(stop_criterium<10^(-4), arr.ind = FALSE)
# Comp. Decrease/Increase between iteration of default, random and ones model.
comp_rel_change_ELBO_MSE <- function(model,seq){
   no_col <- length(seq)
   mat <-  matrix(0,5, no_col)
   mat[1,]  <- model$bounds[seq]
   mat[2, 2:dim(mat)[2]] <-    model$bounds[seq[2:no_col]] / model$bounds[seq[1:(no_col-1)]] - 1
   mat[3, 2:dim(mat)[2]] <-    model$MSE_insample[seq[2:no_col]] / model$MSE_insample[seq[1:(no_col-1)]] -1
   mat[4, 2:dim(mat)[2]] <-    model$MSE_insample_rescaled[seq[2:no_col]] / model$MSE_insample_rescaled[seq[1:(no_col-1)]] -1
   mat[5, 2:dim(mat)[2]] <-    model$MSE_outofsample[seq[2:no_col]] / model$MSE_outofsample[seq[1:(no_col-1)]] -1
   colnames(mat) <- seq
   rownames(mat) <- c("ELBO", "Delta ELBO", "Delta MSE", "Delta MSE res.", "Delta MSE o.f.s.")
   mat
}
#
seq <- seq(5,100,10)
#
currenttable <- comp_rel_change_ELBO_MSE(allones, seq)
xtable(currenttable, digits=3)

#######################################################################################################################
#######################################################################################################################
# compare ranking of predictions of different models:
###################################################################################################
# ranking # rank largerst value with 1!
myrank <- function(x) {return(rank(-x, ties.method = "first"))}  # largest values have to be ranked first (i.e. largest= 1)
#! Sorting of cell-lines is changing ranking for 
# 1.) NA values: Ranked in order of appearance after all non-NA values
# 2.) Equal values: Ranked in order of appearance. 
###################################################################################################
install.packages("magrittr")
library(magrittr)
###################################################################################################
# reload test data including drug 26: 
filepath_testdataoutcome= "submission data/DREAM7_DrugSensitivity1_test_data.txt"
drugResponsetest=as.matrix(read.table(filepath_testdataoutcome, header=T, row.names=1)) 
drugResponsetest=preprocessData(drugResponsetest, cellline_names_test)
testingCelllineNames=rownames(drugResponsetest)
View(drugResponsetest)
cat("mean of drug 26 is constant: ",mean(drugResponsetest[, "Drug26"], na.rm=TRUE))
summary(drugResponsetest[, "Drug26"])
constantvalue_drug26 <- mean(drugResponsetest[, "Drug26"], na.rm=TRUE)
#
drugReponseALLTRUE <-rbind (traindata$DrugResponse[1:35,] , 
                            as.matrix(read.table(filepath_testdataoutcome, header=T, row.names=1)) 
)
###################################################################################################
# rank test data for comparison
#   a) Comparison Alltrue data ranked and only test data ranked
#   b) Comparision betw. alphabetially sorted test cell lines and order of challenge.
cellline_names_test_alphabetically <- sort(rownames(drugResponsetest))  # alphabetically ordered cell lines
drugResponsetest_alphabetically <- preprocessData(drugResponsetest, cellline_names_test_alphabetically)   # First sort (rank first and sort would produce other Ranking.)
rank_true_testdata_alphabetically    <- apply(drugResponsetest_alphabetically, 2, myrank)
rank_true_testdata_alphabetically    <- cbind(as.vector(rownames(rank_true_testdata_alphabetically)), rank_true_testdata_alphabetically)
colnames(rank_true_testdata_alphabetically)[1] <- "DrugAnonID"
View(rank_true_testdata_alphabetically)
write.csv(rank_true_testdata_alphabetically, file="./data/rankings/rank_true_testdata_alphabetically.csv", quote=FALSE, row.names=FALSE)
#
rank_true_testdata <- apply(drugResponsetest, 2, myrank)
rank_true_testdata <- cbind(as.vector(rownames(rank_true_testdata)), rank_true_testdata)
colnames(rank_true_testdata)[1] <- "DrugAnonID"
write.csv(rank_true_testdata ,file="./data/rankings/rank_true_testdata.csv", quote=FALSE, row.names=FALSE)
View(rank_true_testdata)
#
rank_true_all <- drugReponseALLTRUE %>% apply(., 2, myrank)
rank_true_all <- cbind(as.vector(rownames(rank_true_all)), rank_true_all)
colnames(rank_true_all)[1] <- "DrugAnonID"
write.csv(rank_true_all ,file="./data/rankings/rank_true_all.csv", quote=FALSE, row.names=FALSE)
View(rank_true_all)
#
rank_true_all_alphabetically <- drugReponseALLTRUE %>% preprocessData(.,sort(cellline_names_all)) %>% apply(., 2, myrank)
rank_true_all_alphabetically <- cbind(as.vector(rownames(rank_true_all_alphabetically)), rank_true_all_alphabetically)
colnames(rank_true_all_alphabetically)[1] <- "DrugAnonID"
write.csv(rank_true_all_alphabetically, file="./data/rankings/rank_true_all_alphabetically.csv", quote=FALSE, row.names=FALSE)
View(rank_true_all_alphabetically)
###################################################################################################
#
# rankings:
getpredictionstable <- function(predictions){
  predictiontable   <- matrix(0, nrow=18, ncol=T)
  for (t in 1:T) {  
    predictiontable[,t]       <- predictions$y[[t]]$mu
    }
    predictiontable <- cbind(predictiontable[,1:25], constantvalue_drug26, predictiontable[,26:T])
    dimnames(predictiontable) <- list(cellline_names_test, colnames(drugResponsetest))
  return(predictiontable)
}
# save ranking in file-format required in DREAM7 D7C4 challenge:
save_ranks_inchallengeformat <- function(ranks){
  ranks <- cbind(as.vector(rownames(ranks)), ranks)
  colnames(ranks)[1] <- "DrugAnonID"
  filename <- paste0("./data/rankings/", match.call()$ranks)
  print(filename)
  write.csv(ranks, file = filename, quote=FALSE, row.names=FALSE)
}
####
rank_allones <- pred_allones %>% getpredictionstable(.)%>% apply(.,2, myrank)
save_ranks_inchallengeformat(rank_allones)
#
rank_default <- pred_default %>% getpredictionstable(.)%>% apply(.,2, myrank)
save_ranks_inchallengeformat(rank_default)
#
rank_maxELBO <- pred_maxELBO %>% getpredictionstable(.)%>% apply(.,2, myrank)
save_ranks_inchallengeformat(rank_maxELBO)
#
rank_minELBO <- pred_maxELBO %>% getpredictionstable(.)%>% apply(.,2, myrank)
save_ranks_inchallengeformat(rank_minELBO)
#
rank_random <- pred_random %>% getpredictionstable(.)%>% apply(.,2, myrank)
save_ranks_inchallengeformat(rank_random)
#
cbind(rank_true_testdata[,10],rank_allones[,10], rank_default[,10], rank_maxELBO[,10], rank_random[,10])

#unscaled results:
rank_allones_ns <- pred_allones_ns %>% getpredictionstable(.)%>% apply(.,2, myrank)
save_ranks_inchallengeformat(rank_allones_ns)
#
rank_default_ns <- pred_default_ns %>% getpredictionstable(.)%>% apply(.,2, myrank)
save_ranks_inchallengeformat(rank_default_ns)
#
rank_maxELBO_ns <- pred_maxELBO_ns %>% getpredictionstable(.)%>% apply(.,2, myrank)
save_ranks_inchallengeformat(rank_maxELBO_ns)
#
rank_minELBO_ns <- pred_maxELBO_ns %>% getpredictionstable(.)%>% apply(.,2, myrank)
save_ranks_inchallengeformat(rank_minELBO_ns)
#
rank_random_ns <- pred_random_ns %>% getpredictionstable(.) %>% apply(., 2, myrank)
save_ranks_inchallengeformat(rank_random_ns)
###################################################################################################
get_ofs_predictions_of_grid <- function(i, iter, ytrain= ytrain){
  pred <- setparameters_from_grid(i,iter) %>% bmtmkl_train(Ktrain, ytrain, .) %>% bmtmkl_test(Ktest, .) %>% getpredictionstable(.)
  return(pred)
}
###################################################################################################
rank_random_40  <- get_ofs_predictions_of_grid(51688,40, ytrain) %>% apply(.,2, myrank)
save_ranks_inchallengeformat(rank_random_40)
rank_allones_40 <- get_ofs_predictions_of_grid(51668,40, ytrain) %>% apply(.,2, myrank)
save_ranks_inchallengeformat(rank_allones_40)
rank_default_40 <- get_ofs_predictions_of_grid(1,40, ytrain)     %>% apply(.,2, myrank)
save_ranks_inchallengeformat(rank_default_40)
rank_maxELBO_40 <- get_ofs_predictions_of_grid(52452,40, ytrain) %>% apply(.,2, myrank)
save_ranks_inchallengeformat(rank_maxELBO_40)
rank_minELBO_40 <- get_ofs_predictions_of_grid(73,40, ytrain)    %>% apply(.,2, myrank)
save_ranks_inchallengeformat(rank_minELBO_40)
#
set_min_mse_ofs <- c(45278,45359,45440) # out-of-sample min MSE
rank_min_mse_ofs_1 <- get_ofs_predictions_of_grid(45278, 40, ytrain) %>% apply(.,2, myrank)
save_ranks_inchallengeformat(rank_min_mse_ofs_1)
rank_min_mse_ofs_2 <- get_ofs_predictions_of_grid(45359, 40, ytrain) %>% apply(.,2, myrank)
save_ranks_inchallengeformat(rank_min_mse_ofs_2)
rank_min_mse_ofs_3 <- get_ofs_predictions_of_grid(45440, 40, ytrain) %>% apply(.,2, myrank)
save_ranks_inchallengeformat(rank_min_mse_ofs_3)
#
set_mse_min_ins <- c(59041,59042,59043) # in-sample min MSE
rank_min_mse_ins_1 <- get_ofs_predictions_of_grid(59041, 40, ytrain) %>% apply(.,2, myrank)
save_ranks_inchallengeformat(rank_min_mse_ins_1)
rank_min_mse_ins_2 <- get_ofs_predictions_of_grid(59042, 40, ytrain) %>% apply(.,2, myrank)
save_ranks_inchallengeformat(rank_min_mse_ins_2)
rank_min_mse_ins_3 <- get_ofs_predictions_of_grid(59043, 40, ytrain) %>% apply(.,2, myrank)
save_ranks_inchallengeformat(rank_min_mse_ins_3)
####
#### get preditions for optimal number of iterations determined previously:
# default:
rank_default_120 <- get_ofs_predictions_of_grid(1,120, ytrain)     %>% apply(.,2, myrank)
save_ranks_inchallengeformat(rank_default_120)
unlist(getMSE(setparameters_from_grid(1, 120)))
#
rank_default_145 <- get_ofs_predictions_of_grid(1,145, ytrain)     %>% apply(.,2, myrank)
save_ranks_inchallengeformat(rank_default_145)
unlist(getMSE(setparameters_from_grid(1, 145)))
#
rank_default_106 <- get_ofs_predictions_of_grid(1,106, ytrain)     %>% apply(.,2, myrank)
save_ranks_inchallengeformat(rank_default_106)
unlist(getMSE(setparameters_from_grid(1, 106)))
####
## ones:
rank_allones_33 <- get_ofs_predictions_of_grid(51668,33, ytrain) %>% apply(.,2, myrank)
save_ranks_inchallengeformat(rank_allones_33)
unlist(getMSE(setparameters_from_grid(51668, 33)))
#
rank_allones_31 <- get_ofs_predictions_of_grid(51668,31, ytrain) %>% apply(.,2, myrank)
save_ranks_inchallengeformat(rank_allones_31)
unlist(getMSE(setparameters_from_grid(51668, 31)))
#
rank_allones_41 <- get_ofs_predictions_of_grid(51668,41, ytrain) %>% apply(.,2, myrank)
save_ranks_inchallengeformat(rank_allones_41)
unlist(getMSE(setparameters_from_grid(51668, 41)))
####
## random:
rank_random_91  <- get_ofs_predictions_of_grid(51688,91, ytrain) %>% apply(.,2, myrank)
save_ranks_inchallengeformat(rank_random_91)
unlist(getMSE(setparameters_from_grid(51688, 91)))
#
rank_random_87  <- get_ofs_predictions_of_grid(51688,87, ytrain) %>% apply(.,2, myrank)
save_ranks_inchallengeformat(rank_random_87)
unlist(getMSE(setparameters_from_grid(51688, 87)))
#
rank_random_158  <- get_ofs_predictions_of_grid(51688,158, ytrain) %>% apply(.,2, myrank)
save_ranks_inchallengeformat(rank_random_158)
unlist(getMSE(setparameters_from_grid(51688, 158)))
###################################################################################################
rank_random_ns_220  <- get_ofs_predictions_of_grid(51688,220, ic50_train) %>% apply(.,2, myrank)
save_ranks_inchallengeformat(rank_random_ns_220)
#
rank_allones_ns_40 <- get_ofs_predictions_of_grid(51688,40, ic50_train) %>% apply(.,2, myrank)
save_ranks_inchallengeformat(rank_allones_ns_40)
#rank_allones_40 == rank_allones_ns_40
#
rank_default_ns_477 <- get_ofs_predictions_of_grid(1,477, ic50_train)     %>% apply(.,2, myrank)
save_ranks_inchallengeformat(rank_default_ns_477)
#
# rank_maxELBO_40 <- get_ofs_predictions_of_grid(52452,40, ic50_train) %>% apply(.,2, myrank)
# save_ranks_inchallengeformat(rank_maxELBO_40)
# rank_minELBO_40 <- get_ofs_predictions_of_grid(73,40, ic50_train)    %>% apply(.,2, myrank)
# save_ranks_inchallengeformat(rank_minELBO_40)
