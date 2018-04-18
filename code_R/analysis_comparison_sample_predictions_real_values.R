## author: Henry Webel
## content: tables for presentation on BMTMKL
## version: 20/01/2018
###################################################################################################
# Laptop:
setwd("C:/Users/Henry/Documents/1_Statistik/statisticsmaster/cancerprototyp/")
# Linux:
setwd("/home/enryh/Desktop/statisticsmaster/cancerprototyp/")
#
source("code_R/bayesian_multitask_multiple_kernel_learning_train_uncommented.R")
# file.edit("code_R/bayesian_multitask_multiple_kernel_learning_train_uncommented.R")
source("code_R/bayesian_multitask_multiple_kernel_learning_test.R")
# load(".RData")
## load data and compute:
source("code_R/data_master.R")         # load all necessary data and calculate kernels
# file.edit("code_R/data_master.R")    # inspect and understand!
source("code_R/bmtmkl_setparfct.R")
# file.edit("code_R/bmtmkl_setparfct.R")
###### setup ######################################################################################
latexfiles_path <- "/presentation"
#
if(!require("xtable")) install.packages("xtable")  # load package for latex tables
library("xtable")
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
if(!require("pastecs")) install.packages("pastecs") 
library("pastecs")
###
#take a sample of drugs
sample <- c(1,8,23)
cells  <- c(1:11)
###
# summary stats:
stat.desc(y_train[,1:14], basic=FALSE) 
stat.desc(y_train[,1:14], desc =FALSE)
stat.desc(y_train[,15:dim(y_train)[2]], basic=FALSE) 
stat.desc(y_train[,15:dim(y_train)[2]], desc =FALSE)

stat.desc(y_train[,sample])

stat.desc(y_train[,sample]) 
tab_desc <- xtable(stat.desc(y_train[,sample]))
#####################################################################################
parameters <- setgammadisthyperparameter(1,1,1,1,1,1,1,1,1,1,40)
print(unlist(parameters)) # check
state <- bayesian_multitask_multiple_kernel_learning_train(Ktrain, ytrain, parameters)
t(state$be$mu)
str(state)
######################################################################################
#combine out-of-sample predictions and true values:
colnames(Ktest[[1]]) == rownames(drugResponsetest)
#####################################################################################
prediction <- bayesian_multitask_multiple_kernel_learning_test(Ktest, state)
# create matrix of the true and predicted drug outcomes for the selected sample:
outofsample <- NULL
colnames<- vector("list", 2*length(sample))
i = 1
for (t in sample) {
  outofsample <- cbind(outofsample , drugResponsetest[cells,t], (prediction$y[[t]]$mu[cells] * drugscale[t] + druglevel[t] )  )
  colnames[(2*i-1):(2*i)] <- c(sprintf("D%dt" ,t) , sprintf("D%dp" ,t))
  i <- i+1
 # dimnames(outofsample)[2] <- c(sprintf("drug_%d_p" ,t), sprintf("drug_%d_t" ,t))
}
dimnames(outofsample)[[2]] <- as.list(colnames)
outofsample
xtable(outofsample)
currenttable <- xtable(outofsample)
print(currenttable, type="latex", file=paste0(getwd(),latexfiles_path,"/table_outofsample_predictions.tex")) #,NA.string = getOption("xtable.NA.string", "-"))
####################################################################################
#combine in-sample predictions and true values:
fit <- bayesian_multitask_multiple_kernel_learning_test(Ktrain, state)
#####################################################################################
# create matrix of the true and predicted drug outcomes for the selected sample:
trainingCelllineNames=rownames(y_train)
insample <- NULL
colnames<- vector("list", 2*length(sample))
i = 1
for (t in sample) {
  temp <- cbind((ytrain[[t]] * drugscale[t] + druglevel[t]), 
                (fit$y[[t]]$mu * drugscale[t] + druglevel[t]))
  dimnames(temp)[[2]]<- c(sprintf("D%dt" ,t) , sprintf("D%dp" ,t))
  insample <- cbind(insample, preprocessData(temp, trainingCelllineNames))
  # dimnames(insample)[2] <- c(sprintf("drug_%d_p" ,t), sprintf("drug_%d_t" ,t))
}
insample
currenttable <- xtable(insample[cells,])
print(currenttable, type="latex", file=paste0(getwd(),latexfiles_path,"/table_insample_predictions.tex"))#,NA.string = getOption("xtable.NA.string", "NA"))
####################################################################################
# table with hyperparameters:
currenttable <- (as.matrix(state$parameters))
colnames(currenttable)<- "set to "
print(xtable(currenttable), type="latex", file=paste0(getwd(),latexfiles_path,"/table_state_hyperparameters_set.tex"))
####################################################################################
# table with ELBO and relative change between iterations for selected model. 
last <- length(state$bounds)
currenttable <- (as.matrix(state$bounds[(last-9):last]))
colnames(currenttable)<- "Value of ELBO"
rownames(currenttable) <- (last-9):last
currenttable
ELBO_rel_change <- 1- state$bounds[2:last]/state$bound[1:(last-1)]  ## change in ELBO:
currenttable <- cbind(currenttable, 0)
colnames(currenttable)[2]  <- "rel. change"
currenttable[2:10,2] <- ELBO_rel_change[(last-9):(last-1)]
currenttable
print(xtable(currenttable), type="latex", file=paste0(latexfiles_path,"/table_ELBO.tex"))
#
drugno = 8
#
colnames(state$G[[drugno]]$mu) <- colnames(Ktrain[[drugno]])
rownames(state$G[[drugno]]$mu) <- c(paste0("mu_g_{10,",c(1:22),"}"))
#
expectation_upsilon <- function(model, drugno){print(model$upsilon$alpha[drugno]*model$upsilon$beta[drugno])}
#
xtable(t(state$G[[drugno]]$mu[1:6,])[fit_sample,] , digits=3, caption="Intermediate outputs of ones model")
expectation_upsilon(state, drugno)
library("pastecs")
#
xtable(stat.desc(t(state$G[[drugno]]$mu[1:6,]), basic=TRUE)[c(1,4:6,8,9,13),],
       digits=3, caption="intermediate outputs summary statistics - ones model" )  # allones
xtable(rbind(state$be$mu[31:36] , (state$omega$alpha * state$omega$beta)[1:6]),digits=4)
###########################################################


table_coef <- function(state, i){
  table <- c(state$be$mu[i],     # mean of bias of drug 10
             state$be$mu[(T+1):(T+P)],    # mean of kernel coefficients for 22 views/input kernels
             state$gamma$alpha[i] * state$gamma$beta[i] , # precision of bias
             state$omega$alpha * state$omega$beta , # precisions for kernel coefficients
             state$epsilon$alpha[i] * state$epsilon$beta[i])  # precision for outcome
  #
  names(table)<- c(paste0("bias_drug_", i), paste0("e_", c(1:22)),paste0("E_gamma",i), paste0("E_omega", c(1:22)), paste0("E_epsilon_", i))
  #
  return(table)
}
#
tab_coef_8 <- as.matrix(table_coef(state, 8))
xtable((tab_coef_8), digits=4)
