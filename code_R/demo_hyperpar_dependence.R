#review...
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
T <- 28  ## in Costello et. al. only 28 drugs were mentioned
#set the number of kernels (e.g., the number of views in Nature Biotechnology paper)
P <- 22
## load data and compute:
source("code_R/data_master.R")
source("code_R/bmtmkl_setparfct.R")
parameters_default <- setgammadisthyperparameter(10^(-10),10^(-10),10^(-10),10^(-10),10^(-10),10^(-10),10^(-10),10^(-10),10^(-10),10^(-10),200,1606)
parameters_allone <- setgammadisthyperparameter(1,1,1,1,1,1,1,1,1,1,200,1606)
######################################################################################
#combine out-of-sample predictions and true values:
colnames(Ktest[[1]]) == rownames(drugResponsetest)
###################################################################################################
getfullmodel <- function(parameters){
  start.time <- Sys.time()
  state <- bayesian_multitask_multiple_kernel_learning_train(Ktrain, ytrain, parameters)
  trainingfit <- bayesian_multitask_multiple_kernel_learning_test(Ktrain, state)
  prediction  <- bayesian_multitask_multiple_kernel_learning_test(Ktest, state)
  ### squared error (in-sample-fit)
  evaluationstats = list(NULL)
  MSE_insample_rescaled = 0
  for (t in 1:T){
    # cat(colnames(y_train)[t], "\t has ", sum(!is.na(y_train[,t])), "\t observed values. \n")
    MSE_insample_rescaled <- MSE_insample_rescaled + 
       1/ sum(!is.na(y_train[,t])) * 
       sum((((ytrain[[t]]* drugscale[t] + druglevel[t]) - (trainingfit$y[[t]]$mu * drugscale[t] + druglevel[t]))^2))
  }
  evaluationstats$MSE_insample_rescaled <- MSE_insample_rescaled
  #
  MSE_insample = 0
  for (t in 1:T){
    # cat(colnames(y_train)[t], "\t has ", sum(!is.na(y_train[,t])), "\t observed values. \n")
    MSE_insample <- ( MSE_insample + 
                        1/ sum(!is.na(y_train[,t])) * 
                        sum(((ytrain[[t]] - trainingfit$y[[t]]$mu)^2)))
  } 
  evaluationstats$MSE_insample <- MSE_insample  # reference valus (TSS) is just T, so here 28.
  ### MSE for out-of-sample
  MSE_outofsample_rescaled_pred = 0
  for (t in 1:T){
    # cat(colnames(drugResponsetest)[t], "\t has ", sum(!is.na(drugResponsetest[,t])), "\t observed values. \n")
    MSE_outofsample_rescaled_pred <- MSE_outofsample_rescaled_pred + 
            1/ sum(!is.na(drugResponsetest[,t])) * 
            sum(((drugResponsetest[,t] - (prediction$y[[t]]$mu * drugscale[t] + druglevel[t]))^2) , na.rm=TRUE)
  }
  evaluationstats$MSE_outofsample_rescaled_pred <- MSE_outofsample_rescaled_pred
  #
  MSE_outofsample_scaled_drug = 0
  for (t in 1:T){
    # cat(colnames(drugResponsetest)[t], "\t has ", sum(!is.na(drugResponsetest[,t])), "\t observed values. \n")
    MSE_outofsample_scaled_drug <- MSE_outofsample_scaled_drug + 
      1/ sum(!is.na(drugResponsetest[,t])) * 
      sum(((scale(drugResponsetest[,t]) - (prediction$y[[t]]$mu))^2) , na.rm=TRUE)
  } 
  evaluationstats$MSE_outofsample_scaled_drug <- MSE_outofsample_scaled_drug 
  #
  MSE_outofsample_scaled_drug_2 = 0
  for (t in 1:T){
    # cat(colnames(drugResponsetest)[t], "\t has ", sum(!is.na(drugResponsetest[,t])), "\t observed values. \n")
    MSE_outofsample_scaled_drug_2 <- MSE_outofsample_scaled_drug_2 + 
      1/ sum(!is.na(drugResponsetest[,t])) * 
      sum(((((drugResponsetest[,t]- druglevel[t])/drugscale[t] ) - (prediction$y[[t]]$mu))^2) , na.rm=TRUE)
  }
  evaluationstats$MSE_outofsample_scaled_drug_2 <- MSE_outofsample_scaled_drug_2
  ###
  RSS  <- 0
  for (t in 1:T) {  RSS <- RSS + sum((ytrain[[t]]-trainingfit$y[[t]]$mu)^2) }
  evaluationstats$RSS <- RSS
  evaluationstats$TSS <- sum(unlist(ytrain)^2)  # equals N - T
  evaluationstats$R_sq <- 1 - ( evaluationstats$RSS / evaluationstats$TSS)
  ### squared error (sq err): all observations in test sample of 18 cell lines - out-of-sample
  RSS_rescaled = 0
  for (t in 1:T){
    # cat(colnames(drugResponsetest)[t], "\t has ", sum(!is.na(drugResponsetest[,t])), "\t observed values. \n")
    RSS_rescaled <- RSS_rescaled + 
      sum(((drugResponsetest[,t] - (prediction$y[[t]]$mu * drugscale[t] + druglevel[t]))^2) , na.rm=TRUE)
  }
  evaluationstats$RSS_rescaled <- RSS_rescaled
  ##
  R_Squared_0_test<- 0
  for (t in 1:T) {
    y_hat <- prediction$y[[t]]$mu 
    y     <- y_train[,t]
    s     <- ((sum(((drugResponsetest[,t]-druglevel[t])/drugscale[t]) * (y_hat), na.rm=TRUE))
              / (sum((y_hat)[!is.na(drugResponsetest[,t])])^2))
    R_Squared_0_test <- sum((((drugResponsetest[,t]-druglevel[t])/drugscale[t]) - s * y_hat)^2 , na.rm=TRUE) / sum((((drugResponsetest[,t]-druglevel[t])/drugscale[t]) - 0 )^2, na.rm=TRUE)
  }
  R_Squared_0_test <- 1- R_Squared_0_test
  ###
  ############################
  R_Squared_0_test_rescaled <- 0
  for (t in 1:T) {
    y_hat <- prediction$y[[t]]$mu * drugscale[t] + druglevel[t]
    y     <- drugResponsetest[,t]
    s     <- ((sum(drugResponsetest[,t] * (y_hat), na.rm=TRUE))
              / (sum((y_hat)[!is.na(drugResponsetest[,t])])^2))
    R_Squared_0_test_rescaled <- sum((drugResponsetest[,t] - s * y_hat)^2 , na.rm=TRUE) / sum((drugResponsetest[,t] - mean(y_train[,t], na.rm=TRUE))^2, na.rm=TRUE)
  }
  R_Squared_0_test_rescaled <- 1- R_Squared_0_test_rescaled
  ###
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  ### print results:
  cat("\n################################################################################\n",
      "Time:\t\t ", time.taken , "\n",
      "ELBO after ",parameters$iteration, ":\t",state$bounds[parameters$iteration],"\n",
      "Mean-Squared-Error in-sample:    \t", MSE_insample , "\n",
      "Squared-Error out-of-sample:    \t", RSS_rescaled , "\n",
      "Mean-Squared-Error out-of-sample:\t ", MSE_outofsample_rescaled_pred , "\n",
      "R-Squared 0 on out-of-sample:\t", R_Squared_0_test, "\n",
      "gammapriors: \t")
      cat(unlist(state$parameters[1:10]),"\n", sep=", ") 
  cat("##################################################################################")
  #
  # cat(unlist(parameters$iteration), unlist(state$bounds[parameters$iteration]), sq_err , unlist(state$parameters[1]),unlist(state$parameters[2]),
  #     unlist(state$parameters[3]),unlist(state$parameters[4]),unlist(state$parameters[5]),unlist(state$parameters[6]),
  #     unlist(state$parameters[7]),unlist(state$parameters[8]),unlist(state$parameters[9]),unlist(state$parameters[10]),"\n",
  #     file="gridsearch_v2.csv" , append=TRUE, sep=";")
  return(list("state" = state, "trainfit" = trainingfit, "predictions" = prediction, "evaluationstats" =evaluationstats))
}
###################################################################################################
default <- getfullmodel(parameters_default)
allone <- getfullmodel(parameters_allone)
###################################################################################################
# ### build true IC50 values matrix including NAs.
drugReponseALLTRUE <-rbind (traindata$DrugResponse[1:35,] , 
                             as.matrix(read.table(filepath_testdataoutcome, header=T, row.names=1)) 
                          )
drugResponsetest=preprocessData(drugResponsetest, cellline_names_test)
str(drugReponseALLTRUE)
View(drugReponseALLTRUE)
# rank and save  (all 53 cell lines and 31 drugs):
# rank predicitons
myrank <- function(x) {return(rank(x, ties.method = "first"))}
drug_TrueRanking <- apply(-drugReponseALLTRUE, 2, myrank)
View(drug_TrueRanking)
file_conda<-"/home/enryh/anaconda3/lib/python3.6/site-packages/dreamtools/dream7/D7C4/templates/rankingtrue.csv"
drug_TrueRanking <- cbind(rownames(drugReponseALLTRUE), drug_TrueRanking)
colnames(drug_TrueRanking)=  c("DrugAnonID", colnames(traindata$DrugResponse))
write.csv(drug_TrueRanking, file="./data/rankingtrue.csv", quote=FALSE, row.names=FALSE)
write.csv(drug_TrueRanking, file=file_conda, quote=FALSE, row.names=FALSE)
# produce tables to include in analysis:
if(!require("xtable")) install.packages("xtable") 
library("xtable")
currenttable <- xtable(outofsample)
print(currenttable, type="latex", file=paste0(latexfiles_path,"/table_outofsample_predictions.tex"))
###################################################################################################
model <- getfullmodel(parameters_allone)
predictiontable <- matrix(0, nrow=18, ncol=T)
for (t in 1:T) {
  # predictiontable[,t] <- c(traindata$DrugResponse[1:35,t], (model$predictions$y[[t]]$mu * drugscale[t] + druglevel[t]))
  predictiontable[,t] <- model$predictions$y[[t]]$mu 
}
# dimnames(predictiontable) <- list(cellline_names_all, colnames(drugResponsetest))
dimnames(predictiontable) <- list(cellline_names_test, colnames(drugResponsetest))
View(predictiontable)
#####
# y_train <- traindata$DrugResponse[1:35,c(-5,-24,-26)] # exclude some drugs
y_train_31drugs <- traindata$DrugResponse[1:35,] # exclude some drugs
y_train_31drugs <- scale(y_train_31drugs)
druglevel_31drugs=attr(y_train_31drugs, "scaled:center")
drugscale_31drugs=attr(y_train_31drugs, "scaled:scale")
# ### add previously excluded drug (with their training mean)  # add drugs as they have been?!
# mymean <- function(x) { return(mean(x, na.rm=TRUE))}
# means <- apply(traindata$DrugResponse, 2 , mymean)
#insert truth instead of mean...
# mean(y_train_31drugs[,5], na.rm=TRUE)
testtable <- rbind(y_train_31drugs,
             cbind(predictiontable[,c(1:4)], 0, predictiontable[,c(5:22)], 0, predictiontable[,c(23)], NaN, predictiontable[,c(24:T)])
              )
# insert training data...

# colnames(testtable) <- colnames(traindata$DrugResponse)
View(testtable)
# rank predicitons
myrank <- function(x) {return(rank(-x, ties.method = "first"))}  # largest values have to be ranked first (i.e. largest= 1)
testtable_ranked <- apply(testtable, 2, myrank)
View(testtable_ranked)
#
file_conda <- "/home/enryh/anaconda3/lib/python3.6/site-packages/dreamtools/dream7/D7C4/templates/test.csv"
testtable_ranked <- cbind(rownames(testtable_ranked), testtable_ranked)
# testtable_ranked <- cbind(cellline_names_test, testtable_ranked)
colnames(testtable_ranked)=  c("DrugAnonID", paste0("Drug",1:31))
View(testtable_ranked)
write.csv(testtable_ranked, file="./data/test.csv", quote=FALSE, row.names=FALSE)
write.csv(testtable_ranked, file=file_conda, quote=FALSE, row.names=FALSE)
write.csv(testtable_ranked[order(row.names(testtable)),], file="./data/test_ordered.csv", quote=FALSE, row.names=FALSE)
write.csv(testtable_ranked[order(row.names(testtable)),], 
          file= "/home/enryh/anaconda3/lib/python3.6/site-packages/dreamtools/dream7/D7C4/templates/test_orderd.csv", quote=FALSE, row.names=FALSE)
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
# train model:
parameters <- setgammadisthyperparameter(1E-10,	1,	1E-10,	1E-10,	1E-10,	1E-10,	1E-10,	1E-10,	1,	10000000000, 200, 1606)
state <- bayesian_multitask_multiple_kernel_learning_train(Ktrain, ytrain, parameters)
trainingfit <- bayesian_multitask_multiple_kernel_learning_test(Ktrain, ytrain, parameters)
prediction <- bayesian_multitask_multiple_kernel_learning_test(Ktest, state)
drugs <- traindata
# create matrix of predicted drug outcomes for the selected sample and run python scripts:  
trainingCelllineNames=rownames(y_train)
predictions_realvalues <- NULL
i = 1
for (t in 1:T) {
  temp <- (prediction$y[[t]]$mu * drugscale[t] + druglevel[t])
  colnames(temp)<- c(sprintf("drug%d" ,t))
  predictions_realvalues <- cbind(predictions_realvalues, preprocessData(temp, CelllineNames))
  # dimnames(insample)[2] <- c(sprintf("drug_%d_p" ,t), sprintf("drug_%d_t" ,t))
}
dim(predictions_realvalues)
### rank
# https://stackoverflow.com/questions/17888940/rank-all-columns-in-matrix-and-then-reorder-a-different-matrix-using-rank
predictions_ranks <- apply(predictions_realvalues ,2, function(x){rank(x, ties.method="random")})
### call Python-Script:
# https://www.r-bloggers.com/calling-python-from-r-with-rpython/
###################################################################################################
# select different hyperparameters and get the same squared error: 
source("code_R/bmtmkl_setparfct.R")
parameters_1 <- setgammadisthyperparameter(1E-10,	1,	1E-10,	1E-10,	1E-10,	1E-10,	1E-10,	1,	1,	1E-10, 200, 1606)
parameters_2 <- setgammadisthyperparameter(1E-10,	1,	1E-10,	1E-10,	1E-10,	1E-10,	1E-10,	1,	1,	1E-10, 200, 1506)
parameters_3 <- setgammadisthyperparameter(1E-10,	1,	1E-10,	1E-10,	1E-10,	1E-10,	1E-10,	1,	1,	1E-10, 600, 1606)
model_1 <- runmodel(parameters_1)
runmodel(parameters_2)
runmodel(parameters_3)

###################################################################################################
set <- c(1E-10,	1,	1E-10,	1E-10,	1E-10,	1E-10,	1E-10,	1E-10,	1,	10)
parameters_1 <- setgammadisthyperparameter(1E-10,	1,	1E-10,	1E-10,	1E-10,	1E-10,	1E-10,	1E-10,	1,	10, 200, 1606)
parameters_2 <- setgammadisthyperparameter(1E-10,	1,	1E-10,	1E-10,	1E-10,	1E-10,	1E-10,	1E-10,	1,	10, 200, 1506)
parameters_3 <- setgammadisthyperparameter(1E-10,	1,	1E-10,	1E-10,	1E-10,	1E-10,	1E-10,	1E-10,	1,	10, 600, 1606)
runmodel(parameters_1)
runmodel(parameters_2)
runmodel(parameters_3)

### relative veränderung von 0.01