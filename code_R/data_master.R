# Authors: Henry Webel
# code-Snippets from Ammad-ud-din Muhammad, https://users.ics.aalto.fi/ammad/ included
# data available at: https://www.synapse.org/#!Synapse:syn2785782

#DZA:
# setwd("C:/Users/Webel/Desktop/Versuche/statisticsmaster/cancerprototyp/")
# # Laptop:
# setwd("C:/Users/Henry/Documents/1_Statistik/statisticsmaster/cancerprototyp/")
# setwd("/home/enryh/Desktop/statisticsmaster/cancerprototyp/")
# source('code_R/DREAM7ChallengeDataLoaderfunctions.R')
# source('code_R/data_loading.R')
# # View(loadDream7Data)
# dir("task1/training")
# dataLocation="H:/WS1617/bioinformatics_topic/data/NCI_dream_challenge/task1/training/DrugSensitivity1/"
# list.files(dataLocation) 
# traindata <- loadDream7Data(dataLocation)  # provided by 	Ammad-ud-din Muhammad, https://users.ics.aalto.fi/ammad/
# rm("loadDream7Data", "preprocessData", "processExomeSeqDataset", "dataLocation")
###################################################################################################
load("traindata.RData")  ## in case no files have been provided
###################################################################################################
str(traindata)                      # first look, all have 53 cell-lines 
names_datasets <- names(traindata)  # list of dataset names
cellline_names_all <- rownames(traindata$DrugResponse)
###################################################################################################
# View(t(traindata$GeneExpression))
# View(t(traindata$CNV))  ## some missing values in between 
###################################################################################################
# install.packages("kernlab")  # install kernlab to compute kernel matrices
if(!require("kernlab")) install.packages("kernlab") 
library("kernlab")
# euclidean dist betw. two datasets (rows)
# if(!require("caret")) install.packages("caret") 
if(!require("Rfast")) install.packages("Rfast") 
# if(!require("ggplot2")) install.packages("ggplot2") 
library(Rfast)
# library(caret)
# install.packages("Rfast")
# install.packages("caret", dep=TRUE)
##################################################################################
# # input matrices for kernels have to be imputed. Missing values in columns of 
###################################################################################################
#' @title A function for processing the DREAM7 data sources
#'
#' \code{preprocessData} is a utility function and is only called
#' from \code{\link{loadDream7Data}}.
#'
#' @param inputData the input data source.
#' @param allCelllineNames cell line names such that 
#' 1:35 are training cell lines and 36:53 are test cell lines.
#' @return processedInputData processed data such that the order of the cell lines 
#' matches the order in the "allCelllineNames" variable and 
#' for missing cell lines, NA's are added
preprocessData <- function(inputData, allCelllineNames)
{
  
  processedInputData <- matrix(nrow=length(allCelllineNames), ncol=dim(inputData)[2])
  dim(processedInputData)
  rownames(processedInputData) <- allCelllineNames
  colnames(processedInputData) <- colnames(inputData)
  for(i in 1:length(allCelllineNames))
  {
    if(allCelllineNames[i] %in% rownames(inputData))
    {
      processedInputData[i,1:dim(inputData)[2]] <- inputData[which(allCelllineNames[i]==rownames(inputData)),1:dim(inputData)[2]]
    }
    else
    {
      processedInputData[i,] <- NA
    }
  }
  return(processedInputData=processedInputData)
  
}
##################################################################################
# select real valued views and process
selection <- c(2,3,4,6,7,9) # 6 realvalued views
# select first most varying features, then perform scaling.
i=0  # skip first dataset consisting of drugs sensitivity values
data_used = list()
## define variance function with "na.rm=TRUE"
myvarfct <- function(x) {
  return(var(x, na.rm=TRUE))
}
mymean <- function(x) {
  return(var(x, na.rm=TRUE))
}
## apply to all selected real-valued datasets:
for (dataset in traindata[selection]) {
  i = i + 1 
  data_used[[i]] <- dataset
  # data_used[[i]] <- dataset[ , sort((apply(dataset, 2, myvarfct)) , decreasing=TRUE, index.return=TRUE)$ix[1:5000]]
  # data_used[[i]] <- (data_used[[i]])[,1:5000]  
  data_used[[i]] <- scale(data_used[[i]])
  #data_used is already demeanded.
  data_used[[i]][is.na(data_used[[i]])] <- 0
}
str(data_used)
names(data_used) <- names(traindata)[selection]
(data_used[[3]])[1:10,1:10]
#################################################################################
views = list()
i=0
# help(rbfdot): The Gaussian RBF kernel k(x,x') = \exp(-?? \|x - x'\|^2) 
for (dataset in data_used) {
  i = i + 1
  sigma =  1 /(2* dim(dataset)[2])  # see definition of RBF kernel in "kernlab" package.
  views[[i]] <- kernelMatrix(rbfdot(sigma = sigma) , dataset) # compute kernel
  print(dim(views[[i]]))
  print((views[[i]])[10:15,10:15])
}
names(views) <- names(data_used) # reassign names for datasets
str(views)
# ##################################################################################
# kernelsToTrain <- array(0, dim=c(53,53, length(views)))
# for (i in 1:length(views)) {
#   kernelsToTrain[,,i] <- views[[i]]
# }
# #################################################################################
##################################################################################
### binary data:
binarydata <- c(5,8)
str(traindata[binarydata])
# https://stats.stackexchange.com/questions/176613/jaccard-similarity-in-r
# https://stats.stackexchange.com/questions/49453/calculating-jaccard-or-other-association-coefficient-for-binary-data-using-matri
library(Matrix)
jaccard <- function(m) {
  ## common values:
  A = tcrossprod(m)
  ## indexes for non-zero common values
  im = which(A > 0, arr.ind=TRUE)
  ## counts for each row
  b = rowSums(m)
  
  ## only non-zero values of common
  Aim = A[im]
  
  ## Jacard formula: #common / (#i + #j - #common)
  J = sparseMatrix(
    i = im[,1],
    j = im[,2],
    x = Aim / (b[im[,1]] + b[im[,2]] - Aim),
    dims = dim(A)
  )
  
  return( J )
}
#################################
i= length(views) 
for(dataset in traindata[binarydata]) {
  i = i +1
  views[[i]]<- as.matrix(jaccard(dataset))
  dimnames(views[[i]]) <- list(rownames(dataset), rownames(dataset))
  print(views[[i]])
  # a <- append(a, as.matrix(jaccard(traindata[[8]])))
}
names(views) <- names(traindata[c(selection, binarydata)])
str(views)
# View(views[[7]])
# View(views[[8]])
# rm("binarydata", "i", "selection", "sigma", "myvarfct", "jaccard", "dataset")
#################################################################################
#################################################################################
# add pathways: All real valued # not scaled?!
load("pathways.RData")
#
i= length(views)
for(dataset in pathways) {
  dataset[is.infinite(dataset)] <- 0
  dataset <- scale(dataset)
  dataset[is.na(dataset)] <- 0
  i = i + 1
  sigma =  1 /(2* dim(dataset)[2])  # see definition of RBF kernel in "kernlab" package.
  views[[i]] <- kernelMatrix(rbfdot(sigma = sigma) , dataset) # compute kernel
  print(dim(views[[i]]))
}
names(views) <- c(names(traindata[c(selection, binarydata)]), names(pathways))
for (k in 9:18){
  print(names(views)[[k]])
  print((views[[k]])[10:15,10:15])
}
##################################################################################################################################################################
# # add interactions on DNA based information
print(names(views)[1:6]) # interaction between real valued original datasets.
# #
i= length(views) # Expression * Methylation
i= i+1
views[[i]] <- views[[1]]*views[[3]] # Expression * Copy Number Variation (CNV)
i= i+1
views[[i]] <- views[[1]]*views[[2]] # Copy Number Variation (CNV) * Methylation
i= i+1
views[[i]] <- views[[2]]*views[[3]] # Expression * Copy Number Variation (CNV) * Methylation
i= i+1
views[[i]] <- views[[1]]*views[[2]]*views[[3]] # all three DNA measures
cat("a total of ",i, "views has been calculated\n")
names(views)<- c(names(traindata[c(selection, binarydata)]), names(pathways) , c("Expression_Methylation", "Expression_CNV","CNV_Methylation", "Expres_CNV_Methyl"))
names_views <- names(views)
##################################################################################################################################################################
##################################################################################################################################################################
# create a kernel matrix  
N_train   <- 35
N_test    <- 18
N_pathway <- length(views)
P         <- length(views)
##################################################################################
# createkernelmatrix <- function(views) {  ## has to return results...todo
N_train   <- 35
N_test    <- 18
N_pathway <- length(views)
K_train <- array(0, dim = c(N_train, N_train, N_pathway))
K_test <- array(0, dim = c(N_train, N_test, N_pathway))
###################################################################################
for (k in 1:N_pathway) {
  K_train[,,k] <- views[[k]][1:N_train,1:N_train] 
  K_test[,,k]  <- views[[k]][1:N_train , (N_train+1):(N_train+N_test)] 
}
rm("k")
cat("Dim: \tN_train, N_test, N_Views \nDim K_train:" ,
    dim(K_train) , "\nDim _test:  ", dim(K_test) ,"\n") 
####################################################################################
# name rows and columns
cellline_names_train <- cellline_names_all[1:N_train]
cellline_names_test  <- cellline_names_all[(N_train+1):(N_train+N_test)]
dimnames(K_train) <- list(cellline_names_train , cellline_names_train, names(views))
dimnames(K_test)  <- list(cellline_names_train , cellline_names_test, names(views))
####################################################################################
y_train <- traindata$DrugResponse[1:35,] # exclude some drugs c(-12,26,27) are not scored in wpc
# y_train <- traindata$DrugResponse[1:35,]
T <- dim(y_train)[2]
###################################################################################################
# #initialize the kernels and outputs of each task for training
Ktrain <- vector("list", T)
Ktest <- vector("list", T)
ytrain <- vector("list", T)
###################################################################################################
### y_train <- scale(y_train, center=TRUE, scale=TRUE)
y_train <- scale(y_train) # gives precise results... 
# mystddev <- function(x){
#   mean <- sum(x[!is.na(x)])/sum(!is.na(x))
#   sd <- sqrt(1/ (sum(!is.na(x))-1) * sum((x[!is.na(x)]- mean)^2))
# }
# sd_temp <- apply(y_train, 2, sd , na.rm=TRUE)
# mean_temp <- apply(y_train, 2, mean , na.rm=TRUE)
druglevel=attr(y_train, "scaled:center")
drugscale=attr(y_train, "scaled:scale")
#
test_indices <- c(36:53)
train_indices <- c(1:35)
rm_t <- NULL
for (t in 1:T)
{
  na_index <- which(is.na(y_train[,t])==0) ## not na
  N_na_index <- length(na_index)
  if(N_na_index!=0)
  {
    Ktrain[[t]] = K_train[na_index,na_index,]
    Ktest[[t]] = K_test[na_index,1:length(test_indices),]
    #
    dimnames(Ktrain[[t]]) <- list(cellline_names_train[na_index] , cellline_names_train[na_index], names(views))
    dimnames(Ktest[[t]])  <- list(cellline_names_train[na_index] , cellline_names_test, names(views))
    #
    ytrain[[t]] <- y_train[na_index,t]
  }else{
    rm_t <- c(rm_t,t)
  }
}

### y_train <-  Y[train_indices,]
if(length(rm_t)>0){
  Ktrain <- Ktrain[-rm_t]
  Ktest <- Ktest[-rm_t]
  ytrain <- ytrain[-rm_t]
  #y_test <- y_test[,-rm_t]
  y_train <- y_train[,-rm_t]
  #
  T <- T-length(rm_t)
  druglevel<- druglevel[-rm_t]
  drugscale<- drugscale[-rm_t]
}
cat("drug ",rm_t, " has been removed")
rm("na_index", "N_na_index")
# createkernelmatrix(views)  ## programm function to be able to work with different sets of data
###################################################################################################
# same for non-standardized out-comes... missing has to be set manually
ic50 <- traindata$DrugResponse[1:35,-c(26)] # exclude some drugs
ic50_train <- vector("list", T)
for (t in 1:T){
  na_index <- which(is.na(ic50[,t])==0) ## not na
  N_na_index <- length(na_index)
  ic50_train[[t]] <- ic50[na_index,t]
  }
###################################################################################################
cellline_names_train <- cellline_names_all[1:N_train]
cellline_names_test  <- cellline_names_all[(N_train+1):(N_train+N_test)]
###################################################################################################
filepath_testdataoutcome= "submission data/DREAM7_DrugSensitivity1_test_data.txt"
drugResponsetest=as.matrix(read.table(filepath_testdataoutcome, header=T, row.names=1))[,-rm_t]  
drugResponsetest=preprocessData(drugResponsetest, cellline_names_test)
testingCelllineNames=rownames(drugResponsetest)
###################################################################################################
rm("data_used", "dataset", "names_datasets", "views" )
###################################################################################################