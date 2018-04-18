# Mehmet Gonen (mehmet.gonen@gmail.com)
# Helsinki Institute for Information Technology HIIT
# Department of Information and Computer Science
# Aalto University School of Science
# source: https://github.com/mehmetgonen/bmtmkl/blob/master/bayesian_multitask_multiple_kernel_learning_train.R
###################################################################################################
bmtmkl_test <- function(Km, state) {
  T <- length(Km)
  N <- matrix(0, T, 1)
  for (t in 1:T) {
    N[t] <- dim(Km[[t]])[2]
  }
  P <- dim(Km[[1]])[3]
  
  G <- vector("list", T)
  for (t in 1:T) {
    G[[t]] <- list(mu = matrix(0, P, N[t]), sigma = matrix(0, P, N[t]))
    for (m in 1:P) {
      G[[t]]$mu[m,]   <-     crossprod(state$a[[t]]$mu , Km[[t]][,,m])  #mean
      G[[t]]$sigma[m,] <- 1 / (state$upsilon$alpha[t] * state$upsilon$beta[t]) + diag(crossprod(Km[[t]][,,m], state$a[[t]]$sigma) %*% Km[[t]][,,m])
    } 
  }
  
  y <- vector("list", T)
  for (t in 1:T) {
    y[[t]] <- list(mu = matrix(0, N[t], 1), sigma = matrix(0, N[t], 1))
    y[[t]]$mu <- crossprod(rbind(matrix(1, 1, N[t]), G[[t]]$mu), state$be$mu[c(t, (T + 1):(T + P))])
    y[[t]]$sigma <- 1 / (state$epsilon$alpha[t] * state$epsilon$beta[t]) + diag(crossprod(rbind(matrix(1, 1, N[t]), G[[t]]$mu), state$be$sigma[c(t, (T + 1):(T + P)), c(t, (T + 1):(T + P))]) %*% rbind(matrix(1, 1, N[t]), G[[t]]$mu))
  }
  
  prediction <- list(G = G, y = y)
}
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