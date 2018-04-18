# #calculate determinant ? - log determinant
logdet <- function(Sigma) {
  2 * sum(log(diag(chol(Sigma))))
}

# #replicate matrix - repmat
repmat <- function(M, row, column) {
  kronecker(matrix(1, row, column), M) # row x columsn matrix of ones
}

set.seed(parameters$seed)

T <- length(Km)
D <- matrix(0, T, 1)
N <- matrix(0, T, 1)
for (t in 1:T) {
  D[t] <- dim(Km[[t]])[1]
  N[t] <- dim(Km[[t]])[2]
}
P <- dim(Km[[1]])[3]

log2pi <- log(2 * pi)

lambda <- vector("list", T)
for (t in 1:T) {
  lambda[[t]] <- list(shape =matrix(parameters$alpha_lambda + 0.5, D[t], 1), scale = matrix(parameters$beta_lambda, D[t], 1))
}
upsilon <- list(shape =matrix(parameters$alpha_upsilon + 0.5 * N * P, T, 1), scale = matrix(parameters$beta_upsilon, T, 1))
a <- vector("list", T)
for (t in 1:T) {
  a[[t]] <- list(mean = matrix(rnorm(D[t]), D[t], 1), covariance =diag(1, D[t], D[t]))
}
G <- vector("list", T)
for (t in 1:T) {
  G[[t]] <- list(mean = matrix(rnorm(P * N[t]), P, N[t]), covariance =diag(1, P, P))
}
gamma <- list(shape =matrix(parameters$alpha_gamma + 0.5, T, 1), scale = matrix(parameters$beta_gamma, T, 1))
omega <- list(shape =matrix(parameters$alpha_omega + 0.5, P, 1), scale = matrix(parameters$beta_omega, P, 1))
epsilon <- list(shape =matrix(parameters$alpha_epsilon + 0.5 * N, T, 1), scale = matrix(parameters$beta_epsilon, T, 1))
be <- list(mean = rbind(matrix(0, T, 1), matrix(1, P, 1)), covariance =diag(1, T + P, T + P))

KmKm <- vector("list", T)
for (t in 1:T) {
  KmKm[[t]] <- matrix(0, D[t], D[t])
  for (m in 1:P) {
    KmKm[[t]] <- KmKm[[t]] + tcrossprod(Km[[t]][,,m], Km[[t]][,,m])
  }
  Km[[t]] <- matrix(Km[[t]], D[t], N[t] * P)
}

if (parameters$progress == 1) {
  bounds <- matrix(0, parameters$iteration, 1)
}

atimesaT.mean <- vector("list", T)
for (t in 1:T) {
  atimesaT.mean[[t]] <- tcrossprod(a[[t]]$mean, a[[t]]$mean) + a[[t]]$covariance
}
GtimesGT.mean <- vector("list", T)
for (t in 1:T) {
  GtimesGT.mean[[t]] <- tcrossprod(G[[t]]$mean, G[[t]]$mean) + N[t] * G[[t]]$covariance
}
btimesbT.mean <- tcrossprod(be$mean[1:T], be$mean[1:T]) + be$covariance[1:T, 1:T]
etimeseT.mean <- tcrossprod(be$mean[(T + 1):(T + P)], be$mean[(T + 1):(T + P)]) + be$covariance[(T + 1):(T + P), (T + 1):(T + P)]
etimesb.mean <- matrix(0, P, T)
for (t in 1:T) {
  etimesb.mean[,t] <- be$mean[(T + 1):(T + P)] * be$mean[t] + be$covariance[(T + 1):(T + P), t]
}
KmtimesGT.mean <- vector("list", T)
for (t in 1:T) {
  KmtimesGT.mean[[t]] <- Km[[t]] %*% matrix(t(G[[t]]$mean), N[t] * P, 1)
}
###################################################################################################
# #
# #Input: lists of Kernelmatrices (T times NxNxP), vector of labels, 
# # initial set parameters of DBG (list of different parameters)
# # bayesian_multitask_multiple_kernel_learning_train <- function(Km, y, parameters) {
# #preparations:
# set.seed(parameters$seed) # fixed seed parameter in  demo.R
# 
# #parameter
# T <- length(Km)       # Number of Tasks T (here: # of drugs)
# D <- matrix(0, T, 1)  # initialize vector with 0s of length T (# of drugs=T)
# N <- matrix(0, T, 1)  
# for (o in 1:T) {
#   D[o] <- dim(Km[[o]])[1] #rows    ->task-specific
#   N[o] <- dim(Km[[o]])[2] #columns ->task-specific # since kernel matrices are provided: both dimension are the same... (different for preditions on test data)
# }
# P <- dim(Km[[1]])[3]  # third dimension: gives the number of views-> fixed!!! (imputation necessary?)
# # some overall constants (total number of cell-lines and multiplicatives)
# log2pi <- log(2 * pi) # constant
# cat("overall number of drugs tested on cell lines:", sum(N))
# NP <- sum(N) * P
# NT <- sum(N) * T
# #################################################################################################################
# #### Gamma distribution: lambda is alpha and beta
# # q(lambda)
# lambda <- vector("list", T)  # creates a list of T entries
# #list of list with matrices of parameters (gamma dist):
# #initialize variational parameters of lambda,
# # which is the gamma prior used for the precision (inverse variance) of input kernel weights a. 
# ## ! ## note: variational shape (alpha) parameter is constant (i.e. not updated)!
# for (o in 1:T) {
#   lambda[[o]] <- list(shape = matrix(parameters$alpha_lambda + 0.5, D[o], 1), 
#                       scale = matrix(parameters$beta_lambda,        D[o], 1))  # inital: hyperpar.
# }
# # initialize parameters of variationa gamma dist used for precision of intermediate outputs
# upsilon <- list(shape = matrix(parameters$alpha_upsilon + 0.5 * N * P, T, 1), 
#                 scale = matrix(parameters$beta_upsilon,                T, 1))
# ############# 
# # q(A): mean and covariances for Normal-dist of input kernel weights
# a <- vector("list", T) # creates a list of T, normal 
# for (o in 1:T) {
#   #random initialization:  rnorm(n, mean = 0, sd = 1)
#   # why not take hyperparameters for initialization? todo
#   a[[o]] <- list(mean = matrix(rnorm(D[o]), D[o], 1), covariance = diag(1, D[o], D[o]))
# }
# #############
# # again. 
# G <- vector("list", T)
# for (o in 1:T) {
#   # random initialization:  rnorm(n, mean = 0, sd = 1)
#   G[[o]] <- list(mean = matrix(rnorm(P * N[o]), P , N[o]), covariance = diag(1, P, P))  # mean_o: PxN_o ; COV_o: PxP
# }
# #set the hyperparameters of gamma prior used for bias:
# gamma   <- list( shape = matrix(parameters$alpha_gamma   + 0.5, T, 1), scale = matrix(parameters$beta_gamma, T, 1))
# #set the hyperparameters of gamma prior used for kernel weights:
# omega   <- list( shape = matrix(parameters$alpha_omega   + 0.5, P, 1), scale = matrix(parameters$beta_omega, P, 1))
# #set the hyperparameters of gamma prior used for output noise/bias/intercept:
# epsilon <- list( shape = matrix(parameters$alpha_epsilon + 0.5 * N, T, 1),   scale = matrix(parameters$beta_epsilon, T, 1))
# # bias und e together. (intercept and kernel weights for target output estimation)
# be      <- list( mean  = rbind(matrix(0, T, 1), matrix(1, P, 1)), covariance = diag(1, T + P, T + P))
# ########################################################################################################################
# # list of Kernel-Matrices
# KmKm <- vector("list", T) # initialise list of length T
# for (o in 1:T) {
#   KmKm[[o]] <- matrix(0, D[o], D[o])    # N_o x N_o matrix of 0s at list position o # each task has set of available cell-lines?
#   #use: D[o], N[o], P
#   for (m in 1:P) {
#     # views squared K_t,k^\prime K_t,k
#     KmKm[[o]] <- KmKm[[o]] + tcrossprod(Km[[o]][,,m], Km[[o]][,,m])   ## hier werden die einzelnen Kernel-Matrizen der views quadriert addiert
#   }
#   Km[[o]] <- matrix(Km[[o]], D[o], N[o] * P) # original data is reshaped from 3dim array to 2dim (as in Viewer) # km_o = N_o + (P*N_o)
#   # Km (N_t x (N_t*P))  
# }
# 
# if (parameters$progress == 1) {
#   bounds <- matrix(0, parameters$iteration, 1)
# }
# #####################################################################################################################
# # some constructs for updating -> ToDo: Identify updates
# # Projection matrices: each randomly generated previously
# # Normal covariances (mean):  E[X?]= E[X]*E[X]+VAR[X]
# atimesaT.mean <- vector("list", T) # initialize vector with T  lists 
# for (o in 1:T) {
#   #  x %*% t(x) (tcrossprod).
#   atimesaT.mean[[o]] <- tcrossprod(a[[o]]$mean, a[[o]]$mean) + a[[o]]$covariance  # dimension? N_txN_t
# }
# # Why does N_t appear in this equation? -> product/sum over all cell-lines per drug.
# # recall what is saved in G[[o]]
# G <- vector("list", T)
# for (o in 1:T) {
#   # random initialization:  rnorm(n, mean = 0, sd = 1)
#   G[[o]] <- list(mean = matrix(rnorm(P * N[o]), P, N[o]), covariance = diag(1, P, P))  # mean_o: PxN_o ; COV_o: PxP
# }
# # -> product of normal-distributions with different means and same variance.
# GtimesGT.mean <- vector("list", T)
# for (o in 1:T) {
#   GtimesGT.mean[[o]] <- tcrossprod(G[[o]]$mean, G[[o]]$mean) + N[o] * G[[o]]$covariance  # PxP
# }
# #
# btimesbT.mean <- tcrossprod(be$mean[1:T], be$mean[1:T]) + be$covariance[1:T, 1:T]
# etimeseT.mean <- tcrossprod(be$mean[(T + 1):(T + P)], be$mean[(T + 1):(T + P)]) + be$covariance[(T + 1):(T + P), (T + 1):(T + P)]
# etimesb.mean  <- matrix(0, P, T)
# for (o in 1:T) {
#   etimesb.mean[,o] <- be$mean[(T + 1):(T + P)] * be$mean[o] + be$covariance[(T + 1):(T + P), o]
# }
# KmtimesGT.mean <- vector("list", T)
# for (o in 1:T) {
#   KmtimesGT.mean[[o]] <- Km[[o]] %*% matrix(t(G[[o]]$mean), N[o] * P, 1) # hat dimension N_t
# }
# #