###################################################################################################
# author:  Henry Webel, henry.webel@hu-berlin.de
# content: function to set hyperparameters for Bayesian Multiple Kernel Multiple Task Learning
###################################################################################################
# source("bmtmkl_setparfct.R") # load function in script
###################################################################################################
setgammadisthyperparameter <-function(alpha_lambda, beta_lambda, 
                                      alpha_upsilon, beta_upsilon, 
                                      alpha_gamma, beta_gamma, 
                                      alpha_omega, beta_omega,
                                      alpha_epsilon, beta_epsilon,
                                      iterations= 200, seed= 1606) {
  parameters <- list()
  #set the hyperparameters of gamma prior used for sample weights
  parameters$alpha_lambda <-  alpha_lambda
  parameters$beta_lambda <-   beta_lambda
  
  #set the hyperparameters of gamma prior used for intermediate noise
  parameters$alpha_upsilon <-  alpha_upsilon
  parameters$beta_upsilon <-  beta_upsilon
  
  #set the hyperparameters of gamma prior used for bias
  parameters$alpha_gamma <- alpha_gamma   #1e-10
  parameters$beta_gamma  <- beta_gamma #1e-10
  
  #set the hyperparameters of gamma prior used for kernel weights
  parameters$alpha_omega <- alpha_omega
  parameters$beta_omega <-  beta_omega
  
  #set the hyperparameters of gamma prior used for output noise
  parameters$alpha_epsilon <-  alpha_epsilon
  parameters$beta_epsilon <-  beta_epsilon
  
  ### IMPORTANT ###
  #For gamma priors, you can experiment with three different (alpha, beta) values
  #(1, 1) => default priors
  #(1e-10, 1e+10) => good for obtaining sparsity
  #(1e-10, 1e-10) => good for small sample size problems (like in Nature Biotechnology paper)
  
  #set the number of iterations
  parameters$iteration <- iterations
  
  #determine whether you want to calculate and store the lower bound values
  parameters$progress <- 1
  
  #set the seed for random number generator used to initalize random variables
  parameters$seed <- seed
  return(parameters)
}
###################################################################################################
# former Version in demo.R of Mehmet Gönen ########################################################
# #initalize the parameters of the algorithm
# hyper_alpha <-10^(-10)
# hyper_beta  <-10^(-10)
# setgammadisthyperparameter<- function(alpha, beta) {
# parameters <- list()
# hyper_alpha <- alpha
# hyper_beta  <- beta
# #set the hyperparameters of gamma prior used for sample weights
# parameters$alpha_lambda <-  hyper_alpha
# parameters$beta_lambda <-   hyper_beta
# 
# #set the hyperparameters of gamma prior used for intermediate noise
# parameters$alpha_upsilon <-  hyper_alpha
# parameters$beta_upsilon <-  hyper_beta
# 
# #set the hyperparameters of gamma prior used for bias
# parameters$alpha_gamma <- hyper_alpha   #1e-10
# parameters$beta_gamma  <- hyper_beta #1e-10
# 
# #set the hyperparameters of gamma prior used for kernel weights
# parameters$alpha_omega <- hyper_alpha
# parameters$beta_omega <-  hyper_beta
# 
# #set the hyperparameters of gamma prior used for output noise
# parameters$alpha_epsilon <-  hyper_alpha
# parameters$beta_epsilon <-  hyper_beta
# 
# ### IMPORTANT ###
# #For gamma priors, you can experiment with three different (alpha, beta) values
# #(1, 1) => default priors
# #(1e-10, 1e+10) => good for obtaining sparsity
#   #(1e-10, 1e-10) => good for small sample size problems (like in Nature Biotechnology paper)
# 
#   #set the number of iterations
#   parameters$iteration <- 300
# 
#   #determine whether you want to calculate and store the lower bound values
#   parameters$progress <- 1
# 
#   #set the seed for random number generator used to initalize random variables
#   parameters$seed <- 1606
#   return(parameters)
# }