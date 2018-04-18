# Mehmet Gonen (mehmet.gonen@gmail.com)
# Helsinki Institute for Information Technology HIIT
# Department of Information and Computer Science
# Aalto University School of Science
###################################################################################
# source:  https://github.com/mehmetgonen/kbmtl 
# https://github.com/mehmetgonen/bmtmkl
# article: https://doi.org/10.1093/bioinformatics/btu464
###################################################################################
# in order to play with data-> load manually:
# DZA:
#   setwd("C:/Users/Webel/Desktop/Versuche/statisticsmaster/cancerprototyp/")
# # Laptop:
# setwd("C:/Users/Henry/Documents/1_Statistik/statisticsmaster/cancerprototyp/")
# Linux Latop:
setwd("/home/enryh/Desktop/statisticsmaster/cancerprototyp/")
## load data and compute:
source("code_R/data_master.R")
source("code_R/bmtmkl_setparfct.R")
## set three parameters which the function takes by hand. 
###################################################################################
Km <- Ktrain
y  <- ytrain
parameters <- setgammadisthyperparameter(1,1,1,1,1,1,1,1,1,1,200,1606)
# parameters <- setgammadisthyperparameter(10^{-10},10^{-10},10^{-10},10^{-10},10^{-10},10^{-10},10^{-10},10^{-10},10^{-10},10^{-10},200,1606)
source("ch_2_initialization_random_variable_parameters.R", echo=TRUE)
###############################################################################################
# train 10 times:
source("ch_2.4_updates_variational_parameters.R")
# todo: get uncommented algorithm and write 10fold update.  check: where is "mu" from , where "mean".
###############################################################################################
source("ch_2.5.2_gradients.R")

source("ch2_5.2_ELBO_for_gradient_simple.R")
source("code_R/ch_2.5.2_gradients_simple.R")
source("ch2_5.2_hyperpar_gradients_emprical.R")
#
# check gradients empirically:
  # gradient lambda shape/alpha hyperparameter: 415.577
  # gradient lambda scale/theta hyperparameter: 204.7362
  # gradient upsilon shape/alpha hyperparameter: 29.53556
  # gradient upsilon scale/theta: -28.81287
# differential quotient: (ELBO(alpha+h)- ELBO(alpha))/h 
source("ch_2.3_ELBO.R")
# differencequotient <- function(hyperparameter, interval){
#     h  <- interval    # assign given interval for difference
#     assign(deparse(substitute(hyperparameter)), hyperparameter + h)
#     lb_h <- ELBO_for_gradient_calc((hyperparameter))
#     print(lb_h)
#     assign(deparse(substitute(hyperparameter)), hyperparameter - h)
#     lb <- ELBO_for_gradient_calc(hyperparameter)
#     print(lb)
#     print(hyperparameter+h)
#     # assign(deparse(substitute(hyperparameter)), (hyperparameter - h), env=globalenv())
#     print ((lb_h - lb)/ h)
#     return((lb_h - lb)/ h)
# }
# differencequotient(parameters$alpha_lambda, 10^(-1))
# # ToDo: Only PointerSet? How to give an alias to the original parameter?
# gradient lambda shape/alpha hyperparameter: 257.138 
# gradient lambda scale/theta hyperparameter: 52.97952 
source("ch2_5.2_hyperpar_gradients_emprical.R")
#
####
# gradient lambda shape/alpha hyperparameter: 415.577
h <- 10^(-6)
parameters$alpha_lambda <- parameters$alpha_lambda + h
lb_h <- ELBO()
parameters$alpha_lambda <- parameters$alpha_lambda - h
lb <- ELBO()
dq <- (lb_h- lb)/ h
cat("approx. gradient lambda shape/alpha hyperparameter:", dq, "\n")
# gradient lambda scale/theta hyperparameter: 275.742 
# approx. gradient lambda shape/alpha hyperparameter: 341.5846 
h <- 10^(-7)
parameters$beta_lambda <- parameters$beta_lambda + h
lb_h <- ELBO()
parameters$beta_lambda <- parameters$beta_lambda - h
lb <- ELBO()
dq <- (lb_h- lb)/ h
cat("approx. gradient lambda scale/beta hyperparameter:", dq ,"\n")
# approx. gradient lambda scale/beta hyperparameter: 209.8016
#####################################################################
# gradient gamma shape/alpha hyperparameter: 29.53556
h <- 10^(-7)
parameters$alpha_gamma <- parameters$alpha_gamma + h
lb_h <- ELBO()
parameters$alpha_gamma <- parameters$alpha_gamma - h
lb <- ELBO()
dq <- (lb_h- lb)/ h
cat("approx. gradient gamma shape/alpha hyperparameter:", dq, "\n")
# gradient gamma scale/theta: 17.97984
h <- 10^(-7)
lb <- ELBO()
parameters$beta_gamma <- parameters$beta_gamma + h
lb_h <- ELBO()
parameters$beta_gamma <- parameters$beta_gamma - h
dq <- (lb_h- lb)/ h
cat("approx. gradient lambda scale/beta hyperparameter:", dq, "\n")
#####################################################################
# gradient omega shape/alpha hyperparameter: 29.53556
h <- 10^(-7)
parameters$alpha_omega <- parameters$alpha_omega + h
lb_h <- ELBO()
parameters$alpha_omega <- parameters$alpha_omega - h
lb <- ELBO()
dq <- (lb_h- lb)/ h
cat("approx. gradient omega shape/alpha hyperparameter:", dq, "\n")
# gradient omega scale/theta: 17.97984
h <- 10^(-7)
lb <- ELBO()
parameters$beta_omega <- parameters$beta_omega + h
lb_h <- ELBO()
parameters$beta_omega <- parameters$beta_omega - h
dq <- (lb_h- lb)/ h
cat("approx. gradient lambda scale/beta hyperparameter:", dq, "\n")
#####################################################################
# gradient upsilon shape/alpha hyperparameter: 29.53556
h <- 10^(-7)
parameters$alpha_upsilon <- parameters$alpha_upsilon + h
lb_h <- ELBO()
parameters$alpha_upsilon <- parameters$alpha_upsilon - h
lb <- ELBO()
dq <- (lb_h- lb)/ h
cat("approx. gradient upsilon shape/alpha hyperparameter:", dq, "\n")
# gradient upsilon scale/theta: 17.97984
h <- 10^(-7)
lb <- ELBO()
parameters$beta_upsilon <- parameters$beta_upsilon + h
lb_h <- ELBO()
parameters$beta_upsilon <- parameters$beta_upsilon - h
dq <- (lb_h- lb)/ h
cat("approx. gradient lambda scale/beta hyperparameter:", dq, "\n")
#####################################################################
# gradient epsilon shape/alpha hyperparameter: 29.53556
h <- 10^(-7)
parameters$alpha_epsilon <- parameters$alpha_epsilon + h
lb_h <- ELBO()
parameters$alpha_epsilon <- parameters$alpha_epsilon - h
lb <- ELBO()
dq <- (lb_h- lb)/ h
cat("approx. gradient epsilon shape/alpha hyperparameter:", dq, "\n")
# gradient epsilon scale/theta: 17.97984
h <- 10^(-7)
lb <- ELBO()
parameters$beta_epsilon <- parameters$beta_epsilon + h
lb_h <- ELBO()
parameters$beta_epsilon <- parameters$beta_epsilon - h
dq <- (lb_h- lb)/ h
cat("approx. gradient lambda scale/beta hyperparameter:", dq, "\n")
###############################################################################################
  # Hyperparameters: ELBO depends on Hyperparameters!
    b <- 0.5
    c <- 10
    kappa <- 0.9
    learningrate <-  b * (c+ l)^(-kappa)    # reset learning rate
    function(parameter, gradient){
      return(parameter + l^learnrate * gradient)
    }
    # l=1 # fix in the start
    # new = old + (l)^hyperlernrate * gradient 
    # update hyperparameters every 10th iteration:
    # if (iter%%10 == 0) {
    # gradient lambda alpha
    # parameters$alpha_lambda <- parameters$alpha_lambda + l^learnrate * gradient_lambda_alpha
    # parameters$beta_lambda  <- 
    # parameters$alpha_gamma  <- 
    # parameters$beta_gamma   <- 
    # parameters$alpha_omega  <- 
    # parameters$beta_omega   <- 
    
    parameters$beta_upsilon  <- parameters$beta_upsilon  + learningrate  * gradient_beta_upsilon
    parameters$alpha_upsilon <- parameters$alpha_upsilon + learningrate  * gradient_alpha_upsilon
    print(gradient_beta_upsilon)
    # parameters$alpha_epsilon <- 
    # parameters$beta_epsilon  <- 

    ###############################################################################################
    ###############################################################################################
    # ELBO, if set
    # if (iter%%5== 0) { 
    # if (parameters$progress == 1) {

    lb <- ELBO()
    
      # bounds[iter] <- lb
    # }}  # two condition for lower bound calculation
  #} # closes iter loop
  
  ## outputs ####
  # if (parameters$progress == 1) {
  #   # includes bounds just calculated
  #   state <- list(lambda = lambda, upsilon = upsilon, a = a, gamma = gamma, omega = omega, epsilon = epsilon, 
  #                 be = be, bounds = bounds, parameters = parameters)
  # }
  # else {
  #   # no bounds
  #   state <- list(lambda = lambda, upsilon = upsilon, a = a, gamma = gamma, omega = omega, epsilon = epsilon, 
  #                 be = be, parameters = parameters)
  # }
# } # end of function
