lb <- ELBO_for_gradient_calc()   # current ELBO
h <- 10^(-7)   # set stepsize for quotient differential
#########################################################################################
# gradient lambda shape/alpha hyperparameter: 415.577
parameters$alpha_lambda <- parameters$alpha_lambda + h
lb_h <- ELBO_for_gradient_calc()
parameters$alpha_lambda <- parameters$alpha_lambda - h
gradient_lambda_alpha <- (lb_h- lb)/ h
cat("approx. gradient lambda shape/alpha hyperparameter:", gradient_lambda_alpha, "\n")
# gradient lambda scale/theta hyperparameter: 275.742 
# approx. gradient lambda shape/alpha hyperparameter: 341.5846 
parameters$beta_lambda <- parameters$beta_lambda + h
lb_h <- ELBO_for_gradient_calc()
parameters$beta_lambda <- parameters$beta_lambda - h
gradient_lambda_beta <- (lb_h- lb)/ h
cat("approx. gradient lambda scale/beta hyperparameter:", gradient_lambda_beta ,"\n")
# approx. gradient lambda scale/beta hyperparameter: 209.8016
#####################################################################
# gradient gamma shape/alpha hyperparameter: 29.53556
parameters$alpha_gamma <- parameters$alpha_gamma + h
lb_h <- ELBO_for_gradient_calc()
parameters$alpha_gamma <- parameters$alpha_gamma - h
gradient_gamma_alpha <- (lb_h- lb)/ h
cat("approx. gradient gamma shape/alpha hyperparameter:", gradient_gamma_alpha, "\n")
# gradient gamma scale/theta: 17.97984
parameters$beta_gamma <- parameters$beta_gamma + h
lb_h <- ELBO_for_gradient_calc()
parameters$beta_gamma <- parameters$beta_gamma - h
gradient_gamma_beta <- (lb_h- lb)/ h
cat("approx. gradient gamma scale/beta hyperparameter:", gradient_gamma_beta, "\n")
####################################################################################
# gradient omega shape/alpha hyperparameter: 29.53556
parameters$alpha_omega <- parameters$alpha_omega + h
lb_h <- ELBO_for_gradient_calc()
parameters$alpha_omega <- parameters$alpha_omega - h
gradient_omega_alpha <- (lb_h- lb)/ h
cat("approx. gradient omega shape/alpha hyperparameter:", gradient_omega_alpha, "\n")
# gradient omega scale/theta: 17.97984
parameters$beta_omega <- parameters$beta_omega + h
lb_h <- ELBO_for_gradient_calc()
parameters$beta_omega <- parameters$beta_omega - h
gradient_gamma_beta <- (lb_h- lb)/ h
cat("approx. gradient omega scale/beta hyperparameter:", gradient_gamma_beta, "\n")
####################################################################################
# gradient upsilon shape/alpha hyperparameter: 29.53556
parameters$alpha_upsilon <- parameters$alpha_upsilon + h
lb_h <- ELBO_for_gradient_calc()
parameters$alpha_upsilon <- parameters$alpha_upsilon - h
gradient_upsilon_alpha <- (lb_h- lb)/ h
cat("approx. gradient upsilon shape/alpha hyperparameter:", gradient_upsilon_alpha, "\n")
# gradient upsilon scale/theta: 17.97984
parameters$beta_upsilon <- parameters$beta_upsilon + h
lb_h <- ELBO_for_gradient_calc()
parameters$beta_upsilon <- parameters$beta_upsilon - h
gradient_upsilon_beta <- (lb_h- lb)/ h
cat("approx. gradient upsilon scale/beta hyperparameter:", gradient_upsilon_beta, "\n")
####################################################################################
# gradient epsilon shape/alpha hyperparameter: 29.53556
parameters$alpha_epsilon <- parameters$alpha_epsilon + h
lb_h <- ELBO_for_gradient_calc()
parameters$alpha_epsilon <- parameters$alpha_epsilon - h
gradient_epsilon_alpha <- (lb_h- lb)/ h
cat("approx. gradient epsilon shape/alpha hyperparameter:", gradient_epsilon_alpha, "\n")
# gradient epsilon scale/theta: 17.97984
parameters$beta_epsilon <- parameters$beta_epsilon + h
lb_h <- ELBO_for_gradient_calc()
parameters$beta_epsilon <- parameters$beta_epsilon - h
gradient_epsilon_beta <- (lb_h- lb)/ h
cat("approx. gradient epsilon scale/beta hyperparameter:", gradient_epsilon_beta, "\n")
####################################################################################