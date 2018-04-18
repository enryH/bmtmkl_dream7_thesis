### main function to get derivative of prior and entropy w.r.t to hyperparameter of random variable of interest
#
# for alpha:
GradientAlpha <- function(rv.gamma.shape, rv.gamma.scale, hyperalpha, hyperbeta, cardinality=1) {
  gradient <- sum(- digamma(hyperalpha) 
                  - log(hyperbeta) 
                  - (1 - hyperalpha) * trigamma(rv.gamma.shape)
                  + (digamma(rv.gamma.shape))
                  + (log(rv.gamma.scale)) 
                  - (rv.gamma.scale)/ hyperbeta
                  - (rv.gamma.shape - 1) * trigamma(rv.gamma.shape)
                  + 1
  )
  return(gradient)
}
GradientAlpha(upsilon$shape, upsilon$scale, parameters$alpha_upsilon, parameters$beta_upsilon, T)
#
###
# for theta which is the scale - (called beta for hyperparameter (parameters$beta_?)) 
# GradientTheta <-  function(rv.gamma.shape, rv.gamma.scale, hyperalpha, hyperbeta) {
#   gradient <- sum(  - hyperalpha / hyperbeta 
#                     + hyperalpha / hyperbeta^2 * rv.gamma.scale 
#                     + 1/ hyperbeta^2 * rv.gamma.shape * rv.gamma.scale^2 
#                     - 1/ hyperbeta^3 * rv.gamma.shape * rv.gamma.scale^2
#   )
#   return(gradient)
# }
GradientTheta <-  function(rv.gamma.shape, rv.gamma.scale, hyperalpha, hyperbeta) {
  gradient <- sum(  - (hyperalpha / hyperbeta) 
                    + (hyperalpha / hyperbeta^2) * rv.gamma.scale 
                    + (1/ hyperbeta^2) * rv.gamma.shape * rv.gamma.scale
                    - (1/ hyperbeta^3) * rv.gamma.shape * rv.gamma.scale^2
  )
  return(gradient)
}
GradientTheta(upsilon$shape, upsilon$scale, parameters$alpha_upsilon, parameters$beta_upsilon)
###############################################################################################
###############################################################################################    
### for testing:  ####
#' rv.gamma.shape <- upsilon$shape
#' rv.gamma.scale <- upsilon$scale
#' hyperalpha     <- parameters$alpha_upsilon
#' hyperbeta      <- parameters$beta_upsilon
#' C              <- T
#' # E[log(prior)] : derivative w.r.t to alpha/shape hyperparameter alpha_rv (rv: random variable)
#' ExpectedLogPriorGradientAlpha <-
#'   function(rv.gamma.shape,
#'            rv.gamma.scale,
#'            hyperalpha,
#'            hyperbeta,
#'            cardinality) {
#'     #' provide shape and scale parameter of gamma dist. random variables as vector
#'     #' hyperalpha and hyperbeta are previously set hyperparameters
#'     C <- cardinality
#'     gradient <-
#'       (
#'         -C * (digamma(hyperalpha) + log(hyperbeta))
#'         + (hyperalpha - 1) * sum(trigamma(rv.gamma.shape))
#'         + sum(digamma(rv.gamma.shape)) + sum(log((rv.gamma.scale)))
#'         - 1 / hyperbeta * sum((rv.gamma.scale))
#'       )
#'     return(gradient)
#'   }
#' #
#' print(
#'   ExpectedLogPriorGradientAlpha(
#'     upsilon$shape,
#'     upsilon$scale,
#'     parameters$alpha_upsilon,
#'     parameters$beta_upsilon,
#'     T
#'   )
#' )
#' # E[log(q)] . derivative w.r.t to alpha/shape hyperparameter alpha_rv (rv: random variable)
#' EntropyGradientAlpha <-
#'   function(rv.gamma.shape,
#'            rv.gamma.scale,
#'            hyperalpha,
#'            hyperbeta,
#'            cardinality) {
#'     gradient <- sum((1 - rv.gamma.shape) * trigamma(rv.gamma.shape) + 1)
#'     # gradient <- sum( (1- rv.gamma.shape) * trigamma(rv.gamma.shape)) + C
#'     return(gradient)
#'   }
#' #
#' print(
#'   EntropyGradientAlpha(
#'     upsilon$shape,
#'     upsilon$scale,
#'     parameters$alpha_upsilon,
#'     parameters$beta_upsilon,
#'     T
#'   )
#' )
#' # E[log(prior)] : derivative w.r.t to theta/scale hyperparameter theta_rv (rv: random variable)
#' ExpectedLogPriorGradientTheta <-
#'   function(rv.gamma.shape,
#'            rv.gamma.scale,
#'            hyperalpha,
#'            hyperbeta,
#'            cardinality) {
#'     C <- cardinality
#'     return(
#'       -C * hyperalpha / hyperbeta
#'       + (hyperalpha - 1) / hyperbeta ^ 2 * sum(rv.gamma.scale)
#'       + 1 / hyperbeta ^ 2 * sum(rv.gamma.shape * rv.gamma.scale)
#'       - 1 / hyperbeta ^ 3 * sum(rv.gamma.shape * rv.gamma.scale ^
#'                                   2)
#'     )
#'   }
#' # return(gradient)
#' print(
#'   ExpectedLogPriorGradientTheta(
#'     upsilon$shape,
#'     upsilon$scale,
#'     parameters$alpha_upsilon,
#'     parameters$beta_upsilon,
#'     T
#'   )
#' )
#' #
#' EntropyGradientThetaScale <-
#'   function(rv.gamma.shape,
#'            rv.gamma.scale,
#'            hyperalpha,
#'            hyperbeta,
#'            cardinality) {
#'     gradient <- 1 / hyperbeta ^ 2 * sum(rv.gamma.scale)
#'     return(gradient)
#'   }
#' #
#' print(
#'   EntropyGradientThetaScale(
#'     upsilon$shape,
#'     upsilon$scale,
#'     parameters$alpha_upsilon,
#'     parameters$beta_upsilon,
#'     T
#'   )
#' )
#' #
#' GradientAlpha <-
#'   function(rv.gamma.shape,
#'            rv.gamma.scale,
#'            hyperalpha,
#'            hyperbeta,
#'            cardinality) {
#'     C <- cardinality
#'     gradient <- sum(
#'       -digamma(hyperalpha)
#'       + (digamma(rv.gamma.shape))
#'       # + sum(log(rv.gamma.scale/hyperbeta)) #+ C * log(hyperbeta) +(log(rv.gamma.scale)) - log(hyperbeta)
#'       - 1 / hyperbeta * (rv.gamma.scale)
#'       - (rv.gamma.shape - hyperalpha) * trigamma(rv.gamma.shape)
#'       + 1
#'     )
#'     return(gradient)
#'   }
#
# GradientAlpha(upsilon$shape, upsilon$scale, parameters$alpha_upsilon, parameters$beta_upsilon, T)
# print(   ExpectedLogPriorGradientAlpha(upsilon$shape, upsilon$scale, parameters$alpha_upsilon, parameters$beta_upsilon, T)
#          +          EntropyGradientAlpha(upsilon$shape, upsilon$scale, parameters$alpha_upsilon, parameters$beta_upsilon, T)  )
# print(       EntropyGradientThetaScale(upsilon$shape, upsilon$scale, parameters$alpha_upsilon, parameters$beta_upsilon, T)
#              + ExpectedLogPriorGradientTheta(upsilon$shape, upsilon$scale, parameters$alpha_upsilon, parameters$beta_upsilon, T) )
##############################################
### gradients for Lambda
gradient_lambda_alpha     <-  0
for (o in 1:T) {
  # a_ti_sq <- diag(atimesaT.mean[[o]])  ## a_{t,i}^2 as a vector for drug t of length N_t
  # C_at    <- tcrossprod(KmtimesGT.mean[[o]], KmtimesGT.mean[[o]]) * upsilon$shape[o] * upsilon$scale[o]
  gradient_lambda_alpha <- (gradient_lambda_alpha 
                            + GradientAlpha(lambda[[o]]$shape, lambda[[o]]$scale, parameters$alpha_lambda, parameters$beta_lambda)
  )
}
cat("gradient lambda shape/alpha hyperparameter:",gradient_lambda_alpha,"\n")
# gradients lambda beta (in scale notation)
gradient_lambda_beta     <-  0
for (o in 1:T) {
  # a_ti_sq <- diag(atimesaT.mean[[o]])
  # C_at    <- tcrossprod(KmtimesGT.mean[[o]], KmtimesGT.mean[[o]]) * upsilon$shape[o] * upsilon$scale[o]
  gradient_lambda_beta <- ( gradient_lambda_beta 
                            + GradientTheta(lambda[[o]]$shape, lambda[[o]]$scale, parameters$alpha_lambda, parameters$beta_lambda)
  )
}
cat("gradient lambda scale/theta hyperparameter:", gradient_lambda_beta,"\n")
#
gradient_gamma_alpha <- ( GradientAlpha(gamma$shape, gamma$scale, parameters$alpha_gamma, parameters$beta_gamma)
                            # + sum(0.5 * trigamma(gamma$shape))
                            # - 0.5 * sum(diag(diag(as.vector(gamma$shape * gamma$scale), T, T) %*% btimesbT.mean))
)           
cat("gradient gamma shape/alpha hyperparameter:",gradient_gamma_alpha,"\n")
gradient_gamma_beta  <- (GradientTheta(gamma$shape, gamma$scale, parameters$alpha_gamma, parameters$beta_gamma)
                         # + sum( 0.5 * gamma$scale / parameters$beta_gamma^2
                         # - 0.5 *  diag(diag(as.vector(gamma$shape * (gamma$scale)^2/(parameters$beta_gamma)^2))  %*% btimesbT.mean) )
)
cat("gradient gamma scale/theta hyperparameter:", gradient_gamma_beta,"\n")
gradient_omega_alpha <- (GradientAlpha(omega$shape, omega$scale, parameters$alpha_omega, parameters$beta_omega)
                         # +sum( 0.5* trigamma(omega$shape)
                         # - 0.5 * diag(diag(as.vector(omega$shape * omega$scale), P, P) %*% etimeseT.mean)  )
)
cat("gradient omega shape/alpha hyperparameter:",gradient_omega_alpha,"\n")
gradient_omega_beta  <- (GradientTheta(omega$shape, omega$scale, parameters$alpha_omega, parameters$beta_omega)
                         # + sum(0.5* omega$scale / parameters$beta_omega^2
                         # - 0.5 * diag(diag(as.vector(omega$shape * omega$scale^2/ parameters$beta_omega^2), P, P) %*% etimeseT.mean) )
)
cat("gradient omega scale/theta hyperparameter:", gradient_omega_beta,"\n")
#
# define gradients for upsilon ####
# calculate gradients:
# ExpectedLogPrior("upsilon")
# Gt.constant <- matrix(0,T,1)
  # for (o in 1:T){     
  #   Gt.constant[o,1] <- ( sum(diag(GtimesGT.mean[[o]])) 
  #     - 2 * crossprod(a[[o]]$mean, KmtimesGT.mean[[o]])
  #     + sum(diag(KmKm[[o]] %*% atimesaT.mean[[o]]))
  #   )
  # }
gradient_upsilon_alpha <- GradientAlpha(upsilon$shape, upsilon$scale, parameters$alpha_upsilon, parameters$beta_upsilon)
#   sum(
#   # log p(upsilon)
#   - digamma(parameters$alpha_upsilon)
#   + digamma(upsilon$shape)
#   + (parameters$alpha_upsilon - 1) * trigamma(upsilon$shape)
#   + log(parameters$beta_upsilon)
#   + log(upsilon$scale)
#   - 1/ parameters$beta_upsilon * upsilon$scale
#   # entropy: log q(upsilon)
#   - (upsilon$shape -1 ) * trigamma(upsilon$shape)
#   + 1
#   # log p(G) contribution:
#   # +  0.5 * P *N* trigamma(upsilon$shape)
#   # -  0.5 *  upsilon$scale  * Gt.constant
# ) 
cat("gradient upsilon shape/alpha hyperparameter:", gradient_upsilon_alpha,"\n")
#
gradient_upsilon_beta  <-  GradientTheta(upsilon$shape, upsilon$scale, parameters$alpha_upsilon, parameters$beta_upsilon) # log p(upsilson)  # log q(upsilon)
  # log p(G)  # ist P*N correct? todo? yes -> N ist vector containing N_1,...,N_t...
  # + sum( 0.5 * P * N * upsilon$scale / parameters$beta_upsilon^2
  #      - 0.5 * upsilon$shape * upsilon$scale^2/parameters$beta_upsilon^2 * Gt.constant )
# )
cat("gradient upsilon scale/theta:", gradient_upsilon_beta,"\n")
# }
#
# yt.constant <- matrix(0, T, 1)
# for (o in 1:T) {
#   yt.constant[o,1]<- ( -0.5 * crossprod(y[[o]], y[[o]])
#   + crossprod(y[[o]], crossprod(G[[o]]$mean, be$mean[(T + 1):(T + P)]))
#   + sum(be$mean[o] * y[[o]])
#   - 0.5 * sum(diag(etimeseT.mean %*% GtimesGT.mean[[o]]))
#   - sum(crossprod(G[[o]]$mean, etimesb.mean[, o]))
#   - 0.5 * N[o] * btimesbT.mean[o, o]
#   )
# }
gradient_epsilon_alpha   <- (GradientAlpha(epsilon$shape, epsilon$scale, parameters$alpha_epsilon, parameters$beta_epsilon)
                            # + sum( 0.5 * P *N* trigamma(epsilon$shape)
                            # -  0.5 *  epsilon$scale  * yt.constant)
                            )
cat("gradient epsilon shape/alpha hyperparameter:", gradient_epsilon_alpha,"\n")
gradient_epsilon_beta    <- (GradientTheta(epsilon$shape, epsilon$scale, parameters$alpha_epsilon, parameters$beta_epsilon)
                           # + sum( 0.5 * P * N * epsilon$scale / parameters$beta_epsilon^2
                           # - 0.5 * epsilon$shape * epsilon$scale^2/parameters$beta_epsilon^2 * yt.constant )
                           )
cat("gradient epsilon scale/theta:", gradient_epsilon_beta,"\n")