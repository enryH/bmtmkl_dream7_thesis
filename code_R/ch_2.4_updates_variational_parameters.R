##Optimization ###########################################################################################################
# :
### todo: Implement stopping rule (relative improvement smaller than 0.001)?
## how fast?
# for (iter in 1:parameters$iteration) {
# runteniterations <- function(){
for (iter in 1:10) {
  # update lambda
  for (t in 1:T) {
    lambda[[t]]$scale <- 1 / (1 / parameters$beta_lambda + 0.5 * diag(atimesaT.mean[[t]]))
  }
  # update upsilon
  for (t in 1:T) {
    upsilon$scale[t] <- 1 / (1 / parameters$beta_upsilon + 0.5 * (sum(diag(GtimesGT.mean[[t]])) - 2 * sum(matrix(crossprod(a[[t]]$mean, Km[[t]]), N[t], P) * t(G[[t]]$mean)) + sum(diag(KmKm[[t]] %*% atimesaT.mean[[t]]))))
  }
  # update a
  for (t in 1:T) {
    a[[t]]$covariance <- chol2inv(chol(diag(as.vector(lambda[[t]]$shape * lambda[[t]]$scale), D[t], D[t]) + upsilon$shape[t] * upsilon$scale[t] * KmKm[[t]]))
    a[[t]]$mean <- a[[t]]$covariance %*% (upsilon$shape[t] * upsilon$scale[t] * KmtimesGT.mean[[t]])
    atimesaT.mean[[t]] <- tcrossprod(a[[t]]$mean, a[[t]]$mean) + a[[t]]$covariance
  }
  # update G
  for (t in 1:T) {
    G[[t]]$covariance <- chol2inv(chol(diag(upsilon$shape[t] * upsilon$scale[t], P, P) + epsilon$shape[t] * epsilon$scale[t] * etimeseT.mean))
    G[[t]]$mean <- G[[t]]$covariance %*% (upsilon$shape[t] * upsilon$scale[t] * t(matrix(crossprod(a[[t]]$mean, Km[[t]]), N[t], P)) + epsilon$shape[t] * epsilon$scale[t] * (tcrossprod(be$mean[(T + 1):(T + P)], y[[t]]) - matrix(etimesb.mean[,t], P, N[t], byrow = FALSE)))
    GtimesGT.mean[[t]] <- tcrossprod(G[[t]]$mean, G[[t]]$mean) + N[t] * G[[t]]$covariance
    KmtimesGT.mean[[t]] <- Km[[t]] %*% matrix(t(G[[t]]$mean), N[t] * P, 1)
  }
  # update gamma
  gamma$scale<- 1 / (1 / parameters$beta_gamma + 0.5 * diag(btimesbT.mean))
  # update omega
  omega$scale<- 1 / (1 / parameters$beta_omega + 0.5 * diag(etimeseT.mean))
  # update epsilon
  for (t in 1:T) {
    epsilon$scale[t] <- 1 / (1 / parameters$beta_epsilon + 0.5 * as.double(crossprod(y[[t]], y[[t]]) - 2 * crossprod(y[[t]], crossprod(rbind(matrix(1, 1, N[t]), G[[t]]$mean), be$mean[c(t, (T + 1):(T + P))])) + N[t] * btimesbT.mean[t, t] + sum(diag(GtimesGT.mean[[t]] %*% etimeseT.mean)) + 2 * sum(diag(crossprod(rowSums(G[[t]]$mean), etimesb.mean[,t])))))
  }
  # update b and e
  be$covariance <- rbind(cbind(diag(as.vector(gamma$shape * gamma$scale), T, T) + diag(as.vector(N * epsilon$shape * epsilon$scale), T, T), matrix(0, T, P)), cbind(matrix(0, P, T), diag(as.vector(omega$shape * omega$scale), P, P)))
  for (t in 1:T) {
    be$covariance[(T + 1):(T + P), t] <- epsilon$shape[t] * epsilon$scale[t] * rowSums(G[[t]]$mean)
    be$covariance[t, (T + 1):(T + P)] <- epsilon$shape[t] * epsilon$scale[t] * t(rowSums(G[[t]]$mean))
    be$covariance[(T + 1):(T + P), (T + 1):(T + P)] <- be$covariance[(T + 1):(T + P), (T + 1):(T + P)] + epsilon$shape[t] * epsilon$scale[t] * GtimesGT.mean[[t]]
  }
  be$covariance <- chol2inv(chol(be$covariance))
  be$mean <- matrix(0, T + P, 1)
  for (t in 1:T) {
    be$mean[t] <- epsilon$shape[t] * epsilon$scale[t] * sum(y[[t]])
    be$mean[(T + 1):(T + P)] <- be$mean[(T + 1):(T + P)] + epsilon$shape[t] * epsilon$scale[t] * G[[t]]$mean %*% y[[t]]
  }
  be$mean <- be$covariance %*% be$mean
  btimesbT.mean <- tcrossprod(be$mean[1:T], be$mean[1:T]) + be$covariance[1:T, 1:T]
  etimeseT.mean <- tcrossprod(be$mean[(T + 1):(T + P)], be$mean[(T + 1):(T + P)]) + be$covariance[(T + 1):(T + P), (T + 1):(T + P)]
  for (t in 1:T) {
    etimesb.mean[,t] <- be$mean[(T + 1):(T + P)] * be$mean[t] + be$covariance[(T + 1):(T + P), t]
  }
}
#   # update lambda 
#   for (o in 1:T) {  
#     #extracts diagonal elements
#     lambda[[o]]$scale <- 1 / ((1 / parameters$beta_lambda) + 0.5 * diag(atimesaT.mean[[o]]))  # in 1x N_o
#   }
#   # update upsilon
#   for (o in 1:T) {
#     # scalar:
#     c_vt <- (sum(diag(GtimesGT.mean[[o]])) 
#              - 2 * sum(matrix(crossprod(a[[o]]$mean, Km[[o]]), N[o], P) * t(G[[o]]$mean)) 
#              + sum(diag(KmKm[[o]] %*% atimesaT.mean[[o]])))
#     upsilon$scale[o] <- 1 / (1 / parameters$beta_upsilon 
#                              + 0.5 * c_vt)
#   }
#   # update weight a
#   for (o in 1:T) {
#     #todo: what is KmKm-> check in detail!
#     a[[o]]$covariance <- chol2inv(chol(diag(as.vector(lambda[[o]]$shape * lambda[[o]]$scale), D[o], D[o]) 
#                                        + upsilon$shape[o] * upsilon$scale[o] * KmKm[[o]]))  ## KmKm matrix is multiplied with scalars
#     a[[o]]$mean <- a[[o]]$covariance %*% (upsilon$shape[o] * upsilon$scale[o] * KmtimesGT.mean[[o]])
#     atimesaT.mean[[o]] <- tcrossprod(a[[o]]$mean, a[[o]]$mean) + a[[o]]$covariance  # Erwartungswerte addiert? Varianzen addiert?
#   }
#   # update intermediate outputs G
#   for (o in 1:T) {
#     G[[o]]$covariance <- chol2inv(chol(diag(upsilon$shape[o] * upsilon$scale[o], P, P) 
#                                        + epsilon$shape[o] * epsilon$scale[o] * etimeseT.mean))
#     #G_t_mean has P x N_t observations... just all values next to each other stored...
#     G[[o]]$mean <- (G[[o]]$covariance %*% (upsilon$shape[o] * upsilon$scale[o] * t(matrix(crossprod(a[[o]]$mean, Km[[o]]), N[o], P)) 
#                                            + epsilon$shape[o] * epsilon$scale[o] * (tcrossprod(be$mean[(T + 1):(T + P)], y[[o]]) 
#                                                                                     - repmat(etimesb.mean[,o], 1, N[o]))) )  
#     #! dimension of kernels (N[o]) and tasks do not necessarily match 
#     #-> drugresponse has to be imputed as well.
#     GtimesGT.mean[[o]] <- tcrossprod(G[[o]]$mean, G[[o]]$mean) + N[o] * G[[o]]$covariance # PXP 
#     KmtimesGT.mean[[o]]<- Km[[o]] %*% matrix(t(G[[o]]$mean), N[o] * P, 1) # dimension? N_t*P x N_t*P
#   }
#   # update gamma
#   gamma$scale <- 1 / (1 / parameters$beta_gamma + 0.5 * diag(btimesbT.mean)) 
#   # update omega
#   omega$scale <- 1 / (1 / parameters$beta_omega + 0.5 * diag(etimeseT.mean))  # auf Diagonalelementen steht Erwartungswert von e_k Quadrat
#   # update epsilon
#   for (o in 1:T) {
#     epsilon$scale[o] <- 1 / (1 / parameters$beta_epsilon 
#                              + 0.5 * as.double(crossprod(y[[o]], y[[o]]) 
#                                                - 2 * crossprod(y[[o]], crossprod(rbind(matrix(1, 1, N[o]), G[[o]]$mean), be$mean[c(o, (T + 1):(T + P))])) 
#                                                + N[o] * btimesbT.mean[o, o] + sum(diag(GtimesGT.mean[[o]] %*% etimeseT.mean)) 
#                                                + 2 * sum(diag(crossprod(rowSums(G[[o]]$mean), etimesb.mean[,o])))))
#   }
#   # update b and e
#   be$covariance <- rbind(cbind(diag(as.vector(gamma$shape * gamma$scale), T, T) 
#                                + diag(as.vector(N * epsilon$shape * epsilon$scale), T, T), matrix(0, T, P)), 
#                          cbind(matrix(0, P, T), diag(as.vector(omega$shape * omega$scale), P, P)))
#   for (o in 1:T) {
#     be$covariance[(T + 1):(T + P), o] <- epsilon$shape[o] * epsilon$scale[o] * rowSums(G[[o]]$mean)      ##lower left
#     be$covariance[o, (T + 1):(T + P)] <- epsilon$shape[o] * epsilon$scale[o] * t(rowSums(G[[o]]$mean))   ## upper right
#     be$covariance[(T + 1):(T + P), (T + 1):(T + P)] <- be$covariance[(T + 1):(T + P), (T + 1):(T + P)] + epsilon$shape[o] * epsilon$scale[o] * GtimesGT.mean[[o]]
#   }
#   be$covariance <- chol2inv(chol(be$covariance))
#   be$mean <- matrix(0, T + P, 1)
#   for (o in 1:T) {
#     be$mean[o] <- epsilon$shape[o] * epsilon$scale[o] * sum(y[[o]])
#     be$mean[(T + 1):(T + P)] <- be$mean[(T + 1):(T + P)] + epsilon$shape[o] * epsilon$scale[o] * G[[o]]$mean %*% y[[o]]
#   }
#   be$mean <- be$covariance %*% be$mean
#   btimesbT.mean <- tcrossprod(be$mean[1:T], be$mean[1:T]) + be$covariance[1:T, 1:T]
#   etimeseT.mean <- tcrossprod(be$mean[(T + 1):(T + P)], be$mean[(T + 1):(T + P)]) + be$covariance[(T + 1):(T + P), (T + 1):(T + P)]
#   for (o in 1:T) {
#     etimesb.mean[,o] <- be$mean[(T + 1):(T + P)] * be$mean[o] + be$covariance[(T + 1):(T + P), o]
#   }
# } 
# }