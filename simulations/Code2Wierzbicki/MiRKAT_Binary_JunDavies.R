# n = 100
# y = (runif(n) < 0.5)^2
# Z = matrix((runif(n*5) < 0.6)^2, n, 5)
# K = Z %*% t(Z)
# K2 = (Z %*% t(Z) +1)^2
# Ks = list(K, K2)
# 
# MiRKAT_omniKernel_logistic_davies(y,X = NULL, Ks)
# MiRKAT_omniKernel_logistic_davies(y,X = NULL, Ks = list(K))
MiRKAT_singleKernel_logistic_davies = function(y, X = NULL, K, family= "binomial"){
  library(CompQuadForm)
  family= "binomial"
  # try(library(BiasedUrn))
  n = length(y)
  
  if (is.null(X)) {
    px   = 1
    X1  = matrix(1, nrow=n)
    mod = glm(y~1 ,family = family)
    mu    = rep(mean(y), n)
    res   = y - mu
    
  } else {
    X1  = cbind(1, X)
    px  = NCOL(X1)
    mod = glm(y~X, family = family)
    mu  = mod$fitted.values
    res = y - mu
  }
  w   = mu*(1-mu)
  D0  =  sqrt(w)
  
  
  DX12 <- D0 * X1
  P0 <- diag(n) - DX12 %*% solve(t(DX12) %*% (DX12)) %*% t(DX12)
 #  PKP = P0 %*% (D0 %*% K %*% D0) %*% P0  
  PKP = P0 %*% (D0*t(D0 * K)) %*% P0
  Q <- as.numeric(res %*% K %*% res)

#   ee = eigen(PKP - Q * P0/n, symmetric = T)  
#   lambda0 = ee$values[abs(ee$values) >= 1e-10]
#   lambda0 = eigen(PKP)$values -Q*eigen(P0)$values/n
#   lambda0 = lambda0[lambda01 > 1e-10] 
#   davies(0, lambda=lambda01)$Qq
  
  eP0 = c(rep(1, n-px), rep(0, px))
  lambda0 = eigen(PKP, symmetric = T)$values - Q*eP0/n
  lambda0 = lambda0[abs(lambda0) >= 1e-10]
  davies(0, lambda=lambda0,acc=0.01)$Qq
}

# Now consider the omnibusTest


MiRKAT_omniKernel_logistic_davies = function(y, X = NULL, Ks, family= "binomial", n.sim = 1000){
  # try(library(BiasedUrn, lib ="/netscr/nzhao/Dissertation/Rlib"))
  try(library(BiasedUrn))
  library(CompQuadForm)
  n = length(y)
  family= "binomial"
  if (is.null(X)) {
    px   = 1
    X1  = matrix(1, nrow=n)
    mod = glm(y~1 ,family = family)
    mu    = rep(mean(y), n)
    res   = y - mu
    
  } else {
    X1  = cbind(1, X)
    px  = NCOL(X1)
    mod = glm(y~X, family = family)
    mu  = mod$fitted.values
    res = y - mu
  }
  w   = mu*(1-mu)
  D0  =  sqrt(w)  
  DX12 <- D0 * X1
  P0 <- diag(n) - DX12 %*% solve(t(DX12) %*% (DX12)) %*% t(DX12)
   
  getIndivP = function(K, res, D0, px){
    n = length(res)
    Q <- as.numeric(res %*% K %*% res) 
    PKP = P0 %*% (D0*t(D0 * K)) %*% P0 # For different p, this is different. 
    eP0 = c(rep(1, n-px), rep(0, px))
    ePKP = eigen(PKP, symmetric = T)$values
    lambda0 = ePKP - Q*eP0/n
    
    lambda0 = lambda0[abs(lambda0) >= 1e-10]
    p = davies(0, lambda=lambda0)$Qq
    return(list(Q = Q, ePKP = ePKP, p = p ))
  }
 
  S = sapply(Ks, getIndivP, res,  D0, px)
 
  odds    = exp(mod$linear.predictors)
  n.case  = sum(y)
  m1      = c(rep(1, length(y)))
  y.sim   = rMFNCHypergeo(n.sim, m1, n.case, odds)
  # Obtaining residuals
  if (is.null(X)){ # no need to rerun the logistic regression, thus faster.
    mu.sim = mean(y)
  }else{
    mu.sim  = sapply(1:n.sim, function(i){
      mu  = glm(y.sim[,i]~X, family = family)$fitted.values
    })
  }
  res.sim = y.sim - mu.sim
  Q.sim   = sapply(1:length(Ks), function(j){
    # j represents the kernels and i represents the simulated datas
    sapply(1:(n.sim), function(i){res.sim[,i] %*% Ks[[j]] %*% res.sim[,i]})
  })
  
  p_sim = sapply(1:length(Ks), function(j){
    ePKP = S[,j]$ePKP; eP0 = c(rep(1, n-px), rep(0, px))
    sapply(1:n.sim, function(i){
      lambda0 = ePKP - Q.sim[i,j]*eP0/n
      lambda0 = lambda0[abs(lambda0) >= 1e-10]
      davies(0, lambda=lambda0)$Qq
    })
  })
  #    p_sim1 = matrix(NA, n.sim, length(Ks))   
  # for (j in 1:length(Ks)){
  #   ePKP = S[,j]$ePKP; eP0 = c(rep(1, n-px), rep(0, px))
  #   p_sim1[,j] = sapply(1:n.sim, function(i){
  #     lambda0 = ePKP - Q.sim[i,j]*eP0/n
  #     lambda0 = lambda0[abs(lambda0) >= 1e-10]
  #     davies(0, lambda=lambda0)$Qq
  #   })}
  P_perm = apply(p_sim ,1, min)
  minP   = min(unlist(S[3,]))
  finalP= mean(P_perm <= minP)
  return(list(indP = unlist(S[3,]), obj = S, finalP = finalP)) 

}

