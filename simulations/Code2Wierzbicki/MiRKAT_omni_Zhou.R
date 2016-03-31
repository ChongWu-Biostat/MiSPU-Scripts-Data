

MiRKAT_omni_resid_Zhou = function(y, Ks, X = NULL, nperm = 1000) {
  try(library(CompQuadForm))
  try(library(CompQuadForm, lib.loc = "~/Rlib"))
  
  #Helper function
  permuted.index = function (n){
    out = sample.int(n, n)
    return(out)
  }
  getQ = function(K, res, s2){    
    Q = 1/s2*res %*% K %*% res
  }
  getLambda = function(K, P0){
    PKP = P0 %*% K %*% P0
    ee = eigen(PKP, symmetric = T)         
    lambda0 = ee$values[ee$values > 1e-10]
    return(lambda0)    
  }
  
  getindivP = function(Q, lambda0, n, px){
    
    if (length(lambda0) >= n-px){ 
      # In fact, it can not go bigger than n-p because P0 has rank n-p
      lambda = c(lambda0 - Q/(n-px))
      k = length(lambda)
      h = rep(1,k)
    }else{
      lambda = c(lambda0 - Q/(n-px), -Q/(n-px))
      k = length(lambda0)
      h = c(rep(1, k), n-px - k)
    }
    
    p_davies = davies(q = 0, lambda = lambda, h = h, acc = 0.00005)$Qq
    # p_davies = davies(q = Q,lambda =lambda0, acc = 0.00005)$Qq
    p_davies = ifelse(p_davies < 0, 0, p_davies) 
    return(p_davies)  
  }
  ############
  n = length(y)
  if (is.null(X)) {
    p    = 1
    X1   = matrix(1, nrow=n)
    mod  = lm(y~1)
  } else {
    p    = ncol(X)
    X1   = cbind(1, X)
    mod  = lm(y~X)
  }
  s2  = summary(mod)$s**2
  D0  = diag(n)#diag(mu*(1-mu)) for logistic
  res = resid(mod)
  P0  = D0 - X1%*%solve(t(X1)%*%X1)%*%t(X1)
  px  = NCOL(X1)
  # P-VALUES 
  Qs = lapply(Ks, getQ, res, s2)
  lambda0 = lapply(Ks, getLambda, P0)
  ps = rep(NA, length(Ks))
  for (i in 1:length(Ks)){
    ps[i] = getindivP(Qs[[i]], lambda0[[i]], n, px)    
  }
  minP    = min(ps)
  # permutations for the final result
  perm    = sapply(1:nperm, function(x) permuted.index(n))
  y_star  = mod$fitted + matrix(res[perm],n,nperm)
  
  res_sim = qr.resid(qr(X1), y_star)     
  modelVar= function(x, px){sum((x - mean(x))^2)/(n - px)}  
  
  sigma2_sim= apply(res_sim, 2, modelVar, px) # Already the variance
  
  Q_sim   = sapply(1:length(Ks), function(j){
    sapply(1:nperm, function(i){
      res_sim[,i] %*% Ks[[j]] %*%res_sim[,i]/sigma2_sim[i]})
  })

  p_sim = sapply(1:length(Ks), function(j){
   sapply(1:nperm, function(i){
     getindivP(Q_sim[i,j], lambda0[[j]], n, px)   
   })
 })

  minP_sim= apply(p_sim,1, min)
  p_final = mean(minP >= minP_sim) # The first one is from the data and all others are simualted ones
  
  return(list(indivP = ps , omnibus_p = p_final))
}


