D2K     = function(D){
    n = NROW(D)
    centerM = diag(n) - 1/n
    K = -0.5*centerM %*% (D*D) %*% centerM
    eK= eigen(K, symmetric = T)
    # K = eK$vector %*% diag(abs(eK$values)) %*% t(eK$vector)
    K = eK$vector  %*% (t(eK$vector) * abs(eK$values))
    return(K)
}


swak<- function(Y, X,cov,tree, B=1000,pow=c(0.5,1,1.5,2),out_type = "D"){
    
    n <- length(Y)
    if (is.null(X) && length(X)>0) X=as.matrix(X, ncol=1)
    k <- ncol(X)
    
    if (class(X)!= "matrix") stop("X is not a matrix")
    if (nrow(X)!=n) stop("Dimensions of y and X do not match")
    
    m = ncol(X)
    if (out_type == "D") {
        family = "binomial"
        if (is.null(cov)) {
            mu    = rep(mean(Y), n)
            residual   = Y - mu
        } else {
            mod = glm(Y~cov, family = family)
            mu  = mod$fitted.values
            residual = Y - mu
        }
    } else {
        if (is.null(cov)) {
            mu    = rep(mean(Y), n)
            residual   = Y - mu
        } else {
            mod  = lm(Y~cov)
            mu  = mod$fitted.values
            residual = Y - mu
        }
    }
    
    ##observed statistics
    out.i = 1
    GuniF.cum = GUniFrac_cum(X,tree)
    cum = GuniF.cum$cum
    
    
    br.len = GuniF.cum$br.len
    br.len = as.matrix(br.len)
    U = residual %*% t(cum)
    U = abs(U)
    
    score.weight = U^0
    GUniF = GUniFracCpp(cum,br.len,score.weight)
    
    tmp = GUniFrac(X, tree, alpha = c(0))$unifracs
    duw = tmp[, , "d_UW"]
    d0 = GUniF$d0
    dw = GUniF$d1  # Weighted UniFrac
    d5 = GUniF$d5    # GUniFrac with alpha 0.5
    dbc = BCdist(cum,t(score.weight))

    ## transform the distance matrix to the kernel
    duw = D2K(duw)
    d0 = D2K(d0)
    dw = D2K(dw)
    d5 = D2K(d5)
    dbc = D2K(dbc)

    ks = list(dbc,duw,d0,d5,dw)
    
    output = matrix(NA,9,10)
    
    if(out_type  == "C"){
        
        res = MiRKAT(y = Y,X = cov,Ks = ks,out_type = "C",nperm = B,method= "permutation")
        output[1,1:6] = c(res$indivP,res$omnibus_p)
        
        res = MiRKAT(y = Y,X = cov,Ks = ks,out_type = "C",nperm = B,method= "davies")
        output[2,1:6] = c(res$indivP,res$omnibus_p)
        
        cum2 = cum *matrix(rep(br.len,each = dim(cum)[2]),nrow = dim(cum)[1],ncol = dim(cum)[2],byrow = TRUE)
        aSPU.res <- aSPU(Y, t(cum2), cov = cov, resample = "perm",model = "gaussian", pow = c(2:8, Inf), n.perm = B)
        output[3,1:9] = aSPU.res$pvs
        aSPU.res <- aSPU(Y, X, cov = cov, resample = "perm",model = "gaussian", pow = c(1:8, Inf), n.perm = 1000)
        output[4,1:10] = aSPU.res$pvs
        
        aSPU.res <- aSPU(Y, t(cum), cov = cov, resample = "perm",model = "gaussian", pow = c(1:8, Inf), n.perm = B)
        output[5,1:10] = aSPU.res$pvs
        
        tmp.cum = cum
        tmp.cum[tmp.cum != 0] = 1
        cum3 = tmp.cum *matrix(rep(br.len,each = dim(cum)[2]),nrow = dim(cum)[1],ncol = dim(cum)[2],byrow = TRUE)
        aSPU.res <- aSPU(Y, t(cum3), cov = cov, resample = "perm",model = "gaussian", pow = c(2:8, Inf), n.perm = B)
        output[6,1:9] = aSPU.res$pvs
        
        
        MiSPU.res <- MiSPUR(Y,t(cum2),t(cum3), cov = cov, model = "gaussian",pow = c(2:8, Inf), n.perm = B)
        
        output[7,1:9] = MiSPU.res$Unweighted$pvs
        
        output[8,1:9] = MiSPU.res$Weighted$pvs
        output[9,1] = MiSPU.res$Final$pvalue
        
        out = list(output1= output)
    }
    
    if(out_type  == "D"){
        res = MiRKAT(y = Y,X = cov,Ks = ks,out_type = "D",nperm = B,method= "permutation")
        res2 = MiRKAT(y = Y,X = cov,Ks = ks,out_type = "D",nperm = B, method = "davies")
        output[1,1:6] = c(res$indivP,res$omnibus_p)
        output[2,1:6] = c(res2$indivP,res2$omnibus_p)
        
        cum2 = cum *matrix(rep(br.len,each = dim(cum)[2]),nrow = dim(cum)[1],ncol = dim(cum)[2],byrow = TRUE)
        aSPU.res <- aSPU(Y, t(cum2), cov = cov, resample = "perm",model = "binomial", pow = c(1:8, Inf), n.perm = B)
        output[3,1:10] = aSPU.res$pvs
        
        aSPU.res <- aSPU(Y, X, cov = cov, resample = "perm",model = "binomial", pow = c(1:8, Inf), n.perm = B)
        output[4,1:10] = aSPU.res$pvs
        
        aSPU.res <- aSPU(Y, t(cum), cov = cov, resample = "perm",model = "binomial", pow = c(1:8, Inf), n.perm = B)
        output[5,1:10] = aSPU.res$pvs
        
        tmp.cum = cum
        tmp.cum[tmp.cum != 0] = 1
        cum3 = tmp.cum *matrix(rep(br.len,each = dim(cum)[2]),nrow = dim(cum)[1],ncol = dim(cum)[2],byrow = TRUE)
        aSPU.res <- aSPU(Y, t(cum3), cov = cov, resample = "perm",model = "binomial", pow = c(2:8, Inf), n.perm = B)
        output[6,1:9] = aSPU.res$pvs
        
        MiSPU.res <- MiSPUR(Y,t(cum2),t(cum3), cov = cov, model = "binomial",pow = c(1:8, Inf), n.perm = B)
        
        output[7,1:10] = MiSPU.res$Unweighted$pvs
        
        output[8,1:10] = MiSPU.res$Weighted$pvs
        output[9,1] = MiSPU.res$Final$pvalue
        
        out = list(output1= output)
    }
    
    
    return(out)
}






