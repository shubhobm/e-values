# misc functions
## function to generate from multivariate normal
my.mvrnorm = function(n, mu, Sigma){
    p = length(mu)
    # compute square root of covariance matrix
    eo=eigen(Sigma, symmetric=TRUE)
    sigma.sqrt=eo$vec%*%diag(eo$val^0.5)%*%t(eo$vec)
    
    # generate random normals from runif by box-muller transform
    rnorm.vec = sqrt(-2*log(runif(n*p)))*cos(2*pi*runif(n*p))
    
    # generate sample matrix
    sample.matrix = matrix(rep(mu, n), nrow=n, byrow=T) +
        matrix(rnorm.vec, nrow=n, ncol=p)%*%sigma.sqrt
    return(sample.matrix)
}

ecmeml1 = function (y, subj, pred, xcol, zcol, vmax, occ, start, maxits = 1000, 
                    eps = 1e-04) 
{
    tmp <- table(subj)
    m <- length(tmp)
    nmax <- max(tmp)
    ntot <- length(y)
    pcol <- ncol(pred)
    q <- length(zcol)
    p <- length(xcol)
    {
        if (missing(vmax)) {
            vmax <- diag(rep(1, nmax))
            occ <- integer(ntot)
            iflag <- as.integer(1)
        }
        else iflag <- as.integer(0)
    }
    storage.mode(vmax) <- "double"
    {
        if (!missing(start)) {
            beta <- start$beta
            sigma2 <- start$sigma2
            xi <- start$psi/start$sigma2
            storage.mode(beta) <- "double"
            storage.mode(xi) <- "double"
            storage.mode(sigma2) <- "double"
            sflag <- as.integer(1)
        }
        else {
            beta <- numeric(p)
            xi <- matrix(0, q, q)
            sigma2 <- 0
            sflag <- as.integer(0)
        }
    }
    # if(trace==F){
    #   cat("Performing ECME...")
    # }
    now <- proc.time()
    err <- 0
    tmp <- .Fortran("ecmeml", ntot, as.integer(subj), m, ist = integer(m), 
                    ifin = integer(m), as.integer(occ), nmax, vmax, w = array(0, 
                                                                              c(nmax, nmax, m)), vinv = array(0, c(nmax, nmax, 
                                                                                                                   m)), pcol, as.double(pred), q, as.integer(zcol), 
                    ztvinv = array(0, c(q, nmax, m)), ztvinvz = array(0, 
                                                                      c(q, q, m)), iflag = iflag, err = as.integer(err), 
                    msg = integer(1), u = array(0, c(q, q, m)), iter = integer(1), 
                    sflag, sigma2 = sigma2, p, as.integer(xcol), beta = beta, 
                    as.double(y), delta = rep(0, ntot), xtw = matrix(0, p, 
                                                                     nmax), xtwx = matrix(0, p, p), xtwy = numeric(p), 
                    xtwxinv = matrix(0, p, p), wkqq1 = matrix(0, q, q), wkqq2 = matrix(0, 
                                                                                       q, q), xi = xi, wkqnm = array(0, c(q, nmax, m)), 
                    b = matrix(0, q, m), cvgd = integer(1), obeta = rep(0, 
                                                                        p), oxi = matrix(0, q, q), maxits = as.integer(maxits), 
                    llvec = numeric(as.integer(maxits)), eps = as.double(eps), 
                    PACKAGE = "lmm")
    clock <- proc.time() - now
    # if(trace==F){
    #   cat("\n")
    # }
    iter <- tmp$iter
    msg <- tmp$msg
    {
        if (msg == 1) 
            warning("Supplied V <- i matrix is not positive definite")
        else if (msg == 2) 
            warning("GLS failed for start vals, t(X)%*%inv(V)%*%X not full rank")
        else if (msg == 3) 
            warning("Inadequate information to obtain starting value of psi")
        else if (msg == 4) 
            warning("Value of psi became non-pos.def. during iterations")
        else if (msg == 5) 
            warning("t(X)%*%W%*%X became non-pos.def. during iterations")
    }
    llvec <- tmp$llvec[1:iter]
    converged <- tmp$cvgd == as.integer(1)
    cov.beta <- tmp$xtwxinv * tmp$sigma2
    b.hat <- tmp$b
    cov.b <- tmp$u * tmp$sigma2
    psi <- tmp$xi * tmp$sigma2
    beta <- tmp$beta
    if (!is.null(dimnames(pred)[[2]])) {
        colnames <- dimnames(pred)[[2]]
        names(beta) <- colnames[xcol]
        dimnames(psi) <- list(colnames[zcol], colnames[zcol])
    }
    list(beta = beta, sigma2 = tmp$sigma2, psi = psi, converged = converged, 
         iter = iter, loglik = llvec, cov.beta = cov.beta, b.hat = b.hat, 
         cov.b = cov.b)
}

get_outputs = function(beta, sig, D, ni, m, sdm, nboot, th, n.iter, ncores=1){
    
    ## set up data
    set.seed(1e3)
    n = ni*m
    p = length(beta)
    subj = rep(1:m, rep(ni,m))
    X = cbind(1,matrix(runif(n*p, -2, 2), ncol=p))
    Z = X[,1:pZ]
    W.true = matrix(0, ncol=n, nrow=n)
    for(im in 1:m){
        iind = ((im-1)*ni+1):(im*ni)
        Zi = Z[iind,]
        W.true[iind,iind] = sig * diag(ni) + Zi %*% D %*% t(Zi)
    }
    W.true = solve(W.true)
    
    ## depth model selection function
    loopfun = function(seed){
        ## Generate samples
        set.seed(1e3*seed)
        y = rep(0,n)
        for(im in 1:m){
            iind = ((im-1)*ni+1):(im*ni)
            y[iind] = my.mvrnorm(1, X[iind,] %*% c(0,beta), W.true[iind,iind])
        }
        
        ## full model constants
        lmmfull = ecmeml1(y=y, subj=subj, pred=X, xcol=1:(p+1), zcol=1:pZ)
        W = matrix(0, nrow=n, ncol=n)
        for(im in 1:m){
            iind = ((im-1)*ni+1):(im*ni)
            Zi = Z[iind,]
            W[iind,iind] = with(lmmfull, solve(sigma2 * diag(ni) + Zi %*% psi %*% t(Zi)))
        }
        beta.hat = lmmfull$beta
        sigma2 = lmmfull$sigma2
        H = lmmfull$cov.beta %*% t(X) %*% W
        r = y - X %*% beta.hat; r = r - mean(r)
        
        nsd = length(sdm)
        beta_hat_list = vector("list",nsd)
        b_hat_list = vector("list",nsd)
        for(ind in 1:nsd){
            # initialize
            beta_hat_list[[ind]] = rep(0,p)
            sd = sdm[ind]            
            
            # matrix of full model estimates: Fn is approx by Fn^b1
            beta.mat = matrix(0 ,nrow=nboot, ncol=p+1)
            SSPmat.d = matrix(0, nrow=nboot, ncol=p+1)
            for(i in 1:nboot){
                iresid = as.matrix(sd * rep(rgamma(m,1,1)-1, rep(ni,m)) * r, ncol=1)
                beta.mat[i,] = as.numeric(beta.hat) + as.numeric(H %*% iresid)
            }
            
            beta.mat1 = matrix(0 ,nrow=nboot, ncol=p+1)
            for(i in 1:nboot){
                iresid = as.matrix(sd * rep(rgamma(m,1,1)-1, rep(ni,m)) * r, ncol=1)
                beta.mat1[i,] = as.numeric(beta.hat) + as.numeric(H %*% iresid)
            }
            SSPmat.d[,p+1] = depth.Mahalanobis(beta.mat1[,-1], beta.mat[,-1])
            
            ## now marginal models
            for(j in 1:p){
                jbeta.mat1 = beta.mat1
                jbeta.mat1[,j+1] = 0
                SSPmat.d[,j] = depth.Mahalanobis(jbeta.mat1[,-1], beta.mat[,-1])
            }
            
            Cn.vec = apply(SSPmat.d, 2, mean)
            which.sel = sort(which(Cn.vec[-(p+1)] < th*Cn.vec[p+1]))
            Cn.vec
            
            ## selected model
            nsel = length(which.sel)
            if(nsel>0){
                X.sel = X[,c(1,which.sel+1)]
                lmmsel = ecmeml1(y=y, subj=subj, pred=X.sel, xcol=1:(nsel+1), zcol=1:pZ)
                beta_hat_list[[ind]][which.sel] = lmmsel$beta[-1]
                b_hat_list[[ind]] = lmmsel$b.hat
            }
        }
        bic_vec = sapply(1:nsd, function(ind){
            betahat = beta_hat_list[[ind]]
            r = y - X %*% c(0,betahat)
            for(im in 1:m){ # subtract random effect part
                iind = ((im-1)*ni+1):(im*ni)
                r[iind] = r[iind] - Z[iind,] %*% b_hat_list[[ind]][,im]
            }
            crossprod(r)/sigma2 + sqrt(n)*sd*sum(betahat!=0)/2 # calculate bic with sd
        })
        beta_hat_list[[which.min(bic_vec)]]
    }
    set.seed(11092015)
    if(ncores==1){
        all_beta_hat = lapply(1:1e2, function(x) try(loopfun(x)))
    } else{
        all_beta_hat = mclapply(1:1e2, function(x) try(loopfun(x)), mc.cores=ncores)
    }
    all_beta_hat = all_beta_hat[which(sapply(all_beta_hat,class)=="numeric")]
    
    # get metrics
    nobj = length(all_beta_hat)
    metric.mat = matrix(0, ncol=6, nrow=nobj)
    for(obi in 1:nobj){
        beta.i = all_beta_hat[[obi]]
        metric.mat[obi,] = c(
            sum(beta.i != 0 & beta != 0)/sum(beta.i != 0), # TPR
            sum(beta.i != 0 & beta == 0)/sum(beta.i != 0), # FPR
            sum(beta.i == 0 & beta == 0)/(p-sum(beta.i != 0)), # TNR
            sum(beta.i == 0 & beta != 0)/(p-sum(beta.i != 0)), # FNR
            paste(which(beta.i != 0), collapse="|") == "1|2", # accuracy
            sum(beta.i != 0)) # size
    }
    rbind(apply(metric.mat,2,mean),apply(metric.mat,2,sd))
}
