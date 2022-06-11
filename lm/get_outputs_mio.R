# master function
get_outputs_mio = function(rho, sigma, k, eq){

    # simulation setup
    Sigma.X = diag(p)
    for(i in 2:p){
        for(j in 1:(i-1)){
            Sigma.X[i,j] = rho^abs(i-j)
            Sigma.X[j,i] = Sigma.X[i,j]
        }
    }
    
    ## Set up coef vector
    if(eq==1){
        beta.sig = rep(1,k*5)
    } else{
        beta.sig = rep(c(1.5, 0.5, 1, 1.5, 1), k)
    }        
    n.sig = length(beta.sig)
    beta = c(beta.sig, rep(0,p - n.sig))
    
    ## iterator function
    run_mio = function(seed){
        
        # generate data
        set.seed(1e3*seed)
        X = mvrnorm(n, mu=rep(0,p), Sigma=Sigma.X)
        y = X %*% beta + sigma*rnorm(n)
        beta.hat = solve(crossprod(X)) %*% t(X) %*% y
        sigma2 = sum(crossprod(y - X %*% beta.hat))/(n-1)
        
        ## MIO forward selection
        runtime = system.time({
        bs.obj = bs(X, y, intercept=F, k=seq(1,29,by=2), time.limit=10)
        beta.mat = coef(bs.obj)
        bicvals = apply(beta.mat, 2, function(b) 
          sum(crossprod(y - X %*% b))/sigma2 + log(n)*sum(b!=0)/2)
        })
        cat("Iteration",seed,"done\n")
        c(runtime[3], beta.mat[, which.min(bicvals)])
    }
    
    # get outputs
    cl = makeCluster(7)
    registerDoParallel(cl) 
    clusterExport(cl,list('bs','run_mio','mvrnorm','n','p','Sigma.X','beta','sigma'))
    out_obj = parLapply(cl, 1:niter, function(x) try(run_mio(x)))
    stopCluster(cl)
    out_obj1 = out_obj[which(sapply(out_obj,class)=="numeric")] # exclude error iterations
    runtimes = sapply(out_obj1, function(x) x[1])
    betamat = matrix(unlist(lapply(out_obj1, function(x) x[-1])), ncol=p, byrow=T)

    # compute metrics
    metric.mat = matrix(0, ncol=6, nrow=niter)
    for(i in 1:niter){
      set.seed(8312020*i)
      beta.i = betamat[i,]
      Xt = mvrnorm(n, mu=rep(0,p), Sigma=Sigma.X)
      yt = Xt %*% beta + sigma*rnorm(n)
      metric.mat[i,] = c(
        runtimes[i],
        sum(beta.i != 0 & beta != 0),
        sum(beta.i != 0 & beta == 0),
        sum(beta.i != 0),
        crossprod(Xt %*% (beta.i - beta))/crossprod(Xt %*% beta),
        mean(crossprod(yt - Xt %*% beta.i)))
    }
    c(apply(metric.mat,2,mean), apply(metric.mat,2,sd))
}