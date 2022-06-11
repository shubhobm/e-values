# master function
get_outputs = function(rho, sigma, k, eq){

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
    run_all_methods = function(seed){
        
        # generate data
        set.seed(1e3*seed)
        X = mvrnorm(n, mu=rep(0,p), Sigma=Sigma.X)
        y = X %*% beta + sigma*rnorm(n)
        
        # initialize
        all_outputs = matrix(0, nrow=5, ncol=p+1)
        
        ## The e-value method
        # compute full model constraints
        XtY = t(X) %*% y
        H.eig = eigen(crossprod(X))
        H2 = with(H.eig, vectors %*% diag(1/sqrt(values)) %*% t(vectors)) %*% t(X)
        beta.hat = with(H.eig, vectors %*% diag(1/values) %*% t(vectors)) %*% XtY
        r = y - X %*% beta.hat; r = r - mean(r)
        sigma2 = sum(r^2)/(n-1)
        
        # draw bootstrap weights and calculate weighted residuals
        runtime = system.time({
        H2_w_r = lapply(1:nboot, function(x) H2 %*% ((rgamma(n,1,1) - 1)*r))

        # repeat for multiple tau and select one with highest bic
        SSPmat.d = matrix(0, nrow=nboot, ncol=p+1)
        beta.mat.m = matrix(0 ,nrow=nboot, ncol=p)
        
        ntau = length(tau_vec)
        beta_hat_list = vector("list",ntau)
        for(ind in 1:ntau){
            # initialize
            beta_hat_list[[ind]] = rep(0,p)
            tau = tau_vec[ind]
            
            # Cn for full model estimates: Fn is approx by Fn^b1
            for(i in 1:nboot){
                # w = rgamma(n,1,1) - 1
                # beta.mat.m[i,] = beta.hat + tau * H2 %*% (w*r)
                beta.mat.m[i,] = beta.hat + tau * H2_w_r[[i]]
            }
            SSPmat.d[,p+1] = depth.Mahalanobis(beta.mat.m, beta.mat.m)
            
            # marginal models
            for(j in 1:p){
                jbeta.mat.m = beta.mat.m
                jbeta.mat.m[,j] = 0
                SSPmat.d[,j] = depth.Mahalanobis(jbeta.mat.m, beta.mat.m)
            }
            mean.vec = apply(SSPmat.d, 2, mean)
            which.sel = which(mean.vec[-(p+1)] < th*mean.vec[p+1])
            
            ## selected model
            X.sel = X[,which.sel]
            if(length(which.sel)>0){
                beta.hat.sel = solve(crossprod(X.sel)) %*% t(X.sel) %*% y
                beta_hat_list[[ind]][which.sel] = beta.hat.sel
            }
        }
        bic_vec = sapply(beta_hat_list, function(b)
            crossprod(y - X %*% b)/sigma2 + sqrt(n)*tau*sum(b!=0)/2) # calculate bic with tau
        })
        all_outputs[1,] = c(runtime[3], beta_hat_list[[which.min(bic_vec)]])
        
        # ## MIO forward selection
        # runtime = system.time({
        # fs.obj = bestsubset::fs(X1, y, intercept=F)
        # beta.mat = coef(fs.obj)
        # bicvals = apply(beta.mat, 2, function(b) n*log(mean((y - X1 %*% b)^2))+sum(b!=0)*log(n))
        # })
        # all_outputs[2,] = c(runtime[3], beta.mat[, which.min(bicvals)])
        
        ## lasso with BIC selection
        runtime = system.time({
        lasso.obj = ncvreg(X, y, penalty="lasso")
        bic.vec = with(lasso.obj, (loss/sigma2 + log(n)*apply(beta[-1,]!=0,2,sum)/2))
        })
        all_outputs[2,] = c(runtime[3], lasso.obj$beta[-1,which.min(bic.vec)])
        
        ## scad with BIC selection
        runtime = system.time({
        scad.obj = ncvreg(X, y, penalty="SCAD")
        bic.vec = with(scad.obj, (loss/sigma2 + log(n)*apply(beta[-1,]!=0,2,sum)/2))
        })
        all_outputs[3,] = c(runtime[3], scad.obj$beta[-1,which.min(bic.vec)])
        
        ## stepwise
        beta.step = rep(0,p)
        runtime = system.time({
        df = data.frame(cbind(y,X))
        names(df) = c("y",names(df)[-(p+1)])
        lmod = lm(y~., data=df)
        #stepmod = stepAIC(lmod, direction="backward", trace=F, k=log(n))
        stepmod = stepAIC(lmod, direction="backward", trace=F)
        #lmod0 = lm(y~1, data=df)
        #scope=paste(c("y~",paste0("X",1:p, collapse="+")), collapse=" ")
        #stepmod = stepAIC(lmod0, direction="forward", trace=F, scope=as.formula(scope), k=log(n))
        which.sel = which(colnames(df) %in% names(coef(stepmod))) - 1
        if(length(which.sel)>0){
            beta.step[which.sel] = coef(lm.fit(X[,which.sel],y))
        }
        })
        all_outputs[4,] = c(runtime[3], beta.step)
        
        ## knockoffs
        beta.ko = rep(0,p)
        runtime = system.time({
        suppressWarnings(kmod <- knockoff.filter(X,y,offset=0))
        which.sel = kmod$selected
        if(length(which.sel)>0){
            beta.ko[which.sel] = coef(lm.fit(X[,which.sel],y))
        }
        })
        all_outputs[5,] = c(runtime[3], beta.ko)
        
        # return
        row.names(all_outputs) = c("evalue","lasso","scad","step","ko")
        t(all_outputs)
    }
    
    # evaluator function
    evaluate = function(betamat, beta){
        niter = nrow(betamat)
        metric.mat = matrix(0, ncol=5, nrow=niter)
        
        for(i in 1:niter){
            set.seed(8312020*i)
            beta.i = betamat[i,]
            Xt = mvrnorm(n, mu=rep(0,p), Sigma=Sigma.X)
            yt = Xt %*% beta + sigma*rnorm(n)
            metric.mat[i,] = c(
                sum(beta.i != 0 & beta != 0),
                sum(beta.i != 0 & beta == 0),
                sum(beta.i != 0),
                crossprod(Xt %*% (beta.i - beta))/crossprod(Xt %*% beta),
                mean(crossprod(yt - Xt %*% beta.i)))
        }
        c(apply(metric.mat,2,mean), apply(metric.mat,2,sd))
    }
    
    # get outputs
    out_obj = mclapply(1:niter, function(x) try(run_all_methods(x)), mc.cores=50)
    out_obj = out_obj[which(sapply(out_obj,class)=="matrix")] # exclude error iterations
    runtimes = sapply(out_obj, function(x) x[1,])
    nobj = length(out_obj)
    coef.array = array(0, c(p,5,nobj))
    for(it in 1:nobj){
        coef.array[,,it] = out_obj[[it]][-1,]
    }
    
    # compute metrics
    metric.mat = matrix(0, nrow=5, ncol=11)
    row.names(metric.mat) = c("evalue","lasso","scad","step","ko")
    colnames(metric.mat) = c(paste0(c("runtime","tp","fp","sp","pe","mspe"),"_mean"),
                             paste0(c("tp","fp","sp","pe","mspe"),"_sd"))
    metric.mat[,1] = rowMeans(runtimes)
    for(i in 1:5){
        metric.mat[i,-1] = evaluate(t(coef.array[,i,]),beta)
    }
    metric.mat
}
