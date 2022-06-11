## function to do depth-based model selection in linear model
# AR(4) linear fit
lm.fit.ar4 = function(x,y){
  n = length(y)
  ar.y = cbind(y[1:(n-4)],y[2:(n-3)],y[3:(n-2)],y[4:(n-1)])
  total.x = cbind(ar.y, x[5:n,])
  scale.x = cbind(1,scale(total.x[,-7]))
  scale.y = scale(y[5:n])
  list(fit=lm.fit(scale.x, scale.y), X=scale.x, y=scale.y)
}

lm.fit.spatial = function(x,coords){
  ci = coords[1]; cj = coords[2]; ck = coords[3]
  if(ci==1 | cj==1 | ck==1 | ci==64 | cj==64 | ck==33){
    neighbor.y = matrix(NA, nrow=nrow(YData), ncol=26)
  } else{
    grid = expand.grid((ci-1):(ci+1),(cj-1):(cj+1),(ck-1):(ck+1))
    a = as.numeric((grid[,1] - 1)*Dim[2]*Dim[3] + (grid[,2] -1)*Dim[3] + grid[,3])
    neighbor.y = YData[,a]
  }
  ## control for all-zero responses
  if(min(diag(var(neighbor.y)))>0){
    # cat(coords)
    # cat("\n")
    # cat(a)
    # cat("\n")
    total.x = cbind(neighbor.y[,-14],x)
    scale.x = scale(total.x[,-29])
    scale.y = scale(neighbor.y[,14])
    out.list = list(fit=lm.fit(scale.x, scale.y), X=scale.x, y=scale.y)
  } else{
    out.list = list(fit=NA, X=NA, y=NA)
  }
  
  # return
  out.list
}

# depth-based selection
lm.select.depth = function(X,y,method="lm", threshold=.9, proportion=.9,
                           sd.vec = seq(.05, .25, by=.01),
                           nboot=1e3){
  list0 = ls()
  # train full model
  set.seed(10252016)
  # train = sample(1:n, ceiling(n*.8), replace=F)
  # ntrain = length(train)
  
  if(method=="lm.ar4"){
    fitObj = lm.fit.ar4(X,y)
    X = fitObj$X
    y = fitObj$y
    fitObj = fitObj$fit
  } else{
    fitObj = lm.fit(X,y)
  }
  beta.hat = fitObj$coef
  beta.cov = var(fitObj$residuals) * solve(crossprod(X))
  H = beta.cov %*% t(X)
  r = fitObj$residuals
  
  # initialize quantities
  n = nrow(X)
  p = ncol(X)
  nboot = 1e3
  sdn.vec = n^sd.vec
  
  ## matrix of full model bootstrap betas
  beta.mat = matrix(0 ,nrow=nboot, ncol=p)
  resid.mat = matrix(rnorm(n*nboot), ncol=nboot) * matrix(r, nrow=n, ncol=nboot, byrow=F)
  score.mat = t(H %*% resid.mat)
  depth.full = mdepth.RP(score.mat, score.mat)$dep
  Cn.full = mean(depth.full)
  ## we are going to recycle this score matrix for different values of bootstrap standard deviation
  
  ## loop over the parameters to save on memory:
  ## recycle memory by not storing H, covariance matrix etc for truncated models
  loopfun = function(j){
    
    ## calculate quantities for truncated model
    Hj = H[-j,] - outer(beta.cov[j,-j], as.numeric(X[,j]))
    rj = r + X[,j] * beta.hat[j]
    jresid.mat = matrix(rnorm(n*nboot), ncol=nboot) * matrix(rj, nrow=n, ncol=nboot, byrow=F)
    jscore.mat = t(Hj %*% jresid.mat)
    
    ## calculate Cn for truncated model, for a range of bootstrap variances
    jdepth.mat = matrix(0, ncol=length(sdn.vec), nrow=nboot)
    for(i in 1:length(sdn.vec)){
      sdn = sdn.vec[i]
      jbeta.mat = matrix(0, ncol=p, nrow=nboot)
      jbeta.mat[,-j] = matrix(beta.hat[-j], nrow=nboot, ncol=p-1, byrow=T) + sdn*jscore.mat
      jdepth.mat[,i] = mdepth.RP(matrix(beta.hat, nrow=nboot, ncol=p, byrow=T) + sdn*score.mat,
                                 jbeta.mat)$dep
    }
    
    # return mean depth of truncated model, for all values of bootstrap sd
    Cn.j = apply(jdepth.mat, 2, mean)
    Cn.j
  }
  
  # Cn.mat1 = lapply(1:p, loopfun)
  # system.time(
  Cn.mat1 <- lapply(1:p, loopfun)
  # )
  Cn.mat = matrix(unlist(Cn.mat1), ncol=length(sdn.vec), byrow=T)
  
  ## save best index list
  index.list = list()
  for(i in 1:length(sdn.vec)){
    index.list[[i]] = which(Cn.mat[,i] < proportion*Cn.full)
  }
  
  # find best model
  cnt.list = list()
  for(j in 1:p){
    cnt.list[[j]] = as.numeric(lapply(index.list, function(x) j %in% x))
  }
  cnt.mat = matrix(unlist(cnt.list),nrow=p, byrow=TRUE)
  best.index = which(rowSums(cnt.mat) > threshold*length(sdn.vec))
  
  # clear memory for unneeded elements
  list1 = ls()
  keep.index = which(list1 %in% list0 | list1=="best.index" | list1=="beta.hat")
  rm(list=list1[-keep.index])
  
  # return
  list(best.index=best.index, beta=beta.hat)
}

# depth-based selection.. spatial version
lm.select.spatial = function(X, coords, threshold=.9,
                           sd.vec = seq(.05, .25, by=.01),
                           nboot=1e3){
  list0 = ls()
  # train full model
  set.seed(10252016)

  fitObj = lm.fit.spatial(X,coords)
  X = fitObj$X
  y = fitObj$y
  fitObj = fitObj$fit
  n = nrow(X)
  p = ncol(X)
  
  # do selection only when there is a model
  if(length(fitObj)==1){
    out.list = list(position=coords, best.index=NA, beta=NA)
  } else{
    beta.hat = fitObj$coef
    beta.cov = var(fitObj$residuals) * solve(crossprod(X))
    H = beta.cov %*% t(X)
    r = fitObj$residuals
    
    # initialize quantities
    nboot = 1e3
    sdn.vec = n^sd.vec
    
    ## matrix of full model bootstrap betas
    beta.mat = matrix(0 ,nrow=nboot, ncol=p)
    resid.mat = matrix(rnorm(n*nboot), ncol=nboot) * matrix(r, nrow=n, ncol=nboot, byrow=F)
    score.mat = t(H %*% resid.mat)
    depth.full = mdepth.RP(score.mat, score.mat)$dep
    Cn.full = mean(depth.full)
    ## we are going to recycle this score matrix for different values of bootstrap standard deviation
    
    ## loop over the parameters to save on memory:
    ## recycle memory by not storing H, covariance matrix etc for truncated models
    loopfun = function(j){
      
      ## calculate quantities for truncated model
      Hj = H[-j,] - outer(beta.cov[j,-j], as.numeric(X[,j]))
      rj = r + X[,j] * beta.hat[j]
      jresid.mat = matrix(rnorm(n*nboot), ncol=nboot) * matrix(rj, nrow=n, ncol=nboot, byrow=F)
      jscore.mat = t(Hj %*% jresid.mat)
      
      ## calculate Cn for truncated model, for a range of bootstrap variances
      jdepth.mat = matrix(0, ncol=length(sdn.vec), nrow=nboot)
      for(i in 1:length(sdn.vec)){
        sdn = sdn.vec[i]
        jbeta.mat = matrix(0, ncol=p, nrow=nboot)
        jbeta.mat[,-j] = matrix(beta.hat[-j], nrow=nboot, ncol=p-1, byrow=T) + sdn*jscore.mat
        jdepth.mat[,i] = mdepth.RP(matrix(beta.hat, nrow=nboot, ncol=p, byrow=T) + sdn*score.mat,
                                   jbeta.mat)$dep
      }
      
      # return mean depth of truncated model, for all values of bootstrap sd
      Cn.j = apply(jdepth.mat, 2, mean)
      Cn.j
    }
    
    # Cn.mat1 = lapply(1:p, loopfun)
    # system.time(
    Cn.mat1 <- lapply(1:p, loopfun)
    # )
    Cn.mat = matrix(unlist(Cn.mat1), ncol=length(sdn.vec), byrow=T)
    
    ## save best index list
    index.list = list()
    for(i in 1:length(sdn.vec)){
      index.list[[i]] = which(Cn.mat[,i] < Cn.full)
    }
    
    # find best model
    cnt.list = list()
    for(j in 1:p){
      cnt.list[[j]] = as.numeric(lapply(index.list, function(x) j %in% x))
    }
    cnt.mat = matrix(unlist(cnt.list),nrow=p, byrow=TRUE)
    best.index = which(rowSums(cnt.mat) > threshold*length(sdn.vec))
    
    # flush memory for unneeded elements
    list1 = ls()
    keep.index = which(list1 %in% list0 | list1=="best.index" | list1=="beta.hat")
    rm(list=list1[-keep.index])
    
    out.list = list(position=coords, best.index=best.index, beta=beta.hat)
  }
  
  # return
  out.list
}