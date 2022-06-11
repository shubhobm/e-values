# ----------------------------------------------------------------------
# load libraries
# ======================================================================
Required.Packages <- c("R.utils","data.table","ddalpha","parallel","rrcov","lme4","tidyverse")
sapply(Required.Packages, FUN = function(x) {suppressMessages(require(x, character.only = TRUE))})

# read in data
rainsmall = data.frame(fread("rainsmall.csv"))

# train full model
rainsmall[,-(1:3)] = scale(rainsmall[,-(1:3)])
rainsmall = data.table(rainsmall)[order(year)]
trainset = which(rainsmall[,year]<2003)
testset = which(rainsmall[,year]>=2003)
varnames = names(rainsmall)[-(1:3)]
formula = paste(varnames, collapse="+")
random_terms = "+ (1|year)"
formula = as.formula(paste("log(PRCP+1) ~", formula, random_terms))
mod.full = lmer(formula, data=rainsmall, subset=trainset)
y = getME(mod.full, 'y')
x = getME(mod.full, 'X')
n = nrow(x)
p = ncol(x)-1
nboot = 1e3

# full model quantities
ni_vec = rainsmall[trainset,.N,by=year][,N]
m = length(ni_vec)
sigma2 = getME(mod.full, 'sigma')^2
Z = getME(mod.full, 'Z')
Psi = getME(mod.full, 'Lambda')
W = sigma2 + Z %*% Psi %*% t(Z)
beta.hat = getME(mod.full, 'fixef')
H = vcov(mod.full) %*% t(x) %*% W
r = y - x %*% beta.hat; r = r - mean(r)

set.seed(09042020)
sdm = seq(.05,1,by=.05)
nsd = length(sdm)
mod_list = vector("list",nsd)
cn_list = vector("list",nsd)
for(ind in 1:nsd){
    # initialize
    beta_hat_list[[ind]] = rep(0,p)
    mu_list[[ind]] = 0
    sd = sdm[ind]            
    
    # matrix of full model estimates: Fn is approx by Fn^b1
    beta.mat = matrix(0 ,nrow=nboot, ncol=p+1)
    SSPmat.d = matrix(0, nrow=nboot, ncol=p+1)
    for(i in 1:nboot){
        iresid = as.matrix(sd * rep(rgamma(m,1,1)-1, ni_vec) * r, ncol=1)
        beta.mat[i,] = as.numeric(beta.hat) + as.numeric(H %*% iresid)
    }
    
    beta.mat1 = matrix(0 ,nrow=nboot, ncol=p+1)
    for(i in 1:nboot){
        iresid = as.matrix(sd * rep(rgamma(m,1,1)-1, ni_vec) * r, ncol=1)
        beta.mat1[i,] = as.numeric(beta.hat) + as.numeric(H %*% iresid)
    }
    SSPmat.d[,p+1] = depth.projection(beta.mat1[,-1], beta.mat[,-1], seed=1e3)
    
    ## now marginal models
    for(j in 1:p){
        jbeta.mat1 = beta.mat1
        jbeta.mat1[,j+1] = 0
        SSPmat.d[,j] = depth.projection(jbeta.mat1[,-1], beta.mat[,-1], seed=1e3)
    }
    Cn.vec = apply(SSPmat.d, 2, mean)
    which.sel = which(Cn.vec[-(p+1)] < th*Cn.vec[p+1])
    cn_list[[ind]] = Cn.vec

    ## selected model
    nsel = length(which.sel)
    if(nsel>0){
        formula.sel = as.formula(
            paste("log(PRCP+1) ~", paste(varnames[which.sel], collapse="+"), random_terms))
        mod_list[[ind]] = lmer(formula.sel, data=rainsmall, subset=trainset)
    }
}

# calculate bic and prediction errors
ytest = rainsmall[testset,log(PRCP+1)]
xtest = as.matrix(data.frame(rainsmall[testset])[,-(1:3)])
bic_vec = rep(crossprod(y)/sigma2, nsd)
pe_vec = rep(crossprod(ytest)/sigma2, nsd)
for(ind in which(sapply(mod_list,class)=="lmerMod")){
    mod.ind = mod_list[[ind]]
    betahat.ind = getME(mod.ind, "beta")[-1]
    which.sel = which(cn_list[[ind]][-(p+1)] < th*cn_list[[ind]][p+1])
    bic_vec[[ind]] = crossprod(
        y - getME(mod.ind, "mu"))/sigma2 +
        sqrt(n)*sdm[ind]*sum(betahat.ind!=0)/2 # calculate bic with sd
    pe_vec[[ind]] = crossprod(ytest - xtest[,which.sel] %*% betahat.ind)
}

# select using bic
best_bic = which.min(bic_vec)
which.sel1 = which(cn_list[[best_bic]][-(p+1)] < th*cn_list[[best_bic]][p+1])
beta.hat.best1 = rep(0,p)
beta.hat.best1[which.sel1] = getME(mod_list[[best_bic]], "beta")[-1]
names(beta.hat.best1) = names(beta.hat)[-1]
# t statistics
tstat1 = fixef(mod_list[[best_bic]])[-1]/sqrt(diag(vcov(mod_list[[best_bic]]))[-1])

# select using holdout prediction
best_pe = which.min(pe_vec)
which.sel2 = which(cn_list[[best_pe]][-(p+1)] < th*cn_list[[best_pe]][p+1])
beta.hat.best2 = rep(0,p)
beta.hat.best2[which.sel2] = getME(mod_list[[best_pe]], "beta")[-1]
names(beta.hat.best2) = names(beta.hat)[-1]
# t statistics
tstat2 = fixef(mod_list[[best_pe]])[-1]/sqrt(diag(vcov(mod_list[[best_pe]]))[-1])

# make graph
dt = merge(data.table(tstat1, var=names(tstat1)),
                  data.table(tstat2, var=names(tstat2)),
                  by="var",all=T)
saveRDS(list(cbind(beta.hat.best1, beta.hat.best2), dt),
        file="output.rds")
dt[is.na(tstat2), tstat2 := 0]
dt[abs(tstat2)>=10, tstat2 := 10*sign(tstat2)]
dt[is.na(tstat2), tstat2 := 0]
dt[abs(tstat1)>=10, tstat1 := 10*sign(tstat1)]
dt
names(dt)[2:3] = c("GBIC","Test Error")
dt_melt = melt(dt, id.vars="var", measure.vars=names(dt)[2:3])
names(dt_melt)[2] = "Method"

pdf("monsoon.pdf", width=7, height=4)
ggplot(dt_melt, aes(x=var, y=value, fill=Method)) +
    geom_bar(stat = "identity", position=position_dodge()) +
    theme_bw() +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    geom_hline(yintercept=qnorm(.975), linetype="dashed") +
    geom_hline(yintercept=-qnorm(.975), linetype="dashed") +
    xlab("Covariate") +
    ylab("t-statistic")
dev.off()
