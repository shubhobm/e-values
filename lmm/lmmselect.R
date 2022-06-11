# ----------------------------------------------------------------------
# load libraries
# ======================================================================
Required.Packages <- c("R.utils","ddalpha","parallel","lmm")
sapply(Required.Packages, FUN = function(x) {suppressMessages(require(x, character.only = TRUE))})
source("get_outputs.R")

## Run it!!
set.seed(11092015)
p = 9
sig = 1
nboot = 1e3
n.iter = 1e2
sdm = seq(1,5,by=.2)
pZ = 4
beta = c(rep(1,2),rep(0,p-2))
D = matrix(c(9,   4.8, 0.6, 0,
             4.8, 4,   1,   0,
             0.6, 1,   1,   0,
             0,   0,   0,   0),
           ncol=4, byrow=T)
th_vec = c(1,.99,.95,.9,.85)
nth = length(th.vec)

out1 = vector("list",nth)
out2 = vector("list",nth)
for(ith in 1:nth){
    out1[[ith]] = get_outputs(beta=beta, sig=sig, D=D, ni=5, m=30, sdm=sdm, nboot=nboot, th=th_vec[[ith]], n.iter=n.iter, ncores=50)
    saveRDS(out1, file="s1_ni5m30.rds")
    out2[[ith]] = get_outputs(beta=beta, sig=sig, D=D, ni=10, m=60, sdm=sdm, nboot=nboot, th=th_vec[[ith]], n.iter=n.iter, ncores=50)
    saveRDS(out2, file="s2_ni10m60.rds")
}
