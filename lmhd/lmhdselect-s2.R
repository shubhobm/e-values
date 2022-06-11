# ----------------------------------------------------------------------
# load libraries
# ======================================================================
Required.Packages <- c("R.utils","ddalpha","parallel","ncvreg","knockoff","SIS","optparse")
sapply(Required.Packages, FUN = function(x) {suppressMessages(require(x, character.only = TRUE))})
source("get_outputs_hd.R")

## parse arguments *****************************************************
# **********************************************************************
option_list = list(
    make_option(c("-n", "--numsamp"), type="double", default=100, help="sample size [default = %default]", metavar="double"),
    make_option(c("-p", "--params"), type="double", default=500, help="params size [default = %default]", metavar="double"),
    make_option(c("-b", "--bootsamp"), type="double", default=1e3, help="resample size [default = %default]", metavar="double"),
    make_option(c("-s", "--sigma"), type="double", default=1, help="error sd [default = %default]", metavar="double"),
    make_option(c("-e", "--equal"), type="double", default=1, help="equal coef indicator [default = %default]", metavar="double")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
n = opt$n
p = opt$p
nboot = opt$b
sigma = opt$s
eq = opt$e

# grid of parameters
tau_vec = log(n)*seq(.2,2,by=.4)*sqrt(log(p)/n)
th = 1
rho = 0.5
niter = 1e2
k_vec = 1:5
nk = length(k_vec)
out_list = vector("list",nk)
for(rk in 1:nk){
    cat("Setting 2: Doing k =", k_vec[rk],"-> ")
    out_list[[rk]] = get_outputs(rho=rho, sigma=sigma, k=k_vec[rk], eq=eq)
    saveRDS(out_list, file=paste0("s2hd_n",n,"p",p,"s",10*sigma,"b",nboot,"e",as.numeric(eq),".rds")) # save after each rho
    cat("done\n")
}