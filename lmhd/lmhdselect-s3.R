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
    make_option(c("-e", "--equal"), type="double", default=1, help="equal coef indicator [default = %default]", metavar="double")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
n = opt$n
p = opt$p
nboot = opt$b
eq = opt$e

# grid of parameters
tau_vec = log(n)*seq(.2,2,by=.4)*sqrt(log(p)/n)
th = 1
k = 1
niter = 1e2
rho = 0.5
sigma_vec = seq(.3, 2.5, by=.2)
nsigma = length(sigma_vec)
out_list = vector("list",nsigma)
for(rs in 1:nsigma){
    cat("Setting 3: Doing sigma =", sigma_vec[rs],"-> ")
    out_list[[rs]] = get_outputs(rho=rho, sigma=sigma_vec[rs], k=k, eq=eq)
    saveRDS(out_list, file=paste0("s3hd_n",n,"p",p,"b",nboot,"e",as.numeric(eq),".rds")) # save after each rho
    cat("done\n")
}