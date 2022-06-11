# ----------------------------------------------------------------------
# load libraries
# ======================================================================
Required.Packages <- c("R.utils","ddalpha","parallel","doParallel","bestsubset")
sapply(Required.Packages, FUN = function(x) {suppressMessages(require(x, character.only = TRUE))})
source("get_outputs_mio.R")

# grid of parameters
n = 500
p = 100
niter = 7

# setting 1
sigma = 1
eq = 0
k = 1
rho_vec = (1:9)/10
nrho = length(rho_vec)
out_list = vector("list",nrho)
for(r in 1:nrho){
    cat("Setting 1: Doing rho =", rho_vec[r],"-> ")
    out_list[[r]] = get_outputs_mio(rho=rho_vec[r], sigma=sigma, k=k, eq=eq)
    saveRDS(out_list, file=paste0("s1mio_n",n,"p",p,"s",10*sigma,"e",as.numeric(eq),".rds")) # save after each rho
    cat("done\n")
}

# setting 2
rho = 0.5
eq = 1
k_vec = 1:5
nk = length(k_vec)
out_list = vector("list",nk)
for(r in 1:nk){
  cat("Setting 2: Doing k =", k_vec[r],"-> ")
  out_list[[r]] = get_outputs_mio(rho=rho, sigma=sigma, k=k_vec[r], eq=eq)
  saveRDS(out_list, file=paste0("s2mio_n",n,"p",p,"s",10*sigma,"e",as.numeric(eq),".rds")) # save after each
  cat("done\n")
}

# setting 3
k = 1
eq = 1
rho = 0.5
sigma_vec = seq(.3, 2.5, by=.2)
nsigma = length(sigma_vec)
out_list = vector("list",nsigma)
for(r in 1:nsigma){
  cat("Setting 3: Doing sigma =", sigma_vec[r],"-> ")
  sigma = sigma_vec[r]
  out_list[[r]] = get_outputs_mio(rho=rho, sigma=sigma, k=k, eq=eq)
  saveRDS(out_list, file=paste0("s3mio_n",n,"p",p,"e",as.numeric(eq),".rds")) # save after each
  cat("done\n")
}