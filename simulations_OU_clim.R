## ------------------------------------------------------------------------- ##
## Simulations to assess model selection with the OU-climatic model.         ##
##                                                							             ##	
##  Codes sent on a cluster.  				                                       ##
## ------------------------------------------------------------------------- ##

library(mvMORPH)# >v1.1.9
library(RPANDA) # >v2.3

# define parameters and functions for simulations ####
nsim=100               # number of reps to perform
root_age = 60          # age (Mya) of the root node in tree
ntax=164               # number of tips in simulated trees
nfoss= c(0, 0.05, 0.1, 0.25, 0.5)           # what proportion of tips should be fossil?
delta = 0.001          # time step for the simulations (My)
sim_sigma2 = 0.025     # variance of the random noise = sigma^2
sim_alpha = c((log(2)/(root_age/0.5)),(log(2)/(root_age)),(log(2)/(root_age/3)),(log(2)/(root_age/5)),(log(2)/(root_age/10)))    # alpha of the OU process, expressed as phylogenetic 1/2-life
sim_theta = 4          # Theta value at the root
sim_beta = c(-5,-1,0,1,5)        # Beta, relationship between climate and optimum
st_val1<-seq(-6, 6, 0.75)        # start values for parameter estimation, here a range to avoid local optima in likelihood space
start_val = list(par1=st_val1)

# function for fitting the model
fun_temp <- function(t, env, param, theta0) theta0 + param*env(t)

# load climate data, scaled 0-1 
curb<-read.table("Climate_curve/Cramer2011trend.csv", sep=";", dec=".", header=T)
cur1=scale(curb[,3],min(curb[,3],na.rm=T),max(curb[,3], na.rm=T)-min(curb[,3],na.rm=T))
d18=splinefun(curb[,1],cur1)

# Define the list or argument we will use for each simulations (i.e. each jobs will be a possible combination of both parameters)
parameters_jobs <- expand.grid(sim_alpha, sim_beta, nfoss)

# retrieve job "id" using commandArgs (jobs numbering from 0)
args = commandArgs(TRUE)
i = as.numeric(args[1])# + 1L
param <- as.matrix(parameters_jobs[i , ]) # param[1] is alpha and param[2] is beta, param[3] is the proportion of fossils
fosstips<-param[[3]]


results <- matrix(NA, nrow=nsim, ncol=18)
results_param <- matrix(NA, nrow = nsim, ncol=5)

# start iterate across nsim simulations
for(k in 1:nsim){
  # simulate tree with a given proportion of fossils and scaled to "root_age"
  N_tips<-1
  while(N_tips!=ntax){
    tree <- pbtree(b=1, d=fosstips, n=((1-fosstips)*ntax), scale=root_age, t=NULL, nsim=1, type=c("continuous"))
    N_tips<-Ntip(tree)
  }
  
  # Simulate traits
  message(c("Simulating data rep ", k))
  trait <- sim_t_env_ou(phylo=tree, sigma=sqrt(sim_sigma2), alpha=param[1],
               theta0=sim_theta, param=param[2], model=fun_temp,
               env_data=d18, step=delta, plot=FALSE)
  
  # clim-OU [with previous implementation - on github - if there's issues in parameter search, I switch to Nelder-Mead algorithm].
  fit_ou_clim <- fit_t_env_ou(phylo = tree, data = trait, model=fun_temp, env_data=d18, param = start_val, method = "Nelder-Mead", control=list(maxit=10000))
  # clim model rates
  fit_env <- fit_t_env(phylo = tree, data=trait, env_data=d18, model="EnvLin")
  # Bm
  fit_bm <- mvBM(tree = tree, data = trait, model = "BM1", echo = FALSE, diagnostic=FALSE)
  # OU
  fit_ou <- mvOU(tree = tree, data = trait, model = "OU1", echo = FALSE, diagnostic=FALSE, method="sparse")
  # EB
  fit_eb <- mvEB(tree = tree, data = trait, echo = FALSE, diagnostic=FALSE)
  # BM trend
  fit_bm_trend <- mvBM(tree=tree, data=trait, model="BM1", param=list(trend=TRUE), echo = FALSE, diagnostic=FALSE)
  
  # retrieve results
  results[k,] <- c(fit_ou_clim$aicc, fit_ou_clim$convergence, fit_ou_clim$hess.value,
                   fit_env$aicc, fit_env$convergence, fit_env$hess.value,
                   fit_bm$AICc, fit_bm$convergence, fit_bm$hess.value,
                   fit_ou$AICc, fit_ou$convergence, fit_ou$hess.value,
                   fit_eb$AICc, fit_eb$convergence, fit_eb$hess.value,
                   fit_bm_trend$AICc, fit_bm_trend$convergence, fit_bm_trend$hess.value)
  # save parameter estimates
  results_param[k,] <- unlist(fit_ou_clim$param)
}

colnames(results) = c("env_ou","convergence","hessian",
                      "env_bm","convergence","hessian",
                      "bm","convergence","hessian",
                      "ou","convergence","hessian",
                      "eb","convergence","hessian",
                      "trend","convergence","hessian")

colnames(results_param) = c("sigma2", "alpha", "theta_0", "beta", "root_implied")

filename <- paste("model_selection_",fosstips,"_alpha_",round(param[1], digits = 3),"_beta_",param[2],".csv",sep="")
write.csv(results, file=filename)

filename_par <- paste("model_param_",fosstips,"_alpha_",round(param[1], digits = 3),"_beta_",param[2],".csv",sep="")
write.csv(results_param, file=filename_par)

# Note: with low alpha value, the model fit may converge on unlikely values. Model fit with outlier or that did not converged are not included in the comparisons
