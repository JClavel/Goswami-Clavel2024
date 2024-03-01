## ------------------------------------------------------------------- ##
##                                                                     ##  
## Various functions and wrappers                                      ##
## J. Clavel                                                           ##
## ------------------------------------------------------------------- ##


# Function to gen the joint estimates of variance and mu, in pruning like fashion [see Felsenstein 1973; function from Silvestro et al. 2019]
get_joint_mu_s2 <- function (mu_f,s2_f,mu_g,s2_g){
  s2_fg = (s2_f*s2_g)/(s2_f+s2_g)
  mu_fg = (mu_f/s2_f + mu_g/s2_g) * s2_fg
  return(c(mu_fg,s2_fg))
}

# Sampling of ancestral states (and tips values) from the conditional gaussian distributions using Felesentein pruning like algorithm
# assumes that the branch lengths of the tree are scaled by the variance term "sigma^2"
conditional_sampling <- function(tree, pheno){
  
  # check the order of the tree for the traversal
  if(attr(tree,"order")!="postorder") stop("the tree provided should be in post-order; see ?reorder.phylo")
    
  # parameters
  N <- Ntip(tree)
  pair = seq(from=1, to=(2*N-4), by=2) # to (2*N-2) but (2*N-4) avoid updating the ancestral state? TO CHECK
  edge1 = tree$edge[,1]
  edge2 = tree$edge[,2]
  
  # loop over pair
  for(i in pair){
    j=i+1
    anc = edge1[i]                        # index ancestor of 1 & 2
    index_stem = which(edge2==anc)        # index for ancestor of ancestor
    anc_stem = edge1[index_stem]          # index stem to ancestor of 1 & 2
    v1 = tree$edge.length[i]              # br length of descendant 1
    v2 = tree$edge.length[j]              # br length of descendant 2
    vstem = tree$edge.length[index_stem]  # br length of the stem to ancestor of 1 & 2
    d1 = edge2[i]                         # index of descendant 1
    d2 = edge2[j]                         # index of descendant 2
    p1 = pheno[d1]                        # phenotype of descendant 1
    p2 = pheno[d2]                        # phenotype of descendant 2
    panc_stem = pheno[anc_stem]           # phenotype of ancestor stem
    
    # get conditional for descendents
    val_mu_s2 = get_joint_mu_s2(p1, v1, p2, v2)
    
    # get conditional for stem
    val_mu_s2_stem = get_joint_mu_s2(panc_stem, vstem, val_mu_s2[1], val_mu_s2[2])
    
    # Sample from the normal distribution
    pheno[anc] = rnorm(1, mean=val_mu_s2_stem[1], sd=sqrt(val_mu_s2_stem[2]))
    
  }# end loop
  
  # update the tips states => sample from the normal centered on current node values
  index <- tree$edge[,2]%in%1:Ntip(tree)
  br_length <- tree$edge.length[index]
  pheno[tree$edge[index,2]] = rnorm(Ntip(tree), mean=pheno[tree$edge[index,1]], sd=sqrt(br_length))
  return(pheno)
}

# wrappers
get.descendants <- function( phy , node ) { extract.clade( phy , node )$tip.label }
subset.rows <- function( data , rownames ) { data[ rownames , ] }

# Compute the Mahalanobis distance
cholMaha <- function(X, covX) {
  dec <- chol( covX )
  tmp <- forwardsolve(t(dec), t(X) )
  dist(t(tmp))
}

# Modified code from Navalon et al. 2022
# same as dimensionality.phylo but with disparity 
disparity.phylo <- function( data , phy , fit=NULL, dist="euclidean") {
  nodes.temp <- unique( phy$edge[,1] )
  descendants.temp <- lapply( nodes.temp , FUN = get.descendants , phy = phy )
  names( descendants.temp ) <- nodes.temp			
  node.data.list <- lapply( descendants.temp , FUN = subset.rows , data = data )
  disparity.temp<-list() 
  for (i in 1: length( node.data.list ) ) {
    if(dist=="euclidean"){
      distances<-dist(node.data.list[[i]], method = "euclidean")
    }else{
      distances<-cholMaha(node.data.list[[i]], fit$sigma$Pinv)
    } 
    
    disparity.temp[i]<- mean(distances)  # loop does the same process as Procrustes variance:
  }
  
  names(disparity.temp)<-nodes.temp
  disparity.temp<-unlist(disparity.temp)
  disparity.temp
  disparity.temp.scaled<-disparity.temp / disparity.temp[1] #max(disparity.temp) 
  disparity.temp.scaled
}

# compute branching times - works with fossil data (contrary to branching.times in ape)
branching.times.foss <- function(phylo){
  # Use "phytools" called by mvMORPH
  times<-max(nodeHeights(phylo))-nodeHeights(phylo)[match(1:phylo$Nnode+Ntip(phylo),phylo$edge[,1]),1]
  names(times)<-1:phylo$Nnode+Ntip(phylo)
  return(times)
}

# compute the average pairwise distance on a time bin (modified from Navalon et al. 2022)
compute_avg_bin <- function(disparity, bin=1, tree, reverse=FALSE){
  times <- branching.times.foss(tree)[paste(unique(tree$edge[,1]))] 
  max_time = max(nodeHeights(tree))
  time_bin = seq(from=1, to=max_time, by=bin)
  bin_mat = cbind(time_bin-bin,time_bin)
  list_time_bin = sapply(times, function(v.element) bin_mat[v.element >= bin_mat[,1] & v.element <= bin_mat[,2],2])
  avg_bin_disp <- sapply(time_bin, function(x) mean(disparity[which(list_time_bin==x)]))
  if(any(!is.finite(avg_bin_disp))) avg_bin_disp2 = approx(seq_along(avg_bin_disp),avg_bin_disp,n=length(avg_bin_disp))$y else avg_bin_disp2=avg_bin_disp
  
  if(reverse) avg_bin_disp2 = rev(avg_bin_disp2)
  return(list(avgDisp=avg_bin_disp2, time=time_bin))
}


# compute disparity through time and confidence interval given a model fit
plot_dtt_fit <- function(fit, tree, nsim, dist="euclidean", bin=2, col="red", add=FALSE, reverse=FALSE, col1="black",...){
  
  # Empirical DTT - if a regressor was added to the model we focus on the residuals
  if(fit$dims$m>1) data=fit$residuals else data=fit$variables$Y
  
  # compute the disparity
  empirical_disparity <- disparity.phylo(data=data,
                                         tree, fit=fit, dist=dist) 
  
  # if euclidean distance used we should simulate the traits as independent
  if(dist=="euclidean") fit$sigma$Pinv <- diag(diag(fit$sigma$Pinv))
  
  # simulate data according to model fit
  sim_data <- simulate(fit, nsim=nsim)
  
  # if model fit with covariate, we should remove the covariate from the simulations (e.g., size-free)
  if(fit$dims$m>1){
    effects <- fitted(fit)
    sim_data <- lapply(sim_data, function(sim_data_full) sim_data_full - effects) 
  }
  
  # run the estimation of disparity across the simulations
  results_sim <- sapply(sim_data, function(simulated_dataset){
    disparity.phylo(data=simulated_dataset, tree, fit=fit, dist=dist) 
  })
  
  # retrieve the average and quantiles?
  avg_disparity_simulations <- rowMeans(results_sim)
  lower_disparity_simulations <- apply(results_sim,1,quantile, probs=0.025)
  upper_disparity_simulations <- apply(results_sim,1,quantile, probs=0.975)
  
  
  # plot the disparity through time
  avg_disp = compute_avg_bin(empirical_disparity, bin=bin, tree, reverse=reverse)
  if(!add) plot(avg_disp$avgDisp~avg_disp$time, type="l", xlab="Time (Ma)", ylab="relative subclade disparity",
                col=col1, las=1, ...)
  else
    lines(avg_disp$avgDisp~avg_disp$time, col=col1, ...)
  
  # confidence interval
  avg_sim = compute_avg_bin(avg_disparity_simulations, bin=bin, tree, reverse=reverse)
  lower_sim = compute_avg_bin(lower_disparity_simulations, bin=bin, tree, reverse=reverse)
  upper_sim = compute_avg_bin(upper_disparity_simulations, bin=bin, tree, reverse=reverse)
  
  lines(avg_sim$avgDisp~avg_sim$time, col=col, lty=3)
  lines(lower_sim$avgDisp~lower_sim$time, col=col, lty=2)
  lines(upper_sim$avgDisp~upper_sim$time, col=col, lty=2)
  
  # return dtt measures?
  results <- list(avg_sim=avg_sim, lower_sim=lower_sim, upper_sim=upper_sim, empirical=avg_disp, time=avg_sim$time)
  invisible(results)
}
