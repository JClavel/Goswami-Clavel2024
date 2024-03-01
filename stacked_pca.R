## ------------------------------------------------------------------------- ##
## Reconstruct a stacked PCA for the mammals dataset and with expectations   ##
## from a BM model										                            					 ##	
##  																	                                    	 ##	
## ------------------------------------------------------------------------- ##
source("fun.R")
library(mvMORPH)

# phylogenetic tree
tree <- read.nexus("Mammals_whole_dataset/tree80_85_85.nex")
tree <- reorder.phylo(tree,"postorder")

# tree with branch lengths scaled by the rates
treeScaled = read.nexus("Mammals_whole_dataset/treeScaleRate.nex")
treeScaled <- reorder.phylo(treeScaled,"postorder")

# Morphology (example focused on the landmarks)
data <- read.csv("Mammals_whole_dataset/shape.data.322.csv", row.names = 1)[,1:198]

# tip age (and reordering)
tip_age <- read.csv("Mammals_whole_dataset/tip_ages.csv", row.names = 1, header = TRUE)
tip_age <- tip_age[tree$tip.label,,drop=FALSE]

# number of simulations used to compute CI
nsim = 100L


## ---------- Models fit by penalized likelihood ------------- ##
# Perform a PCA to get the stacked PCA then...
pca <- prcomp(data[tree$tip.label,])

# fit a (multivariate) BM model on the PCA axes
fit_bm<- mvgls(pca$x~1, tree=tree, model="BM", method = "H&L") # replace tree by treeScaled for variable rates

# ancestral states for the fitted process - we use the ML estimate to start the sampling
ancestors <- ancestral(fit_bm)

## Simulate the datasets by sampling from the process
sim_data_bm <- lapply(1:nsim, function(i){
    uncorr_data <- sapply(1:ncol(ancestors), function(i){
      treeScaledbis = tree                                                              # replace by treeScaled for the variable rates model
      treeScaledbis$edge.length <- treeScaledbis$edge.length*diag(fit_bm$sigma$Pinv)[i] # we scale the branch length by the marginal variance
      conditional_sampling(treeScaledbis, c(pca$x[,i],ancestors[,i]))[1:Ntip(tree)]     # we sample the ancestral states and tips values
  })
})


# Use "phytools" called by mvMORPH to obtain the branching times
times<-max(nodeHeights(tree))-nodeHeights(tree)[match(1:tree$Nnode+Ntip(tree),tree$edge[,1]),1]

# time bins for the ancestral states
bin<-findInterval(tip_age$Tip_Date, rev(c(56,33.9,23.03,5.333,0.012,0)))
bin<-as.factor(bin)
levels(bin)<-c(80,60,40,20,0)
bin<-as.numeric(as.character(bin))

# for tips
bin2<-findInterval(tip_age$Tip_Date, rev(c(56,33.9,23.03,5.333,0.012,0)))
bin2<-as.factor(bin2)
levels(bin2)<-c(80,60,40,20,0)
bin2<-as.numeric(as.character(bin2))

# Generate a stacked PCA
library(scatterplot3d)
plot3d1<-scatterplot3d(x = pca$x[,1], y = pca$x[,2], z = bin2, pch =19, box=FALSE, label.tick.marks=F, angle=20,
                       xlab = "PC1", ylab="PC2", zlab = "Times")
plot3d1$plane3d(20,0,0,lty="solid",col="grey")
plot3d1$plane3d(40,0,0,lty="solid",col="grey")
plot3d1$plane3d(60,0,0,lty="solid",col="grey")
plot3d1$plane3d(80,0,0,lty="solid",col="grey")

# plot the simulated datasets on the same PCA
for(i in 1:nsim) plot3d1$points3d(x = sim_data_bm[[i]][,1], y = sim_data_bm[[i]][,2], z = bin, pch=19, col="red")
plot3d1$points3d(x = pca$x[,1], y = pca$x[,2], z = bin2, pch=19, col="black") # superimpose onto the original?







## ------------------------------------------------------------------------- ##
## Reconstruct a stacked PCA for the mammals dataset and with expectations   ##
## from a lambda model with variable rates model 		  											 ##	
## (like in Goswami et al. 2022)                                           	 ##	
## ------------------------------------------------------------------------- ##

# fit a (multivariate) lambda model on the variable rates tree and on the PCA axes
fit_bmm<- mvgls(pca$x~1, tree=treeScaled, model="lambda", method = "H&L") # replace tree by treeScaled for variable rates

# ancestral states for the fitted process - we use the ML estimate to start the sampling
ancestors <- ancestral(fit_bmm)

## Simulate the datasets by sampling from the process
sim_data_bm <- lapply(1:nsim, function(i){
    uncorr_data <- sapply(1:ncol(ancestors), function(i){
      treeScaledbis = treeScaled                                                         # replace by treeScaled for the variable rates model
      treeScaledbis$edge.length <- treeScaledbis$edge.length*diag(fit_bmm$sigma$Pinv)[i] # we scale the branch length by the marginal variance
      conditional_sampling(treeScaledbis, c(pca$x[,i],ancestors[,i]))[1:Ntip(tree)]      # we sample the ancestral states and tips values
  })
})


# Use "phytools" called by mvMORPH to obtain the branching times
times<-max(nodeHeights(tree))-nodeHeights(tree)[match(1:tree$Nnode+Ntip(tree),tree$edge[,1]),1]

# time bins for the ancestral states
bin<-findInterval(tip_age$Tip_Date, rev(c(56,33.9,23.03,5.333,0.012,0)))
bin<-as.factor(bin)
levels(bin)<-c(80,60,40,20,0)
bin<-as.numeric(as.character(bin))

# for tips
bin2<-findInterval(tip_age$Tip_Date, rev(c(56,33.9,23.03,5.333,0.012,0)))
bin2<-as.factor(bin2)
levels(bin2)<-c(80,60,40,20,0)
bin2<-as.numeric(as.character(bin2))

# Generate a stacked PCA
library(scatterplot3d)
plot3d1<-scatterplot3d(x = pca$x[,1], y = pca$x[,2], z = bin2, pch =19, box=FALSE, label.tick.marks=F, angle=20,
                       xlab = "PC1", ylab="PC2", zlab = "Times")
plot3d1$plane3d(20,0,0,lty="solid",col="grey")
plot3d1$plane3d(40,0,0,lty="solid",col="grey")
plot3d1$plane3d(60,0,0,lty="solid",col="grey")
plot3d1$plane3d(80,0,0,lty="solid",col="grey")

# plot the simulated datasets on the same PCA
for(i in 1:nsim) plot3d1$points3d(x = sim_data_bm[[i]][,1], y = sim_data_bm[[i]][,2], z = bin, pch=19, col="green")
plot3d1$points3d(x = pca$x[,1], y = pca$x[,2], z = bin2, pch=19, col="black") # superimpose onto the original?

