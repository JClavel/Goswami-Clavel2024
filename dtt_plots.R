## ------------------------------------------------------------------------- ##
## Reconstruct the DTT for the mammal dataset and compare it to expectations ##
## from a BM, EB, OU, and variable rate model		  						 ##
##  																         ##
## ------------------------------------------------------------------------- ##
source("fun.R")
library(mvMORPH)

# phylogenetic tree
tree <- read.nexus("Mammals_whole_dataset/tree80_85_85.nex")

# tree with branch lengths scaled by the rates estimlated using bayesTraits
treeScaled = read.nexus("Mammals_whole_dataset/treeScaleRate.nex")

# Morphology - the shape data for mammals
data <- read.csv("Mammals_whole_dataset/shape.data.322.csv", row.names = 1)

# compute pPC scores as used in the original study (Goswami et al. 2022)
pPCA <- phyl.pca(tree, data)
dat_pc = pPCA$S[,1:67] # the 67 first PCs ~95% of the variance explained

# prepare the data for the 'mvgls' function and check that all tables are well ordered
dat <- list(Y=as.matrix(dat_pc[tree$tip.label,]) )

# time bin spacing used for the DTT plots
tbin = 1

# number of simulations used to compute CI
nsim = 100L

# Root age of the tree
max_time = max(nodeHeights(tree))

## ----------------- Models fit by penalized likelihood ------------------ ##
# Brownian motion
fit_bm<- mvgls(Y~1, data=dat, tree=tree, model="BM", method = "H&L")
# Ornstein-Uhlenbeck
fit_ou<- mvgls(Y~1, data=dat, tree=tree, model="OU", method = "PL")
# Early-Burst
fit_eb<- mvgls(Y~1, data=dat, tree=tree, model="EB", method = "H&L")
# Lambda with branch lengths scaled by rates - replicate the model fit in Goswami et al. 2022
fit_bmr<- mvgls(Y~1, data=dat, tree=treeScaled, model="lambda", method = "H&L")


## ---------- Compute the DTT on empirical and simulated data ------------ ##
par(mfrow=c(2,2))
plot_dtt_fit(fit_bm, tree, nsim=nsim, dist="euclidean", bin=tbin, col="red", add=FALSE,
             main="Brownian motion", reverse=TRUE, ylim=c(0,1.2), xaxt='n')
axis(1, seq(from=0, to=max_time, by=20), labels = seq(from=80, to=0, by=-20)) 
plot_dtt_fit(fit_ou, tree, nsim=nsim, dist="euclidean", bin=tbin, col="blue", add=FALSE,
             main="Ornstein-Uhlenbeck", reverse=TRUE, ylim=c(0,1.2), xaxt='n')
axis(1, seq(from=0, to=max_time, by=20), labels = seq(from=80, to=0, by=-20)) 
plot_dtt_fit(fit_eb, tree, nsim=nsim, dist="euclidean", bin=tbin, col="orange", add=FALSE,
             main="Early Burst", reverse=TRUE, ylim=c(0,1.2), xaxt='n')
axis(1, seq(from=0, to=max_time, by=20), labels = seq(from=80, to=0, by=-20)) 
plot_dtt_fit(fit_bmr, tree, nsim=nsim, dist="euclidean", bin=tbin, col="green", add=FALSE, 
             main="BM variable rates", reverse=TRUE,ylim=c(0,1.2), xaxt='n')
axis(1, seq(from=0, to=max_time, by=20), labels = seq(from=80, to=0, by=-20)) 

# # add legend to the plot
# legend("bottomleft", col=c("red","blue","orange","green"),
#        lty=1, legend=c("Brownian motion","Ornstein-Uhlenbeck", "Early Burst", "BM variable rates"))

# tweak EB to show the effect of various levels of Early burst
par(mfrow=c(1,1))
library(geiger)
plot_dtt_fit(fit_eb, tree, nsim=nsim, dist="euclidean", bin=tbin, col="orange",
             main="Early Burst", reverse=TRUE, ylim=c(0,1), xaxt='n')
axis(1, seq(from=0, to=max_time, by=20), labels = seq(from=80, to=0, by=-20)) 

# transform the tree to represent various half-lifes
fit_eb$corrSt$phy <- rescale(tree, model="EB", a=-log(2)/(max_time/2)) # 2 half-lives

plot_dtt_fit(fit_eb, tree, nsim=nsim, dist="euclidean", bin=tbin, col="pink", add=TRUE, reverse=TRUE)

fit_eb$corrSt$phy <- rescale(tree, model="EB", a=-log(2)/(max_time/5)) # 5 half-lives

plot_dtt_fit(fit_eb, tree, nsim=nsim, dist="euclidean", bin=tbin, col="purple", add=TRUE, reverse=TRUE)

fit_eb$corrSt$phy <- rescale(tree, model="EB", a=-log(2)/(max_time/10)) # 10 half-lives

plot_dtt_fit(fit_eb, tree, nsim=nsim, dist="euclidean", bin=tbin, col="grey", add=TRUE, reverse=TRUE)


legend("bottomleft", col=c("orange","pink","purple","grey"),
       lty=1, legend=c("BM/EB r=0","EB with 2 half-lives","EB with 5 half-lives", "EB with 10 half-lives"))
