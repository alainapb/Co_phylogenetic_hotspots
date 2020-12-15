## This script runs the cophylogenetic methods of Hadfield et al. 2014 on the 
## host and parasite phylogenies for the GMPD data

# library(asreml) ## this actually requires that I buy a piece of software; however, 
## Hadfield et al. 2014 used this package to verify the results obtained via MCMCglmm,
## so it is not strictly necessary
library(MCMCglmm)
library(gdata)
library(igraph)
library(phytools)
library(tidyverse)

###############################
##      Load the trees       ##
###############################
mam_tree <- read.tree("mammal_tree_clean.tre")
para_tree <- read.tree("helminth_tree_clean.tre")

#####################################################
## Form S^{-1} for each term using the phylogenies ##
#####################################################

## S^{-1} is a phylogenetic covariance matrix with ancestral nodes retained 
paraA <-inverseA(para_tree)$Ainv                                      # parasite main effect
mamA <-inverseA(mam_tree)$Ainv                                        # host main effect
mam.paraA <-as(kronecker(mamA, paraA), "dgCMatrix")                   # coevolutionary effect
mam.paraAS <-as(kronecker(mamA, Diagonal(nrow(paraA))), "dgCMatrix")  # host evolutionary effect
mam.paraSA <-as(kronecker(Diagonal(nrow(mamA)), paraA), "dgCMatrix")  # parasite evolutionary effect

rownames(mam.paraA)<-apply(expand.grid(rownames(paraA), rownames(mamA)), 1, function(x){paste(x[2],x[1], sep=".")})
rownames(mam.paraAS)<-rownames(mam.paraSA)<-rownames(mam.paraA)

###############################
##  Format the GMPD Data     ##
###############################
data <- read.csv("GMPD_clean.csv")

# ## verify that we have concordance between the host phylogeny and the host data in the GMPD
# unique(gsub(" ","_",data$HostCorrectedName))%in%mam_tree$tip.label %>% all
# mam_tree$tip.label%in%unique(gsub(" ","_",data$HostCorrectedName)) %>% all
# 
# ## verify that we have concordance between the parasite phylogeny and the parasite data in the GMPD
# unique(gsub(" ","_",data$ParasiteCorrectedName))%in%para_tree$tip.label %>% all
# para_tree$tip.label%in%unique(gsub(" ","_",data$ParasiteCorrectedName)) %>% all

## create a new data.frame with all possible host-parasite combinations
expand.grid(Host.species=(data$HostCorrectedName %>% unique),
            Parasite.species=(data$ParasiteCorrectedName %>% unique)) -> ndata
## columns for each of the phylogenetic and non-phylogenetic effects
ndata$Host.species <- gsub(" ", "_", ndata$Host.species)                            # phylogenetic main effect for hosts
ndata$Parasite.species <- gsub(" ", "_", ndata$Parasite.species)                    # phylogenetic main effect for parasites
ndata$Parasite.species.ide<-ndata$Parasite.species                                  # non-phylogenetic main effect for parasites
ndata$Host.species.ide<-ndata$Host.species                                          # non-phylogenetic main effect for hosts
ndata$Host.Parasite<-paste(ndata$Host.species, ndata$Parasite.species, sep=".")      # phylogenetic coevolutionary effect
##ndata$Host.Parasite.ide<-paste(ndata$Host.species, ndata$Parasite.species, sep=".")  # non-phylogenetic interaction effect
## We actually cannot measure a non-phylogenetic interaction effect because we don't have replicate incidence matrices
ndata$Host.Parasite.ide2<-paste(ndata$Host.species, ndata$Parasite.species, sep=".") # phylogenetic host evolutionary effect
ndata$Host.Parasite.ide3<-paste(ndata$Host.species, ndata$Parasite.species, sep=".") # phylogenetic parasite evolutionary effect

## Based on the code below, there are 892 unique host-parasite associations in the GMPD 
## paste(gsub(" ", "_", data$HostCorrectedName), gsub(" ", "_", data$ParasiteCorrectedName), sep=".") %>% unique %>% length
##
## create a column called 'presence' that is '1' if the host-parasite combination
## exists in the GMPD, and '0' otherwise
## create a column called 'ncount' that counts the number of times each association 
## is recorded in the GMPD - this could potentially allow us to control for sampling later?
## we can also account for sampling by counting the number of times each host and parasite shows up in the dataset
## create a column called 'nhosts.sampled' that counts how many times each *host* species show up in the dataset
## create a column called 'nparas.sampled' that counts how many times each *parasite* species shows up in the dataset
ndata$presence <- ndata$ncount <- ndata$nhosts.sampled <- ndata$nparas.sampled <- 0
for (i in 1:nrow(data)) {
  this.host <- gsub(" ", "_", data$HostCorrectedName[i])
  this.parasite <- gsub(" ", "_", data$ParasiteCorrectedName[i])
  ndata$presence[which(ndata$Host.species==this.host & ndata$Parasite.species==this.parasite)] <- 1
  ndata$ncount[which(ndata$Host.species==this.host & ndata$Parasite.species==this.parasite)] <- ndata$ncount[which(ndata$Host.species==this.host & ndata$Parasite.species==this.parasite)] + 1
}
ndata$nhosts.sampled <- sapply(ndata$Host.species, function(h) sum(gsub(" ", "_", data$HostCorrectedName)==h))
ndata$nparas.sampled <- sapply(ndata$Parasite.species, function(p) sum(gsub(" ", "_", data$ParasiteCorrectedName)==p))


################
## Run models ## 
################

priorI=list(R=list(V=1, fix=1))

## create parameter expanded priors for each variance component
## There are 7 such components here: 
## Parasite.species -> the phylogenetic main effect of parasites
## Host.species -> the phylogenetic main effect of hosts
## Parasite.species.ide -> the nonphylogenetic main effect of parasites
## Host.species.ide -> the nonphylogenetic main effect of hosts
## Host.parasite -> the phylogenetic coevolutionary effect
## Host.parasite.ide2 -> the phylogenetic host evolutionary effect
## Host.parasite.ide3 -> the phylogenetic parasite evolutionary effect
priorI$G<-lapply(1:7, function(x){list(V=1, nu=1, alpha.mu=0, alpha.V=1000)})
names(priorI$G)<-paste("G", 1:7, sep="")

## main model: includes controls for sampling effort
mI.MCMCa<-MCMCglmm(presence~log(nhosts.sampled)+log(nparas.sampled), random=~Parasite.species+Host.species+Parasite.species.ide+Host.species.ide+Host.Parasite+Host.Parasite.ide2+Host.Parasite.ide3,family="categorical", data=ndata, ginverse=list(Parasite.species=paraA, Host.species=mamA, Host.Parasite=mam.paraA, Host.Parasite.ide2=mam.paraAS, Host.Parasite.ide3=mam.paraSA), prior=priorI, slice=T, nitt=400000, thin=400, burnin=100000)
save(mI.MCMCa, file=paste(paste(results.filepath, "mI.MCMCa", sep=.Platform$file.sep), ptree, format(Sys.time(), "%d-%m-%Y"), "R", sep="."))

## model b: does not control for sampling effort
mI.MCMCb<-MCMCglmm(presence~1, random=~Parasite.species+Host.species+Parasite.species.ide+Host.species.ide+Host.Parasite+Host.Parasite.ide2+Host.Parasite.ide3,family="categorical", data=ndata, ginverse=list(Parasite.species=paraA, Host.species=mamA, Host.Parasite=mam.paraA, Host.Parasite.ide2=mam.paraAS, Host.Parasite.ide3=mam.paraSA), prior=priorI, slice=T, nitt=400000, thin=400, burnin=100000)
save(mI.MCMCb, file=paste(paste(results.filepath, "mI.MCMCb", sep=.Platform$file.sep), ptree, format(Sys.time(), "%d-%m-%Y"), "R", sep="."))

## model c: does not control for sampling effort, and does not try to account for non-phylogenetic effects (trying to speed things up a bit)
priorI2=list(R=list(V=1, fix=1))
priorI2$G<-lapply(1:5, function(x){list(V=1, nu=1, alpha.mu=0, alpha.V=1000)})
names(priorI2$G)<-paste("G", 1:5, sep="")
mI.MCMCc<-MCMCglmm(presence~1, random=~Parasite.species+Host.species+Host.Parasite+Host.Parasite.ide2+Host.Parasite.ide3,family="categorical", data=ndata, ginverse=list(Parasite.species=paraA, Host.species=mamA, Host.Parasite=mam.paraA, Host.Parasite.ide2=mam.paraAS, Host.Parasite.ide3=mam.paraSA), prior=priorI2, slice=T, nitt=400000, thin=400, burnin=100000)


###################################################
## Run module on data (from Krasnov et al. 2012) ##
###################################################

# A matrices for host and parasite
ht<-vcv(mam_tree, corr=T)
pt<-vcv(para_tree, corr=T)

# incidence data
Y<-table(ndata$Host.species, ndata$Parasite.species, ndata$presence)[,,2]>0

hi<-match(rownames(Y), rownames(ht))
pi<-match(colnames(Y), rownames(pt))

# reorder A matrices so match row/columns of Y
ht<-ht[hi,hi]
pt<-pt[pi,pi]

# assign species to modules
G<-graph.incidence(Y)
MU<-walktrap.community(G)
MU<-membership(MU)

lth<-which(lower.tri(ht))

# correlation between comembership and phylogenetic distance of hosts
chk<-cor(outer(MU[1:length(hi)], MU[1:length(hi)], "==")[lth], c(1-ht)[lth])

ltp<-lower.tri(pt)

cpk<-cor(outer(MU[length(hi)+1:length(pi)], MU[length(hi)+1:length(pi)], "==")[ltp], c(1-pt)[ltp])
# correlation between comembership and phylogenetic distance of parasites

shk.l<-1:1000
spk.l<-1:1000
# storing host & parasite membership/phylogenetic-distance correlations after Legendre Permutations 
shk.h<-1:1000
spk.h<-1:1000
# storing host & parasite membership/phylogenetic-distance correlations after Hommola Permutations 

for(i in 1:1000){
  
  Y2<-apply(Y, 2, sample)                     # Legendre sampling
  
  G2<-graph.incidence(Y2)
  MU2<-walktrap.community(G2)
  MU2<-membership(MU2)
  
  shk.l[i]<-cor(outer(MU2[1:length(hi)], MU2[1:length(hi)], "==")[lth], c(1-ht)[lth])
  spk.l[i]<-cor(outer(MU[length(hi)+1:length(pi)], MU[length(hi)+1:length(pi)], "==")[ltp], c(1-pt)[ltp])
  
  Y2<-Y[sample(1:nrow(Y)),sample(1:ncol(Y))]  # Hommola sampling
  
  G2<-graph.incidence(Y2)
  MU2<-walktrap.community(G2)
  MU2<-membership(MU2)
  
  shk.h[i]<-cor(outer(MU2[1:length(hi)], MU2[1:length(hi)], "==")[lth], c(1-ht)[lth])
  spk.h[i]<-cor(outer(MU[length(hi)+1:length(pi)], MU[length(hi)+1:length(pi)], "==")[ltp], c(1-pt)[ltp])
  print(i)
}

k.tests<-cbind(c(chk, cpk), c(sum(chk>shk.l)/1000, sum(cpk>spk.l)/1000), c(sum(chk>shk.h)/1000, sum(cpk>spk.h)/1000))
# store metrics, and the proportion of times the metrics under permutation were greater 

rownames(k.tests)<-c("H", "P")
colnames(k.tests)<-c("statistic", "L-pval", "H-pval")

saveRDS(k.tests, file="k.tests.RDS")

######################################
## Run parafit on consolidated data ##
######################################

## make sure you have a compiled Hommola.so available in this directory before running hommola()
source("Hommola.R")

# A matrices for host and parasite
ht<-vcv(mam_tree, corr=T)
pt<-vcv(para_tree, corr=T)

Y<-table(ndata$Host.species, ndata$Parasite.species, ndata$presence)[,,2]>0
# incidence data

hi<-match(rownames(Y), rownames(ht))
pi<-match(colnames(Y), rownames(pt))
  
ht<-ht[hi,hi]
pt<-pt[pi,pi]
# reorder A matrices so match row/columns of Y

htp<-pcoa(1-ht)$vectors
ptp<-pcoa(1-pt)$vectors
# get principal coordinate of phylogenetic distance matrix

hte<-t(t(eigen(solve(ht))$vectors)*sqrt(eigen(solve(ht))$values))
pte<-t(t(eigen(solve(pt))$vectors)*sqrt(eigen(solve(pt))$values))

# get unnormalised eigenvectors of A
lD<-t(htp)%*%Y%*%ptp
iD<-t(hte)%*%(Y-mean(Y))%*%pte
# D matrices (see Appendix) for Legendre and Ives methods
  
cl<-sum(diag(t(lD)%*%lD))
ci<-sum(diag(t(iD)%*%iD))
ch<-hommola(Y, ht,pt)
# metrics (see Appendix) of Legendre, Ives & Hommola

sl.l<-1:1000
si.l<-1:1000
sh.l<-1:1000
# stoing metrics of Legendre, Ives & Hommola after Legendre permutations  
sl.h<-1:1000
si.h<-1:1000
sh.h<-1:1000
# stoing metrics of Legendre, Ives & Hommola after Hommola permutations  

##################
## Permutations ##
##################
## Compute the Legendre, Ives, and Hommola metrics for randomly permuted 
## incidence data. We do this twice, using different methods for performing
## the permutation. The Legendre permutation simply generates a random permutation 
## of the incidence data columnwise (e.g., since columns represent parasites, in the 
## permuted data, each parasite infects the same number of hosts, but the identity
## of those hosts has been randomly determined). The Hommola permutation
for(j in 1:1000){
  print(j)
  
  Y2<-apply(Y, 2, sample)                     # Legendre permtations
  
  lD<-t(htp)%*%Y2%*%ptp                       # Legendre D matrix
  iD<-t(hte)%*%(Y2-mean(Y2))%*%pte            # Ives D matrix
  
  sl.l[j]<-sum(diag(t(lD)%*%lD))              # Legendre metric (ParafitGlobal)
  si.l[j]<-sum(diag(t(iD)%*%iD))              # Ives metric (MSEb)
  sh.l[j]<-hommola(Y2, ht,pt)                 # Hommola metric
  
  Y2<-Y[sample(1:nrow(Y)),sample(1:ncol(Y))]  # Hommola sampling
  
  lD<-t(htp)%*%Y2%*%ptp
  iD<-t(hte)%*%(Y2-mean(Y2))%*%pte
  
  sl.h[j]<-sum(diag(t(lD)%*%lD)) 
  si.h[j]<-sum(diag(t(iD)%*%iD)) 
  sh.h[j]<-hommola(Y2, ht,pt)    
}

other.tests<-cbind(c(cl, ci, ch), c(sum(cl<sl.l)/1000, sum(ci>si.l)/1000, sum(ch<sh.l)/1000), c(sum(cl<sl.h)/1000, sum(ci>si.h)/1000, sum(ch<sh.h)/1000))
# store metrics, and the proportion of times the metrics under permutation were greater 

rownames(other.tests)<-c("L", "I", "H")
colnames(other.tests)<-c("statistic", "L-pval", "H-pval")

saveRDS(other.tests, file="other.tests.RDS")

## Another analysis to potentially run: we can compare the host ranges of two parasites
## by computing the phylogenetic distance between the hosts of two parasites. The smaller
## this distance, the more similar the host ranges of the two parasites. To identify 
## putative cases of "niche partitioning", we find parasites with very low PD between
## host ranges, but whose hosts don't overlap at all (e.g., the parasitize related, but
## not identical, hosts). You can do a similar analysis, comparing the PD between the 
## parasites that infect different hosts