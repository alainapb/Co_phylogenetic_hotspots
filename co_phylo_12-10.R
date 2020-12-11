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

##load trees
mammal_tree<-read.tree("mammal_w_humans.trees") #191 taxa (includes homo sapiens)

mam_tree<-read.tree("mam_tree.trees") #190 taxa

para_tree<-read.tree("helminth_tree.trees") #249 taxa
## need to make this tree ultrametric to use the method
## however, you can see that it is *very* close to ultrametric as is
library(ouch)
as(ape2ouch(para_tree),"data.frame") %>% tail(.,100)
## so while the code from Hadfield et al provides methods for making 
## a tree ultrametric, having tried those methods out, there are far
## too clunky to handle this situation, where I just need to round branch
## lengths up by the tiniest amount
## Liam Revell has a nice function for doing described on his website
## http://blog.phytools.org/2017/03/forceultrametric-method-for-ultrametric.html
## that is incorporated into the phytools packages
para_tree2 <- force.ultrametric(para_tree, method='extend')
## you can check to see how little this changes anything
is.ultrametric(para_tree2)
plot(para_tree, cex=0.2)
plot(para_tree2, cex=0.2)
para_tree <- para_tree2

###############################
## Form S^{-1} for each term ##
###############################

## S^{-1} is a phylogenetic covariance matrix with ancestral nodes retained 
paraA <-inverseA(para_tree)$Ainv                                        # parasite main effect
mamA <-inverseA(mam_tree)$Ainv                                                # host main effect
mam.paraA <-as(kronecker(mamA, paraA), "dgCMatrix")                   # coevolutionary effect
mam.paraAS <-as(kronecker(mamA, Diagonal(nrow(paraA))), "dgCMatrix")  # host evolutionary effect
mam.paraSA <-as(kronecker(Diagonal(nrow(mamA)), paraA), "dgCMatrix")  # parasite evolutionary effect

rownames(mam.paraA)<-apply(expand.grid(rownames(paraA), rownames(mamA)), 1, function(x){paste(x[2],x[1], sep=".")})
rownames(mam.paraAS)<-rownames(mam.paraSA)<-rownames(mam.paraA)

###############################
##          GMPD Data        ##
###############################
data <- read.csv("GMPDsmall.csv")


## Another analysis to potentially run: we can compare the host ranges of two parasites
## by computing the phylogenetic distance between the hosts of two parasites. The smaller
## this distance, the more similar the host ranges of the two parasites. To identify 
## putative cases of "niche partitioning", we find parasites with very low PD between
## host ranges, but whose hosts don't overlap at all (e.g., the parasitize related, but
## not identical, hosts). You can do a similar analysis, comparing the PD between the 
## parasites that infect different hosts