## Cleaning up the data
library(phytools)
library(tidyverse)

###############################
##      Load the trees       ##
###############################
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
##          GMPD Data        ##
###############################
data <- read.csv("GMPDsmall.csv")

## Are all parasites in the GMPD found in the parasite phylogeny?
## No. There are 251 parasites in the GMPD. There are 249 species
## in the parasite phylogeny. The phylogeny has removed two species
## that were misplaced in the tree.
unique(data$ParasiteCorrectedName)[which(sapply(gsub(" ", "_", data$ParasiteCorrectedName) %>% unique, 
                                                function(s) !(s%in% para_tree$tip.label)))]
## remove these two parasites from the GMPD
data[-which(data$ParasiteCorrectedName%in%c("Prosthenorchis elegans","Filaria martis")),] -> data2

## Are all hosts in the GMPD found in the mammal phylogeny?
## No. There are 206 hosts in the GMPD (only 200 after removing all of the hosts of P. elegans and F. martis.) 
## There are only 190 species in the mammal phylogeny.
unique(data$HostCorrectedName) %>% length ## 206
unique(data2$HostCorrectedName) %>% length ## 200
mam_tree$tip.label %>% length ## 190
unique(data$HostCorrectedName)[which(sapply(gsub(" ", "_", data$HostCorrectedName) %>% unique,
                                            function(s) !(s%in%mam_tree$tip.label)))]
unique(data2$HostCorrectedName)[which(sapply(gsub(" ", "_", data2$HostCorrectedName) %>% unique,
                                             function(s) !(s%in%mam_tree$tip.label)))]
## however, it is the same set of hosts that are missing from both data and data2, meaning that data2 will be further reduced

## Because I'm not sure where this mammal tree came from, compare against some other mammal supertrees
## mammal tree from https://data.vertlife.org/ - note that there are actually 10,000 of these trees available (interested in the distribution of branch lengths, etc., across giant trees)
mtree <- read.tree("MamPhy_BDvr_Completed_5911sp_topoCons_FBDasZhouEtAl_v2_tree0000.tre")
## mammal tree from Faurby
mtree2 <- read.tree("faurby_tree_1.tre")
## let's see how many of the species in the GMPD show up in these three phylogenies
unique(gsub(" ","_",data2$HostCorrectedName))%in%mam_tree$tip.label %>% sum ## 184
unique(gsub(" ","_",data2$HostCorrectedName))%in%mtree$tip.label %>% sum ## 186
unique(gsub(" ","_",data2$HostCorrectedName))%in%mtree2$tip.label %>% sum ## 187

## so I have three different mammal phylogenies, each with a different number of taxa that can be found in the GMPD
## presumably there is a lot of overlap?
mam_tree$tip.label%in%mtree$tip.label ## not complete overlap here - 6 species in mam_tree that are not found in mtree
mam_tree$tip.label%in%mtree2$tip.label ## all of the species in mam_tree are found in the Faurby tree, but the Faurby tree contains 3 extra species

## So let's work with the Faurby tree
## prune the Faurby tree to the 187 mammals that are also found in the GMPD
keep.tip(mtree2, mtree2$tip.label[which(mtree2$tip.label%in%gsub(" ","_",data2$HostCorrectedName))]) -> mam_tree2

## prune the GMPD to only include these 187 mammals as well
data2[which(gsub(" ","_",data2$HostCorrectedName)%in%mam_tree2$tip.label),] -> data3

## verify that we have concordance between the host phylogeny and the host data in the GMPD
unique(gsub(" ","_",data3$HostCorrectedName))%in%mam_tree2$tip.label %>% all
mam_tree2$tip.label%in%unique(gsub(" ","_",data3$HostCorrectedName)) %>% all

## How many parasites remain in this pruned GMPD?
unique(data3$ParasiteCorrectedName) %>% length ## 229

## prune the parasite tree to only include these 229 parasite species
keep.tip(para_tree, gsub(" ","_",unique(data3$ParasiteCorrectedName))) -> para_tree2

## verify that we have concordance between the parasite phylogeny and the parasite data in the GMPD
unique(gsub(" ","_",data3$ParasiteCorrectedName))%in%para_tree2$tip.label %>% all
para_tree2$tip.label%in%unique(gsub(" ","_",data3$ParasiteCorrectedName)) %>% all

## save all of these new things
write.tree(mam_tree2, file="mammal_tree_clean.tre")
write.tree(para_tree2, file="helminth_tree_clean.tre")
write.csv(data3, "GMPD_clean.csv")
