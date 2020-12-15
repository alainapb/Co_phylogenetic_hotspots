## Cleaning up the data
library(phytools)
library(tidyverse)

###############################
##      Load the trees       ##
###############################
## mammal tree from Faurby
mam_tree <- read.tree("faurby_tree_1.tre")

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

## verify that we have concordance between the parasite phylogeny and the parasite data in the GMPD
unique(gsub(" ","_",data$ParasiteCorrectedName))%in%para_tree$tip.label %>% all
para_tree$tip.label%in%unique(gsub(" ","_",data$ParasiteCorrectedName)) %>% all
## what's in the tree that's not in the database?
para_tree$tip.label[which(!(para_tree$tip.label%in%unique(gsub(" ","_",data$ParasiteCorrectedName))))] ## Schistosoma bovis
## drop this tip
para_tree <- drop.tip(para_tree, "Schistosoma_bovis")

## verify that we have concordance between the parasite phylogeny and the parasite data in the GMPD
unique(gsub(" ","_",data$ParasiteCorrectedName))%in%para_tree$tip.label %>% all
para_tree$tip.label%in%unique(gsub(" ","_",data$ParasiteCorrectedName)) %>% all

## Are all hosts in the GMPD found in the mammal phylogeny?
## No. 
unique(data$HostCorrectedName)[which(sapply(gsub(" ", "_", data$HostCorrectedName) %>% unique,
                                            function(s) !(s%in%mam_tree$tip.label)))]
## Do some renaming
gsub("Lama glama", "Lama guanicoe", data$HostCorrectedName) -> data$HostCorrectedName 
gsub("Equus asinus", "Equus hemionus", data$HostCorrectedName) -> data$HostCorrectedName 
gsub("Bos frontalis", "Bos gaurus", data$HostCorrectedName) -> data$HostCorrectedName 
gsub("Monachus schauinslandi", "Neomonachus schauinslandi", data$HostCorrectedName) -> data$HostCorrectedName 
gsub("Neotragus moschatus", "Neotragus pygmaeus", data$HostCorrectedName) -> data$HostCorrectedName 
gsub("Ovis aries", "Ovis orientalis", data$HostCorrectedName) -> data$HostCorrectedName 

## Drop this one remaining species from the GMPD
data[-which(data$HostCorrectedName=="Taurotragus oryx"),] -> data

## verify that we STILL have concordance between the parasite phylogeny and the parasite data in the GMPD
unique(gsub(" ","_",data$ParasiteCorrectedName))%in%para_tree$tip.label %>% all
para_tree$tip.label%in%unique(gsub(" ","_",data$ParasiteCorrectedName)) %>% all

## prune the tree to include only the species from the GMPD
mam_tree <- keep.tip(mam_tree, gsub(" ","_",unique(data$HostCorrectedName)))

## verify that we have concordance between the host phylogeny and the host data in the GMPD
unique(gsub(" ","_",data$HostCorrectedName))%in%mam_tree$tip.label %>% all
mam_tree$tip.label%in%unique(gsub(" ","_",data$HostCorrectedName)) %>% all

## how many hosts and parasites in this new dataset?
length(unique(data$HostCorrectedName)) ## 197 hosts
length(unique(data$ParasiteCorrectedName)) ## 248 parasites

## save all of these new things
write.tree(mam_tree, file="mammal_tree_clean.tre")
write.tree(para_tree, file="helminth_tree_clean.tre")
write.csv(data, "GMPD_clean.csv")
