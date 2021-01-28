## Cleaning up the data
library(phytools)
library(tidyverse)

###############################
##      Load the trees       ##
###############################
## mammal tree from Faurby
mam_tree <- read.tree("faurby_tree_1.tre")

## helminth tree
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
## But this is actually due to some naming issues that are easily sorted
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

## prune the mammal tree to include only the species from the GMPD
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

###################################
##          Nearctic Data        ##
###################################

## Reload the trees so we start with the full phylogenies and then trim to the species in the dataset
## mammal tree from Faurby
mam_tree <- read.tree("faurby_tree_1.tre")
## helminth tree
para_tree<-read.tree("helminth_tree.trees") #249 taxa

#host-parasite association data
data<-read.csv("Host_Range_Nearctic_Mammals_list.csv", header=T)
#change to characters instead of factors
data$Host<-as.character(nearctic_assoc_data$Host)
data$Parasite<-as.character(nearctic_assoc_data$Parasite)

## which parasites in the Nearctic dataset are not found in our phylogeny?
setdiff(gsub(" ", "_", data$Parasite),para_tree$tip.label)
## remove this species from the data
data <- data[-which(data$Parasite=="Ascaris lumbricoides"),]

## which parasites in the phylogeny are not found in the Nearctic dataset?
setdiff(para_tree$tip.label, gsub(" ","_",data$Parasite))
## remove these species from the tree
para_tree <- drop.tip(para_tree, setdiff(para_tree$tip.label, gsub(" ","_",data$Parasite)))

## verify that we have concordance between the parasite phylogeny and the parasite data
unique(gsub(" ","_",data$Parasite))%in%para_tree$tip.label %>% all
para_tree$tip.label%in%unique(gsub(" ","_",data$Parasite)) %>% all

## which mammals in the Nearctic dataset are not found in the phylogeny?
setdiff(gsub(" ", "_", data$Host),mam_tree$tip.label)
## rename American mink in the data
gsub("Mustela vison", "Neovison vison", data$Host) -> data$Host
## rename Alces alces (Eurasian elk) to Alces americanus (American moose) - these are sister species (or even subspecies) so any phylogenetic conclusions are unaltered in an analysis containing only one of the two species
gsub("Alces americanus", "Alces alces", data$Host) -> data$Host
## rename Arctic fox in the data
gsub("Alopex lagopus", "Vulpes lagopus", data$Host) -> data$Host
## rename American porcupine in the data
gsub("Erethizon dorsatus", "Erethizon dorsatum", data$Host) -> data$Host

## now prune the phylogeny to only include the species that are found the data
mam_tree <- keep.tip(mam_tree, unique(gsub(" ", "_", data$Host)))

## verify that we have concordance between the host phylogeny and the host data
unique(gsub(" ","_",data$Host))%in%mam_tree$tip.label %>% all
mam_tree$tip.label%in%unique(gsub(" ","_",data$Host)) %>% all

## how many hosts and parasites in this new dataset?
length(unique(data$Host)) ## 69 hosts
length(unique(data$Parasite)) ## 90 parasites

## save all of these new things
write.tree(mam_tree, file="nearctic_mammal_tree_clean.tre")
write.tree(para_tree, file="nearctic_helminth_tree_clean.tre")
write.csv(data, "nearctic_data_clean.csv")

