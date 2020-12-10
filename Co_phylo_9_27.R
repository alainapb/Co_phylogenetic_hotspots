library(ape)
library(plyr)
library(tidyverse)
library(ouch)
library(dplyr)
library(phytools)
library(geiger)
library(caper)
library(picante)
library(igraph)

##load trees
mammal_tree<-read.tree("mammal_w_humans.trees") #191 taxa (includes homo sapiens)

mam_tree<-read.tree("mam_tree.trees") #190 taxa

para_tree<-read.tree("helminth_tree.trees") #249 taxa

##load association data
GMPDsmall<-read.csv("GMPDsmall.csv")

#interspecific distance matrix of parasites
dis_para<-cophenetic.phylo(para_tree) 

#intespecific distance matrix of mammals 
dis_mam<-cophenetic.phylo(mam_tree)

##make the pipeline to get the absense/presence of 
##each parasite (rows) that infect a given host (column)

mamlist<-as.vector(mam_tree$tip.label) #list of mammals in tip label order

paralist<-as.vector(para_tree$tip.label) #list of parasites in tip label order

#make a dataframe with host names as columns and parasite names as rows
df <- data.frame(matrix(ncol=length(mamlist),nrow=length(paralist)), row.names = paralist)
colnames(df)<-mamlist

# complete dataframe with association data, 
# where 0 = not found in host and 1 = found in host

for (i in 1:nrow(df))
 df[i,] <- (colnames(df) %in% unique(subset(GMPDsmall, ParasiteCorrectedName==rownames(df)[i])$HostCorrectedName)) %>% as.numeric

#save matrix
write.csv(df,file = "association_matrix_co_phylo.csv", row.names = TRUE, col.names = TRUE)

## Create MPD scores for every parasite using the Phylacine tree
para.mpd <- mpd(df, dis_mam)
para.mpd[is.na(para.mpd)] <- 0 ## one host parasites have an mpd of 0
names(para.mpd)<-paralist


## Create MPD scores for every host using the helminth tree
mam.mpd <- mpd(t(df), dis_para)
mam.mpd[is.na(mam.mpd)] <- 0 ## hosts infected by only one parasite have MPD of 0
names(mam.mpd)<-mamlist


## Calculate PD for each mammal to humans

#calculate the pairwise distances
mam_dis<-cophenetic.phylo(mammal_tree)

#pull out the PD for each mammal to humans
mam_dis = as.data.frame(mam_dis)
distances<-mam_dis[,"Homo sapiens"]
names(distances)<-rownames(mam_dis)
view(distances)

#remove humans from distances
distances<-distances[names(distances) != "Homo sapiens"]




##11/30/20 JD: phylo diversity of pathogens vs. host breadth of pathogens

##need to calculate the average para.mpd for each host

ggplot(aes(x=mam.mpd, y=))+
  geom_point()



##### Figure Options ####

#this makes a network association---too much data for this to be used
#ww<-graph_from_incidence_matrix(df, mode="out") #this makes a network association
#plot(ww, size=.2)

## MAIN FIGURE ##
##need a matrix of mam:para with the joined mpd to make this heatmap

heatmap(global_dis_mam) # just the heat map

heatmap(global_dis_mam, distance) #can add on the layer of host distance to humans as a second layering 



### plots trees against each other ###
cophyloplot(mam_tree,para_tree, assoc=df,length.line = 4, space = 28, gap = 3, fsize=0.3)

##not sure what this does...
phylo.heatmap(mam_tree, df, standardize = TRUE, fsize=0.03)

