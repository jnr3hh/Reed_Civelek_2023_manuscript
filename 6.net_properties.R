#6. Network properties



######################## Read in networks, basic properties ########################
setwd("~/Civelek Lab/RIMBANET/Final Net Data/REDO MET SEQ/STARandGTExNetsFinal/Newest")

library(splineTimeR)
library(bc3net)
library(igraph)
library(dplyr)
library(qgraph)

#prepare the RIMBANET output file by 
#1. replacing the the arrow -> with " " whitespace
#2. removing the "digraph {}"
#3. removing unconnected genes
#4. saving as a textfile

BN = read.table("BN_digraph_pruned_SFV.txt", sep = " ", fill = T)
BN = BN %>% dplyr::select(!V2)
IG = graph_from_edgelist(as.matrix(BN), directed = T)

# distance matrix
dist = distances(IG, v = V(IG), to = V(IG), mode = "out")


# scale free test
networkProperties(IG)
#small world lest
smallworldIndex(IG)











