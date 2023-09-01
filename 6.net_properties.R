#6. Network properties



######################## Read in networks, basic properties ########################
setwd("~/Civelek Lab/RIMBANET/Final Net Data/REDO MET SEQ/STARandGTExNetsFinal/Newest")

library(splineTimeR)
library(bc3net)
library(igraph)
library(dplyr)
library(qgraph)
BN = read.table("BN_digraph_pruned_SFV.txt", sep = " ", fill = T)
IG = graph_from_edgelist(as.matrix(BN[1:nrow(BN),]), directed = T)

# distance matrix
dist = distances(IG, v = V(IG), to = V(IG), mode = "out")


# scale free test
networkProperties(IG)
#small world lest
smallworldIndex(IG)











