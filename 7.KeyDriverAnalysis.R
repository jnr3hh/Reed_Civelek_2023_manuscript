# 7. Key driver analysis
library(splineTimeR)
library(bc3net)
library(igraph)
library(dplyr)



# run once for each network, switch the file names
BN = read.table("BN_digraph_pruned_SFV.txt", sep = " ", fill = T)
BN = BN %>% dplyr::select(!V2)
IG = graph_from_edgelist(as.matrix(BN), directed = T)
# distance matrix
dist = distances(IG, v = V(IG), to = V(IG), mode = "out")

######################

# downstream gene matrix
nHopCounts = data.frame(matrix(,nrow =nrow(dist),ncol =10))
s = 1
for(gene in 1:nrow(dist)){
  vect = dist[s,]
  vect = vect[!is.infinite(vect)]
  t = tabulate(vect, nbins = 10)
  nHopCounts[s,] = t
  s = s + 1
}

## what is the average number of downstream connections at each dist 1:10
nHopCounts = t(nHopCounts)
colnames(nHopCounts) = row.names(dist)
mean = as.vector(rowMeans(nHopCounts))
## what is one standard deviation away from the mean at each dist 1:10
sd = c(1:10)
for(r in c(1:10)){
  sd[r] = sd(nHopCounts[r,])
}

# weights nearby genes as more important
dropoff = c(1,1,1,.75,.5,.25,.125,.0625,.03125,0.015625)
#plot(dropoff)

nHopLogic = nHopCounts
for(n in c(1:nrow(nHopCounts))){
  for(m in c(1:ncol(nHopCounts))){
    #print(mean[n]+sd[n])
    if(nHopCounts[n,m]>(mean[n])){
      #print("m")
      #print(m)
      #print('n')
      #print(n)
      #print((nHopCounts[n,m]-mean[n])/sd[n])
      #print(avConperHop[n,m])
      #print(dropoff[n])
      #nHopLogic[n,m] = ((nHopCounts[n,m]-mean[n])/sd[n])*avConperHop[n,m]*dropoff[n]
      nHopLogic[n,m] = ((nHopCounts[n,m]-mean[n])/sd[n])*dropoff[n]
    }
    else
      nHopLogic[n,m] = 0
  }
}

# Find top ranked KDs
rank = rowSums(t(nHopLogic))
rank[is.na(rank)] = 0
gene = colnames(nHopCounts)
kds = as.character()
scores = as.character()
for(r in c(1:length(rank))){
  if(rank[r]>=0.0){
    kds = cbind(kds, gene[r])
    scores = c(scores, rank[r])
    
  }
}
l = names(scores)
scores = as.numeric(scores)

names(scores) = l
KDs = as.data.frame(t(kds))
KDs$score = scores



KDs = arrange(KDs, desc(score))
SFV_kd1 = KDs[1:(length(scores)/10),]
View(as.data.frame(SFV_kd1))
# change file names
write.table(SFV_kd1, 'SFVstructKDs.txt', quote = F, row.names = F)
write.table(SFV_kd1$V1, 'structSFV_topKDS_sub.txt', quote = F, row.names = F, col.names = F)

#################

KDs$highQual1 = 1
KDs$highQual2 = 1
KDs$highQual3 = 1
KDs$highQual4 = 1

#colocBFD = read.table("495GWASgenes.txt")

#colnames(colocBFD) = c('gene_symbol')

#nearestBFD = read.table("PCEMTsig.txt", header = T)
#nearestBFD = colocBFD
#colnames(nearestBFD) = 'gene_symbol'

#new = read.table('newshort.txt', header = F)


#over1 = new[!is.na(match(new$V1, colocBFD$V1)),]
#over2 = new[!is.na(match(new$V1, nearestBFD$V1)),]


s=1
for(kd in KDs$V1){
  downstream = dist[rownames(dist)== kd]
  genes = rownames(dist)
  downstreamGenes = genes[downstream <= 1 & downstream!= 0]
  GWASgenes = unique(c((as.character(nearestBFD$genes)), as.character(colocBFD$gene_symbol)))
  GWAS <- list(GWASg = GWASgenes)
  enrich = enrichment(downstreamGenes, genes, GWAS, adj = "fdr", verbose = FALSE)
  KDs$highQual1[KDs$V1 == kd] = enrich$padj
  s=s+1
  print(s)
}
for(kd in KDs$V1){
  downstream = dist[rownames(dist)== kd]
  genes = rownames(dist)
  downstreamGenes = genes[downstream <= 2 & downstream!= 0]
  GWASgenes = unique(c((as.character(nearestBFD$genes)), as.character(colocBFD$gene_symbol)))
  GWAS <- list(GWASg = GWASgenes)
  enrich = enrichment(downstreamGenes, genes, GWAS, adj = "fdr", verbose = FALSE)
  KDs$highQual2[KDs$V1 == kd] = enrich$padj
}
for(kd in KDs$V1){
  downstream = dist[rownames(dist)== kd]
  genes = rownames(dist)
  downstreamGenes = genes[downstream <= 3 & downstream!= 0]
  GWASgenes = unique(c((as.character(nearestBFD$genes)), as.character(colocBFD$gene_symbol)))
  GWAS <- list(GWASg = GWASgenes)
  enrich = enrichment(downstreamGenes, genes, GWAS, adj = "fdr", verbose = FALSE)
  KDs$highQual3[KDs$V1 == kd] = enrich$padj
}
for(kd in KDs$V1){
  downstream = dist[rownames(dist)== kd]
  genes = rownames(dist)
  downstreamGenes = genes[downstream <= 4 & downstream!= 0]
  GWASgenes = unique(c((as.character(nearestBFD$genes)), as.character(colocBFD$gene_symbol)))
  GWAS <- list(GWASg = GWASgenes)
  enrich = enrichment(downstreamGenes, genes, GWAS, adj = "fdr", verbose = FALSE)
  KDs$highQual4[KDs$V1 == kd] = enrich$padj
}


enrich1 = KDs$V1[KDs$highQual1<0.05]
enrich2 = KDs$V1[KDs$highQual2<0.05]
enrich3 = KDs$V1[KDs$highQual3<0.05]
enrich4 = KDs$V1[KDs$highQual4<0.05]

SFVenrichedKDs = unique(c(enrich1,enrich2,enrich3,enrich4))
#change variable names
write.table(SFVenrichedKDs, 'enrichedSFV.txt', quote = F, row.names = F)

















################## run this part once

GMSQenrichedKDs = read.table("enrichedGMSQ.txt", header = T)
GMSQKDshigh = read.table("structGMSQ_topKDS_sub.txt", header = F)
SMSQenrichedKDs = read.table("enrichedSMSQ.txt", header = T)
SMSQKDshigh = read.table("structSMSQ_topKDS_sub.txt", header = F)
GMVenrichedKDs = read.table("enrichedGMV.txt", header = T)
GMVKDshigh = read.table("structGMV_topKDS_sub.txt", header = F)
SMVenrichedKDs = read.table("enrichedSMV.txt", header = T)
SMVKDshigh = read.table("structSMV_topKDS_sub.txt", header = F)

GFSQenrichedKDs = read.table("enrichedGFSQ.txt", header = T)
GFSQKDshigh = read.table("structGFSQ_topKDS_sub.txt", header = F)
SFSQenrichedKDs = read.table("enrichedSFSQ.txt", header = T)
SFSQKDshigh = read.table("structSFSQ_topKDS_sub.txt", header = F)
GFVenrichedKDs = read.table("enrichedGFV.txt", header = T)
GFVKDshigh = read.table("structGFV_topKDS_sub.txt", header = F)
SFVenrichedKDs = read.table("enrichedSFV.txt", header = T)
SFVKDshigh = read.table("structSFV_topKDS_sub.txt", header = F)


GMSQ_total_KDs = unique(c(GMSQenrichedKDs$x,GMSQKDshigh$V1))
SMSQ_total_KDs = unique(c(SMSQenrichedKDs$x,SMSQKDshigh$V1))
GMV_total_KDs = unique(c(GMVenrichedKDs$x,GMVKDshigh$V1))
SMV_total_KDs = unique(c(SMVenrichedKDs$x,SMVKDshigh$V1))

GFSQ_total_KDs = unique(c(GFSQenrichedKDs$x,GFSQKDshigh$V1))
SFSQ_total_KDs = unique(c(SFSQenrichedKDs$x,SFSQKDshigh$V1))
GFV_total_KDs = unique(c(GFVenrichedKDs$x,GFVKDshigh$V1))
SFV_total_KDs = unique(c(SFVenrichedKDs$x,SFVKDshigh$V1))


#KD overlap

overlapMV = GMV_total_KDs[!is.na(match(GMV_total_KDs,SMV_total_KDs))]
overlapMSQ = GMSQ_total_KDs[!is.na(match(GMSQ_total_KDs,SMSQ_total_KDs))]
overlapFV = GFV_total_KDs[!is.na(match(GFV_total_KDs,SFV_total_KDs))]
overlapFSQ = GFSQ_total_KDs[!is.na(match(GFSQ_total_KDs,SFSQ_total_KDs))]




