# 8. pathway enrichments within downstream geneset
library(topGO)
library(org.Hs.eg.db)

keydrivers = read.table("53kds.txt", header = F)
keydrivers = keydrivers$V1

sigRes <- data.frame(GO.ID=character(),
                    Term=character(), 
                    Annotated=numeric(),
                    Significant=numeric(),
                    Expected=numeric(),
                    classicFisher=numeric(),
                    classicFisherA = numeric(),
                    mod=character(),
                    stringsAsFactors=FALSE) 


# run once for each network, switch the file names
BN = read.table("BN_digraph_pruned_GMV.txt", sep = " ", fill = T)
BN = BN %>% dplyr::select(!V2)
IG = graph_from_edgelist(as.matrix(BN), directed = T)

# distance matrix
dist = distances(IG, v = V(IG), to = V(IG), mode = "out")

geneNames = unionBFD

#downstreamGenes = Moccur3
#module = 'M in 3/4 nets'
intModules= keydrivers
i = 1
print(i)
for (module in intModules)
{
  print(i)
  print(module)
  i = i+1
  downstream = dist[rownames(dist)== module]
  genes = rownames(dist)
  downstreamGenes = genes[downstream <= 4 & downstream != 0]
  if(length(downstreamGenes != 0)){
    
    genesInTargetModule = downstreamGenes
    myInterestingGenes <- (t(genesInTargetModule)) # int genes ( all the genes in each modules)
    geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
    names(geneList) <- geneNames # gene names and 1 or 0 in their cell s
    
    GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "symbol")
    #change ID = "ensembl" if using ens ids
    
    resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    resultFisher
    allGO = usedGO(object = GOdata) 
    allRes <- GenTable(GOdata, classicFisher = resultFisher,ranksOf = "classicFisher", topNodes = length(allGO))
    a <- p.adjust(allRes$classicFisher, method = "fdr", n = length(allRes$classicFisher))
    allRes$classicFisherA <- a 
    # GO term result file 
    
    sig = allRes[allRes$classicFisherA <= 0.05,]
    if(nrow(sig) != 0){
      
      sig$mod = module
      sigRes = rbind(sigRes, sig)
    }
    
    
  }
}
GMVGO = sigRes
write.table(GMVGO, 'GMV_dsGO_221.txt', quote = F, row.names = F, sep = '\t')


