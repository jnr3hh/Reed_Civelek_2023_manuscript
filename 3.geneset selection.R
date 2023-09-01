#.3 geneset selection

setwd("~/Civelek Lab/RIMBANET/Final Net Data/REDO MET SEQ/STARandGTExNetsFinal/Newest")
library(WGCNA)

######################## read in iterWGCNA modules#################

modsGFV = read.table("GFV_merged-0.05-membership.txt", header = T)
eigengeneGFV = read.table("GFV_merged-0.05-eigengenes.txt", header = T)

modsGMV = read.table("GMV_merged-0.05-membership.txt", header = T)
eigengeneGMV = read.table("GMV_merged-0.05-eigengenes.txt", header = T)

modsSMV = read.table("SMV_merged-0.05-membership.txt", header = T)
eigengeneSMV = read.table("SMV_merged-0.05-eigengenes.txt", header = T)

modsSFV = read.table("SFV_merged-0.05-membership.txt", header = T)
eigengeneSFV = read.table("SFV_merged-0.05-eigengenes.txt", header = T)

modsGMSQ = read.table("GMSQ_merged-0.05-membership.txt", header = T)
eigengeneGMSQ = read.table("GMSQ_merged-0.05-eigengenes.txt", header = T)

modsGFSQ = read.table("GFSQ_merged-0.05-membership.txt", header = T)
eigengeneGFSQ = read.table("GFSQ_merged-0.05-eigengenes.txt", header = T)

modsSMSQ = read.table("SMSQ_merged-0.05-membership.txt", header = T)
eigengeneSMSQ = read.table("SMSQ_merged-0.05-eigengenes.txt", header = T)

modsSFSQ = read.table("SFSQ_merged-0.05-membership.txt", header = T)
eigengeneSFSQ = read.table("SFSQ_merged-0.05-eigengenes.txt", header = T)



######################## read in preprocessed gene expression files #################

inGFV = read.table('iterGFV.txt', header = T)
inGMV = read.table('iterGMV.txt', header = T)
inSFV = read.table('iterSFV.txt', header = T)
inSMV = read.table('iterSMV.txt', header = T)
inGFSQ = read.table('iterGFSQ.txt', header = T)
inGMSQ = read.table('iterGMSQ.txt', header = T)
inSFSQ = read.table('iterSFSQ.txt', header = T)
inSMSQ = read.table('iterSMSQ.txt', header = T)


######################## determine bayes input geneset ###########################

inModGMV = modsGMV$Gene[modsGMV$Module != "UNCLASSIFIED"]
inModSMV = modsSMV$Gene[modsSMV$Module != "UNCLASSIFIED"]
inModGFV = modsGFV$Gene[modsGFV$Module != "UNCLASSIFIED"]
inModSFV = modsSFV$Gene[modsSFV$Module != "UNCLASSIFIED"]
inModGMSQ = modsGMSQ$Gene[modsGMSQ$Module != "UNCLASSIFIED"]
inModSMSQ = modsSMSQ$Gene[modsSMSQ$Module != "UNCLASSIFIED"]
inModGFSQ = modsGFSQ$Gene[modsGFSQ$Module != "UNCLASSIFIED"]
inModSFSQ = modsSFSQ$Gene[modsSFSQ$Module != "UNCLASSIFIED"]

inModGMVs = modsGMV[modsGMV$Module != "UNCLASSIFIED",]
inModGMVs$Group = 'GMV'
inModSMVs = modsSMV[modsSMV$Module != "UNCLASSIFIED",]
inModSMVs$Group = 'SMV'
inModGFVs = modsGFV[modsGFV$Module != "UNCLASSIFIED",]
inModGFVs$Group = 'GFV'
inModSFVs = modsSFV[modsSFV$Module != "UNCLASSIFIED",]
inModSFVs$Group = 'SFV'
inModGMSQs = modsGMSQ[modsGMSQ$Module != "UNCLASSIFIED",]
inModGMSQs$Group = 'GMSQ'
inModSMSQs = modsSMSQ[modsSMSQ$Module != "UNCLASSIFIED",]
inModSMSQs$Group = 'SMSQ'
inModGFSQs = modsGFSQ[modsGFSQ$Module != "UNCLASSIFIED",]
inModGFSQs$Group = 'GFSQ'
inModSFSQs = modsSFSQ[modsSFSQ$Module != "UNCLASSIFIED",]
inModSFSQs$Group = 'SFSQ'


fSQmembership = rbind(inModSFSQs,inModGFSQs)
fSQccurences = as.data.frame(table(fSQmembership$Gene))
fSQgenes = fSQccurences$Var1[fSQccurences$Freq==2]
mSQmembership = rbind(inModSMSQs,inModGMSQs)
mSQccurences = as.data.frame(table(mSQmembership$Gene))
mSQgenes = mSQccurences$Var1[mSQccurences$Freq==2]
fVmembership = rbind(inModSFVs, inModGFVs)
fVccurences = as.data.frame(table(fVmembership$Gene))
fVgenes = fVccurences$Var1[fVccurences$Freq==2]
mVmembership = rbind(inModGMVs, inModSMVs )
mVccurences = as.data.frame(table(mVmembership$Gene))
mVgenes = mVccurences$Var1[mVccurences$Freq==2]


library(dplyr)
library(tidyr)
##
mVGeness = as.data.frame(mVgenes)
colnames(mVGeness) = "Gene"
inModSMVgene = inModSMVs[,1:2]
MVGenes = left_join(mVGeness, inModSMVgene, by = "Gene")
colnames(MVGenes) = c("Gene", "Starnet Module")
inModGMVgene = inModGMVs[,1:2]
MVGenes = left_join(MVGenes, inModGMVgene, by = "Gene")
colnames(MVGenes) = c("Gene", "Starnet Module", "GTEx Module")

mSQGeness = as.data.frame(mSQgenes)
colnames(mSQGeness) = "Gene"
inModSMSQgene = inModSMSQs[,1:2]
MSQGenes = left_join(mSQGeness, inModSMSQgene, by = "Gene")
colnames(MSQGenes) = c("Gene", "Starnet Module")
inModGMSQgene = inModGMSQs[,1:2]
MSQGenes = left_join(MSQGenes, inModGMSQgene, by = "Gene")
colnames(MSQGenes) = c("Gene", "Starnet Module", "GTEx Module")

fVGeness = as.data.frame(fVgenes)
colnames(fVGeness) = "Gene"
inModSFVgene = inModSFVs[,1:2]
FVGenes = left_join(fVGeness, inModSFVgene, by = "Gene")
colnames(FVGenes) = c("Gene", "Starnet Module")
inModGFVgene = inModGFVs[,1:2]
FVGenes = left_join(FVGenes, inModGFVgene, by = "Gene")
colnames(FVGenes) = c("Gene", "Starnet Module", "GTEx Module")

fSQGeness = as.data.frame(fSQgenes)
colnames(fSQGeness) = "Gene"
inModSFSQgene = inModSFSQs[,1:2]
FSQGenes = left_join(fSQGeness, inModSFSQgene, by = "Gene")
colnames(FSQGenes) = c("Gene", "Starnet Module")
inModGFSQgene = inModGFSQs[,1:2]
FSQGenes = left_join(FSQGenes, inModGFSQgene, by = "Gene")
colnames(FSQGenes) = c("Gene", "Starnet Module", "GTEx Module")

write.table(FSQGenes, "FSQ Shared Module Genes.txt",quote  = F, row.names = F, sep = "\t")
write.table(MSQGenes, "MSQ Shared Module Genes.txt",quote  = F, row.names = F, sep = "\t")
write.table(FVGenes, "FV Shared Module Genes.txt",quote  = F, row.names = F, sep = "\t")
write.table(MVGenes, "MV Shared Module Genes.txt",quote  = F, row.names = F, sep = "\t")




conditionSpecific = unique(c(as.character(mVgenes),as.character(fVgenes),as.character(mSQgenes),as.character(fSQgenes)))

#other genesets
klf14net = read.table("klf14network.txt")
colocBFD = read.table('raulieColoc.txt')
nearestBFD = read.table("allWHRGWASgenes.txt")
overlap = nearestBFD[!is.na(match(nearestBFD$V1,colocBFD$V1)),]
GWAS = unique(c(as.character(colocBFD$V1), as.character(nearestBFD$V1)))

unionBFD = unique(c(as.character(conditionSpecific), as.character(nearestBFD$V1),as.character(t(klf14net$V1)), as.character(colocBFD$gene_symbol)))
length(conditionSpecific[!is.na(match(conditionSpecific, GWAS))])


write.table(unionBFD, "allInputNetGenes.txt", quote = F, row.names = F, sep = "\t")
unionBFD = read.table("allInputNetGenes.txt", header = T, sep = '\t')






