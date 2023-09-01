# 4. write input files

######################## find geneset in input files, discretize, write Bayes input files#######

contGFV = inGFV[!is.na(match(row.names(inGFV), as.character(unionBFD))),]
contGMV = inGMV[!is.na(match(row.names(inGMV), as.character(unionBFD))),]
contGFSQ = inGFSQ[!is.na(match(row.names(inGFSQ), as.character(unionBFD))),]
contGMSQ = inGMSQ[!is.na(match(row.names(inGMSQ), as.character(unionBFD))),]

library(arules)
discGFSQ =  discretizeDF(contGFSQ, default = list(method = "cluster", breaks = 3, labels = c("0","1", "2")))
discGFV =  discretizeDF(contGFV, default = list(method = "cluster", breaks = 3, labels = c("0","1", "2")))
discGMSQ =  discretizeDF(contGMSQ, default = list(method = "cluster", breaks = 3, labels = c("0","1", "2")))
discGMV =  discretizeDF(contGMV, default = list(method = "cluster", breaks = 3, labels = c("0","1", "2")))

#discGTEx =  discretizeDF(adipose7, default = list(method = "cluster", breaks = 3, labels = c("0","1", "2")))
gtexcisqtls = read.table("Adipose_Subcutaneous.v8.egenes.txt", header = T, sep = "\t")
sigeqtl = gtexcisqtls[gtexcisqtls$qval<0.05,]
egenesSQ = sigeqtl$gene_name
gtexviscciseqtls = read.table(gzfile('Adipose_Visceral_Omentum.v8.egenes.txt.gz'),sep = "\t", header = T)
sigeqtlvisc = gtexcisqtls[gtexviscciseqtls$qval<0.05,]
egenesV = sigeqtlvisc$gene_name

qtlsinSetFSQ = egenesSQ[!is.na(match(egenesSQ,rownames(contGFSQ)))]
qtlsinSetFV = egenesV[!is.na(match(egenesV,rownames(contGFV)))]
qtlsinSetMSQ = egenesSQ[!is.na(match(egenesSQ,rownames(contGMSQ)))]
qtlsinSetMV = egenesV[!is.na(match(egenesV,rownames(contGMV)))]




contSFV = inSFV[!is.na(match(row.names(inSFV), as.character(unionBFD))),]
contSMV = inSMV[!is.na(match(row.names(inSMV), as.character(unionBFD))),]
contSFSQ = inSFSQ[!is.na(match(row.names(inSFSQ), as.character(unionBFD))),]
contSMSQ = inSMSQ[!is.na(match(row.names(inSMSQ), as.character(unionBFD))),]

library(arules)
#discGTEx =  discretizeDF(adipose7, default = list(method = "cluster", breaks = 3, labels = c("0","1", "2")))
discSFSQ =  discretizeDF(contSFSQ, default = list(method = "cluster", breaks = 3, labels = c("0","1", "2")))
discSFV =  discretizeDF(contSFV, default = list(method = "cluster", breaks = 3, labels = c("0","1", "2")))
discSMSQ =  discretizeDF(contSMSQ, default = list(method = "cluster", breaks = 3, labels = c("0","1", "2")))
discSMV =  discretizeDF(contSMV, default = list(method = "cluster", breaks = 3, labels = c("0","1", "2")))


setwd("~/Civelek Lab/RIMBANET/Final Net Data/REDO MET SEQ/STARandGTExNetsFinal")

starcisqtlsSQ = read.table("starnetSFeQTLs.txt", header = T, sep = "\t")
sigeqtlSQ = starcisqtlsSQ$genename[starcisqtlsSQ$SFegenes=='yes']
egenesSQ = sigeqtlSQ

starcisqtlsV = read.table("starnetVAFeQTLs.txt", header = T, sep = "\t")
sigeqtlV = starcisqtlsV$genename[starcisqtlsV$VAFegenes=='yes']
egenesV = sigeqtlV

qtlsinSetFSQ = egenesSQ[!is.na(match(egenesSQ,rownames(contSFSQ)))]
qtlsinSetFV = egenesV[!is.na(match(egenesV,rownames(contSFV)))]
qtlsinSetMSQ = egenesSQ[!is.na(match(egenesSQ,rownames(contSMSQ)))]
qtlsinSetMV = egenesV[!is.na(match(egenesV,rownames(contSMV)))]






setwd("~/Civelek Lab/RIMBANET/Final Net Data/REDO MET SEQ/STARandGTExNetsFinal/Newest")

write.table(discSFSQ, "discSFSQ.txt", col.names = F, quote = F)
write.table(contSFSQ, "contSFSQ.txt", col.names = F, quote = F)
write.table(qtlsinSetFSQ, "qtlsSFSQ.txt", col.names = F, quote = F, row.names = F)

write.table(discSFV, "discSFV.txt", col.names = F, quote = F)
write.table(contSFV, "contSFV.txt", col.names = F, quote = F)
write.table(qtlsinSetFV, "qtlsSFV.txt", col.names = F, quote = F, row.names = F)

write.table(discSMSQ, "discSMSQ.txt", col.names = F, quote = F)
write.table(contSMSQ, "contSMSQ.txt", col.names = F, quote = F)
write.table(qtlsinSetMSQ, "qtlsSMSQ.txt", col.names = F, quote = F, row.names = F)

write.table(discSMV, "discSMV.txt", col.names = F, quote = F)
write.table(contSMV, "contSMV.txt", col.names = F, quote = F)
write.table(qtlsinSetMV, "qtlsSMV.txt", col.names = F, quote = F, row.names = F)




setwd("~/Civelek Lab/RIMBANET/Final Net Data/REDO MET SEQ/STARandGTExNetsFinal/Newest")

write.table(discGFSQ, "discGFSQ.txt", col.names = F, quote = F)
write.table(contGFSQ, "contGFSQ.txt", col.names = F, quote = F)
write.table(qtlsinSetFSQ, "qtlsGFSQ.txt", col.names = F, quote = F, row.names = F)

write.table(discGFV, "discGFV.txt", col.names = F, quote = F)
write.table(contGFV, "contGFV.txt", col.names = F, quote = F)
write.table(qtlsinSetFV, "qtlsGFV.txt", col.names = F, quote = F, row.names = F)

write.table(discGMSQ, "discGMSQ.txt", col.names = F, quote = F)
write.table(contGMSQ, "contGMSQ.txt", col.names = F, quote = F)
write.table(qtlsinSetMSQ, "qtlsGMSQ.txt", col.names = F, quote = F, row.names = F)

write.table(discGMV, "discGMV.txt", col.names = F, quote = F)
write.table(contGMV, "contGMV.txt", col.names = F, quote = F)
write.table(qtlsinSetMV, "qtlsGMV.txt", col.names = F, quote = F, row.names = F)



