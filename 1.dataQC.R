#1. Data Preprocessing

setwd("~/Civelek Lab/RIMBANET/Final Net Data/REDO MET SEQ")
library(WGCNA)
#library(CePa)
library(dplyr)
library(stringr)
library(readxl)

#read in data (input gene exp and annotation data)
# do on cluster
annotgenes = read_excel("gtexannot.xlsx")
gtex = read.delim(file = "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", skip = 2)

#select adipose tissue sample names
gtexidsS = annotgenes$SAMPID[annotgenes$SMTSD == "Adipose - Subcutaneous"]
gtexidsV = annotgenes$SAMPID[annotgenes$SMTSD == "Adipose - Visceral (Omentum)"]

#same format
gtexidsS = str_replace_all(gtexidsS, "[-]", ".")
gtexidsV = str_replace_all(gtexidsV, "[-]", ".")

#select adipose tissue gene expression
gtexVisc = gtex[,!is.na(match(colnames(gtex),as.character(gtexidsV)))]
gtexSubq = gtex[,!is.na(match(colnames(gtex),as.character(gtexidsS)))]
#rownames as gene names
ensNames = gtex$Names
geneNames = gtex$Description
rownames(gtexVisc) = ensNames
rownames(gtexSubq) = ensNames
#remove object 
rm(gtex)


# separate males and females
expgtexVisc = as.numeric(as.character(gtexVisc[rownames(gtexVisc) == 'ENSG00000229807',]))
hist(expgtexVisc, breaks = 100)
abline(v = 0.5)
malesVisc = gtexVisc[,expgtexVisc<0.5]
femalesVisc = gtexVisc[,expgtexVisc>=0.5]

expgtexSubq = as.numeric(as.character(gtexSubq[rownames(gtexSubq) == 'ENSG00000229807',]))
hist(expgtexSubq, breaks = 100)
abline(v = 0.5)
malesSubq = gtexSubq[,expgtexSubq<0.5]
femalesSubq = gtexSubq[,expgtexSubq>=0.5]


#which genes are protein coding
library(biomaRt)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
MVids = rownames(maleVisc)
MSids = rownames(maleSubq)
MVidss =  gsub("\\..*","",MVids)
MSidss =  gsub("\\..*","",MSids)
FVids = rownames(femaleVisc)
FSids = rownames(femaleSubq)
FVidss =  gsub("\\..*","",FVids)
FSidss =  gsub("\\..*","",FSids)


annotFV = getBM(attributes = c("ensembl_gene_id", 'external_gene_name', 'gene_biotype'), filters = ("ensembl_gene_id"), values = FVidss, mart = ensembl)
proteincodingFV = annotFV[annotFV$gene_biotype == 'protein_coding',]
annotFS = getBM(attributes = c("ensembl_gene_id", 'external_gene_name', 'gene_biotype'), filters = ("ensembl_gene_id"), values = FSidss, mart = ensembl)
proteincodingFS = annotFS[annotFS$gene_biotype == 'protein_coding',]

annotMV = getBM(attributes = c("ensembl_gene_id", 'external_gene_name', 'gene_biotype'), filters = ("ensembl_gene_id"), values = MVidss, mart = ensembl)
proteincodingMV = annotMV[annotMV$gene_biotype == 'protein_coding',]
annotMS = getBM(attributes = c("ensembl_gene_id", 'external_gene_name', 'gene_biotype'), filters = ("ensembl_gene_id"), values = MSidss, mart = ensembl)
proteincodingMS = annotMS[annotMS$gene_biotype == 'protein_coding',]

#protein coding gene expression
maleSubq_PC = maleSubq[!is.na(match(MSidss,proteincodingMS$ensembl_gene_id)),]
maleVisc_PC = maleVisc[!is.na(match(MVidss,proteincodingMV$ensembl_gene_id)),]
femaleSubq_PC = femaleSubq[!is.na(match(FSidss,proteincodingFS$ensembl_gene_id)),]
femaleVisc_PC = femaleVisc[!is.na(match(FVidss,proteincodingFV$ensembl_gene_id)),]

MVidss_PC = MVidss[!is.na(match(MSidss,proteincodingMS$ensembl_gene_id))]
MSidss_PC = MSidss[!is.na(match(MVidss,proteincodingMV$ensembl_gene_id))]
MVidss_PC_short = gsub("\\..*","",Vidss_PC)
MSidss_PC_short = gsub("\\..*","",Sidss_PC)
FVidss_PC = FVidss[!is.na(match(FSidss,proteincodingFS$ensembl_gene_id))]
FSidss_PC = FSidss[!is.na(match(FVidss,proteincodingFV$ensembl_gene_id))]
FVidss_PC_short = gsub("\\..*","",FVidss_PC)
FSidss_PC_short = gsub("\\..*","",FSidss_PC)


#collapse on single gene
gtex_FSubq_PC_collapsed = collapseRows(as.matrix(femaleSubq_PC), FSidss_PC_short, rownames(femaleSubq_PC))
gtex_FSubq_PC_collapsed_df = as.data.frame(gtex_FSubq_PC_collapsed$datETcollapsed)

gtex_FVisc_PC_collapsed = collapseRows(as.matrix(femaleVisc_PC), FVidss_PC_short, rownames(femaleVisc_PC))
gtex_FVisc_PC_collapsed_df = as.data.frame(gtex_FVisc_PC_collapsed$datETcollapsed)

gtex_MSubq_PC_collapsed = collapseRows(as.matrix(maleSubq_PC), MSidss_PC_short, rownames(maleSubq_PC))
gtex_MSubq_PC_collapsed_df = as.data.frame(gtex_MSubq_PC_collapsed$datETcollapsed)

gtex_MVisc_PC_collapsed = collapseRows(as.matrix(maleVisc_PC), MVidss_PC_short, rownames(maleVisc_PC))
gtex_MVisc_PC_collapsed_df = as.data.frame(gtex_MVisc_PC_collapsed$datETcollapsed)


#replace low TPMs with 0
gtex_fsubq_PC_coll_zero = replace(gtex_FSubq_PC_collapsed_df, gtex_FSubq_PC_collapsed_df<=0.1, 0)
gtex_fvisc_PC_coll_zero = replace(gtex_FVisc_PC_collapsed_df, gtex_FVisc_PC_collapsed_df<=0.1, 0)
gtex_msubq_PC_coll_zero = replace(gtex_MSubq_PC_collapsed_df, gtex_MSubq_PC_collapsed_df<=0.1, 0)
gtex_mvisc_PC_coll_zero = replace(gtex_MVisc_PC_collapsed_df, gtex_MVisc_PC_collapsed_df<=0.1, 0)

MVsums = rep(0, nrow(gtex_mvisc_PC_coll_zero))
MSQsums = rep(0,nrow(gtex_msubq_PC_coll_zero))
FVsums = rep(0, nrow(gtex_fvisc_PC_coll_zero))
FSQsums = rep(0,nrow(gtex_fsubq_PC_coll_zero))

#calculate number of zeros
for(k in 1:nrow(gtex_fsubq_PC_coll_zero)){
  FVsums[k] = sum(gtex_fvisc_PC_coll_zero[k,]==0)
  FSQsums[k]= sum(gtex_fsubq_PC_coll_zero[k,]==0)
  MVsums[k] = sum(gtex_mvisc_PC_coll_zero[k,]==0)
  MSQsums[k]= sum(gtex_msubq_PC_coll_zero[k,]==0)
  print(k)
}

#remove rows with more than 80% zero
gtex_fsubq_PC_coll_nozero = gtex_fsubq_PC_coll_zero[(FSQsums< ncol(femaleSubq)*0.8),]
gtex_fvisc_PC_coll_nozero = gtex_fvisc_PC_coll_zero[(FVsums < ncol(femaleVisc)*0.8),]
gtex_msubq_PC_coll_nozero = gtex_msubq_PC_coll_zero[(MSQsums< ncol(maleSubq)*0.8),]
gtex_mvisc_PC_coll_nozero = gtex_mvisc_PC_coll_zero[(MVsums < ncol(maleVisc)*0.8),]


#log 
gtex_FSQ_gene_exp_processed = log2(gtex_fsubq_PC_coll_nozero+1)
gtex_FV_gene_exp_processed = log2(gtex_fvisc_PC_coll_nozero+1)
gtex_MSQ_gene_exp_processed = log2(gtex_msubq_PC_coll_nozero+1)
gtex_MV_gene_exp_processed = log2(gtex_mvisc_PC_coll_nozero+1)

write.table(gtex_FSQ_gene_exp_processed, "iterGFSQ.txt", quote = F, sep = "\t")
write.table(gtex_FV_gene_exp_processed, "iterGFV.txt", quote = F, sep = "\t")
write.table(gtex_MSQ_gene_exp_processed, "iterGMSQ.txt", quote = F, sep = "\t")
write.table(gtex_MV_gene_exp_processed, "iterGMV.txt", quote = F, sep = "\t")

#Starnet


setwd("~/Civelek Lab/RIMBANET/Final Net Data/REDO MET SEQ")

subqinfo = read.table("SF_sampleinfo.txt",header = T)
viscinfo = read.table("VAF_sampleinfo.txt", header = T)
Starsubq = read.csv('SFtpm.csv', header = T)
Starvisc = read.csv('VAFtpm.csv', header = T)

fullStarsubq = Starsubq[,2:574]
fullStarvisc = Starvisc[,2:532]
rownames(fullStarvisc) = Vids
rownames(fullStarsubq) = Sids
StarsubqS = fullStarsubq %>% select(!subq_missing_ids[2:25])
StarviscS = fullStarvisc %>% select(!visc_missing_ids[2:23])

malesSTARvisc = viscinfo$sampleid[viscinfo$Sex == 'male']
mStarvisc = StarviscS[,!is.na(match(colnames(StarviscS),malesSTARvisc))]

femalesSTARvisc = viscinfo$sampleid[viscinfo$Sex == 'female']
fStarvisc = StarviscS[,!is.na(match(colnames(StarviscS),femalesSTARvisc))]

malesSTARsubq = subqinfo$sampleid[subqinfo$Sex == 'male']
mStarsubq = StarsubqS[,!is.na(match(colnames(StarsubqS),malesSTARsubq))]

femalesSTARsubq = subqinfo$sampleid[subqinfo$Sex == 'female']
fStarsubq = StarsubqS[,!is.na(match(colnames(StarsubqS),femalesSTARsubq))]



#which genes are protein coding
library(biomaRt)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
MVids = rownames(mStarvisc)
MSids = rownames(mStarsubq)
FVids = rownames(fStarvisc)
FSids = rownames(fStarsubq)

annotMV = getBM(attributes = c("ensembl_gene_id", 'external_gene_name', 'gene_biotype'), filters = ("ensembl_gene_id"), values = MVids, mart = ensembl)
proteincodingMV = annotMV[annotMV$gene_biotype == 'protein_coding',]
annotMS = getBM(attributes = c("ensembl_gene_id", 'external_gene_name', 'gene_biotype'), filters = ("ensembl_gene_id"), values = MSids, mart = ensembl)
proteincodingMS = annotMS[annotMS$gene_biotype == 'protein_coding',]

annotFV = getBM(attributes = c("ensembl_gene_id", 'external_gene_name', 'gene_biotype'), filters = ("ensembl_gene_id"), values = FVids, mart = ensembl)
proteincodingFV = annotFV[annotFV$gene_biotype == 'protein_coding',]
annotFS = getBM(attributes = c("ensembl_gene_id", 'external_gene_name', 'gene_biotype'), filters = ("ensembl_gene_id"), values = FSids, mart = ensembl)
proteincodingFS = annotFS[annotFS$gene_biotype == 'protein_coding',]


fStarsubq_PC = fStarsubq[!is.na(match(rownames(fStarsubq), proteincodingFS$ensembl_gene_id)),]
fStarvisc_PC = fStarvisc[!is.na(match(rownames(fStarvisc), proteincodingFV$ensembl_gene_id)),]

mStarsubq_PC = mStarsubq[!is.na(match(rownames(mStarsubq), proteincodingMS$ensembl_gene_id)),]
mStarvisc_PC = mStarvisc[!is.na(match(rownames(mStarvisc), proteincodingMV$ensembl_gene_id)),]

fstar_subq_PC_zero = replace(fStarsubq_PC, fStarsubq_PC<=0.1, 0)
fstar_visc_PC_zero = replace(fStarvisc_PC, fStarvisc_PC<=0.1, 0)
mstar_subq_PC_zero = replace(mStarsubq_PC, mStarsubq_PC<=0.1, 0)
mstar_visc_PC_zero = replace(mStarvisc_PC, mStarvisc_PC<=0.1, 0)

fVsums = rep(0, nrow(fstar_visc_PC_zero))
fSQsums = rep(0,nrow(fstar_subq_PC_zero))
mVsums = rep(0, nrow(mstar_visc_PC_zero))
mSQsums = rep(0,nrow(mstar_subq_PC_zero))

#calculate number of zeros
for(k in 1:nrow(fstar_subq_PC_zero)){
  FVsums[k] = sum(fstar_visc_PC_zero[k,]==0)
  FSQsums[k]= sum(fstar_subq_PC_zero[k,]==0)
  MVsums[k] = sum(mstar_visc_PC_zero[k,]==0)
  MSQsums[k]= sum(mstar_subq_PC_zero[k,]==0)
  print(k)
}

#remove rows with more than 80% zero
fstar_subq_PC_nozero = fstar_subq_PC_zero[(FSQsums< ncol(fStarsubq)*0.8),]
fstar_visc_PC_nozero = fstar_visc_PC_zero[(FVsums < ncol(fStarvisc)*0.8),]
mstar_subq_PC_nozero = mstar_subq_PC_zero[(MSQsums< ncol(mStarsubq)*0.8),]
mstar_visc_PC_nozero = mstar_visc_PC_zero[(MVsums < ncol(mStarvisc)*0.8),]


#log 
fstar_SQ_gene_exp_processed = log2(fstar_subq_PC_nozero+1)
fstar_V_gene_exp_processed = log2(fstar_visc_PC_nozero+1)
mstar_SQ_gene_exp_processed = log2(mstar_subq_PC_nozero+1)
mstar_V_gene_exp_processed = log2(mstar_visc_PC_nozero+1)

write.table(fstar_SQ_gene_exp_processed, "iterSFSQ.txt", quote = F, sep = "\t")
write.table(fstar_V_gene_exp_processed, "iterSFV.txt", quote = F, sep = "\t")
write.table(mstar_SQ_gene_exp_processed, "iterSMSQ.txt", quote = F, sep = "\t")
write.table(mstar_V_gene_exp_processed, "iterSMV.txt", quote = F, sep = "\t")

