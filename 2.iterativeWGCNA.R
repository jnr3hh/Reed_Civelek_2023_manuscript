#2. iterativeWGCNA setup

# determine WGCNA power that optimizes scale-free network

setwd("~/") ## change to your working directory

library(Biobase)
library(GEOquery)
library(limma)
library(dplyr)

library(WGCNA)
library(sva)
library(colorspace)
library(qvalue)
library(cluster)
library(flashClust)
options(stringsAsFactors = FALSE);
enableWGCNAThreads()
allowWCGNAThreads()


#change input, file should be genes on rows, samples on columns
fstar_SQ_gene_exp_processed = read.table("iterSFSQ.txt", header = T, sep = "\t", fill= T)
input = t(fstar_SQ_gene_exp_processed)

# check for genes and samples with too many missing values in the null data
# if gsg$allOK returns true, all genes have passed the cuts. If not,
# we remove the offending genes and samples from the data set
gsg <- goodSamplesGenes(input, verbose=3)
dim0 = dim(input)
if(!gsg$allOK) {  
  if(sum(!gsg$goodGenes) > 0)
    printFlush(paste('Removing genes:', paste(names(input)[!gsg$goodGenes], collapse = ', ')))
  if(sum(!gsg$goodSamples) > 0)
    printFlush(paste('Removing samples:', paste(rownames(input)[!gsg$goodSamples], collapse = ', ')))
  # Remove offending genes and samples from the data
  input <- input[gsg$goodSamples, gsg$goodGenes]
}
if(dim0[1]!=nrow(input)){paste0("Null samples were removed")}

# Cluster the Null samples to ID outliers (Euclidian distance)
sampleTree = flashClust(dist(input), method = "average")
plot(sampleTree, main = "sample clustering to detect outliers", sub="", xlab="", 
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

# apply a cut after reviewing this plot (where?)
cut = 200 #maybe change this??????
# Plot a line to show the cut
pdf(file= "prePro_Null_sample_clustering_to_detect_outliers.pdf", width=12, height=9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "sample clustering to detect outliers", sub="", xlab="", 
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
abline(h = cut, col = "red");
dev.off()
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = cut, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
input = input[keepSamples, ]
length(which(keepSamples==F))

#rewrite table if outlier samples were removed
write.table(t(input), "iterSFSQ.txt", quote = F, sep = "\t")






###############################################################
###### initial analysis: identify exponential power for sparsification
# used to create adjacency matrix in WGCNA

# range of powers for consideration
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# select power(s) that maximizes scale independence (left) and minimizes mean connectivity (right) 
# ideal powers are typically 6-12 for these datasets
sft = pickSoftThreshold(input, powerVector = powers, verbose = 5)

par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


# specify power for further network analysis
pow = 10





###################################################################
#RUN ITERATIVE WGCNA

#caution, these commands must be run in python


#use the python3 in command line or anaconda window
#move to folder containing files
cd C:\Users\jnr3hh\Documents\iterative\Nutri_indiv

# install the module with:
python3 -m pip install iterativeWGCNA
# make sure you have a working version of rpy2

# confirm that your iterativeWGCNA was installed and loaded correctly
python3 -m iterativeWGCNA --help

# run iterative on one dataset:
python3 -m iterativeWGCNA -i C:\Users\jnr3hh\Documents\iterative\Nutri_indiv\star_SQ_processed_geneExp.txt --enableWGCNAThreads --wgcnaParameters maxBlockSize=30000,power=10

#make sure you edit maxBlockSize to be just bigger than the geneset size
#edit the power
# can use -o to specify an outfile
# --enableWGCNAThreads will paralellize the processes over n-1 cores on your machine

#if you start with 20000 or more genes, it will take 12-36+ hours to run on a normal machine with 8 cores (can move to a computing cluster, but it may be more difficult to get the correct versions of R, python, rpy2, etc)

# keep the output files "merged-0.05-eigengenes.txt" and "merged-0.05-membership.txt", rename with a dataset identifier such as "SFSQ_merged-0.05-eigengenes.txt"







