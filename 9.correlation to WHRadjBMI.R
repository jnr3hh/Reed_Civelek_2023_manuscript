#9. Correlations to WHRadjBMI in STARNET


setwd("~/Civelek Lab/RIMBANET/Final Net Data/REDO MET SEQ/STARandGTExNetsFinal/Newest")
library(ggplot2)
library(dplyr)
library(tidyr)


inGFV = read.table('iterGFV.txt', header = T)
inGMV = read.table('iterGMV.txt', header = T)
inSFV = read.table('iterSFV.txt', header = T)
inSMV = read.table('iterSMV.txt', header = T)
inGFSQ = read.table('iterGFSQ.txt', header = T)
inGMSQ = read.table('iterGMSQ.txt', header = T)
inSFSQ = read.table('iterSFSQ.txt', header = T)
inSMSQ = read.table('iterSMSQ.txt', header = T)

WntKDs = c("TYRO3","ANTXR1","ANAPC2","PSME3","RSPO1","BAZ1B","HELZ","MTMR9",'ARFGEF2',"ARMCX3","BNIP2","KIAA1522","ZNF148")
mitoKDs = c("TRIP12","A4GALT","BAD","MIGA1","C1QTNF3","PSME3","UBR1","ZNF148","INO80D","SPART","ARMCX3","NMT1","YME1L1")

genes = unique(c(WntKDs, mitoKDs))

SFVgenes = as.data.frame(t(inSFV[!is.na(match(rownames(inSFV), genes)),]))
SMVgenes = as.data.frame(t(inSMV[!is.na(match(rownames(inSMV), genes)),]))
GFVgenes = as.data.frame(t(inGFV[!is.na(match(rownames(inGFV), genes)),]))
GMVgenes = as.data.frame(t(inGMV[!is.na(match(rownames(inGMV), genes)),]))
SFSQgenes = as.data.frame(t(inSFSQ[!is.na(match(rownames(inSFSQ), genes)),]))
SMSQgenes = as.data.frame(t(inSMSQ[!is.na(match(rownames(inSMSQ), genes)),]))
GFSQgenes = as.data.frame(t(inGFSQ[!is.na(match(rownames(inGFSQ), genes)),]))
GMSQgenes = as.data.frame(t(inGMSQ[!is.na(match(rownames(inGMSQ), genes)),]))
GMSQgenes$net = 'GMSQ'
GFSQgenes$net = 'GFSQ'
SMSQgenes$net = 'SMSQ'
SFSQgenes$net = 'SFSQ'
GMVgenes$net = 'GMV'
GFVgenes$net = 'GFV'
SMVgenes$net = 'SMV'
SFVgenes$net = 'SFV'
GMSQgenes$net2 = 'GTEx'
GFSQgenes$net2 = 'GTEx'
SMSQgenes$net2 = 'STARNET'
SFSQgenes$net2 = 'STARNET'
GMVgenes$net2 = 'GTEx'
GFVgenes$net2 = 'GTEx'
SMVgenes$net2 = 'STARNET'
SFVgenes$net2 = 'STARNET'
GMSQgenes$sex = 'male'
GFSQgenes$sex = 'female'
SMSQgenes$sex = 'male'
SFSQgenes$sex = 'female'
GMVgenes$sex = 'male'
GFVgenes$sex = 'female'
SMVgenes$sex = 'male'
SFVgenes$sex = 'female'
GMSQgenes$depot = 'SubQ'
GFSQgenes$depot = 'SubQ'
SMSQgenes$depot = 'SubQ'
SFSQgenes$depot = 'SubQ'
GMVgenes$depot = 'Visc'
GFVgenes$depot = 'Visc'
SMVgenes$depot = 'Visc'
SFVgenes$depot = 'Visc'


full= rbind(SFVgenes, SMVgenes, SFSQgenes, SMSQgenes, GMVgenes, GFVgenes, GFSQgenes, GMSQgenes)


# plot expression in depot and sex

full %>% 
dplyr::select(ANAPC2,PSME3,RSPO1,TYRO3, net2, depot, sex) %>% 
  pivot_longer(names_to = "Gene",values_to = "Expression",cols = 1:4) %>% 
  filter(sex == 'female') %>%
  filter(net2 == "GTEx") %>% 
  ggplot(aes(depot, Expression))+ geom_boxplot(color = 'darkred', fatten = 3)+
  geom_jitter(width = .2)+ylab('Expression in Females')+theme_classic()+xlab('Depot')+
  ggtitle('GTEx')+
  theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=12),axis.title=element_text(size=12))+facet_wrap(~Gene)

#annotation data
viscStarphen = read.table("VAF_sampleinfo.txt", header = T)
subqStarphen = read.table("SF_sampleinfo.txt", header = T)


subqStarphenF = subqStarphen %>% filter(Sex == "female") %>% arrange(sampleid)
subqStarphenM = subqStarphen %>% filter(Sex == "male") %>% arrange(sampleid)

viscStarphenF = viscStarphen %>% filter(Sex == "female") %>% arrange(sampleid)
viscStarphenM = viscStarphen %>% filter(Sex == "male") %>% arrange(sampleid)



datTraits= viscStarphenF
WHR = as.data.frame(datTraits$"Waist.Hip.ratio");
names(WHR) = "WHR"
BMI = as.data.frame(datTraits$"BMI");
names(BMI) = "BMI"

WHR<-as.numeric(unlist(WHR))
BMI<-as.numeric(unlist(BMI))
whr.bmi<-lm(WHR~BMI)
viscStarphenF$WHR_BMIadj<-resid(whr.bmi)

datTraits= viscStarphenM
WHR = as.data.frame(datTraits$"Waist.Hip.ratio");
names(WHR) = "WHR"
BMI = as.data.frame(datTraits$"BMI");
names(BMI) = "BMI"

WHR<-as.numeric(unlist(WHR))
BMI<-as.numeric(unlist(BMI))
whr.bmi<-lm(WHR~BMI)
viscStarphenM$WHR_BMIadj<-resid(whr.bmi)

datTraits= subqStarphenF
WHR = as.data.frame(datTraits$"Waist.Hip.ratio");
names(WHR) = "WHR"
BMI = as.data.frame(datTraits$"BMI");
names(BMI) = "BMI"

WHR<-as.numeric(unlist(WHR))
BMI<-as.numeric(unlist(BMI))
whr.bmi<-lm(WHR~BMI)
subqStarphenF$WHR_BMIadj<-resid(whr.bmi)

datTraits= subqStarphenM
WHR = as.data.frame(datTraits$"Waist.Hip.ratio");
names(WHR) = "WHR"
BMI = as.data.frame(datTraits$"BMI");
names(BMI) = "BMI"

WHR<-as.numeric(unlist(WHR))
BMI<-as.numeric(unlist(BMI))
whr.bmi<-lm(WHR~BMI)
subqStarphenM$WHR_BMIadj<-resid(whr.bmi)




wntMV = as.data.frame(t(inSMV[!is.na(match(rownames(inSMV), WntKDs)),]))
wntMV = wntMV %>% arrange(rownames(wntMV))
mitoMV= as.data.frame(t(inSMV[!is.na(match(rownames(inSMV), mitoKDs)),])) 
mitoMV = mitoMV %>% arrange(rownames(mitoMV))
wntFV= as.data.frame(t(inSFV[!is.na(match(rownames(inSFV), WntKDs)),]))
wntFV = wntFV %>% arrange(rownames(wntFV))
mitoFV= as.data.frame(t(inSFV[!is.na(match(rownames(inSFV), mitoKDs)),]))
mitoFV = mitoFV %>% arrange(rownames(mitoFV))
wntMSQ= as.data.frame(t(inSMSQ[!is.na(match(rownames(inSMSQ), WntKDs)),]))
wntMSQ = wntMSQ %>% arrange(rownames(wntMSQ))
mitoMSQ= as.data.frame(t(inSMSQ[!is.na(match(rownames(inSMSQ), mitoKDs)),]))
mitoMSQ = mitoMSQ %>% arrange(rownames(mitoMSQ))
wntFSQ= as.data.frame(t(inSFSQ[!is.na(match(rownames(inSFSQ), WntKDs)),]))
wntFSQ = wntFSQ %>% arrange(rownames(wntFSQ))
mitoFSQ= as.data.frame(t(inSFSQ[!is.na(match(rownames(inSFSQ), mitoKDs)),]))
mitoFSQ = mitoFSQ%>% arrange(rownames(mitoFSQ))


wntFSQcor = as.data.frame(t(wntFSQ[1,]))
wntFSQcor$cor = 0
wntFSQcor$pval = 0
for (i in 1:ncol(wntFSQ)){
  l = cor.test(subqStarphenF$WHR_BMIadj, wntFSQ[,i])
  wntFSQcor$cor[i] = l$estimate 
  wntFSQcor$pval[i] = l$p.value 
}
wntFSQcor$padj = p.adjust(wntFSQcor$pval, method = "fdr", n = length(wntFSQcor$pval))

mitoFSQcor = as.data.frame(t(mitoFSQ[1,]))
mitoFSQcor$cor = 0
mitoFSQcor$pval = 0
for (i in 1:ncol(mitoFSQ)){
  l = cor.test(subqStarphenF$WHR_BMIadj, mitoFSQ[,i])
  mitoFSQcor$cor[i] = l$estimate 
  mitoFSQcor$pval[i] = l$p.value 
}
mitoFSQcor$padj = p.adjust(mitoFSQcor$pval, method = "fdr", n = length(mitoFSQcor$pval))


wntFVcor = as.data.frame(t(wntFV[1,]))
wntFVcor$cor = 0
wntFVcor$pval = 0
for (i in 1:ncol(wntFV)){
  l = cor.test(viscStarphenF$WHR_BMIadj, wntFV[,i])
  wntFVcor$cor[i] = l$estimate 
  wntFVcor$pval[i] = l$p.value 
}
wntFVcor$padj = p.adjust(wntFVcor$pval, method = "fdr", n = length(wntFVcor$pval))

mitoFVcor = as.data.frame(t(mitoFV[1,]))
mitoFVcor$cor = 0
mitoFVcor$pval = 0
for (i in 1:ncol(mitoFV)){
  l = cor.test(viscStarphenF$WHR_BMIadj, mitoFV[,i])
  mitoFVcor$cor[i] = l$estimate 
  mitoFVcor$pval[i] = l$p.value 
}
mitoFVcor$padj = p.adjust(mitoFVcor$pval, method = "fdr", n = length(mitoFVcor$pval))



wntMVcor = as.data.frame(t(wntMV[1,]))
wntMVcor$cor = 0
wntMVcor$pval = 0
for (i in 1:ncol(wntMV)){
  l = cor.test(viscStarphenM$WHR_BMIadj, wntMV[,i])
  wntMVcor$cor[i] = l$estimate 
  wntMVcor$pval[i] = l$p.value 
}
wntMVcor$padj = p.adjust(wntMVcor$pval, method = "fdr", n = length(wntMVcor$pval))

mitoMVcor = as.data.frame(t(mitoMV[1,]))
mitoMVcor$cor = 0
mitoMVcor$pval = 0
for (i in 1:ncol(mitoMV)){
  l = cor.test(viscStarphenM$WHR_BMIadj, mitoMV[,i])
  mitoMVcor$cor[i] = l$estimate 
  mitoMVcor$pval[i] = l$p.value 
}
mitoMVcor$padj = p.adjust(mitoMVcor$pval, method = "fdr", n = length(mitoMVcor$pval))



wntMSQcor = as.data.frame(t(wntMSQ[1,]))
wntMSQcor$cor = 0
wntMSQcor$pval = 0
for (i in 1:ncol(wntMSQ)){
  l = cor.test(subqStarphenM$WHR_BMIadj, wntMSQ[,i])
  wntMSQcor$cor[i] = l$estimate 
  wntMSQcor$pval[i] = l$p.value 
}
wntMSQcor$padj = p.adjust(wntMSQcor$pval, method = "fdr", n = length(wntMSQcor$pval))

mitoMSQcor = as.data.frame(t(mitoMSQ[1,]))
mitoMSQcor$cor = 0
mitoMSQcor$pval = 0
for (i in 1:ncol(mitoMSQ)){
  l = cor.test(subqStarphenM$WHR_BMIadj, mitoMSQ[,i])
  mitoMSQcor$cor[i] = l$estimate 
  mitoMSQcor$pval[i] = l$p.value 
}
mitoMSQcor$padj = p.adjust(mitoMSQcor$pval, method = "fdr", n = length(mitoMSQcor$pval))















wntcor = as.data.frame(cbind(wntMSQcor$cor,wntFSQcor$cor,wntMVcor$cor,wntFVcor$cor))
rownames(wntcor) = colnames(wntFSQ)
wntcor$gene = colnames(wntFSQ)
colnames(wntcor) = c("MSQ","FSQ","MV","FV", "gene")
wntcor$gene <- factor(wntcor$gene, levels=c("ANAPC2","PSME3","RSPO1","TYRO3","ANTXR1",'ARFGEF2',"ARMCX3","BAZ1B","BNIP2","HELZ","KIAA1522","MTMR9","ZNF148"))
wntcor = wntcor %>% arrange(gene)

wntpadj = as.data.frame(cbind(wntMSQcor$padj,wntFSQcor$padj,wntMVcor$padj,wntFVcor$padj))
rownames(wntpadj) = colnames(wntFSQ)
wntpadj$gene = colnames(wntFSQ)
colnames(wntpadj) = c("MSQ","FSQ","MV","FV", "gene")
wntpadj$gene <- factor(wntpadj$gene, levels=c("ANAPC2","PSME3","RSPO1","TYRO3","ANTXR1",'ARFGEF2',"ARMCX3","BAZ1B","BNIP2","HELZ","KIAA1522","MTMR9","ZNF148"))
wntpadj = wntpadj %>% arrange(gene)
wntpadj_ast = wntpadj[,1:4]
wntpadj_ast$FSQ = ""
wntpadj_ast$FV = ""
wntpadj_ast$MSQ = ""
wntpadj_ast$MV = ""
for (i in 1:nrow(wntpadj)){
  for(j in 1:ncol(wntpadj_ast)){
    print(wntpadj[i,j])
    print(i)
    print(j)
    if(wntpadj[i,j] <= 0.05){
      wntpadj_ast[i,j] = "*"
    }
    if(wntpadj[i,j] <= 0.01){
      wntpadj_ast[i,j] = "**"
    }
    if(wntpadj[i,j] <= 0.001){
      wntpadj_ast[i,j] = "***"
    }
  }
}

newx = str_wrap(colnames(wntcor), width = 6)
colnames(wntcor) = newx 

library(stringr)
library(pheatmap)
WntHet = pheatmap((as.matrix(wntcor[,1:4])), cluster_rows = F, cluster_cols = F, fontsize = 8, gaps_row = 4, display_numbers = (as.matrix(wntpadj_ast)))




mitocor = as.data.frame(cbind(mitoMSQcor$cor,mitoFSQcor$cor,mitoMVcor$cor,mitoFVcor$cor))
rownames(mitocor) = colnames(mitoFSQ)
colnames(mitocor) = c("Subq Male","Subq Female","Visc Male","Visc Female")
mitocor$gene = rownames(mitocor)
mitocor$gene = factor(mitocor$gene, levels = c("C1QTNF3","MIGA1","PSME3","UBR1","A4GALT","BAD","TRIP12","ZNF148","ARMCX3","INO80D","NMT1","SPART","YME1L1"))
mitocor = mitocor %>% arrange(gene)

mitopadj = as.data.frame(cbind(mitoMSQcor$padj,mitoFSQcor$padj,mitoMVcor$padj,mitoFVcor$padj))
rownames(mitopadj) = colnames(mitoFSQ)
mitopadj$gene = colnames(mitoFSQ)
colnames(mitopadj) = c("MSQ","FSQ","MV","FV", "gene")
mitopadj$gene = factor(mitopadj$gene, levels = c("C1QTNF3","MIGA1","PSME3","UBR1","A4GALT","BAD","TRIP12","ZNF148","ARMCX3","INO80D","NMT1","SPART","YME1L1"))
mitopadj = mitopadj %>% arrange(gene)
mitopadj_ast = mitopadj[,1:4]
mitopadj_ast$FSQ = ""
mitopadj_ast$FV = ""
mitopadj_ast$MSQ = ""
mitopadj_ast$MV = ""
for (i in 1:nrow(mitopadj)){
  for(j in 1:ncol(mitopadj_ast)){
    print(mitopadj[i,j])
    print(i)
    print(j)
    if(mitopadj[i,j] <= 0.05){
      mitopadj_ast[i,j] = "*"
    }
    if(mitopadj[i,j] <= 0.01){
      mitopadj_ast[i,j] = "**"
    }
    if(mitopadj[i,j] <= 0.001){
      mitopadj_ast[i,j] = "***"
    }
  }
}

library(stringr)
library(pheatmap)
mitoHet = pheatmap((as.matrix(mitocor[,1:4])), cluster_rows = F, cluster_cols = F, fontsize = 8, gaps_row = c(4,8), display_numbers = (as.matrix(mitopadj_ast)))

