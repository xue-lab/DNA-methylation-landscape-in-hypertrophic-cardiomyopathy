# DNA-methylation-landscape-in-hypertrophic-cardiomyopathy

library("ChAMP")

myLoad <- champ.load("./Human_HCM_850k",arraytype = "EPIC")
myLoad$pd

CpG.GUI(CpG=rownames(myLoad$beta), arraytype="EPIC")
QC.GUI(arraytype="EPIC")

myNorm <- champ.norm(plotBMIQ = T, arraytype="EPIC", cores=30)
QC.GUI(myNorm,arraytype="EPIC")
champ.QC(beta = myNorm,Feature.sel="SVD")

myCombat <- champ.runCombat(beta=myNorm,pd=myLoad$pd,variablename="Sample_Group",batchname=c("Array_Plate"))
champ.SVD(beta = myCombat, pd = myLoad$pd)

myDMP <- champ.DMP(beta = myCombat,pheno=myLoad$pd$Sample_Group,arraytype="EPIC")
DMP.GUI()

myDMR <- champ.DMR(beta = myCombat, pheno=myLoad$pd$Sample_Group, arraytype = "EPIC",method="DMRcate", cores=20)
DMR.GUI(arraytype="EPIC")

myGSEA <- champ.GSEA(arraytype = "EPIC", beta = myNorm, DMP = myDMP[[1]], DMR = myDMR, adjPval = 0.05, method = "fisher")
save(myGSEA,file="myGSEA.rda")
head(myGSEA$DMP)
head(myGSEA$DMR)

sig.DMP <- subset(myDMP[["Control_to_Case"]], abs(myDMP[["Control_to_Case"]][["deltaBeta"]]) > 0.1 & myDMP[["Control_to_Case"]][["adj.P.Val"]] < 0.05)
sig.cpg <- rownames(sig.DMP)
index <- which(rownames(myCombat) %in% sig.cpg)
sig.myCombat <- myCombat[index,]
head(sig.myCombat,5)

library(Rtsne)
df <- t(sig.myCombat)  
df[1:4,1:4] 

set.seed(30) 
tsne_out <- Rtsne(df,pca=FALSE,perplexity=20,theta=0, pca_scale = F, normalize = F)

tsnes=tsne_out$Y
colnames(tsnes) <- c("tSNE1", "tSNE2") 
tsnes=as.data.frame(tsnes)

group <- read.csv("./Sample_HCM_850k.csv")
groups <- as.factor(group$Sample_Group)

tsnes$group=groups
ggplot(tsnes, aes(x = tSNE1, y = tSNE2)) +
  scale_colour_manual(name="",values = c("#E7B800","#FF9289")) +
  geom_point(aes(colour=groups))




