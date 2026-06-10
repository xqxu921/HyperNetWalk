#DEGs For all layer

#Install packages
install.packages(c("GEOquery", "limma", "pheatmap", "ggplot2", "gplots", "reshape2", "plyr","Biobase"))
install.packages("Rgraphviz")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("minet")
#install.packages("igraph",dependencies = TRUE)

install.packages(c("igraph","graphlayouts","ggraph","threejs"))
install.packages("Matrix")
install.packages("tidyverse")
####  Libraries
library(GEOquery)
library(limma)
library(Biobase)
library(pheatmap)
library(gplots)
library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr) 
library(Rgraphviz)
library(minet)
library(igraph)
library(graphlayouts) 
library(ggraph)
library(threejs)
library(Matrix)

#set directory
setwd("/home/zcx/Fatemeh/")
#memory limit
# memory.limit(size=25000000)
####Install Packages
setRepositories()
# 1 2
#-- enter 1 2 to select BioC Software

#### Insert Genes Data
series<- "GSE89377"
gset1 <- getGEO(series, GSEMatrix =TRUE, AnnotGPL = TRUE, destdir ="/home/zcx/Fatemeh/Data")
gset<- gset1[[1]]
# Gene Expression
ex<- exprs(gset)

# pdf("Plots/ex_boxplot.pdf",width=107)
# boxplot(ex)
# dev.off()

gr<- c(rep("Normal",13),rep("chHLG",8),rep("chHHG",12), rep("Cirr",12),
       rep("DyNLG",11), rep("DyNHG",11), rep("earlyhc",5), rep("TG1",9),
       rep("TG2",12), rep("TG3",14))


# pdf("Plots/cor_heatmap_ex.pdf",width=15,height=15)
# pheatmap(cor(ex),labels_row = gr,labels_col = gr,color = bluered(256),border_color = NA,cluster_rows = F,cluster_cols = F)
# pheatmap(cor(ex),labels_row = gr,labels_col = gr,color = bluered(256),border_color = NA)
# dev.off()


# pc<-prcomp(ex)
# pdf("Plots/PC_ex.pdf",width=15,height=15)
# plot(pc)
# plot(pc$x[,1:2])
# dev.off()

#ex.scale<- t(scale(t(ex),scale=F))
#pc<- prcomp(ex.scale)
#pdf("Plots/PC_scale.pdf")
#plot(pc)
#plot(pc$x[,1:2])
#dev.off()

#pcr<-data.frame(pc$r[,1:3], Group=gr)
#pdf("Plots/PCA_samples.pdf")
#ggplot(pcr,aes(PC1,PC2,color=Group))+geom_point(size=3)+theme_bw()
#dev.off()

#BY Using LIMMA and Bayes on All Levels
gr<- factor(gr)
gset$description<- gr
design<- model.matrix(~description+0,gset)
colnames(design)<-levels(gr)
fit_All<- lmFit(gset, design)

cont.matrix_All<- makeContrasts(chHLG-Normal,chHHG-Normal,Cirr-Normal,DyNLG-Normal,
                                DyNHG-Normal,earlyhc-Normal,TG1-Normal,TG2-Normal,
                                TG3-Normal,levels=design)


fit_All<- contrasts.fit(fit_All, cont.matrix_All)
fit_All1<- eBayes(fit_All,0.01)
tT_All<- topTable(fit_All1,adjust="fdr", sort.by = "B", number=Inf)


subsetList<- c("ID","Gene.symbol","Gene.ID","chHLG...Normal","chHHG...Normal",
               "Cirr...Normal","DyNLG...Normal","DyNHG...Normal","earlyhc...Normal",
               "TG1...Normal","TG2...Normal","TG3...Normal","adj.P.Val")
tT_All1s<-subset(tT_All,select=subsetList)


aml.UpandDown_All<- subset(tT_All1s, adj.P.Val <0.05)

aml.UpandDown_AllC<- aml.UpandDown_All %>%
  filter(Gene.ID != "") %>%
  group_by(Gene.ID) %>%
  dplyr::slice(1) %>%
  ungroup()


aml.UpandDown_AllC<- data.frame(aml.UpandDown_AllC)

GI_ID_All<-unique(aml.UpandDown_AllC$Gene.ID)
tT_All_sN<- aml.UpandDown_AllC
dim(tT_All_sN)

#write.table(tT_All_sN[,c(2,3)],"DifferentStrategies/GeneIDSymbol.txt",quote=F,sep="\t",row.names=FALSE,col.names = FALSE)
###
LFC_Layers<-data.frame(tT_All_sN[,c(3:12)])
#write.table(LFC_Layers, "DifferentStrategies/LFC_Layers/LFC_Layers.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)

# LFC_L1<-data.frame(tT_All_sN[abs(tT_All_sN[,4])>1,3])
# LFC_L2<-data.frame(tT_All_sN[abs(tT_All_sN[,5])>1,3])
# LFC_L3<-data.frame(tT_All_sN[abs(tT_All_sN[,6])>1,3])
# LFC_L4<-data.frame(tT_All_sN[abs(tT_All_sN[,7])>1,3])
# LFC_L5<-data.frame(tT_All_sN[abs(tT_All_sN[,8])>1,3])
# LFC_L6<-data.frame(tT_All_sN[abs(tT_All_sN[,9])>1,3])
# LFC_L7<-data.frame(tT_All_sN[abs(tT_All_sN[,10])>1,3])
# LFC_L8<-data.frame(tT_All_sN[abs(tT_All_sN[,11])>1,3])
# LFC_L9<-data.frame(tT_All_sN[abs(tT_All_sN[,12])>1,3])

# write.table(LFC_L1, "DifferentStrategies/LFC_Layers/LFC_L1.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)
# write.table(LFC_L2, "DifferentStrategies/LFC_Layers/LFC_L2.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)
# write.table(LFC_L3, "DifferentStrategies/LFC_Layers/LFC_L3.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)
# write.table(LFC_L4, "DifferentStrategies/LFC_Layers/LFC_L4.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)
# write.table(LFC_L5, "DifferentStrategies/LFC_Layers/LFC_L5.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)
# write.table(LFC_L6, "DifferentStrategies/LFC_Layers/LFC_L6.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)
# write.table(LFC_L7, "DifferentStrategies/LFC_Layers/LFC_L7.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)
# write.table(LFC_L8, "DifferentStrategies/LFC_Layers/LFC_L8.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)
# write.table(LFC_L9, "DifferentStrategies/LFC_Layers/LFC_L9.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)



#####Input GGI-RegNet
fileGGI="Data/human-GGI-RegNet.source"
GGIRN<- read.delim(fileGGI,header = FALSE ,sep="\t")

GGIRN<-GGIRN[,c(2,4)]
colnames(GGIRN)<- c("SourceNode","TargetNode")


#Selection of asked genes from among the entire database connections
GGIRN1<- distinct(GGIRN)
GGIRN2<- GGIRN1[GGIRN1$SourceNode %in% tT_All_sN$Gene.ID,]
GGIRN3<- GGIRN2[GGIRN2$TargetNode %in% tT_All_sN$Gene.ID,]

Network1<-GGIRN3
#write.table(Network1, "DifferentStrategies/DEG_Layers/Network1_GGI_common.txt", quote=F,sep="\t",row.names=FALSE)
GNNodes1<-distinct(data.frame(Network1$SourceNode))
GNNodes2<-distinct(data.frame(Network1$TargetNode))
colnames(GNNodes1)<-"Genes"
colnames(GNNodes2)<-"Genes"
GNNodes<- rbind(GNNodes1,GNNodes2)
NetNodes<- distinct(GNNodes)
dim(NetNodes)
#c(NetNodes$Genes)
##############################

C_df1<-tT_All_sN[tT_All_sN$Gene.ID %in% (NetNodes$Genes),]
geneListID<- C_df1
cID <-c(geneListID$ID)
cID
GeneIDL<- c(geneListID$Gene.ID)
GeneIDL


####Export Data to compute MI for each Layer####
#######Layers#######
####################
##Normal
gset_N<- gset[cID, gset$description=="Normal"]
gsetN_df<- as.data.frame(gset_N)
gset_N_Clear<- gsetN_df[1:nrow(gsetN_df),1:nrow(geneListID)]
colnames(gset_N_Clear)<- GeneIDL
Normal_EXp<-t(gset_N_Clear)

########  L1  #######
gset1<- gset[cID, gset$description=="chHLG"]
gset1_df<- as.data.frame(gset1)
gset1Clear<- gset1_df[1:nrow(gset1_df),1:nrow(geneListID)]
colnames(gset1Clear)<- GeneIDL

#Output for MI Function
L1_Exp<-t(gset1Clear)
write.table(L1_Exp, "DifferentStrategies/DEG_Layers/L1_Exp.txt", quote=F,sep="\t",row.names=TRUE)

ex_L1_N<-cbind(Genes=rownames(Normal_EXp),L1_Exp,Normal_EXp)

write.table(ex_L1_N, "DifferentStrategies/DEG_Layers/ex_L1_N.txt", quote=F,sep="\t",row.names=FALSE)


###################
########  L2  #######
gset2<- gset[cID, gset$description=="chHHG"]
gset2_df<- as.data.frame(gset2)
gset2Clear<- gset2_df[1:nrow(gset2_df),1:nrow(geneListID)]
colnames(gset2Clear)<- GeneIDL

#Output for MI Function
L2_Exp<-t(gset2Clear)
write.table(L2_Exp, "DifferentStrategies/DEG_Layers/L2_Exp.txt", quote=F,sep="\t",row.names=TRUE)

ex_L2_N<-cbind(Genes=rownames(Normal_EXp),L2_Exp,Normal_EXp)
write.table(ex_L2_N, "DifferentStrategies/DEG_Layers/ex_L2_N.txt", quote=F,sep="\t",row.names=FALSE)


########  L3  #######
gset3<- gset[cID, gset$description=="Cirr"]
gset3_df<- as.data.frame(gset3)
gset3Clear<- gset3_df[1:nrow(gset3_df),1:nrow(geneListID)]
colnames(gset3Clear)<- GeneIDL

#Output for MI Function
L3_Exp<-t(gset3Clear)
write.table(L3_Exp, "DifferentStrategies/DEG_Layers/L3_Exp.txt", quote=F,sep="\t",row.names=TRUE)

ex_L3_N<-cbind(Genes=rownames(Normal_EXp),L3_Exp,Normal_EXp)
write.table(ex_L3_N, "DifferentStrategies/DEG_Layers/ex_L3_N.txt", quote=F,sep="\t",row.names=FALSE)

###################
########  L4  #######
gset4<- gset[cID, gset$description=="DyNLG"]
gset4_df<- as.data.frame(gset4)
gset4Clear<- gset4_df[1:nrow(gset4_df),1:nrow(geneListID)]
colnames(gset4Clear)<- GeneIDL

#Output for MI Function
L4_Exp<-t(gset4Clear)
write.table(L4_Exp, "DifferentStrategies/DEG_Layers/L4_Exp.txt", quote=F,sep="\t",row.names=TRUE)

ex_L4_N<-cbind(Genes=rownames(Normal_EXp),L4_Exp,Normal_EXp)
write.table(ex_L4_N, "DifferentStrategies/DEG_Layers/ex_L4_N.txt", quote=F,sep="\t",row.names=FALSE)

###################
########  L5  #######
gset5<- gset[cID, gset$description=="DyNHG"]
gset5_df<- as.data.frame(gset5)
gset5Clear<- gset5_df[1:nrow(gset5_df),1:nrow(geneListID)]
colnames(gset5Clear)<- GeneIDL


#Output for MI Function
L5_Exp<-t(gset5Clear)
write.table(L5_Exp, "DifferentStrategies/DEG_Layers/L5_Exp.txt", quote=F,sep="\t",row.names=TRUE)

ex_L5_N<-cbind(Genes=rownames(Normal_EXp),L5_Exp,Normal_EXp)
write.table(ex_L5_N, "DifferentStrategies/DEG_Layers/ex_L5_N.txt", quote=F,sep="\t",row.names=FALSE)

###################
########  L6  #######
gset6<- gset[cID, gset$description=="earlyhc"]
gset6_df<- as.data.frame(gset6)
gset6Clear<- gset6_df[1:nrow(gset6_df),1:nrow(geneListID)]
colnames(gset6Clear)<- GeneIDL


#Output for MI Function
L6_Exp<-t(gset6Clear)
write.table(L6_Exp, "DifferentStrategies/DEG_Layers/L6_Exp.txt", quote=F,sep="\t",row.names=TRUE)

ex_L6_N<-cbind(Genes=rownames(Normal_EXp),L6_Exp,Normal_EXp)
write.table(ex_L6_N, "DifferentStrategies/DEG_Layers/ex_L6_N.txt", quote=F,sep="\t",row.names=FALSE)

###################
########  L7  #######
gset7<- gset[cID, gset$description=="TG1"]
gset7_df<- as.data.frame(gset7)
gset7Clear<- gset7_df[1:nrow(gset7_df),1:nrow(geneListID)]
colnames(gset7Clear)<- GeneIDL

#Output for MI Function
L7_Exp<-t(gset7Clear)
write.table(L7_Exp, "DifferentStrategies/DEG_Layers/L7_Exp.txt", quote=F,sep="\t",row.names=TRUE)

ex_L7_N<-cbind(Genes=rownames(Normal_EXp),L7_Exp,Normal_EXp)
write.table(ex_L7_N, "DifferentStrategies/DEG_Layers/ex_L7_N.txt", quote=F,sep="\t",row.names=FALSE)

###################
########  L8  #######
gset8<- gset[cID, gset$description=="TG2"]
gset8_df<- as.data.frame(gset8)
gset8Clear<- gset8_df[1:nrow(gset8_df),1:nrow(geneListID)]
colnames(gset8Clear)<- GeneIDL

#Output for MI Function
L8_Exp<-t(gset8Clear)
write.table(L8_Exp, "DifferentStrategies/DEG_Layers/L8_Exp.txt", quote=F,sep="\t",row.names=TRUE)

ex_L8_N<-cbind(Genes=rownames(Normal_EXp),L8_Exp,Normal_EXp)
write.table(ex_L8_N, "DifferentStrategies/DEG_Layers/ex_L8_N.txt", quote=F,sep="\t",row.names=FALSE)

###################
########  L9  #######
gset9<- gset[cID, gset$description=="TG3"]
gset9_df<- as.data.frame(gset9)

gset9Clear<- gset9_df[1:nrow(gset9_df),1:nrow(geneListID)]
colnames(gset9Clear)<- GeneIDL

#Output for MI Function
L9_Exp<-t(gset9Clear)
write.table(L9_Exp, "DifferentStrategies/DEG_Layers/L9_Exp.txt", quote=F,sep="\t",row.names=TRUE)

ex_L9_N<-cbind(Genes=rownames(Normal_EXp),L9_Exp,Normal_EXp)
write.table(ex_L9_N, "DifferentStrategies/DEG_Layers/ex_L9_N.txt", quote=F,sep="\t",row.names=FALSE)


############################################################################