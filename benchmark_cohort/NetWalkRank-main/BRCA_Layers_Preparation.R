##Read MI values of each layer edges
library(ggplot2)
edgeMI1File<- "./Data/Results/Results_MI/BRCA_MI_1N_LA_Normalized.txt"
edgeList1<- read.delim(edgeMI1File,header = TRUE ,sep="\t")
colnames(edgeList1)<-c("SourceNode","TargetNode","Weight")
Ledges1<-edgeList1

edgeMI2File<- "./Data/Results/Results_MI/BRCA_MI_2N_LA_Normalized.txt"
edgeList2<- read.delim(edgeMI2File,header = TRUE ,sep="\t")
colnames(edgeList2)<-c("SourceNode","TargetNode","Weight")
Ledges2<-edgeList2

edgeMI3File<- "./Data/Results/Results_MI/BRCA_MI_3N_LA_Normalized.txt"
edgeList3<- read.delim(edgeMI3File,header = TRUE ,sep="\t")
colnames(edgeList3)<-c("SourceNode","TargetNode","Weight")
Ledges3<-edgeList3

edgeMI4File<- "./Data/Results/Results_MI/BRCA_MI_4N_LA_Normalized.txt"
edgeList4<- read.delim(edgeMI4File,header = TRUE ,sep="\t")
colnames(edgeList4)<-c("SourceNode","TargetNode","Weight")
Ledges4<-edgeList4
#######
#Plot MI of Layers
LedgeWeiLs<- data.frame(matrix(0,nrow=0,ncol=4))
colnames(LedgeWeiLs)<- c("SourceNode","TargetNode" ,"Layer","Weight")

LedgeWei<- function(y,tc){
  dfLedWe<- data.frame(SourceNode=y$SourceNode,TargetNode=y$TargetNode,Layer=c(tc),Weight=y$Weight)
  
  LedgeWeiLs<<- rbind(LedgeWeiLs,dfLedWe)
}
LedgeWei(Ledges1,"Stage I")
LedgeWei(Ledges2,"Stage II")
LedgeWei(Ledges3,"Stage III")
LedgeWei(Ledges4,"Stage IV")

pdf("BRCA_MIPlot_LimmaAll.pdf")
ggplot(data =LedgeWeiLs, aes(x=Weight, group=Layer, fill=Layer)) + geom_density(adjust=1.5, alpha=.4)#+ theme_ipsum()
dev.off()

######
##Z Score P-value
LedgeWeiLs$ZScore <- (LedgeWeiLs$Weight-mean(LedgeWeiLs$Weight))/sd(LedgeWeiLs$Weight)
pdf("BRCA_MIPlot_Zscore.pdf")
ggplot(data =LedgeWeiLs, aes(x=ZScore, group=Layer, fill=Layer)) + geom_density(adjust=1.5, alpha=.4)#+ theme_ipsum()
dev.off()


LedgeWeiLs$PValue<- (pnorm(LedgeWeiLs$ZScore,lower.tail=FALSE))
pdf("BRCA_MIPlot_ZscorePval.pdf")
ggplot(data =LedgeWeiLs, aes(x=PValue, group=Layer, fill=Layer)) + geom_density(adjust=1.5, alpha=.4)#+ theme_ipsum()
dev.off()

LedgeWeiLs1<- LedgeWeiLs[LedgeWeiLs$PValue<0.05,]

pdf("BRCA_MIPlot_filter.pdf")
ggplot(data =LedgeWeiLs1, aes(x=Weight, group=Layer, fill=Layer)) + geom_density(adjust=1.5, alpha=.4)#+ theme_ipsum()
dev.off()


pdf("./Data/Results/Results_LeyerPreparation/BRCA_MIPlot_ZscorePval1.pdf")
ggplot(data =LedgeWeiLs, aes(x=Weight, group=Layer, fill=Layer)) + geom_density(adjust=1.5, alpha=.4)+
  geom_vline(xintercept = 0.77, linetype="dashed",color = "red", size=1.5)+theme_bw()
dev.off()

Ledges1_0<- LedgeWeiLs1[LedgeWeiLs1$Layer=="Stage I",c(1,2,4)]
Ledges2_0<- LedgeWeiLs1[LedgeWeiLs1$Layer=="Stage II",c(1,2,4)]
Ledges3_0<- LedgeWeiLs1[LedgeWeiLs1$Layer=="Stage III",c(1,2,4)]
Ledges4_0<- LedgeWeiLs1[LedgeWeiLs1$Layer=="Stage IV",c(1,2,4)]

write.table(Ledges1_0, "./Data/Results/Results_LeyerPreparation/BRCA_edgeList_MI_L1.txt", quote=F,sep="\t",row.names=FALSE)
write.table(Ledges2_0, "./Data/Results/Results_LeyerPreparation/BRCA_edgeList_MI_L2.txt", quote=F,sep="\t",row.names=FALSE)
write.table(Ledges3_0, "./Data/Results/Results_LeyerPreparation/BRCA_edgeList_MI_L3.txt", quote=F,sep="\t",row.names=FALSE)
write.table(Ledges4_0, "./Data/Results/Results_LeyerPreparation/BRCA_edgeList_MI_L4.txt", quote=F,sep="\t",row.names=FALSE)

############################################################################

