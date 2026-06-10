##Read MI values of each layer edges

edgeMI1File<- "DifferentStrategies/Results_MI/MI_1N_LA_Normalized.txt"
edgeList1<- read.delim(edgeMI1File,header = TRUE ,sep="\t")
colnames(edgeList1)<-c("SourceNode","TargetNode","Weight")
Ledges1<-edgeList1

edgeMI2File<- "DifferentStrategies/Results_MI/MI_2N_LA_Normalized.txt"
edgeList2<- read.delim(edgeMI2File,header = TRUE ,sep="\t")
colnames(edgeList2)<-c("SourceNode","TargetNode","Weight")
Ledges2<-edgeList2

edgeMI3File<- "DifferentStrategies/Results_MI/MI_3N_LA_Normalized.txt" 
edgeList3<- read.delim(edgeMI3File,header = TRUE ,sep="\t")
colnames(edgeList3)<-c("SourceNode","TargetNode","Weight")
Ledges3<-edgeList3

edgeMI4File<- "DifferentStrategies/Results_MI/MI_4N_LA_Normalized.txt" 
edgeList4<- read.delim(edgeMI4File,header = TRUE ,sep="\t")
colnames(edgeList4)<-c("SourceNode","TargetNode","Weight")
Ledges4<-edgeList4

edgeMI5File<- "DifferentStrategies/Results_MI/MI_5N_LA_Normalized.txt" 
edgeList5<- read.delim(edgeMI5File,header = TRUE ,sep="\t")
colnames(edgeList5)<-c("SourceNode","TargetNode","Weight")
Ledges5<-edgeList5

edgeMI6File<- "DifferentStrategies/Results_MI/MI_6N_LA_Normalized.txt" 
edgeList6<- read.delim(edgeMI6File,header = TRUE ,sep="\t")
colnames(edgeList6)<-c("SourceNode","TargetNode","Weight")
Ledges6<-edgeList6

edgeMI7File<- "DifferentStrategies/Results_MI/MI_7N_LA_Normalized.txt" 
edgeList7<- read.delim(edgeMI7File,header = TRUE ,sep="\t")
colnames(edgeList7)<-c("SourceNode","TargetNode","Weight")
Ledges7<-edgeList7

edgeMI8File<- "DifferentStrategies/Results_MI/MI_8N_LA_Normalized.txt" 
edgeList8<- read.delim(edgeMI8File,header = TRUE ,sep="\t")
colnames(edgeList8)<-c("SourceNode","TargetNode","Weight")
Ledges8<-edgeList8

edgeMI9File<- "DifferentStrategies/Results_MI/MI_9N_LA_Normalized.txt" 
edgeList9<- read.delim(edgeMI9File,header = TRUE ,sep="\t")
colnames(edgeList9)<-c("SourceNode","TargetNode","Weight")
Ledges9<-edgeList9


#######
#Plot MI of Layers
LedgeWeiLs<- data.frame(matrix(0,nrow=0,ncol=4))
colnames(LedgeWeiLs)<- c("SourceNode","TargetNode" ,"Layer","Weight")

LedgeWei<- function(y,tc){
  dfLedWe<- data.frame(SourceNode=y$SourceNode,TargetNode=y$TargetNode,Layer=c(tc),Weight=y$Weight)
  
  LedgeWeiLs<<- rbind(LedgeWeiLs,dfLedWe)
}
LedgeWei(Ledges1,"Choronic Hepatitis with low grade")
LedgeWei(Ledges2,"Choronic Hepatitis with high grade")
LedgeWei(Ledges3,"Cirrhosis")
LedgeWei(Ledges4,"Dysplastic nodules with low grade")
LedgeWei(Ledges5,"Dysplastic nodules with high grade")
LedgeWei(Ledges6,"Early hepatocellular carcinoma")
LedgeWei(Ledges7,"hepatocellular carcinoma (TG1)")
LedgeWei(Ledges8,"hepatocellular carcinoma (TG2)")
LedgeWei(Ledges9,"hepatocellular carcinoma (TG3)")

# pdf("MIPlot_LimmaAll.pdf")
ggplot(data =LedgeWeiLs, aes(x=Weight, group=Layer, fill=Layer)) + geom_density(adjust=1.5, alpha=.4)#+ theme_ipsum()
# dev.off()

######
##Z Score P-value
LedgeWeiLs$ZScore <- (LedgeWeiLs$Weight-mean(LedgeWeiLs$Weight))/sd(LedgeWeiLs$Weight)
# pdf("MIPlot_Zscore.pdf")
ggplot(data =LedgeWeiLs, aes(x=ZScore, group=Layer, fill=Layer)) + geom_density(adjust=1.5, alpha=.4)#+ theme_ipsum()
# dev.off()


LedgeWeiLs$PValue<- (pnorm(-abs(LedgeWeiLs$ZScore),lower.tail=TRUE))
# pdf("MIPlot_ZscorePval.pdf")
ggplot(data =LedgeWeiLs, aes(x=PValue, group=Layer, fill=Layer)) + geom_density(adjust=1.5, alpha=.4)#+ theme_ipsum()
# dev.off()

LedgeWeiLs1<- LedgeWeiLs[LedgeWeiLs$PValue<0.05,]

# pdf("MIPlot_filter.pdf")
ggplot(data =LedgeWeiLs1, aes(x=Weight, group=Layer, fill=Layer)) + geom_density(adjust=1.5, alpha=.4)#+ theme_ipsum()
# dev.off()


pdf("./Data/Results/Results_LeyerPreparation/MIPlot_ZscorePval1.pdf")
ggplot(data =LedgeWeiLs, aes(x=Weight, group=Layer, fill=Layer)) + geom_density(adjust=1.5, alpha=.4)+
  geom_vline(xintercept = 0.83, linetype="dashed",color = "red", size=1.5)+theme_bw()
dev.off()

Ledges1_0<- LedgeWeiLs1[LedgeWeiLs1$Layer=="Choronic Hepatitis with low grade",c(1,2,4)]
Ledges2_0<- LedgeWeiLs1[LedgeWeiLs1$Layer=="Choronic Hepatitis with high grade",c(1,2,4)]
Ledges3_0<- LedgeWeiLs1[LedgeWeiLs1$Layer=="Cirrhosis",c(1,2,4)]
Ledges4_0<- LedgeWeiLs1[LedgeWeiLs1$Layer=="Dysplastic nodules with low grade",c(1,2,4)]
Ledges5_0<- LedgeWeiLs1[LedgeWeiLs1$Layer=="Dysplastic nodules with high grade",c(1,2,4)]
Ledges6_0<- LedgeWeiLs1[LedgeWeiLs1$Layer=="Early hepatocellular carcinoma",c(1,2,4)]
Ledges7_0<- LedgeWeiLs1[LedgeWeiLs1$Layer=="hepatocellular carcinoma (TG1)",c(1,2,4)]
Ledges8_0<- LedgeWeiLs1[LedgeWeiLs1$Layer=="hepatocellular carcinoma (TG2)",c(1,2,4)]
Ledges9_0<- LedgeWeiLs1[LedgeWeiLs1$Layer=="hepatocellular carcinoma (TG3)",c(1,2,4)]

write.table(Ledges1_0, "./Data/Results/Results_LeyerPreparation/edgeList_MI_L1.txt", quote=F,sep="\t",row.names=FALSE)
write.table(Ledges2_0, "./Data/Results/Results_LeyerPreparation/edgeList_MI_L2.txt", quote=F,sep="\t",row.names=FALSE)
write.table(Ledges3_0, "./Data/Results/Results_LeyerPreparation/edgeList_MI_L3.txt", quote=F,sep="\t",row.names=FALSE)
write.table(Ledges4_0, "./Data/Results/Results_LeyerPreparation/edgeList_MI_L4.txt", quote=F,sep="\t",row.names=FALSE)
write.table(Ledges5_0, "./Data/Results/Results_LeyerPreparation/edgeList_MI_L5.txt", quote=F,sep="\t",row.names=FALSE)
write.table(Ledges6_0, "./Data/Results/Results_LeyerPreparation/edgeList_MI_L6.txt", quote=F,sep="\t",row.names=FALSE)
write.table(Ledges7_0, "./Data/Results/Results_LeyerPreparation/edgeList_MI_L7.txt", quote=F,sep="\t",row.names=FALSE)
write.table(Ledges8_0, "./Data/Results/Results_LeyerPreparation/edgeList_MI_L8.txt", quote=F,sep="\t",row.names=FALSE)
write.table(Ledges9_0, "./Data/Results/Results_LeyerPreparation/edgeList_MI_L9.txt", quote=F,sep="\t",row.names=FALSE)

############################################################################