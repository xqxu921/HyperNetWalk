
setwd("/home/zcx/Fatemeh/")


ResultsTable<- as.matrix(read.delim("DifferentStrategies/REGml/r70/REGml_r70.txt",header = TRUE ,sep="\t")) 

genes_Ranks<-ResultsTable[,c(1,11,12,13,14,15,16,17,18,19,20,21,22,23)]

Ranks_Driver<- data.frame(genes_Ranks)


#############
DriverGenesT<-read.delim("DifferentStrategies/c20.txt")

colnames(DriverGenesT)<-c("genes")
DriverGenesT1<-data.frame(Genes=unique(DriverGenesT$genes))
dim(DriverGenesT1)

Ranks_Driver$DrNo<- "N"
Ranks_Driver[Ranks_Driver$genes %in% DriverGenesT1$Genes,]$DrNo<- "D"

Ranks_Driver[Ranks_Driver$DrNo=="D",]
dim(Ranks_Driver[Ranks_Driver$DrNo=="D",])

########
gmeltDN<-melt(Ranks_Driver[,c(-1,-11,-13,-14)])
DiLength<-dim(Ranks_Driver[Ranks_Driver$DrNo=="D",])
DiLength
#pdf("DifferentStrategies/Plots/DriverNormalBoxPlot_REGmlwpFC07.pdf",width = 15,height =15)
ggplot(gmeltDN,aes(variable,value,fill=DrNo))+geom_boxplot(position = "dodge")+theme_bw()
#dev.off()

# png("DifferentStrategies/Plots/DriverNormalBoxPlot_REGmlwpFC07.png",width = 1000,height =800)
# ggplot(gmeltDN,aes(variable,value,fill=DrNo))+geom_boxplot(position = "dodge")+theme_bw()
# dev.off()


#pdf("Results/DriverNormalViolinPlot_IG07.pdf",width = 15,height =15)
#ggplot(gmeltDN,aes(variable,value,fill=DrNo))+geom_violin(scale="width")+theme_bw()
#dev.off()
##############


##Correlation
#pdf("Results/PheatmapFeatures_IG07.pdf",width = 15,height =15)
#pheatmap(cor(Ranks_Driver[,c(2:10)]))
#dev.off()

Ranks_Driver_C<-Ranks_Driver[,-1]


#PCA
pc<- prcomp(Ranks_Driver[,c(2:10)],scale.=T)
pcx<- data.frame(pc$x,DrNo=Ranks_Driver[,c("DrNo")])
ggplot(pcx,aes(PC1,PC2,color=DrNo))+geom_point()



#### LEARNING AND TEST

pred.summary<- function(m,predict,golden,positive){
  true.positive=sum(m[,predict]==positive&m[,golden]==positive)
  true.negative=sum(m[,predict]!=positive&m[,golden]!=positive)
  false.positive=sum(m[,predict]==positive&m[,golden]!=positive)
  false.negative=sum(m[,predict]!=positive&m[,golden]==positive)
  ppv=round(true.positive/sum(m[,predict]==positive),4)
  npv=round(true.negative/sum(m[,predict]!=positive),4)
  accuracy=round((true.positive+true.negative)/nrow(m),4)
  sensivity=round(true.positive/sum(m[,golden]==positive),4)
  specificity=round(true.negative/sum(m[,golden]!=positive),4)
  precision=true.positive/(true.positive+false.positive)
  Recall= true.positive/(true.positive+false.negative)
  F1=(2* precision*Recall)/(precision+Recall)
  roc_C <- roc(response = as.numeric(m$Golden), predictor =as.numeric(m$Predict))
  
  auc<-roc_C$auc
  list(true.positive=true.positive,true.negative=true.negative,false.positive=false.positive,
       false.negative=false.negative,ppv=ppv,npv=npv,accuracy=accuracy,sensivity=sensivity,
       specificity=specificity,precision=precision,Recall=Recall,F1=F1,auc=auc)
}


find.accuracy<- function(train, test,algorithm){
  model<-algorithm(DrNo~ .,train)
  # model<-algorithm(DrNo~ .,train,type = 'C-classification',kernel = 'radial')
  pr<-predict(model,test)
  if(class(pr)=="list") pr<-pr$class
  x<- data.frame(Golden=test$DrNo, Predict=pr)
  
  pred.summary(x,"Predict","Golden","N")
}


data.partition<- function(data,algorithm){
  
  dataF<-data
  data<-data[,c(-1)]
  dataD<-data[data$DrNo=="D",c(-10)]
  NormalCandidated<-data[data$DrNo=="N",]
  N=nrow(NormalCandidated)
  Rand_loc_p<- (sample(seq(N),nrow(dataD)))
  normalD<-dataF[Rand_loc_p,]
  
  #Save them
  NormalDG<<-cbind(NormalDG,normalD$genes)
  NormalDPR<<-cbind(NormalDPR,normalD$WPR)
  NormalDRank<<-cbind(NormalDRank,normalD$RankWPR)
  
  prepareData<- rbind(dataD,normalD[,c(-1,-11)])
  
  Nn<-nrow(prepareData)
  trainset<-sapply(unique(prepareData$DrNo),function(s) sample(which(prepareData$DrNo==s),Nn/2,replace=T))
  trainset<- as.numeric(trainset)
  testset<-setdiff(seq(Nn),trainset)
  
  unlist(find.accuracy(prepareData[trainset,],prepareData[testset,],algorithm))
  
}



showResultsPred<- function(Re){
  svm_me<-data.frame(Re[[1]])
  lda_me<-data.frame(Re[[2]])
  RF_me<-data.frame(Re[[3]])
  
  svm_mean<- data.frame(rowMeans(svm_me))
  lda_mean<- data.frame(rowMeans(lda_me))
  RF_mean<- data.frame(rowMeans(RF_me))
  
  SLR<-cbind(svm_mean,lda_mean,RF_mean)
  colnames(SLR)<- c("svm","lda","RandomForest")
  
  svm_min<-apply(svm_me, 1, FUN = min)
  lda_min<- apply(lda_me, 1, FUN = min)
  RF_min<- apply(RF_me, 1, FUN = min)
  
  svm_max<- apply(svm_me, 1, FUN = max)
  lda_max<- apply(lda_me, 1, FUN = max)
  RF_max<- apply(RF_me, 1, FUN = max)
  
  svmDom<- cbind(data.frame(svm_min),data.frame(svm_max))
  svmDom$Do<-apply(svmDom, 1, FUN = max)
  ldaDom<- cbind(data.frame(lda_min),data.frame(lda_max))
  ldaDom$Do<-apply(ldaDom, 1, FUN = max)
  RFDom<- cbind(data.frame(RF_min),data.frame(RF_max))
  RFDom$Do<-apply(RFDom, 1, FUN = max)
  
  SLR$svmDo<-abs(SLR$svm-svmDom$Do)
  SLR$ldaDo<-abs(SLR$lda-ldaDom$Do)
  SLR$RFDo<-abs(SLR$RandomForest-RFDom$Do)
  
  return(SLR)
}

showResultsPred_One<- function(Re,colnm){
  
  Re_me<-data.frame(Re)
  Re_mean<- data.frame(rowMeans(Re_me))
  colnames(Re_mean)<- colnm
  
  Re_min<- apply(Re_me, 1, FUN = min)
  Re_max<- apply(Re_me, 1, FUN = max)
  
  ReDom<- cbind(data.frame(Re_min),data.frame(Re_max))
  ReDom$Do<-apply(ReDom, 1, FUN = max)
  
  SLR<- cbind(Re_mean,ReDom$Do)
  colnames(SLR)<- c(colnm,paste(colnm,"Dom",sep="-"))
  return(SLR)
}

#########################################################
RR<-NULL
Ranks_Driver_C$DrNo<- factor(Ranks_Driver_C$DrNo)
Ranks_Driver$DrNo<- factor(Ranks_Driver$DrNo)
#dataD<-Ranks_Driver[Ranks_Driver$DrNo=="D",]
dataD<-Ranks_Driver[Ranks_Driver$DrNo=="D",c(1,11,12)]

NormalCandidated<-Ranks_Driver[Ranks_Driver$DrNo=="N",c(1,11,12)]

NormalDG<<-data.frame(matrix(0,nrow=nrow(dataD)))
NormalDPR<<-data.frame(matrix(0,nrow=nrow(dataD)))
NormalDRank<<-data.frame(matrix(0,nrow=nrow(dataD)))


#######All layers

LayerListNameF<-c("genes","Rank_L1","Rank_L2","Rank_L3","Rank_L4","Rank_L5",
                  "Rank_L6","Rank_L7","Rank_L8","Rank_L9","WPR","RankWPR","DrNo")



# 
ac_L<-mclapply(c(randomForest), function(i) sapply(1:1000, function(c) data.partition(Ranks_Driver[,LayerListNameF],i)))

showResultsPred_One(ac_L,"RF")

Re_Sum<-data.frame(t(data.frame(ac_L)))
Re_Sum_C<-Re_Sum[,c("sensivity","specificity","accuracy","F1","auc")]

NormalDG1<-NormalDG[,-1]
write.table(NormalDG1, "DifferentStrategies/REGml/r70/NormalDG.txt", quote=F,sep="\t",row.names=FALSE)
NormalDPR1<-NormalDPR[,-1]
write.table(NormalDPR1, "DifferentStrategies/REGml/r70/NormalDPR.txt", quote=F,sep="\t",row.names=FALSE)
NormalDRank1<-NormalDRank[,-1]
write.table(NormalDRank1, "DifferentStrategies/REGml/r70/NormalDRank.txt", quote=F,sep="\t",row.names=FALSE)
write.table(dataD, "DifferentStrategies/REGml/r70/dataD.txt", quote=F,sep="\t",row.names=FALSE)
write.table(Re_Sum_C, "DifferentStrategies/REGml/r70/Re_Sum_C.txt", quote=F,sep="\t",row.names=FALSE)

