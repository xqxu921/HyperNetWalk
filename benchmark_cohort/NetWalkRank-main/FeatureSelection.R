
setwd("/home/zcx/Fatemeh/")


ResultsTable<- as.matrix(read.delim("DifferentStrategies/REGmlwpFC/r70/REGmlwpFC_r70.txt",header = TRUE ,sep="\t")) 

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



Ranks_Driver_C<-Ranks_Driver[,-1]



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
  #data<-data[,c(-1)]
  dataD<-data[data$DrNo=="D",]
  NormalCandidated<-data[data$DrNo=="N",]
  N=nrow(NormalCandidated)
  Rand_loc_p<- (sample(seq(N),nrow(dataD)))
  normalD<-dataF[Rand_loc_p,]
  
  
  prepareData<- rbind(dataD,normalD)
  
  Nn<-nrow(prepareData)
  trainset<-sapply(unique(prepareData$DrNo),function(s) sample(which(prepareData$DrNo==s),Nn/2,replace=T))
  trainset<- as.numeric(trainset)
  testset<-setdiff(seq(Nn),trainset)
  
  unlist(find.accuracy(prepareData[trainset,],prepareData[testset,],algorithm))
  
}


############################################################################################
#####SELECT FEATURES#####
check.variables<- function(data,algorithm){
  vars<- setdiff(colnames(data), "DrNo")
  sapply(vars, function(i) data.partition(data[,c(i,"DrNo")],algorithm))
}

pheatmapFUN<- function(x,FS){
  res<-x
  res[is.na(res)]<- 0
  res<-res[-1:-6,]
  pheatmap(data.frame(res),border_color=NA ,fontsize = FS,show_colnames=T,cluster_rows = F)
}

#####Check Combinations
check.comb_N<- function(data,algorithm,n){
  vars<- setdiff(colnames(data), "DrNo")
  co<- combn(vars,n)
  res<- sapply(seq(ncol(co)), function(i) data.partition(data[,c(co[,i],"DrNo")],algorithm))
  colnames(res)<- apply( co, 2, paste,collapse="_")
  res
}

####COMBINATION OF METRICS AND LAYERS
check.comb_N_plus<- function(data,algorithm,n){
  
  LayerPlus<-c("Rank_L1","Rank_L2","Rank_L3","Rank_L4","Rank_L5",
               "Rank_L6","Rank_L7","Rank_L8","Rank_L9")
  LD<-c(LayerPlus,"DrNo")
  
  vars<- setdiff(colnames(data),LD)
  co<- combn(vars,n)
  res<- sapply(seq(ncol(co)), function(i) data.partition(data[,c(co[,i],LayerPlus,"DrNo")],algorithm))
  colnames(res)<- apply( co, 2, paste,collapse="L1...9__")
  res
}

DF_result<-data.frame(matrix(ncol = 0, nrow = 7))
rownames(DF_result)<- c("accuracy", "sensivity", "specificity","precision","Recall","F1","auc")

Find.Best.Comb<- function(resc,DF_result){
  
  RMt<-rowMaxs(resc[-1:-6,])
  DF_result<-cbind(DF_result,RMt)
  return(DF_result)
}



#####Test variables

RepeadSummary1000Fun<- function(resA,n){
  res_T<-do.call(cbind,resA)
  colnL<-colnames(resA[[1]])
  repcolnL<- replicate(n,colnL)
  res_T<-data.frame(res_T)
  colnames(res_T)<-repcolnL
  res_Tmean<-t(apply(res_T,1, function(x) tapply(x,colnames(res_T),mean)))
  return(res_Tmean)
}
########## OUTPUT


colnames(Ranks_Driver_C)<- c("L1","L2","L3","L4","L5","L6","L7","L8","L9",
                             "WPR","RankWPR","BestRank","BestPR","DrNo")

res_RF_1<-lapply(1:1000, function(i) check.variables(Ranks_Driver_C[,c(-10)],randomForest))
res_RF_1_All<-RepeadSummary1000Fun(res_RF_1,1000)
pheatmapFUN(res_RF_1_All,12)


# pdf("1--DifferentStrategies/FeatureSelectionResults/pheatmap_FS_Var1.pdf",width = 15,height =15)
# pheatmapFUN(res_RF_1_All,12)
# dev.off()
# png("1--DifferentStrategies/FeatureSelectionResults/pheatmap_FS_Var1.png")
# pheatmapFUN(res_RF_1_All,15)
# dev.off()
##############################################################################
####COMBINATIONS 

###RandomForest

################  OUTPUT

res_randomF2 <-lapply(1:1000, function(i) check.comb_N(Ranks_Driver_C[,c(-10,-11,-12,-13)],randomForest,2))
res_randomF2_S <- RepeadSummary1000Fun(res_randomF2,1000)
pheatmapFUN(res_randomF2_S ,10)

pdf("1--DifferentStrategies/FeatureSelectionResults/pheatmap_FS_com2.pdf")
pheatmapFUN(res_randomF2_S ,10)
dev.off()
png("1--DifferentStrategies/FeatureSelectionResults/pheatmap_FS_com2.png")
pheatmapFUN(res_randomF2_S ,10)
dev.off()


res_randomF2_1 <-lapply(1:1000, function(i) check.comb_N(Ranks_Driver_C[,c(-10)],randomForest,2))
res_randomF2_1_S <- RepeadSummary1000Fun(res_randomF2_1,1000)
pheatmapFUN(res_randomF2_1_S ,12)

pdf("1--DifferentStrategies/FeatureSelectionResults/pheatmap_FS_com2_1.pdf",width = 20,height =15)
pheatmapFUN(res_randomF2_1_S ,12)
dev.off()
png("1--DifferentStrategies/FeatureSelectionResults/pheatmap_FS_com2_1.png")
pheatmapFUN(res_randomF2_1_S ,12)
dev.off()

########
res_randomF2_1 <-lapply(1:1000, function(i) check.comb_N(Ranks_Driver_C[,c(-10,-12,-13)],randomForest,2))
res_randomF2_1_S <- RepeadSummary1000Fun(res_randomF2_1,1000)
pheatmapFUN(res_randomF2_1_S ,12)

pdf("1--DifferentStrategies/FeatureSelectionResults/pheatmap_FS_com2_1.pdf",width = 18,height =15)
pheatmapFUN(res_randomF2_1_S ,12)
dev.off()
png("1--DifferentStrategies/FeatureSelectionResults/pheatmap_FS_com2_1.png")
pheatmapFUN(res_randomF2_1_S ,10)
dev.off()
#################

RD_C<- Ranks_Driver_C[,c(-10,-11,-12,-13)]
colnames(RD_C)<-c("L1","L2","L3","L4","L5","L6","L7","L8","L9","DrNo")
res_randomF4 <-lapply(1:1000, function(i) check.comb_N(RD_C,randomForest,4))
res_randomF4_S <- RepeadSummary1000Fun(res_randomF4,1000)
pheatmapFUN(res_randomF4_S ,12)

pdf("1--DifferentStrategies/FeatureSelectionResults/pheatmap_FS_com4.pdf",width = 25,height =12)
pheatmapFUN(res_randomF4_S ,12)
dev.off()
png("1--DifferentStrategies/FeatureSelectionResults/pheatmap_FS_com4.png")
pheatmapFUN(res_randomF4_S ,10)
dev.off()

#############
res_randomF4_1 <-lapply(1:1000, function(i) check.comb_N(Ranks_Driver_C[,c(-10,-12,-13)],randomForest,4))
res_randomF4_1_S <- RepeadSummary1000Fun(res_randomF4_1,1000)
pheatmapFUN(res_randomF4_1_S ,9)

pdf("1--DifferentStrategies/FeatureSelectionResults/pheatmap_FS_com4_1.pdf",width = 25,height =12)
pheatmapFUN(res_randomF4_1_S ,12)
dev.off()
png("1--DifferentStrategies/FeatureSelectionResults/pheatmap_FS_com4_1.png",width = 1000,height =800)
pheatmapFUN(res_randomF4_1_S ,12)
dev.off()
################

res_randomF8 <-lapply(1:1000, function(i) check.comb_N(RD_C,randomForest,8))
res_randomF8_S <- RepeadSummary1000Fun(res_randomF8,1000)
pheatmapFUN(res_randomF8_S ,12)

pdf("1--DifferentStrategies/FeatureSelectionResults/pheatmap_FS_com8.pdf",width = 15,height =15)
pheatmapFUN(res_randomF8_S ,12)
dev.off()
png("1--DifferentStrategies/FeatureSelectionResults/pheatmap_FS_com8.png")
pheatmapFUN(res_randomF8_S ,10)
dev.off()



res_randomF8 <-lapply(1:1000, function(i) check.comb_N(Ranks_Driver_C[,c(-10,-12,-13)],randomForest,8))
res_randomF8_S <- RepeadSummary1000Fun(res_randomF8,1000)
pheatmapFUN(res_randomF8_S ,12)

pdf("1--DifferentStrategies/FeatureSelectionResults/pheatmap_FS_com8.pdf",width = 25,height =12)
pheatmapFUN(res_randomF8_S ,12)
dev.off()
png("1--DifferentStrategies/FeatureSelectionResults/pheatmap_FS_com8.png")
pheatmapFUN(res_randomF8_S ,12)
dev.off()


