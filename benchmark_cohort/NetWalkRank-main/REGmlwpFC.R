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
library(e1071)
library(MASS)
library(randomForest)
library(pROC)

#set directory
setwd("/home/zcx/Fatemeh/")
#Set Properties
# setRepositories() #-- enter 1 2 to select BioC Software
# 1 2
options(timeout=500)
#

#########INPUT DATA############
#Insert weighted edges of layers by MI after using threshold 
#L1
edgeFile1<- "./Data/Results/Results_LeyerPreparation/edgeList_MI_L1.txt" 
Ledges1_0<- read.delim(edgeFile1,header = TRUE ,sep="\t")
colnames(Ledges1_0)<-c("SourceNode","TargetNode","Weight")
#L2
edgeFile2<- "./Data/Results/Results_LeyerPreparation/edgeList_MI_L2.txt" 
Ledges2_0<- read.delim(edgeFile2,header = TRUE ,sep="\t")
colnames(Ledges2_0)<-c("SourceNode","TargetNode","Weight")
#L3
edgeFile3<- "./Data/Results/Results_LeyerPreparation/edgeList_MI_L3.txt"
Ledges3_0<- read.delim(edgeFile3,header = TRUE ,sep="\t")
colnames(Ledges3_0)<-c("SourceNode","TargetNode","Weight")
#L4
edgeFile4<- "./Data/Results/Results_LeyerPreparation/edgeList_MI_L4.txt"
Ledges4_0<- read.delim(edgeFile4,header = TRUE ,sep="\t")
colnames(Ledges4_0)<-c("SourceNode","TargetNode","Weight")
#L5
edgeFile5<- "./Data/Results/Results_LeyerPreparation/edgeList_MI_L5.txt"
Ledges5_0<- read.delim(edgeFile5,header = TRUE ,sep="\t")
colnames(Ledges5_0)<-c("SourceNode","TargetNode","Weight")
#L6
edgeFile6<- "./Data/Results/Results_LeyerPreparation/edgeList_MI_L6.txt" 
Ledges6_0<- read.delim(edgeFile6,header = TRUE ,sep="\t")
colnames(Ledges6_0)<-c("SourceNode","TargetNode","Weight")
#L7
edgeFile7<- "./Data/Results/Results_LeyerPreparation/edgeList_MI_L7.txt" 
Ledges7_0<- read.delim(edgeFile7,header = TRUE ,sep="\t")
colnames(Ledges7_0)<-c("SourceNode","TargetNode","Weight")
#L8
edgeFile8<- "./Data/Results/Results_LeyerPreparation/edgeList_MI_L8.txt" 
Ledges8_0<- read.delim(edgeFile8,header = TRUE ,sep="\t")
colnames(Ledges8_0)<-c("SourceNode","TargetNode","Weight")
#L9
edgeFile9<- "./Data/Results/Results_LeyerPreparation/edgeList_MI_L9.txt" 
Ledges9_0<- read.delim(edgeFile9,header = TRUE ,sep="\t")
colnames(Ledges9_0)<-c("SourceNode","TargetNode","Weight")
#

#######
###################
#####Common Genes between Layers#####
Genes_of_Layer<- function(x){
  GL_t1<-distinct(as.data.frame(x[,1]))
  GL_t2<-distinct(as.data.frame(x[,2]))
  colnames(GL_t1)<- "G"
  colnames(GL_t2)<- "G"
  GL_t<-distinct(rbind(GL_t1,GL_t2))
  colnames(GL_t)<-"Genes"
  return(GL_t)
}

##Genes of Layers
GL1<-data.frame(Genes_of_Layer(Ledges1_0))
GL2<-data.frame(Genes_of_Layer(Ledges2_0))
GL3<-data.frame(Genes_of_Layer(Ledges3_0))
GL4<-data.frame(Genes_of_Layer(Ledges4_0))
GL5<-data.frame(Genes_of_Layer(Ledges5_0))
GL6<-data.frame(Genes_of_Layer(Ledges6_0))
GL7<-data.frame(Genes_of_Layer(Ledges7_0))
GL8<-data.frame(Genes_of_Layer(Ledges8_0))
GL9<-data.frame(Genes_of_Layer(Ledges9_0))

commonGenes_1_2 <- merge(GL1,GL2,by.x = "Genes",by.y = "Genes")
commonGenes_2_3 <- merge(GL2,GL3,by.x = "Genes",by.y = "Genes")
commonGenes_3_4 <- merge(GL3,GL4,by.x = "Genes",by.y = "Genes")
commonGenes_4_5 <- merge(GL4,GL5,by.x = "Genes",by.y = "Genes")
commonGenes_5_6 <- merge(GL5,GL6,by.x = "Genes",by.y = "Genes")
commonGenes_6_7 <- merge(GL6,GL7,by.x = "Genes",by.y = "Genes")
commonGenes_7_8 <- merge(GL7,GL8,by.x = "Genes",by.y = "Genes")
commonGenes_8_9 <- merge(GL8,GL9,by.x = "Genes",by.y = "Genes")
#
### L1_L2 ###
n12 <-nrow(commonGenes_1_2)
ct12<- c(commonGenes_1_2$Genes)
Ledges_12<- data.frame(matrix(0,ncol=3,nrow=n12))
colnames(Ledges_12)<- c("SourceNode","TargetNode","Weight")
Ledges_12$SourceNode <- (paste("1", ct12,sep="_" ))
Ledges_12$TargetNode <-  (paste("2", ct12,sep="_" ))
Ledges_12$Weight<- 1
### L2_L3 ###
n23 <-nrow(commonGenes_2_3)
ct23<- c(commonGenes_2_3$Genes)
Ledges_23<- data.frame(matrix(0,ncol=3,nrow=n23))
colnames(Ledges_23)<- c("SourceNode","TargetNode","Weight")
Ledges_23$SourceNode <- (paste("2", ct23,sep="_" ))
Ledges_23$TargetNode <-  (paste("3", ct23,sep="_" ))
Ledges_23$Weight<- 1
### L3_L4 ###
n34 <-nrow(commonGenes_3_4)
ct34<- c(commonGenes_3_4$Genes)
Ledges_34<- data.frame(matrix(0,ncol=3,nrow=n34))
colnames(Ledges_34)<- c("SourceNode","TargetNode","Weight")
Ledges_34$SourceNode <- (paste("3", ct34,sep="_" ))
Ledges_34$TargetNode <-  (paste("4", ct34,sep="_" ))
Ledges_34$Weight<- 1
### L4_L5 ###
n45 <-nrow(commonGenes_4_5)
ct45<- c(commonGenes_4_5$Genes)
Ledges_45<- data.frame(matrix(0,ncol=3,nrow=n45))
colnames(Ledges_45)<- c("SourceNode","TargetNode","Weight")
Ledges_45$SourceNode <- (paste("4", ct45,sep="_" ))
Ledges_45$TargetNode <-  (paste("5", ct45,sep="_" ))
Ledges_45$Weight<- 1
### L5_L6 ###
n56 <-nrow(commonGenes_5_6)
ct56<- c(commonGenes_5_6$Genes)
Ledges_56<- data.frame(matrix(0,ncol=3,nrow=n56))
colnames(Ledges_56)<- c("SourceNode","TargetNode","Weight")
Ledges_56$SourceNode <- (paste("5", ct56,sep="_" ))
Ledges_56$TargetNode <-  (paste("6", ct56,sep="_" ))
Ledges_56$Weight<- 1
### L6_L7 ###
n67 <-nrow(commonGenes_6_7)
ct67<- c(commonGenes_6_7$Genes)
Ledges_67<- data.frame(matrix(0,ncol=3,nrow=n67))
colnames(Ledges_67)<- c("SourceNode","TargetNode","Weight")
Ledges_67$SourceNode <- (paste("6", ct67,sep="_" ))
Ledges_67$TargetNode <-  (paste("7", ct67,sep="_" ))
Ledges_67$Weight<- 1
### L7_L8 ###
n78 <-nrow(commonGenes_7_8)
ct78<- c(commonGenes_7_8$Genes)
Ledges_78<- data.frame(matrix(0,ncol=3,nrow=n78))
colnames(Ledges_78)<- c("SourceNode","TargetNode","Weight")
Ledges_78$SourceNode <- (paste("7", ct78,sep="_" ))
Ledges_78$TargetNode <-  (paste("8", ct78,sep="_" ))
Ledges_78$Weight<- 1
### L8_L9 ###
n89 <-nrow(commonGenes_8_9)
ct89<- c(commonGenes_8_9$Genes)
Ledges_89<- data.frame(matrix(0,ncol=3,nrow=n89))
colnames(Ledges_89)<- c("SourceNode","TargetNode","Weight")
Ledges_89$SourceNode <- (paste("8", ct89,sep="_" ))
Ledges_89$TargetNode <-  (paste("9", ct89,sep="_" ))
Ledges_89$Weight<- 1

#####################################################
#######Make Transition Matrix of each Layer#########

##Genes of all layers##
Genes_Layers<-distinct(rbind(GL1,GL2,GL3,GL4,GL5,GL6,GL7,GL8,GL9))

NR<-nrow(Genes_Layers)
ColRowOrder<-as.character(Genes_Layers$Genes)

#########

fill_RRW<- function(LEdges,ext){
  notInL<- anti_join(Genes_Layers,ext)
  extraG<-data.frame(matrix(0,ncol=3,nrow=nrow(notInL)))
  colnames(extraG)<- c("SourceNode", "TargetNode" ,"Weight")
  extraG[,c(1,2,3)]<-data.frame(notInL$Genes,notInL$Genes,0) 
  
  bind_LG<- rbind(LEdges,extraG)
  Dbind_LG<- distinct(bind_LG)
  Dgraph <- graph.data.frame(Dbind_LG)
  Mgraph<-get.adjacency(Dgraph, sparse = FALSE, attr='Weight')
  Mgraph<- Mgraph[ColRowOrder,]
  Mgraph<- Mgraph[,ColRowOrder]
  return(Mgraph)
}

RW_1<-fill_RRW(Ledges1_0,GL1)
R_1<-RW_1
R_1[R_1 > 0] <- 1 
RCS_1<- R_1

RW_2<-fill_RRW(Ledges2_0,GL2)
R_2<-RW_2
R_2[R_2 > 0] <- 1 
RCS_2<- R_2

RW_3<-fill_RRW(Ledges3_0,GL3)
R_3<-RW_3
R_3[R_3 > 0] <- 1 
RCS_3<- R_3

RW_4<-fill_RRW(Ledges4_0,GL4)
R_4<-RW_4
R_4[R_4 > 0] <- 1 
RCS_4<- R_4

RW_5<-fill_RRW(Ledges5_0,GL5)
R_5<-RW_5
R_5[R_5 > 0] <- 1 
RCS_5<- R_5

RW_6<-fill_RRW(Ledges6_0,GL6)
R_6<-RW_6
R_6[R_6 > 0] <- 1 
RCS_6<- R_6

RW_7<-fill_RRW(Ledges7_0,GL7)
R_7<-RW_7
R_7[R_7 > 0] <- 1 
RCS_7<- R_7

RW_8<-fill_RRW(Ledges8_0,GL8)
R_8<-RW_8
R_8[R_8 > 0] <- 1 
RCS_8<- R_8

RW_9<-fill_RRW(Ledges9_0,GL9)
R_9<-RW_9
R_9[R_9 > 0] <- 1 
RCS_9<- R_9

#
P_edges<- function(y){
  y_0<- y
  y_0$SourceNode<-sub(".*_","",y$SourceNode)
  y_0$TargetNode<-sub(".*_","",y$TargetNode)
  y_0<-distinct(y_0)
  return(y_0)
}

#
RW_12<-fill_RRW(P_edges(Ledges_12),commonGenes_1_2)
R_12<-RW_12
R_12[R_12 > 0] <- 1 
RCS_12<- R_12

RW_23<-fill_RRW(P_edges(Ledges_23),commonGenes_2_3)
R_23<-RW_23
R_23[R_23 > 0] <- 1 
RCS_23<- R_23

RW_34<-fill_RRW(P_edges(Ledges_34),commonGenes_3_4)
R_34<-RW_34
R_34[R_34 > 0] <- 1 
RCS_34<- R_34

RW_45<-fill_RRW(P_edges(Ledges_45),commonGenes_4_5)
R_45<-RW_45
R_45[R_45 > 0] <- 1 
RCS_45<- R_45

RW_56<-fill_RRW(P_edges(Ledges_56),commonGenes_5_6)
R_56<-RW_56
R_56[R_56 > 0] <- 1 
RCS_56<- R_56

RW_67<-fill_RRW(P_edges(Ledges_67),commonGenes_6_7)
R_67<-RW_67
R_67[R_67 > 0] <- 1 
RCS_67<- R_67

RW_78<-fill_RRW(P_edges(Ledges_78),commonGenes_7_8)
R_78<-RW_78
R_78[R_78 > 0] <- 1 
RCS_78<- R_78

RW_89<-fill_RRW(P_edges(Ledges_89),commonGenes_8_9)
R_89<-RW_89
R_89[R_89 > 0] <- 1 
RCS_89<- R_89
#####

Cal_W_L<- function(Rx,RWx){
  R_Sum<-sum(Rx)
  RW_Sum<-sum(RWx)
  
  if(R_Sum>0){
    L_W<- (RW_Sum/R_Sum)  
  }
  else{
    L_W<- 0
  }
  return(L_W)
}

L_W<- c(1:9)*0

#Calculate Weight of Layers by Call Function
L_W[1]<- Cal_W_L(R_1,RW_1)
L_W[2]<- Cal_W_L(R_2,RW_2)
L_W[3]<- Cal_W_L(R_3,RW_3)
L_W[4]<- Cal_W_L(R_4,RW_4)
L_W[5]<- Cal_W_L(R_5,RW_5)
L_W[6]<- Cal_W_L(R_6,RW_6)
L_W[7]<- Cal_W_L(R_7,RW_7)
L_W[8]<- Cal_W_L(R_8,RW_8)
L_W[9]<- Cal_W_L(R_9,RW_9)


sum_LW<-sum(L_W)
#Normalized weight of Layers
Nor_L_W<- ((L_W)/sum_LW)
sumNW<-sum(Nor_L_W)

L_W
Nor_L_W

coefW<- Nor_L_W

########
ColSum2<-function(x){
  #x<-as.matrix(x)
  if(sum(x)==0){
    return(1)
  }
  else{
    return(sum(x))
  }
}

RW_D_R<- function(R_t,RW_t,k){
  sum_value <- apply(R_t,2,ColSum2)
  R_t1<-t(t(R_t)/sum_value)
  RW_t1<- (coefW[k]*(R_t*RW_t))
  return(list(R_t1, RW_t1))
}

RWDR<- RW_D_R(R_1,RW_1,1)
Rd1 <-RWDR[[1]]
RWd1 <-RWDR[[2]]

RWDR<- RW_D_R(R_2,RW_2,2)
Rd2 <-RWDR[[1]]
RWd2 <-RWDR[[2]]

RWDR<- RW_D_R(R_3,RW_3,3)
Rd3 <-RWDR[[1]]
RWd3 <-RWDR[[2]]

RWDR<- RW_D_R(R_4,RW_4,4)
Rd4 <-RWDR[[1]]
RWd4 <-RWDR[[2]]

RWDR<- RW_D_R(R_5,RW_5,5)
Rd5 <-RWDR[[1]]
RWd5 <-RWDR[[2]]

RWDR<- RW_D_R(R_6,RW_6,6)
Rd6 <-RWDR[[1]]
RWd6 <-RWDR[[2]]

RWDR<- RW_D_R(R_7,RW_7,7)
Rd7 <-RWDR[[1]]
RWd7 <-RWDR[[2]]

RWDR<- RW_D_R(R_8,RW_8,8)
Rd8 <-RWDR[[1]]
RWd8 <-RWDR[[2]]

RWDR<- RW_D_R(R_9,RW_9,9)
Rd9 <-RWDR[[1]]
RWd9 <-RWDR[[2]]

# 


RW_D_R_T<- function(R_t,RW_t,k){
  
  sum_value <- apply(R_t,2,ColSum2)
  R_t1<-t(t(R_t)/sum_value)
  RW_t1<- ((((coefW[k]+coefW[k+1]))/2)*(R_t*RW_t))
  
  return(list(R_t1, RW_t1))
}

RWDRT<- RW_D_R_T(R_12,RW_12,1)
RTd12 <-RWDRT[[1]]
RWTd12 <-RWDRT[[2]]

RWDRT<- RW_D_R_T(R_23,RW_23,2)
RTd23 <-RWDRT[[1]]
RWTd23 <-RWDRT[[2]]

RWDRT<- RW_D_R_T(R_34,RW_34,3)
RTd34 <-RWDRT[[1]]
RWTd34 <-RWDRT[[2]]

RWDRT<- RW_D_R_T(R_45,RW_45,4)
RTd45 <-RWDRT[[1]]
RWTd45 <-RWDRT[[2]]

RWDRT<- RW_D_R_T(R_56,RW_56,5)
RTd56 <-RWDRT[[1]]
RWTd56 <-RWDRT[[2]]

RWDRT<- RW_D_R_T(R_67,RW_67,6)
RTd67 <-RWDRT[[1]]
RWTd67 <-RWDRT[[2]]

RWDRT<- RW_D_R_T(R_78,RW_78,7)
RTd78 <-RWDRT[[1]]
RWTd78 <-RWDRT[[2]]

RWDRT<- RW_D_R_T(R_89,RW_89,8)
RTd89 <-RWDRT[[1]]
RWTd89 <-RWDRT[[2]]

####


NG<- 9*NR


#########Personalilzed PageRank
fill_s0_EL<- function(xG){
  nrg<-nrow(xG)
  xtemp<-data.frame(matrix(0,ncol=2,nrow=NR))
  colnames(xtemp)<-c("Genes","Val")
  rownames(xtemp)<- ColRowOrder
  xtemp$Genes<- ColRowOrder
  xtemp[xtemp$Genes %in% xG$Genes,2]<- 1/nrg
  return(xtemp$Val)
}
s0_initial_1 <- fill_s0_EL(GL1)
s0_initial_2 <- fill_s0_EL(GL2)
s0_initial_3 <- fill_s0_EL(GL3)
s0_initial_4 <- fill_s0_EL(GL4)
s0_initial_5 <- fill_s0_EL(GL5)
s0_initial_6 <- fill_s0_EL(GL6)
s0_initial_7 <- fill_s0_EL(GL7)
s0_initial_8 <- fill_s0_EL(GL8)
s0_initial_9 <- fill_s0_EL(GL9)
# 



#All together in one matrix
Rdt_All<-matrix(0,nrow=NG,ncol = NG)
RWdt_All<-matrix(0,nrow=NG,ncol = NG)
RCS_All<-matrix(0,nrow=NG,ncol = NG)

fill_RRW_ALL<- function(Rr,RWrw,rcsx,k,LB){
  if(LB==1){
    sPointR<-((k-1)*NR)+1
    sPointC<-((k-1)*NR)+1
    print(sPointR)
    print((sPointR+NR-1))
  }
  if(LB==2){
    sPointR<-((k-1)*NR)+1
    sPointC<-((k)*NR)+1
    print(sPointR)
    print((sPointR+NR-1))
  }
  if(LB==3){
    sPointR<-((k)*NR)+1
    sPointC<-((k-1)*NR)+1
    print(sPointR)
    print((sPointR+NR-1))
  }
  #Rdt_All[sPointR:(sPointR+NR-1),sPointC:(sPointC+NR-1)] <<- Rr
  RWdt_All[sPointR:(sPointR+NR-1),sPointC:(sPointC+NR-1)] <<- RWrw
  RCS_All[sPointR:(sPointR+NR-1),sPointC:(sPointC+NR-1)] <<- rcsx
}
fill_RRW_ALL(Rd1,RWd1,RCS_1,1,1)
fill_RRW_ALL(Rd2,RWd2,RCS_2,2,1)
fill_RRW_ALL(Rd3,RWd3,RCS_3,3,1)
fill_RRW_ALL(Rd4,RWd4,RCS_4,4,1)
fill_RRW_ALL(Rd5,RWd5,RCS_5,5,1)
fill_RRW_ALL(Rd6,RWd6,RCS_6,6,1)
fill_RRW_ALL(Rd7,RWd7,RCS_7,7,1)
fill_RRW_ALL(Rd8,RWd8,RCS_8,8,1)
fill_RRW_ALL(Rd9,RWd9,RCS_9,9,1)

fill_RRW_ALL(RTd12,RWTd12,RCS_12,1,2)
fill_RRW_ALL(RTd12,RWTd12,RCS_12,1,3)

fill_RRW_ALL(RTd23,RWTd23,RCS_23,2,2)
fill_RRW_ALL(RTd23,RWTd23,RCS_23,2,3)

fill_RRW_ALL(RTd34,RWTd34,RCS_34,3,2)
fill_RRW_ALL(RTd34,RWTd34,RCS_34,3,3)

fill_RRW_ALL(RTd45,RWTd45,RCS_45,4,2)
fill_RRW_ALL(RTd45,RWTd45,RCS_45,4,3)

fill_RRW_ALL(RTd56,RWTd56,RCS_56,5,2)
fill_RRW_ALL(RTd56,RWTd56,RCS_56,5,3)

fill_RRW_ALL(RTd67,RWTd67,RCS_67,6,2)
fill_RRW_ALL(RTd67,RWTd67,RCS_67,6,3)

fill_RRW_ALL(RTd78,RWTd78,RCS_78,7,2)
fill_RRW_ALL(RTd78,RWTd78,RCS_78,7,3)

fill_RRW_ALL(RTd89,RWTd89,RCS_89,8,2)
fill_RRW_ALL(RTd89,RWTd89,RCS_89,8,3)



s0_initial_All<-c(s0_initial_1,s0_initial_2,s0_initial_3,s0_initial_4,s0_initial_5,s0_initial_6,s0_initial_7,
                  s0_initial_8,s0_initial_9)

###PageRank###
##Use All Matrix##
CS<-RCS_All
Trans<- RWdt_All
r=0.7
threshold=10^(-10)
restart =as.matrix(s0_initial_All,ncol=1,nrow=NG)
PR=restart
iter = 1
delta_PR = Inf

while ((delta_PR > threshold) & (iter<200)) {
  prev_PR <- PR
  CST<- CS*Trans
  rCST<-(r*CST)
  PR <<- ((rCST) %*% (PR)) + ((1-r)*(restart))
  delta_PR <- norm((PR-prev_PR),type ="F" )

  print(iter)
  print(delta_PR)

  iter = iter + 1
}


PR_AllLayers_REGmlwpFC70<-PR

PR1<-as.matrix(PR_AllLayers_REGmlwpFC70[,2],ncol=1)
rownames(PR1)<-rownames(PR_AllLayers_REGmlwpFC70)
PR2<-as.matrix(PR_AllLayers_REGmlwpFC70[,3],ncol=1)
rownames(PR2)<-rownames(PR_AllLayers_REGmlwpFC70)
PR3<-as.matrix(PR_AllLayers_REGmlwpFC70[,4],ncol=1)
rownames(PR3)<-rownames(PR_AllLayers_REGmlwpFC70)
PR4<-as.matrix(PR_AllLayers_REGmlwpFC70[,5],ncol=1)
rownames(PR4)<-rownames(PR_AllLayers_REGmlwpFC70)
PR5<-as.matrix(PR_AllLayers_REGmlwpFC70[,6],ncol=1)
rownames(PR5)<-rownames(PR_AllLayers_REGmlwpFC70)
PR6<-as.matrix(PR_AllLayers_REGmlwpFC70[,7],ncol=1)
rownames(PR6)<-rownames(PR_AllLayers_REGmlwpFC70)
PR7<-as.matrix(PR_AllLayers_REGmlwpFC70[,8],ncol=1)
rownames(PR7)<-rownames(PR_AllLayers_REGmlwpFC70)
PR8<-as.matrix(PR_AllLayers_REGmlwpFC70[,9],ncol=1)
rownames(PR8)<-rownames(PR_AllLayers_REGmlwpFC70)
PR9<-as.matrix(PR_AllLayers_REGmlwpFC70[,10],ncol=1)
rownames(PR9)<-rownames(PR_AllLayers_REGmlwpFC70)
###################################

####Input LFC information

LFC_Layers<-read.delim("DifferentStrategies/LFC_Layers/LFC_Layers.txt",header = FALSE)
LFC_info<-LFC_Layers[LFC_Layers$V1 %in% ColRowOrder,]



LFC_C <- function(infox,nx){
  LFCbind<- cbind(infox[,1],abs(as.numeric(abs(infox[,nx]))))
  LFCbind1<- data.frame(LFCbind)
  colnames(LFCbind1)<-c("Genes","LFC")
  LFCbind1$LFC[as.numeric(LFCbind1$LFC)<1]<- 1
  print(dim(LFCbind1[LFCbind1$LFC>1,]))
  ColRowOrderDF<-data.frame(ColRowOrder)
  LFCbind2<-LFCbind1[LFCbind1$Genes %in% ColRowOrder,]
  
  dim(LFCbind2)
  rownames(LFCbind2)<-LFCbind2[,1]
  LFCbind3<- LFCbind2[ColRowOrder,]
  print(dim(LFCbind3))
  return(LFCbind3)
}


LFC_L1<-LFC_C(LFC_info,2)
LFC_L2<-LFC_C(LFC_info,3)
LFC_L3<-LFC_C(LFC_info,4)
LFC_L4<-LFC_C(LFC_info,5)
LFC_L5<-LFC_C(LFC_info,6)
LFC_L6<-LFC_C(LFC_info,7)
LFC_L7<-LFC_C(LFC_info,8)
LFC_L8<-LFC_C(LFC_info,9)
LFC_L9<-LFC_C(LFC_info,10)


head(LFC_L1)
head(PR1)

PRandLFC<- function(PRd,LFCd){
  LFCdata<-as.matrix(LFCd[,2],ncol=1)
  prlfc<-data.frame(cbind(PRd,LFCdata))
  
  prlfc$m<-as.numeric(prlfc[,1])*as.numeric(prlfc[,2])
  print(dim(prlfc[prlfc[,2]>1,]))
  colnames(prlfc)<-c("PR","LFCn","NewPR")
  return (prlfc)
  
}
new_PR1<- PRandLFC(PR1,LFC_L1)
new_PR2<- PRandLFC(PR2,LFC_L2)
new_PR3<- PRandLFC(PR3,LFC_L3)
new_PR4<- PRandLFC(PR4,LFC_L4)
new_PR5<- PRandLFC(PR5,LFC_L5)
new_PR6<- PRandLFC(PR6,LFC_L6)
new_PR7<- PRandLFC(PR7,LFC_L7)
new_PR8<- PRandLFC(PR8,LFC_L8)

PR1<-as.matrix(new_PR1$NewPR,ncol=1)
PR2<-as.matrix(new_PR2$NewPR,ncol=1)
PR3<-as.matrix(new_PR3$NewPR,ncol=1)
PR4<-as.matrix(new_PR4$NewPR,ncol=1)
PR5<-as.matrix(new_PR5$NewPR,ncol=1)
PR6<-as.matrix(new_PR6$NewPR,ncol=1)
PR7<-as.matrix(new_PR7$NewPR,ncol=1)
PR8<-as.matrix(new_PR8$NewPR,ncol=1)



PR_Layers_All<-data.frame(Gene=numeric(),PR_Value=numeric())

fill_PR_Layers_All<- function(k,PR_layer){
  PR_Layers_All[(((k-1)*NR)+1):(k*NR),1]<<-(paste(k,Genes_Layers[,1] ,sep="_" ))
  PR_Layers_All[(((k-1)*NR)+1):(k*NR),2]<<-PR_layer[,1]
}

fill_PR_Layers_All(1,PR1)
fill_PR_Layers_All(2,PR2)
fill_PR_Layers_All(3,PR3)
fill_PR_Layers_All(4,PR4)
fill_PR_Layers_All(5,PR5)
fill_PR_Layers_All(6,PR6)
fill_PR_Layers_All(7,PR7)
fill_PR_Layers_All(8,PR8)
fill_PR_Layers_All(9,PR9)


#######
#Sort Unique vales of PR in whole network
PR_sort_Rank<- data.frame(sort(unique(PR_Layers_All[,2]),decreasing = TRUE))
colnames(PR_sort_Rank)<- "PR"
PR_sort_Rank$Rank<- 1:length(PR_sort_Rank[,1])

##PR values of genes based on all layers 

PR_sort_Ls<- cbind(PR1,PR2,PR3,PR4,PR5,PR6,PR7,PR8,PR9)
sum(PR_sort_Ls ==0)

PR_sort_Ls_Df<-data.frame(PR_sort_Ls)


PR_sort_Ls_Df_LA<-cbind(Genes_Layers[,1],PR_sort_Ls_Df)
colnames(PR_sort_Ls_Df_LA)<-c("Genes","PR","PR","PR","PR","PR","PR","PR","PR","PR")

PR_Values_Round<-round(PR_sort_Ls_Df_LA[,-1],digits = 5)
PR_values_Layers<- cbind(PR_sort_Ls_Df_LA[,1],PR_Values_Round)


##sort and reorder for PR of Layers
Merge_PR_Layer_Rorder<- function(px){
  genes_PR_Rank_x<- merge(px[,c("Genes","PR")],PR_sort_Rank, by.x="PR",by.y="PR")
  genes_PR_Rank_x<- genes_PR_Rank_x[,c(2,3,1)]
  colnames(genes_PR_Rank_x)<-c("genes","Rank","PR")
  rownames(genes_PR_Rank_x)<-genes_PR_Rank_x[,1]
  genes_PR_Rank_x<- genes_PR_Rank_x[ColRowOrder,]
  return(genes_PR_Rank_x)
}
#Layers-Genes-Ranks-PRs
genes_PR_Rank1<-Merge_PR_Layer_Rorder(PR_sort_Ls_Df_LA[,c(1,2)])
genes_PR_Rank2<-Merge_PR_Layer_Rorder(PR_sort_Ls_Df_LA[,c(1,3)])
genes_PR_Rank3<-Merge_PR_Layer_Rorder(PR_sort_Ls_Df_LA[,c(1,4)])
genes_PR_Rank4<-Merge_PR_Layer_Rorder(PR_sort_Ls_Df_LA[,c(1,5)])
genes_PR_Rank5<-Merge_PR_Layer_Rorder(PR_sort_Ls_Df_LA[,c(1,6)])
genes_PR_Rank6<-Merge_PR_Layer_Rorder(PR_sort_Ls_Df_LA[,c(1,7)])
genes_PR_Rank7<-Merge_PR_Layer_Rorder(PR_sort_Ls_Df_LA[,c(1,8)])
genes_PR_Rank8<-Merge_PR_Layer_Rorder(PR_sort_Ls_Df_LA[,c(1,9)])
genes_PR_Rank9<-Merge_PR_Layer_Rorder(PR_sort_Ls_Df_LA[,c(1,10)])


########################################
genes_PR_L_Ranks<- cbind(genes_PR_Rank1[,c(1,2)],genes_PR_Rank2[,2],
                         genes_PR_Rank3[,2],genes_PR_Rank4[,2],
                         genes_PR_Rank5[,2],genes_PR_Rank6[,2],
                         genes_PR_Rank7[,2],genes_PR_Rank8[,2],
                         genes_PR_Rank9[,2])

colnames(genes_PR_L_Ranks)<-c("genes","Rank_L1","Rank_L2","Rank_L3","Rank_L4","Rank_L5",
                              "Rank_L6","Rank_L7","Rank_L8","Rank_L9")



genes_PR_Ranks_D<-cbind(PR_sort_Ls_Df_LA,genes_PR_L_Ranks[,-1])
colnames(genes_PR_Ranks_D)<-c("genes","PR_L1","PR_L2","PR_L3","PR_L4","PR_L5",
                              "PR_L6","PR_L7","PR_L8","PR_L9",
                              "Rank_L1","Rank_L2","Rank_L3","Rank_L4","Rank_L5",
                              "Rank_L6","Rank_L7","Rank_L8","Rank_L9")

genes_Ranks<-genes_PR_L_Ranks

##Calculate Weithed PR and Best
weight_L<- c(Nor_L_W)
library(matrixStats)
########
Sort_PR_Layers<- function(PRL){
  PR_sort_L<- data.frame(sort(unique(PRL[,1]),decreasing = TRUE))
  colnames(PR_sort_L)<- "PR"
  PR_sort_L$Rank<- 1:length(PR_sort_L[,1])
  return(PR_sort_L)
}
########
rowWeMean<- rowWeightedMeans(PR_sort_Ls, w = weight_L)
PR_sort_Ls_Df<-data.frame(PR_sort_Ls)
PR_sort_Ls_Df$WMean<-rowWeMean
PR_sort_Ls_Df_LA<-cbind(Genes_Layers[,1],PR_sort_Ls_Df)
colnames(PR_sort_Ls_Df_LA)<-c("Genes","PR","PR","PR","PR","PR","PR","PR","PR","PR","WPR")
PR_Values_Round<-round(PR_sort_Ls_Df_LA[,-1],digits = 5)
PR_values_Layers<- cbind(PR_sort_Ls_Df_LA[,1],PR_Values_Round)
PR_sort_All<- Sort_PR_Layers(matrix(PR_sort_Ls_Df$WMean,ncol=1))
genes_PR_Rank<- merge(PR_sort_Ls_Df_LA[,c("Genes","WPR")],PR_sort_All, by.x="WPR",by.y="PR")
genes_PR_Rank_WA<- genes_PR_Rank[,c(2,1,3)]
colnames(genes_PR_Rank_WA)<-c("genes","WPR","RankWPR")
rownames(genes_PR_Rank_WA)<-genes_PR_Rank_WA[,1]
genes_PR_Rank_WA<- genes_PR_Rank_WA[ColRowOrder,]

genes_PR_Ranks_D<-cbind(genes_PR_Ranks_D,genes_PR_Rank_WA[,2:3])

##Best
genes_Ranks<-genes_PR_L_Ranks
genes_Ranks$Best<-pmin(genes_Ranks$Rank_L1,genes_Ranks$Rank_L2,
                       genes_Ranks$Rank_L3,genes_Ranks$Rank_L4,
                       genes_Ranks$Rank_L5,genes_Ranks$Rank_L6,
                       genes_Ranks$Rank_L7,genes_Ranks$Rank_L8,
                       genes_Ranks$Rank_L9)

genes_PR_Ranks_D<-cbind(genes_PR_Ranks_D,Best=genes_Ranks$Best)
#####
genes_Ranks$BestPR<-pmin(genes_PR_Ranks_D$PR_L1,genes_PR_Ranks_D$PR_L2,
                         genes_PR_Ranks_D$PR_L3,genes_PR_Ranks_D$PR_L4,
                         genes_PR_Ranks_D$PR_L5,genes_PR_Ranks_D$PR_L6,
                         genes_PR_Ranks_D$PR_L7,genes_PR_Ranks_D$PR_L8,
                         genes_PR_Ranks_D$PR_L9)

BestPR_sort_All<- Sort_PR_Layers(matrix(genes_Ranks$BestPR,ncol=1))
BestPR_sort_All<- merge(genes_Ranks[,c("genes","BestPR")],BestPR_sort_All, by.x="BestPR",by.y="PR")
genes_PR_Rank_BPR<- BestPR_sort_All[,c(2,1,3)]
colnames(genes_PR_Rank_BPR)<-c("genes","BestPR","RankBPR")
rownames(genes_PR_Rank_BPR)<-genes_PR_Rank_BPR[,1]
genes_PR_Rank_BPR<- genes_PR_Rank_BPR[ColRowOrder,]

genes_PR_Ranks_D<-cbind(genes_PR_Ranks_D,BestPR_rank=genes_PR_Rank_BPR$RankBPR)
#####

write.table(genes_PR_Ranks_D, "DifferentStrategies/REGmlwpFC/r70/REGmlwpFC_r70.txt", quote=F,sep="\t",row.names=FALSE)
############################################


