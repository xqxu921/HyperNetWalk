#install.packages("pROC")
library(pROC)
require(pROC)
require(ggpubr)
#install.packages("ggpubr")
require(ggplot2)
require(plyr)
require(gridExtra)
#install.packages("gridExtra")

#############
s_REGmlwpFC <- read.delim("1--DifferentStrategies/REGmlwpFC/r70/Re_Sum_C.txt")


my_fun1 <- function(x){
  mm <- cbind()
  sdd <- cbind()
  for(i in 1:dim(x)[2])
  {
    mm <-cbind(mm,round(mean(x[,i]),4))
    sdd <-cbind(sdd,round(sd(x[,i]),4))
  }
  zz <- rbind(mm,sdd)
  return(zz)
}

info <-c("Se","Sp","acc","F1","AUC")

r_s_REGmlwpFC <- my_fun1(s_REGmlwpFC)
re_all<- function(y,info){
  all <- y
  colnames(all) <- info
  rownames(all) <- c("mean","Std")
  return(all)
}

REGmlwpFC_info_P<-re_all(r_s_REGmlwpFC,info)

my_ROC <- function(N_pr,D_pr)
{
  require(pROC)
  ROC_cur <- NULL
  for(i in 1:dim(N_pr)[2])
  {
    data_n <- as.matrix(N_pr[,i])
    data_t <- D_pr
    label <- matrix(c(rep(1,nrow(data_t)),rep(0,nrow(data_n))))
    ##??????1,??????0 #label
    use <- c(data_t,data_n)
    k1 <- roc(as.factor(label),as.vector(use),plot=F,print.auc=F,print.thres=T,col="blue",lty=1,lwd=2)
    #rk1 <- smooth(k1,method="binormal")
    ROC_cur[[as.character(i)]] <- k1
  }
  return(ROC_cur)
}


########


REfunction<-function (NormalDPR_In,dataD_In,info_P,typ,met){
  
  p_ROC <- my_ROC(as.matrix(NormalDPR_In),as.matrix(dataD_In[,2]))
  d5 <- info_P[,5]-mean(info_P[,5])
  k5 <- which(abs(d5) == min(abs(d5)))[1]
  data_AUC <-as.matrix(info_P[,5])
  
  s5 <- smooth(p_ROC[[as.character(k5)]],method=met)
  #s5 <- smooth(r_60_ROC[[as.character(k5)]],method="logcondens.smooth")
  
  gg <- data.frame(x1=c(1-s5[["specificities"]]),
                   y1=c(s5[["sensitivities"]]),
                   type=c(rep(typ,length(1-s5[["specificities"]]))))
  
  AUC <- c(as.numeric(data_AUC[,1]))
  lab <- c(rep(typ,length(data_AUC[,1])))
  zz <- data.frame(AUC,lab)
  
  
  list(gg,zz)
  
}
#DEGmlwpFC
NormalDPRimlwpFC<-read.delim("1--DifferentStrategies/REGmlwpFC/r70/NormalDPR.txt") 
dataD_imlwpFC<-read.delim("1--DifferentStrategies/REGmlwpFC/r70/dataD.txt")

Re_DEGmlwpFC<-REfunction(NormalDPRimlwpFC,dataD_imlwpFC,s_REGmlwpFC,"REGmlwpFC","logcondens.smooth")
gg_REGmlwpFC1<- Re_DEGmlwpFC[[1]]
zz_REGmlwpFC<- Re_DEGmlwpFC[[2]]

# ggplot(data=gg_REGmlwpFC1,aes(x=x1,y=y1,group = type))+geom_line(aes(x=x1,y=y1,color = type))


##############

D_data<-read.delim( "DawnRank/DR_Ddata.txt")
N_data<-read.delim("DawnRank/DR_Ndata.txt")

my_index <- function(N_pr,D_pr)
{
  require(pROC)
  Index <- rbind()
  for(i in 1:dim(N_pr)[2])
  {
    data_n <- as.matrix(N_pr[,i])
    data_t <- D_pr
    label <- matrix(c(rep(1,nrow(data_t)),rep(0,nrow(data_n))))
    ##??????1,??????0 #label
    use <- c(data_t,data_n)
    k1 <- roc(as.factor(label),as.vector(use),plot=F,print.auc=F,print.thres=T,col="blue",lty=1,lwd=2)
    #rk1 <- smooth(k1,method="binormal")
    zz <- k1[["sensitivities"]]+k1[["specificities"]]
    k <- which(zz==max(zz))
    se_cur <- k1[["sensitivities"]][k[1]]
    sp_cur <- k1[["specificities"]][k[1]]
    auc_cur <- as.numeric(k1[["auc"]])
    acc_cur <- (sp_cur+se_cur)/2
    pre_cur <- se_cur/(se_cur+1-sp_cur)
    f1_cur <- 2*pre_cur*se_cur/(pre_cur+se_cur)
    Index <- rbind(Index,c(se_cur,sp_cur,acc_cur,f1_cur,auc_cur))
  }
  return(Index)
}

# r_085_index <- my_index(as.matrix(r_085[["Normal.pr"]]),as.matrix(r_085[["D"]][,2]))

DawnRank_index <- my_index(as.matrix(N_data),as.matrix(D_data[,2]))


my_fun1 <- function(x){
  mm <- cbind()
  sdd <- cbind()
  for(i in 1:dim(x)[2])
  {
    mm <-cbind(mm,round(mean(x[,i]),4))
    sdd <-cbind(sdd,round(sd(x[,i]),4))
  }
  zz <- rbind(mm,sdd)
  return(zz)
}

info <-c("Se","Sp","acc","F1","AUC")

DawnRank_info <- my_fun1(DawnRank_index)

#median(NetICS_info)
all <- DawnRank_info
colnames(all) <- info
dim(all)
rownames(all) <- c("mean","Std")

#write.table(all,"DawnRank/DawnRank_paramaters.txt",quote = F,sep = " ")



my_ROC <- function(N_pr,D_pr)
{
  require(pROC)
  ROC_cur <- NULL
  for(i in 1:dim(N_pr)[2])
  {
    data_n <- as.matrix(N_pr[,i])
    data_t <- D_pr
    label <- matrix(c(rep(1,nrow(data_t)),rep(0,nrow(data_n))))
    ##??????1,??????0 #label
    use <- c(data_t,data_n)
    k1 <- roc(as.factor(label),as.vector(use),plot=F,print.auc=F,print.thres=T,col="blue",lty=1,lwd=2)
    #rk1 <- smooth(k1,method="binormal")
    ROC_cur[[as.character(i)]] <- k1
  }
  return(ROC_cur)
}


DawnRank_ROC <- my_ROC(as.matrix(N_data),as.matrix(D_data[,2]))

#ɸѡ??????????


d5 <- DawnRank_index[,5]-mean(DawnRank_index[,5])
#min(abs(d4))
k5 <- which(abs(d5) == min(abs(d5)))[1]


colors1 <- c("red")
data_AUC <-as.matrix(DawnRank_index[,5])
library(pROC)
require(pROC)

s5 <- smooth(DawnRank_ROC[[as.character(k5)]],method="binormal")



require(ggpubr)
#install.packages("ggpubr")
require(ggplot2)
require(plyr)
require(gridExtra)
#install.packages("gridExtra")

gg_Dawn <- data.frame(x1=c(1-s5[["specificities"]]),
                      y1=c(s5[["sensitivities"]]),
                      type=c(rep("DawnRank",length(1-s5[["specificities"]]))))

#write.table(gg_Dawn ,"DawnRank/gg_Dawn.txt",quote = F,sep = " ")
####################

D_data<-read.delim( "IMaxDriver/iM_Ddata_Ge.txt")
N_data<-read.delim("IMaxDriver/iM_Ndata_Ge.txt")

# D_data<-read.delim( "IMaxDriver/iM_Ddata.txt")
# N_data<-read.delim("IMaxDriver/iM_Ndata.txt")


my_index <- function(N_pr,D_pr)
{
  require(pROC)
  Index <- rbind()
  for(i in 1:dim(N_pr)[2])
  {
    data_n <- as.matrix(N_pr[,i])
    data_t <- D_pr
    label <- matrix(c(rep(1,nrow(data_t)),rep(0,nrow(data_n))))
    ##??????1,??????0 #label
    use <- c(data_t,data_n)
    k1 <- roc(as.factor(label),as.vector(use),plot=F,print.auc=F,print.thres=T,col="blue",lty=1,lwd=2)
    #rk1 <- smooth(k1,method="binormal")
    zz <- k1[["sensitivities"]]+k1[["specificities"]]
    k <- which(zz==max(zz))
    se_cur <- k1[["sensitivities"]][k[1]]
    sp_cur <- k1[["specificities"]][k[1]]
    auc_cur <- as.numeric(k1[["auc"]])
    acc_cur <- (sp_cur+se_cur)/2
    pre_cur <- se_cur/(se_cur+1-sp_cur)
    f1_cur <- 2*pre_cur*se_cur/(pre_cur+se_cur)
    Index <- rbind(Index,c(se_cur,sp_cur,acc_cur,f1_cur,auc_cur))
  }
  return(Index)
}

#r_085_index <- my_index(as.matrix(r_085[["Normal.pr"]]),as.matrix(r_085[["D"]][,2]))

iMax_index <- my_index(as.matrix(N_data),as.matrix(D_data[,2]))


my_fun1 <- function(x){
  mm <- cbind()
  sdd <- cbind()
  for(i in 1:dim(x)[2])
  {
    mm <-cbind(mm,round(mean(x[,i]),4))
    sdd <-cbind(sdd,round(sd(x[,i]),4))
  }
  zz <- rbind(mm,sdd)
  return(zz)
}

info <-c("Se","Sp","acc","F1","AUC")

iMax_info <- my_fun1(iMax_index)

#median(NetICS_info)
all <- iMax_info
colnames(all) <- info
dim(all)
rownames(all) <- c("mean","Std")

#write.table(all,"IMaxDriver/IMaxDriver_paramaters.txt",quote = F,sep = " ")



my_ROC <- function(N_pr,D_pr)
{
  require(pROC)
  ROC_cur <- NULL
  for(i in 1:dim(N_pr)[2])
  {
    data_n <- as.matrix(N_pr[,i])
    data_t <- D_pr
    label <- matrix(c(rep(1,nrow(data_t)),rep(0,nrow(data_n))))
    ##??????1,??????0 #label
    use <- c(data_t,data_n)
    k1 <- roc(as.factor(label),as.vector(use),plot=F,print.auc=F,print.thres=T,col="blue",lty=1,lwd=2)
    #rk1 <- smooth(k1,method="binormal")
    ROC_cur[[as.character(i)]] <- k1
  }
  return(ROC_cur)
}


iMax_ROC <- my_ROC(as.matrix(N_data),as.matrix(D_data[,2]))

#ɸѡ??????????


d5 <- iMax_index[,5]-mean(iMax_index[,5])
#min(abs(d4))
k5 <- which(abs(d5) == min(abs(d5)))[1]


colors1 <- c("red")
data_AUC <-as.matrix(iMax_index[,5])
library(pROC)
require(pROC)

s5 <- smooth(iMax_ROC[[as.character(k5)]],method="binormal")
# s5 <- smooth(r_085_ROC[[as.character(k5)]],method="density")
# s5 <- smooth(r_085_ROC[[as.character(k5)]],method="logcondens.smooth")

require(ggpubr)
#install.packages("ggpubr")
require(ggplot2)
require(plyr)
require(gridExtra)
#install.packages("gridExtra")

gg_iMax <- data.frame(x1=c(1-s5[["specificities"]]),
                      y1=c(s5[["sensitivities"]]),
                      type=c(rep("iMaxDriver",length(1-s5[["specificities"]]))))

#write.table(gg_iMax,"IMaxDriver/gg_iMax.txt",quote = F,sep = " ")
#################################


D_data<-read.delim( "DriverNet/DN_Ddata_Ge.txt")
N_data<-read.delim("DriverNet/DN_Ndata_Ge.txt")



my_index <- function(N_pr,D_pr)
{
  require(pROC)
  Index <- rbind()
  for(i in 1:dim(N_pr)[2])
  {
    data_n <- as.matrix(N_pr[,i])
    data_t <- D_pr
    label <- matrix(c(rep(1,nrow(data_t)),rep(0,nrow(data_n))))
    ##??????1,??????0 #label
    use <- c(data_t,data_n)
    k1 <- roc(as.factor(label),as.vector(use),plot=F,print.auc=F,print.thres=T,col="blue",lty=1,lwd=2)
    #rk1 <- smooth(k1,method="binormal")
    zz <- k1[["sensitivities"]]+k1[["specificities"]]
    k <- which(zz==max(zz))
    se_cur <- k1[["sensitivities"]][k[1]]
    sp_cur <- k1[["specificities"]][k[1]]
    auc_cur <- as.numeric(k1[["auc"]])
    acc_cur <- (sp_cur+se_cur)/2
    pre_cur <- se_cur/(se_cur+1-sp_cur)
    f1_cur <- 2*pre_cur*se_cur/(pre_cur+se_cur)
    Index <- rbind(Index,c(se_cur,sp_cur,acc_cur,f1_cur,auc_cur))
  }
  return(Index)
}

#r_085_index <- my_index(as.matrix(r_085[["Normal.pr"]]),as.matrix(r_085[["D"]][,2]))

Drivernet_index <- my_index(as.matrix(N_data),as.matrix(D_data[,2]))


my_fun1 <- function(x){
  mm <- cbind()
  sdd <- cbind()
  for(i in 1:dim(x)[2])
  {
    mm <-cbind(mm,round(mean(x[,i]),4))
    sdd <-cbind(sdd,round(sd(x[,i]),4))
  }
  zz <- rbind(mm,sdd)
  return(zz)
}

info <-c("Se","Sp","acc","F1","AUC")

Drivernet_info <- my_fun1(Drivernet_index)

#median(NetICS_info)
all <- Drivernet_info
colnames(all) <- info
dim(all)
rownames(all) <- c("mean","Std")

#write.table(all,"DriverNet/DriverNet_paramaters.txt",quote = F,sep = " ")


my_ROC <- function(N_pr,D_pr)
{
  require(pROC)
  ROC_cur <- NULL
  for(i in 1:dim(N_pr)[2])
  {
    data_n <- as.matrix(N_pr[,i])
    data_t <- D_pr
    label <- matrix(c(rep(1,nrow(data_t)),rep(0,nrow(data_n))))
    ##??????1,??????0 #label
    use <- c(data_t,data_n)
    k1 <- roc(as.factor(label),as.vector(use),plot=F,print.auc=F,print.thres=T,col="blue",lty=1,lwd=2)
    #rk1 <- smooth(k1,method="binormal")
    ROC_cur[[as.character(i)]] <- k1
  }
  return(ROC_cur)
}


Drivernet_ROC <- my_ROC(as.matrix(N_data),as.matrix(D_data[,2]))

#ɸѡ??????????


d5 <- Drivernet_index[,5]-mean(Drivernet_index[,5])
#min(abs(d4))
k5 <- which(abs(d5) == min(abs(d5)))[1]



data_AUC <-as.matrix(Drivernet_index[,5])
library(pROC)
require(pROC)

s5 <- smooth(Drivernet_ROC[[as.character(k5)]],method="binormal")
#s5 <- smooth(r_085_ROC[[as.character(k5)]],method="density")
#s5 <- smooth(r_085_ROC[[as.character(k5)]],method="logcondens.smooth")

require(ggpubr)
#install.packages("ggpubr")
require(ggplot2)
require(plyr)
require(gridExtra)
#install.packages("gridExtra")

gg_DriverNet <- data.frame(x1=c(1-s5[["specificities"]]),
                           y1=c(s5[["sensitivities"]]),
                           type=c(rep("DriverNet",length(1-s5[["specificities"]]))))

#write.table(gg_DriverNet,"DriverNet/gg_DriverNet.txt",quote = F,sep = " ")

######################################

D_data<-read.delim( "MDPFinder/DN_Ddata_Ge.txt")
N_data<-read.delim("MDPFinder/DN_Ndata_Ge.txt")



my_index <- function(N_pr,D_pr)
{
  require(pROC)
  Index <- rbind()
  for(i in 1:dim(N_pr)[2])
  {
    data_n <- as.matrix(N_pr[,i])
    data_t <- D_pr
    label <- matrix(c(rep(1,nrow(data_t)),rep(0,nrow(data_n))))
    ##??????1,??????0 #label
    use <- c(data_t,data_n)
    k1 <- roc(as.factor(label),as.vector(use),plot=F,print.auc=F,print.thres=T,col="blue",lty=1,lwd=2)
    #rk1 <- smooth(k1,method="binormal")
    zz <- k1[["sensitivities"]]+k1[["specificities"]]
    k <- which(zz==max(zz))
    se_cur <- k1[["sensitivities"]][k[1]]
    sp_cur <- k1[["specificities"]][k[1]]
    auc_cur <- as.numeric(k1[["auc"]])
    acc_cur <- (sp_cur+se_cur)/2
    pre_cur <- se_cur/(se_cur+1-sp_cur)
    f1_cur <- 2*pre_cur*se_cur/(pre_cur+se_cur)
    Index <- rbind(Index,c(se_cur,sp_cur,acc_cur,f1_cur,auc_cur))
  }
  return(Index)
}

#r_085_index <- my_index(as.matrix(r_085[["Normal.pr"]]),as.matrix(r_085[["D"]][,2]))

MDPF_index <- my_index(as.matrix(N_data),as.matrix(D_data[,2]))


my_fun1 <- function(x){
  mm <- cbind()
  sdd <- cbind()
  for(i in 1:dim(x)[2])
  {
    mm <-cbind(mm,round(mean(x[,i]),4))
    sdd <-cbind(sdd,round(sd(x[,i]),4))
  }
  zz <- rbind(mm,sdd)
  return(zz)
}

info <-c("Se","Sp","acc","F1","AUC")

MDPF_info <- my_fun1(MDPF_index)

#median(NetICS_info)
all <- MDPF_info
colnames(all) <- info
dim(all)
rownames(all) <- c("mean","Std")
all
#write.table(all,"MDPFinder/MDPFinder_paramaters.txt",quote = F,sep = " ")


my_ROC <- function(N_pr,D_pr)
{
  require(pROC)
  ROC_cur <- NULL
  for(i in 1:dim(N_pr)[2])
  {
    data_n <- as.matrix(N_pr[,i])
    data_t <- D_pr
    label <- matrix(c(rep(1,nrow(data_t)),rep(0,nrow(data_n))))
    ##??????1,??????0 #label
    use <- c(data_t,data_n)
    k1 <- roc(as.factor(label),as.vector(use),plot=F,print.auc=F,print.thres=T,col="blue",lty=1,lwd=2)
    #rk1 <- smooth(k1,method="binormal")
    ROC_cur[[as.character(i)]] <- k1
  }
  return(ROC_cur)
}


MDPF_ROC <- my_ROC(as.matrix(N_data),as.matrix(D_data[,2]))

#ɸѡ??????????


d5 <- MDPF_index[,5]-mean(MDPF_index[,5])
#min(abs(d4))
k5 <- which(abs(d5) == min(abs(d5)))[1]



data_AUC <-as.matrix(MDPF_index[,5])
library(pROC)
require(pROC)

s5 <- smooth(MDPF_ROC[[as.character(k5)]],method="binormal")
#s5 <- smooth(r_085_ROC[[as.character(k5)]],method="density")
#s5 <- smooth(r_085_ROC[[as.character(k5)]],method="logcondens.smooth")

require(ggpubr)
#install.packages("ggpubr")
require(ggplot2)
require(plyr)
require(gridExtra)
#install.packages("gridExtra")

gg_MDPFinder <- data.frame(x1=c(1-s5[["specificities"]]),
                           y1=c(s5[["sensitivities"]]),
                           type=c(rep("MDPFinder",length(1-s5[["specificities"]]))))

#write.table(gg_DriverNet,"DriverNet/gg_DriverNet.txt",quote = F,sep = " ")
###################################################################

gg_OM<-rbind(gg_REGmlwpFC1,gg_Dawn,gg_iMax,gg_DriverNet,gg_MDPFinder)

#gg_OM$type <- factor(gg_OM$Type, levels = c("DawnRank", "iMaxDriver", "DriverNet"))


OMa1 <- ggplot(data=gg_OM,aes(x=x1,y=y1,group = type))+geom_line(aes(x=x1,y=y1,color = type),size = 1,linetype = 1)+labs(y="Sensitivity",x="1-Specificity",
                                                                                                                         title="" ,family = "Times", color="black")+
  theme_bw()+
  theme(legend.justification=c(1,0), legend.position=c(0.98,0.015),legend.key = element_blank())+theme(panel.border = element_blank(),
                                                                                                       panel.grid.major = element_blank(),
                                                                                                       panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))

OMa1<-OMa1+scale_color_brewer(palette="Dark2","         Method              AUC",labels=c(paste("NetWalkRank", as.character(round(REGmlwpFC_info_P[dim(REGmlwpFC_info_P)[1]-1,dim(REGmlwpFC_info_P)[2]],3)),sep="     "),
                                                                                          paste("Dawnrank", as.character(round(DawnRank_info[dim(DawnRank_info)[1]-1,dim(DawnRank_info)[2]],3)),sep="            "),
                                                                                          paste("imaxDriver", as.character(round(iMax_info[dim(iMax_info)[1]-1,dim(iMax_info)[2]],3)),sep="          "),
                                                                                          paste("DriverNet", as.character(round(Drivernet_info[dim(Drivernet_info)[1]-1,dim(Drivernet_info)[2]],3)),sep="             "),
                                                                                          paste("MDPFinder", as.character(round(MDPF_info[dim(MDPF_info)[1]-1,dim(MDPF_info)[2]],3)),sep="          ")))+
  xlim(0,1)+ylim(0,1)+ 
  theme(axis.title.x =element_text( size=10, color="black",family = "Times"),
        axis.title.y=element_text( size=10, color="black",family = "Times"),
        plot.title=element_text( size=12,family = "Times", color="black",hjust = +0.5),
        axis.line = element_line(colour = 'black', size = 0.3))+theme(legend.text = element_text(colour="black", size = 10,family = "Times"),
                                                                      legend.title = element_text( size=10,family = "Times", color="black"),
                                                                      axis.text = element_text(family = "Times",colour="black"),
                                                                      panel.border = element_rect(colour = "black", fill=NA, size=0.5))

OMa1
# pdf("OtherMethods/Plots/methods_AUC.pdf")
# OMa1
# dev.off()
# png("OtherMethods/Plots/methods_AUC.png")
# OMa1
# dev.off()


###########FOR FRAMEWORK
OMa2<- ggplot(data=gg_OM,aes(x=x1,y=y1,group = type))+geom_line(aes(x=x1,y=y1,color = type),size = 1,linetype = 1)+labs(y="",x="",
                                                                                                                        title="" ,family = "Times", color="black")+
  theme_bw()+
  theme(legend.justification=c(1,0), legend.position=c(0.98,0.015),legend.key = element_blank())+theme(panel.border = element_blank(),
                                                                                                       panel.grid.major = element_blank(),
                                                                                                       panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))

OMa2<-OMa2+scale_color_brewer(palette="Dark2","         Method              AUC")+
  xlim(0,1)+ylim(0,1)+  theme(legend.position = "none")+
  theme(axis.title.x =element_text( size=10, color="black",family = "Times"),
        axis.title.y=element_text( size=10, color="black",family = "Times"),
        plot.title=element_text( size=12,family = "Times", color="black",hjust = +0.5),
        axis.line = element_line(colour = 'black', size = 0.3))+
  OMa2

#png("0---PLOTS/methods_AUC1.png")
OMa2
#dev.off()

#######################

zztestOM<-rbind(zz_MDPFinder,zz_DriverNet,zz_imax,zz_Dawn,zz_REGmlwpFC)


zztestOM$lab<-factor(zztestOM$lab,level=c("MDPFinder","DriverNet","iMaxDriver","Dawnrank","NetWalkRank"))



OMp1<-ggplot(zztestOM,aes(x=lab,y=AUC,fill=lab))+
  geom_boxplot(position=position_dodge(0.9), alpha=1)+
  labs(y="AUC",x="",title="") +
  scale_fill_brewer(palette="Dark2",)+theme_bw()

OMp1<-OMp1+guides(fill = FALSE)
OMp1<-OMp1+theme(axis.title.x =element_text( size=10, color="black",family = "Times"),
                 axis.title.y=element_text( size=10, color="black",family = "Times"),
                 plot.title=element_text( color="black",family = "Times"),
                 axis.line = element_line(colour = 'black', size = 0.3))+theme(legend.text = element_text(colour="black", size = 10,family = "Times"),
                                                                               legend.title = element_text( size=10,family = "Times", color="black"),
                                                                               axis.text = element_text(family = "Times",colour="black"),
                                                                               axis.text.x = element_text(angle = 30,hjust =1.2 ),
                                                                               panel.border = element_rect(colour = "black", fill=NA, size=0.5))

# OMp1<-OMp1 + stat_compare_means(method = "t.test",label.y = 0.9,family="Times")
#OMp1<-OMp1 + stat_compare_means(label.y = 0.9,family="Times")
#OMp1<-OMp1+stat_compare_means(comparisons = list(c("MDPFinder","REGmlwpFC")),
#                  label = "p.signif", hide.ns = TRUE)

OMp1<-OMp1+stat_compare_means(comparisons = list(c("MDPFinder","DriverNet"),c("Dawnrank","REGmlwpFC"),c("MDPFinder","REGmlwpFC")))
OMp1<-OMp1+stat_compare_means(comparisons = list(c("MDPFinder","iMaxDriver"),c("MDPFinder","Dawnrank"),c("iMaxDriver","REGmlwpFC")) ,label = "p.signif", hide.ns = TRUE)
OMp1
# pdf("OtherMethods/Plots/methods_BoxPlot_AUC.pdf")
# OMp1
# dev.off()
# png("OtherMethods/Plots/methods_BoxPlot_AUC.png")
# OMp1
# dev.off()




lb1<-paste("A","B", sep = paste(rep(" ", 90), collapse = ""))
lb2<-paste(lb1," ", sep = paste(rep(" ", 50), collapse = ""))


combined_plots <- grid.arrange(
  arrangeGrob(top=lb2,OMa1,OMp1,ncol=2))


# 
# pdf("0---PLOTS/DifferentMethods_grid.pdf",width = 8,height =5)
# combined_plots <- grid.arrange(
#   arrangeGrob(top=lb2,OMa1,OMp1,ncol=2))
# dev.off()

