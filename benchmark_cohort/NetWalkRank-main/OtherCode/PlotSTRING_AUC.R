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

####ADD REG


s_REGm <- read.delim("OtherNetworks/Diff_str/1--STRING_re/m/Re_Sum_C.txt")
s_REGml <- read.delim("OtherNetworks/Diff_str/1--STRING_re/ml/Re_Sum_C.txt")
s_REGmlw <- read.delim("OtherNetworks/Diff_str/1--STRING_re/mlw/Re_Sum_C.txt")
s_REGmlwp <- read.delim("OtherNetworks/Diff_str/1--STRING_re/mlwp/Re_Sum_C.txt")
s_REGmlwpFC <- read.delim("OtherNetworks/Diff_str/1--STRING_re/mlwpFC/Re_Sum_C.txt")


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

r_s_REGm <- my_fun1(s_REGm)
r_s_REGml <- my_fun1(s_REGml)
r_s_REGmlw <- my_fun1(s_REGmlw)
r_s_REGmlwp <- my_fun1(s_REGmlwp)
r_s_REGmlwpFC <- my_fun1(s_REGmlwpFC)

re_all<- function(y,info){
  all <- y
  colnames(all) <- info
  rownames(all) <- c("mean","Std")
  return(all)
}

REGm_info_P<-re_all(r_s_REGm,info)
REGml_info_P<-re_all(r_s_REGml,info)
REGmlw_info_P<-re_all(r_s_REGmlw,info)
REGmlwp_info_P<-re_all(r_s_REGmlwp,info)
REGmlwpFC_info_P<-re_all(r_s_REGmlwpFC,info)

#write.table(REGmlwpFC_info_P,"1--DifferentStrategies/Plots/REGmlwpFC_paramaters.txt",quote = F,sep = " ")



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

#DEGm
NormalDPRim<-read.delim("OtherNetworks/Diff_str/1--STRING_re/m/NormalDPR.txt") 
dataD_im<-read.delim("OtherNetworks/Diff_str/1--STRING_re/m/dataD.txt")

Re_DEGm<-REfunction(NormalDPRim,dataD_im,s_REGm,"REGm","binormal")
gg_REGm1<- Re_DEGm[[1]]
zz_REGm<- Re_DEGm[[2]]

ggplot(data=gg_REGm1,aes(x=x1,y=y1,group = type))+geom_line(aes(x=x1,y=y1,color = type))

#DEGml
NormalDPRil<-read.delim("OtherNetworks/Diff_str/1--STRING_re/ml/NormalDPR.txt") 
dataD_il<-read.delim("OtherNetworks/Diff_str/1--STRING_re/ml/dataD.txt")

Re_DEGml<-REfunction(NormalDPRil,dataD_il,s_REGml,"REGml","binormal")
gg_REGml1<- Re_DEGml[[1]]
zz_REGml<- Re_DEGml[[2]]

ggplot(data=gg_REGml1,aes(x=x1,y=y1,group = type))+geom_line(aes(x=x1,y=y1,color = type))

#DEGmlw
NormalDPRimlw<-read.delim("OtherNetworks/Diff_str/1--STRING_re/mlw/NormalDPR.txt") 
dataD_imlw<-read.delim("OtherNetworks/Diff_str/1--STRING_re/mlw/dataD.txt")

Re_DEGmlw<-REfunction(NormalDPRimlw,dataD_imlw,s_REGmlw,"REGmlw","binormal")
gg_REGmlw1<- Re_DEGmlw[[1]]
zz_REGmlw<- Re_DEGmlw[[2]]

ggplot(data=gg_REGmlw1,aes(x=x1,y=y1,group = type))+geom_line(aes(x=x1,y=y1,color = type))

#DEGmlwp
NormalDPRimlwp<-read.delim("OtherNetworks/Diff_str/1--STRING_re/mlwp/NormalDPR.txt") 
dataD_imlwp<-read.delim("OtherNetworks/Diff_str/1--STRING_re/mlwp/dataD.txt")

Re_DEGmlwp<-REfunction(NormalDPRimlwp,dataD_imlwp,s_REGmlwp,"REGmlwp","logcondens.smooth")
gg_REGmlwp1<- Re_DEGmlwp[[1]]
zz_REGmlwp<- Re_DEGmlwp[[2]]

ggplot(data=gg_REGmlwp1,aes(x=x1,y=y1,group = type))+geom_line(aes(x=x1,y=y1,color = type))

#DEGmlwpFC
NormalDPRimlwpFC<-read.delim("OtherNetworks/Diff_str/1--STRING_re/mlwpFC/NormalDPR.txt") 
dataD_imlwpFC<-read.delim("OtherNetworks/Diff_str/1--STRING_re/mlwpFC/dataD.txt")

Re_DEGmlwpFC<-REfunction(NormalDPRimlwpFC,dataD_imlwpFC,s_REGmlwpFC,"REGmlwpFC","logcondens.smooth")
gg_REGmlwpFC1<- Re_DEGmlwpFC[[1]]
zz_REGmlwpFC<- Re_DEGmlwpFC[[2]]

ggplot(data=gg_REGmlwpFC1,aes(x=x1,y=y1,group = type))+geom_line(aes(x=x1,y=y1,color = type))


ggtest<-rbind(gg_REGml1,gg_REGmlw1,gg_REGmlwp1,gg_REGmlwpFC1)



STRa1 <- ggplot(data=ggtest,aes(x=x1,y=y1,group = type))+geom_line(aes(x=x1,y=y1,color = type),size = 1,linetype = 1)+labs(y="Sensitivity",x="1-Specificity",
                                                                                                                           title="A" ,family = "Times", color="black",face="bold")+
  theme_bw()+
  theme(legend.justification=c(1,0), legend.position=c(1,0),legend.key = element_blank())+theme(panel.border = element_blank(),
                                                                                                panel.grid.major = element_blank(),
                                                                                                panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
  scale_color_manual(values=c('steelblue','hotpink2',"#99cc00","gold1") ,"    Parameter              AUC",labels=c(paste("REGml", as.character(round(REGm_info_P[dim(REGm_info_P)[1]-1,dim(REGm_info_P)[2]],3)),sep="              "),
                                                                                                                   paste("REGmlw", as.character(round(REGml_info_P[dim(REGml_info_P)[1]-1,dim(REGml_info_P)[2]],3)),sep="           "),
                                                                                                                   paste("REGmlwp", as.character(round(REGmlw_info_P[dim(REGmlw_info_P)[1]-1,dim(REGmlw_info_P)[2]],3)),sep="          "),
                                                                                                                   paste("REGmlwpFC", as.character(round(REGmlwp_info_P[dim(REGmlwp_info_P)[1]-1,dim(REGmlwp_info_P)[2]],3)),sep="      ")))+
  xlim(0,1)+ylim(0,1)+ 
  theme(axis.title.x =element_text( size=10, color="black",face="bold"),
        axis.title.y=element_text( size=10, color="black",face="bold"),
        plot.title=element_text( size=12,family = "Times", color="black",face="bold",hjust = -0.007),
        axis.line = element_line(colour = 'black', size = 0.3))+theme(legend.text = element_text(colour="black", size = 10,family = "Times"),
                                                                      legend.title = element_text( size=10,family = "Times", color="black"),
                                                                      axis.text = element_text(family = "Times",colour="black"))

STRa1


#####################
STRa1 <- ggplot(data=ggtest,aes(x=x1,y=y1,group = type))+geom_line(aes(x=x1,y=y1,color = type),size = 1,linetype = 1)+labs(y="Sensitivity",x="1-Specificity",
                                                                                                                           title="STRINGv11" ,family = "Times", color="black")+
  theme_bw()+
  theme(legend.justification=c(1,0), legend.position=c(0.98,0.015),legend.key = element_blank())+theme(panel.border = element_blank(),
                                                                                                       panel.grid.major = element_blank(),
                                                                                                       panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
  scale_color_manual(values=c('steelblue','hotpink2',"#99cc00","gold1") ,"    Parameter               AUC",labels=c(paste("REGml", as.character(round(REGm_info_P[dim(REGm_info_P)[1]-1,dim(REGm_info_P)[2]],3)),sep="              "),
                                                                                                                    paste("REGmlw", as.character(round(REGml_info_P[dim(REGml_info_P)[1]-1,dim(REGml_info_P)[2]],3)),sep="           "),
                                                                                                                    paste("REGmlwp", as.character(round(REGmlw_info_P[dim(REGmlw_info_P)[1]-1,dim(REGmlw_info_P)[2]],3)),sep="        "),
                                                                                                                    paste("REGmlwpFC", as.character(round(REGmlwp_info_P[dim(REGmlwp_info_P)[1]-1,dim(REGmlwp_info_P)[2]],3)),sep="    ")))+
  xlim(0,1)+ylim(0,1)+ 
  theme(axis.title.x =element_text( size=10, color="black",family = "Times"),
        axis.title.y=element_text( size=10, color="black",family = "Times"),
        plot.title=element_text( size=12,family = "Times", color="black",hjust = +0.5),
        axis.line = element_line(colour = 'black', size = 0.3))+theme(legend.text = element_text(colour="black", size = 10,family = "Times"),
                                                                      legend.title = element_text( size=10,family = "Times", color="black"),
                                                                      axis.text = element_text(family = "Times",colour="black"),
                                                                      panel.border = element_rect(colour = "black", fill=NA, size=0.5))

STRa1


# pdf("1--DifferentStrategies/Plots/STRING_StrategiesAUC.pdf")
# STRa1
# dev.off()
# png("1--DifferentStrategies/Plots/STRING_StrategiesAUC.png")
# STRa1
# dev.off()

