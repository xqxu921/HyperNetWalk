#install.packages("pROC")
library(pROC)
require(pROC)
require(ggpubr)
#install.packages("ggpubr")
require(ggplot2)
require(plyr)
require(gridExtra)
#install.packages("gridExtra")

s_REG <- read.delim("9--SEcond HCC/REG/Re_Sum_C.txt")
s_REGm <- read.delim("9--SEcond HCC/REGm/Re_Sum_C.txt")
s_REGml <- read.delim("9--SEcond HCC/REGml/Re_Sum_C.txt")
s_REGmlw <- read.delim("9--SEcond HCC/REGmlw/Re_Sum_C.txt")
s_REGmlwp <- read.delim("9--SEcond HCC/REGmlwp/Re_Sum_C.txt")
s_REGmlwpFC <- read.delim("9--SEcond HCC/REGmlwpFC/Re_Sum_C.txt")

s_REG<-s_REG[!(is.na(s_REG$F1)),]
s_REGm<-s_REGm[!(is.na(s_REGm$F1)),]
s_REGml<-s_REGml[!(is.na(s_REGml$F1)),]

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

r_s_REG <- my_fun1(s_REG)
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

REG_info_P<-re_all(r_s_REG,info)
write.table(REG_info_P,"9--SEcond HCC/secound3/REG_sec_paramaters.txt",quote = F,sep = " ")

REGm_info_P<-re_all(r_s_REGm,info)
write.table(REGm_info_P,"9--SEcond HCC/secound3/REGm_sec_paramaters.txt",quote = F,sep = " ")

REGml_info_P<-re_all(r_s_REGml,info)
write.table(REGml_info_P,"9--SEcond HCC/secound3/REGml_sec_paramaters.txt",quote = F,sep = " ")

REGmlw_info_P<-re_all(r_s_REGmlw,info)
write.table(REGmlw_info_P,"9--SEcond HCC/secound3/REGmlw_sec_paramaters.txt",quote = F,sep = " ")

REGmlwp_info_P<-re_all(r_s_REGmlwp,info)
write.table(REGmlwp_info_P,"9--SEcond HCC/secound3/REGmlwp_sec_paramaters.txt",quote = F,sep = " ")

REGmlwpFC_info_P<-re_all(r_s_REGmlwpFC,info)
write.table(REGmlwpFC_info_P,"9--SEcond HCC/secound3/REGmlwpFC_sec_paramaters.txt",quote = F,sep = " ")




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

#DEG
NormalDPRi<-read.delim("9--SEcond HCC/REG/NormalDPR.txt") 
dataD_i<-read.delim("9--SEcond HCC/REG/dataD.txt")

Re_DEG<-REfunction(NormalDPRi,dataD_i,s_REG,"REG","binormal")
gg_REG1<- Re_DEG[[1]]
zz_REG<- Re_DEG[[2]]

ggplot(data=gg_REG1,aes(x=x1,y=y1,group = type))+geom_line(aes(x=x1,y=y1,color = type))

#DEGm
NormalDPRim<-read.delim("9--SEcond HCC/REGm/NormalDPR.txt") 
dataD_im<-read.delim("9--SEcond HCC/REGm/dataD.txt")

Re_DEGm<-REfunction(NormalDPRim,dataD_im,s_REGm,"REGm","binormal")
gg_REGm1<- Re_DEGm[[1]]
zz_REGm<- Re_DEGm[[2]]

ggplot(data=gg_REGm1,aes(x=x1,y=y1,group = type))+geom_line(aes(x=x1,y=y1,color = type))

#DEGml
NormalDPRil<-read.delim("9--SEcond HCC/REGml/NormalDPR.txt") 
dataD_il<-read.delim("9--SEcond HCC//REGml/dataD.txt")

Re_DEGml<-REfunction(NormalDPRil,dataD_il,s_REGml,"REGml","binormal")
gg_REGml1<- Re_DEGml[[1]]
zz_REGml<- Re_DEGml[[2]]

ggplot(data=gg_REGml1,aes(x=x1,y=y1,group = type))+geom_line(aes(x=x1,y=y1,color = type))

#DEGmlw
NormalDPRimlw<-read.delim("9--SEcond HCC/REGmlw/NormalDPR.txt") 
dataD_imlw<-read.delim("9--SEcond HCC/REGmlw/dataD.txt")

Re_DEGmlw<-REfunction(NormalDPRimlw,dataD_imlw,s_REGmlw,"REGmlw","binormal")
gg_REGmlw1<- Re_DEGmlw[[1]]
zz_REGmlw<- Re_DEGmlw[[2]]

ggplot(data=gg_REGmlw1,aes(x=x1,y=y1,group = type))+geom_line(aes(x=x1,y=y1,color = type))

#DEGmlwp
NormalDPRimlwp<-read.delim("9--SEcond HCC/REGmlwp/NormalDPR.txt") 
dataD_imlwp<-read.delim("9--SEcond HCC/REGmlwp/dataD.txt")

Re_DEGmlwp<-REfunction(NormalDPRimlwp,dataD_imlwp,s_REGmlwp,"REGmlwp","binormal")
gg_REGmlwp1<- Re_DEGmlwp[[1]]
zz_REGmlwp<- Re_DEGmlwp[[2]]

ggplot(data=gg_REGmlwp1,aes(x=x1,y=y1,group = type))+geom_line(aes(x=x1,y=y1,color = type))

#DEGmlwpFC
NormalDPRimlwpFC<-read.delim("9--SEcond HCC/REGmlwpFC/NormalDPR.txt") 
dataD_imlwpFC<-read.delim("9--SEcond HCC/REGmlwpFC/dataD.txt")

Re_DEGmlwpFC<-REfunction(NormalDPRimlwpFC,dataD_imlwpFC,s_REGmlwpFC,"REGmlwpFC","density")
gg_REGmlwpFC1<- Re_DEGmlwpFC[[1]]
zz_REGmlwpFC<- Re_DEGmlwpFC[[2]]

ggplot(data=gg_REGmlwpFC1,aes(x=x1,y=y1,group = type))+geom_line(aes(x=x1,y=y1,color = type))


#ggtest<-rbind(gg_REGm1,gg_REGml1,gg_REGmlw1,gg_REGmlwp1,gg_REGmlwpFC1)
ggtest2<-rbind(gg_REG1,gg_REGm1,gg_REGml1,gg_REGmlw1,gg_REGmlwp1,gg_REGmlwpFC1)



# write.table(ggtest,"1--DifferentStrategies/Plots/AucInfo.txt",quote = F,sep = " ")
# ggtest<-read.table("1--DifferentStrategies/Plots/AucInfo.txt")



a_sec <- ggplot(data=ggtest2,aes(x=x1,y=y1,group = type))+geom_line(aes(x=x1,y=y1,color = type),size = 1,linetype = 1)+labs(y="Sensitivity",x="1-Specificity",
                                                                                                                            title="" ,family = "Times", color="black")+
  theme_bw()+
  theme(legend.justification=c(1,0), legend.position=c(0.98,0.015),legend.key = element_blank())+theme(panel.border = element_blank(),
                                                                                                       panel.grid.major = element_blank(),
                                                                                                       panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
  scale_color_brewer(palette="Set2" ,"    Parameter             AUC",labels=c(paste("REG", as.character(round(REG_info_P[dim(REG_info_P)[1]-1,dim(REG_info_P)[2]],3)),sep="                "),
                                                                              paste("REGm", as.character(round(REGm_info_P[dim(REGm_info_P)[1]-1,dim(REGm_info_P)[2]],3)),sep="             "),
                                                                              paste("REGml", as.character(round(REGml_info_P[dim(REGml_info_P)[1]-1,dim(REGml_info_P)[2]],3)),sep="             "),
                                                                              paste("REGmlw", as.character(round(REGmlw_info_P[dim(REGmlw_info_P)[1]-1,dim(REGmlw_info_P)[2]],3)),sep="          "),
                                                                              paste("REGmlwp", as.character(round(REGmlwp_info_P[dim(REGmlwp_info_P)[1]-1,dim(REGmlwp_info_P)[2]],3)),sep="        "),
                                                                              paste("REGmlwpFC", as.character(round(REGmlwpFC_info_P[dim(REGmlwpFC_info_P)[1]-1,dim(REGmlwpFC_info_P)[2]],3)),sep="    ")))+
  xlim(0,1)+ylim(0,1)+ 
  theme(axis.title.x =element_text( size=10, color="black",family = "Times"),
        axis.title.y=element_text( size=10, color="black",family = "Times"),
        plot.title=element_text( size=12,family = "Times", color="black",hjust = +0.5),
        axis.line = element_line(colour = 'black', size = 0.3))+theme(legend.text = element_text(colour="black", size = 10,family = "Times"),
                                                                      legend.title = element_text( size=10,family = "Times", color="black"),
                                                                      axis.text = element_text(family = "Times",colour="black"),
                                                                      panel.border = element_rect(colour = "black", fill=NA, size=0.5))

a_sec
# pdf("1--DifferentStrategies/Plots/StrategiesAUC_SEC.pdf")
# a_sec
# dev.off()
# png("1--DifferentStrategies/Plots/StrategiesAUC_SEC.png")
# a_sec
# dev.off()

##########AUC 


# zz_REG<-read.delim("1--DifferentStrategies/Plots/REG_AucInfo.txt")

zztest<-rbind(zz_REG,zz_REGm,zz_REGml,zz_REGmlw,zz_REGmlwp,zz_REGmlwpFC)
# zztest<-rbind(zz_REGm,zz_REGml,zz_REGmlw,zz_REGmlwp,zz_REGmlwpFC)

#write.table(zztest,"1--DifferentStrategies/Plots/AucInfoBoxPlot.txt",quote = F,sep = " ")
#zztest<-read.delim("1--DifferentStrategies/Plots/AucInfoBoxPlot.txt")
##zztest$lab<-factor(zztest$lab,level=c("r=0.6","r=0.7","r=0.8","r=0.9"))


#colnames(zz)<-c("auc","l")




p_sec<-ggplot(zztest,aes(x=lab,y=AUC,fill=lab))+
  geom_boxplot(position=position_dodge(0.9), alpha=1)+
  labs(y="AUC",x="",title="") +
  scale_fill_brewer(palette="Set2",)+theme_bw()

p_sec<-p_sec+guides(fill = FALSE)
p_sec<-p_sec+theme(axis.title.x =element_text( size=10, color="black",family = "Times"),
                   axis.title.y=element_text( size=10, color="black",family = "Times"),
                   plot.title=element_text( color="black",family = "Times"),
                   axis.line = element_line(colour = 'black', size = 0.3))+theme(legend.text = element_text(colour="black", size = 10,family = "Times"),
                                                                                 legend.title = element_text( size=10,family = "Times", color="black"),
                                                                                 axis.text = element_text(family = "Times",colour="black"),
                                                                                 axis.text.x = element_text(angle = 30,hjust =1.2 ),
                                                                                 panel.border = element_rect(colour = "black", fill=NA, size=0.5))
# compare_means(AUC ~ lab,  data = zztest,
#               method = "t.test")

p_sec<-p_sec+stat_compare_means(label.y = 0.97,family="Times")

# p_sec<-p_sec+stat_compare_means(comparisons = list(c("REGm", "REGml"),c("REGmlwp", "REGmlwpFC"),c("REG","REGmlwpFC")),
# label = "p.signif", hide.ns = TRUE)


p_sec

# p_sec<-p_sec+stat_compare_means(comparisons = list(c("REG", "REGm"),c("REG","REGmlwpFC")),
#                     label = "p.signif", hide.ns = TRUE)

pdf("0---PLOTS/StrategiesBoxPlot_AUC_SEC_temp.pdf")
p_sec
dev.off()
# png("1--DifferentStrategies/Plots/StrategiesBoxPlot_AUC_SEC.png")
# p_sec
# dev.off()


lb6<-paste("A","B", sep = paste(rep(" ", 75), collapse = ""))
lb7<-paste(lb6," ", sep = paste(rep(" ", 60), collapse = ""))
# grid.arrange(arrangeGrob(top=lb7,p1, nrow = 1))

pdf("0---PLOTS/Sec_HCC_DifferentStrategies_grid1.pdf",width = 8,height =5)
grid.arrange(arrangeGrob(top=lb7,a_sec,p_sec, nrow = 1))
dev.off()
