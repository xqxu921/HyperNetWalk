
setwd("/home/zcx/Fatemeh/")

library(readr) # for reading your sample data
library(ggplot2)
library(ggthemes) # for theme_few()



###############################################################################
REG_D <- read.delim("1--DifferentStrategies/REG/Re1/dataD.txt")
REGm_D <- read.delim("1--DifferentStrategies/REGm/r70/dataD.txt")
REGml_D <- read.delim("1--DifferentStrategies/REGml/r70/dataD.txt")
REGmlw_D <- read.delim("1--DifferentStrategies/REGmlw/r70/dataD.txt")
REGmlwp_D <- read.delim("1--DifferentStrategies/REGmlwp/r70/dataD2.txt")
REGmlwpFC_D <- read.delim("1--DifferentStrategies/REGmlwpFC/r70/dataD3.txt")

REG_D$Type<- "REG"
REGm_D$Type<- "REGm"
REGml_D$Type<- "REGml"
REGmlw_D$Type<- "REGmlw"
REGmlwp_D$Type<- "REGmlwp"
REGmlwpFC_D$Type<- "REGmlwpFC"

colnames(REG_D)<-colnames(REGm_D)

DriverRanks<- rbind(REG_D,REGm_D,REGml_D,REGmlw_D,REGmlwp_D,REGmlwpFC_D)
DriverRanks<-DriverRanks[,3:4]


D1<-ggplot(DriverRanks,aes(x=Type,y=RankWPR,fill=Type,family = "Times"))+
  geom_boxplot(position=position_dodge(1), alpha=1,width=0.9)+
  labs(y="Rank",x="",title="",family = "Times") +
  scale_fill_brewer(palette="Set2")+theme_bw()+guides(fill = FALSE)+
  theme(axis.title.x =element_text( size=10, color="black",family = "Times"),
        axis.title.y=element_text( size=10, color="black",family = "Times"),
        plot.title=element_text( color="black",family = "Times"),
        axis.line = element_line(colour = 'black', size = 0.3))+
  theme(legend.text = element_text(colour="black", size = 10,family = "Times"),
        legend.title = element_text( size=10,family = "Times", color="black"),
        axis.text = element_text(family = "Times",colour="black"),
        axis.text.x = element_text(angle = 30,hjust =1 ),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
  stat_summary(fun=mean, geom="point", shape=19, size=2, color="black", fill="black")+
  theme(plot.title = element_text(hjust = 0.5))
# D1<-D1+stat_summary(fun=mean, geom="point", shape=8, size=4, color="black", fill="black")+
#   stat_summary(fun=mean, geom="point", shape=19, size=4, color="black", fill="black")



D1


# pdf("1--DifferentStrategies/BIGPic/StrategiesBoxPlot_DriversRanks_DEG.pdf")
#  D1
# dev.off()
# 
#  png("0---PLOTS/StrategiesBoxPlot_DriversRanks_DEG.png")
#  D1
#  dev.off()
###############################################################


REGmlwpFC_D3 <- read.delim("2--DifferentLayers/456789/dataD.txt")
REGmlwpFC_D6 <- read.delim("2--DifferentLayers/789/dataD.txt")
REGmlwpFC_D9 <- read.delim("1--DifferentStrategies/REGmlwpFC/r70/dataD3.txt")


REGmlwpFC_D3$Type<- "+3Layer"
REGmlwpFC_D6$Type<- "+6Layer"
REGmlwpFC_D9$Type<- "+9Layer"


DriverRanks<- rbind(REGmlwpFC_D3,REGmlwpFC_D6,REGmlwpFC_D9)
DriverRanks<-DriverRanks[,3:4]


D369<-ggplot(DriverRanks,aes(x=Type,y=RankWPR,fill=Type,family = "Times"))+
  geom_boxplot(position=position_dodge(1), alpha=1,width=0.9)+
  labs(y="Rank",x="",title="",family = "Times") +
  scale_fill_brewer(palette="Set3")+theme_bw()+guides(fill = FALSE)+
  theme(axis.title.x =element_text( size=10, color="black",family = "Times"),
        axis.title.y=element_text( size=10, color="black",family = "Times"),
        plot.title=element_text( color="black",family = "Times"),
        axis.line = element_line(colour = 'black', size = 0.3))+
  theme(legend.text = element_text(colour="black", size = 10,family = "Times"),
        legend.title = element_text( size=10,family = "Times", color="black"),
        axis.text = element_text(family = "Times",colour="black"),
        axis.text.x = element_text(angle = 30,hjust =1 ),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
  stat_summary(fun=mean, geom="point", shape=19, size=2, color="black", fill="black")+
  theme(plot.title = element_text(hjust = 0.5))
# D1<-D1+stat_summary(fun=mean, geom="point", shape=8, size=4, color="black", fill="black")+
#   stat_summary(fun=mean, geom="point", shape=19, size=4, color="black", fill="black")



D369


# 
# pdf("0---PLOTS/Strategies369_DriversRanks_DEG.pdf")
# D369
# dev.off()
# 
# png("0---PLOTS/Strategies369_DriversRanks_DEG.png")
# D369
# dev.off()


################################################################
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

########
gmeltDN<-melt(Ranks_Driver[,c(-1,-11,-12,-13,-14)])
DiLength<-dim(Ranks_Driver[Ranks_Driver$DrNo=="D",])
DiLength


###############PLOT LAYERS
plotData<-Ranks_Driver[,c(-1,-11,-12,-13,-14)]
plotData[plotData$DrNo=="N",10]<-"Non-Driver Genes"
plotData[plotData$DrNo=="D",10]<-"Driver Genes"
colnames(plotData)<-c("Choronic Hepatitis with low grade","Choronic Hepatitis with high grade",
                      "Cirrhosis","Dysplastic nodules with low grade",
                      "Dysplastic nodules with high grade","Early hepatocellular carcinoma",
                      "hepatocellular carcinoma (TG1)","hepatocellular carcinoma (TG2)",
                      "hepatocellular carcinoma (TG3)","DrNo")

gmeltDN1<-melt(plotData)
p_A<-ggplot(gmeltDN1,aes(variable,value,fill=DrNo))+geom_boxplot(position = "dodge")+
  labs(y="Rank",x="",title="")+
  theme_bw()+scale_fill_brewer(palette="Accent",)+
  theme(axis.title.x =element_text( size=10, color="black",family = "Times"),
        axis.title.y=element_text( size=10, color="black",family = "Times"),
        plot.title=element_text( color="black"),
        axis.line = element_line(colour = 'black', size = 0.3))+
  theme(legend.text = element_text(colour="black", size = 10,family = "Times"),
        legend.title = element_blank(),
        axis.text = element_text(family = "Times",colour="black"))+
  theme(legend.position="bottom")+
  scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 2), 
                   labels = function(x) stringr::str_wrap(x, width = 20))


p_A<-ggplot(gmeltDN1,aes(variable,value,fill=DrNo))+geom_boxplot(position = "dodge")+
  labs(y="Rank",x="",title="")+
  theme_bw()+scale_fill_brewer(palette="Accent",)+
  theme(axis.title.x =element_text( size=10, color="black",family = "Times"),
        axis.title.y=element_text( size=10, color="black",family = "Times"),
        plot.title=element_text( color="black"),
        axis.line = element_line(colour = 'black', size = 0.3))+
  theme(legend.text = element_text(colour="black", size = 10,family = "Times"),
        legend.title = element_blank(),
        axis.text = element_text(family = "Times",colour="black",angle = 70,hjust = 1))+
  theme(legend.position="bottom")

p_A
# pdf("1--DifferentStrategies/BIGPic/Layers.pdf")
# p_A
# dev.off()

# png("1--DifferentStrategies/BIGPic/Layers.png")
# p_A
# dev.off()
########################
tPlot<-Ranks_Driver[,c(12,15)]
tPlot[tPlot$DrNo=="D",2]<-"Driver Genes"
tPlot[tPlot$DrNo=="N",2]<-"Non-Driver Genes"
g1<-ggplot(tPlot,aes(x=DrNo,y=RankWPR,fill=DrNo))+
  geom_boxplot(position=position_dodge(0.9), alpha=1)+
  labs(y="Rank",x="",title="") +
  scale_fill_brewer(palette="Accent",)+theme_bw()

g1<-g1+guides(fill = FALSE)
g1<-g1+theme(axis.title.x =element_text( size=10, color="black",family = "Times"),
             axis.title.y=element_text( size=10, color="black",family = "Times"),
             plot.title=element_text( color="black"),
             axis.line = element_line(colour = 'black', size = 0.3))+
  theme(legend.text = element_text(colour="black", size = 10,family = "Times"),
        legend.title = element_text( size=10,family = "Times", color="black"),
        axis.text = element_text(family = "Times",colour="black"),
        axis.text.x = element_text(angle = 30,hjust =1 ))
g1
# pdf("1--DifferentStrategies/BIGPic/RankofDriverGenes.pdf")
# g1
# dev.off()


