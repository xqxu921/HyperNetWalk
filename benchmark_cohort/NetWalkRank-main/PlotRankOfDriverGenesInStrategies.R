setwd("/home/fatemeh/Fatemeh/")

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
# REGm_D$RankWPR<-((REGm_D$RankWPR-min(REGm_D$RankWPR))/(max(REGm_D$RankWPR)-min(REGm_D$RankWPR)))*1000
# REGml_D$RankWPR<-((REGml_D$RankWPR-min(REGml_D$RankWPR))/(max(REGml_D$RankWPR)-min(REGml_D$RankWPR)))*1000
# REGmlw_D$RankWPR<-((REGmlw_D$RankWPR-min(REGmlw_D$RankWPR))/(max(REGmlw_D$RankWPR)-min(REGmlw_D$RankWPR)))*1000
# REGmlwp_D$RankWPR<-((REGmlwp_D$RankWPR-min(REGmlwp_D$RankWPR))/(max(REGmlwp_D$RankWPR)-min(REGmlwp_D$RankWPR)))*1000
# REGmlwpFC_D$RankWPR<-((REGmlwpFC_D$RankWPR-min(REGmlwpFC_D$RankWPR))/(max(REGmlwpFC_D$RankWPR)-min(REGmlwpFC_D$RankWPR)))*1000

DriverRanks<- rbind(REG_D,REGm_D,REGml_D,REGmlw_D,REGmlwp_D,REGmlwpFC_D)
DriverRanks<-DriverRanks[,3:4]


D1<-ggplot(DriverRanks,aes(x=Type,y=RankWPR,fill=Type,family = "Times"))+
  geom_boxplot(position=position_dodge(1), alpha=1,width=0.9)+
  labs(y="Rank",x="",title="Rank of Driver genes in different strategies",family = "Times") +
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
  stat_summary(fun=mean, geom="point", shape=19, size=3, color="black", fill="black")+
  theme(plot.title = element_text(hjust = 0.5))
# D1<-D1+stat_summary(fun=mean, geom="point", shape=8, size=4, color="black", fill="black")+
#   stat_summary(fun=mean, geom="point", shape=19, size=4, color="black", fill="black")



D1


# pdf("0---PLOTS/StrategiesBoxPlot_DriversRanks_DEG.pdf")
#  D1
# dev.off()
# 
#  png("0---PLOTS/StrategiesBoxPlot_DriversRanks_DEG.png")
#  D1
#  dev.off()
