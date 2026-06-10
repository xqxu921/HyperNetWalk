rm(list=ls())
### set your work path

setwd("/home/zcx/Fatemeh/")

### Input your network data
file1 <- "DifferentStrategies/DEG_Layers/Network1_GGI_common.txt"  #your network data
Net <- as.matrix(read.table(file1,sep="\t",header=T))
Net_ID <- Net

### Input your geneExpression data
file2 <- "DifferentStrategies/DEG_Layers/ex_L9_N.txt"
GEData<-read.csv(file2,sep="\t",header =F)
GEData<-GEData[-1,]
Gene_ID <- GEData[,1]
#dim(RNA_seq)
### split data into Disease and Normal

used_GEData <- GEData

library("infotheo")
nc_used<-ncol(used_GEData)
T_data <- used_GEData[,2:(nc_used-13)]
N_data <- used_GEData[,(nc_used-12):nc_used]
dim(T_data)
dim(N_data)
### calculate DMI and store DMI>0
#require(pracma)
#install.packages("pracma")
library(pracma)
selected_net <- rbind()
for(i in 1:dim(Net_ID)[1])
{
  k1 <- which( as.numeric(Net_ID[i,1]) == as.numeric(used_GEData[,1]) )
  k2 <- which( as.numeric(Net_ID[i,2]) == as.numeric(used_GEData[,1]) )
  
  if(isempty(k1)==F)
  {
    if(isempty(k2)==F)
    {
      ### discretize data 
      T1 <- apply(as.matrix(t(T_data[c(k1,k2),])),2,function(x) as.numeric(x))
      N1 <- apply(as.matrix(t(N_data[c(k1,k2),])),2,function(x) as.numeric(x))
      T_data_cur <- discretize(T1,nbins=25)
      N_data_cur <- discretize(N1,nbins=25)
      ### DMI calculation
      DMI_cur <- abs(mutinformation(T_data_cur[,1],T_data_cur[,2],method="emp")-mutinformation(N_data_cur[,1],N_data_cur[,2],method="emp"))
      
      ### store these DMI
      store_net <- c(Net_ID[i,],DMI_cur)
      #store_net <- c(Net_ID[i,],MI_cur)
      selected_net <- rbind(selected_net,store_net)
      print(store_net)
    }
  }
}

dim(selected_net)

##select MI values greater than zero
selected_net1<- as.data.frame(selected_net)
colnames(selected_net1) <- c("ID1","ID2","MI")
selected_net2<- selected_net1[!(selected_net1$MI==0),]
selected_net3<- as.matrix(selected_net2)
dim(selected_net3)
write.table(selected_net3,"DifferentStrategies/Results_MI/MI_9N_LA.txt",sep="\t",row.names = F,quote = F)


##Normalize Data
max_MI <- max(as.numeric(selected_net3[,3]))
min_MI <- min(as.numeric(selected_net3[,3]))

Nor_MI1 <- apply(as.matrix(selected_net2[,3]),2,function(x) (as.numeric(x)-min_MI)/(max_MI-min_MI))

Final_net_store <- cbind(selected_net3[,-3],Nor_MI1)
colnames(Final_net_store) <- c("ID1","ID2","Normalized_MI")

### output results
write.table(Final_net_store,"DifferentStrategies/Results_MI/MI_9N_LA_Normalized.txt",sep="\t",row.names = F,quote = F)
###


