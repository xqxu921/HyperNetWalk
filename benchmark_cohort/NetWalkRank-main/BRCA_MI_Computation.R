rm(list=ls())
library("infotheo")
library(pracma)
### set your work path

### Input your network data
file1 <- "./Data/Results/DEG_Layers/BRCA_Network1_GGI_common.txt"  #your network data
Net <- as.matrix(read.table(file1,sep="\t",header=T))
Net_ID <- Net
file2s <- list(
    "./Data/Results/DEG_Layers/BRCA_ex_L1_N.txt",
    "./Data/Results/DEG_Layers/BRCA_ex_L2_N.txt",
    "./Data/Results/DEG_Layers/BRCA_ex_L3_N.txt",
    "./Data/Results/DEG_Layers/BRCA_ex_L4_N.txt"
)
for(l in 1:length(file2s)){
    file2 <- file2s[[l]]
    ### Input your geneExpression data
    # file2 <- "./Data/Results/DEG_Layers/BRCA_ex_L1_N.txt"
    GEData <- read.csv(file2,sep = "\t",header = T,row.names=1)
    Gene_ID <- rownames(GEData)
    #dim(RNA_seq)
    ### split data into Disease and Normal

    used_GEData <- GEData
    nc_used<-ncol(used_GEData)
    normal_idx <- substr(colnames(used_GEData),14,16) %in% c("11A","11B")
    T_data <- used_GEData[,!normal_idx]
    N_data <- used_GEData[,normal_idx]
    dim(T_data)
    dim(N_data)
    ### calculate DMI and store DMI>0
    selected_net <- rbind()
    for(i in 1:dim(Net_ID)[1])
    {
    k1 <- which(Net_ID[i,1] == rownames(used_GEData))
    k2 <- which(Net_ID[i,2] == rownames(used_GEData))
    
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
    selected_net1 <- as.data.frame(selected_net)
    colnames(selected_net1) <- c("ID1","ID2","MI")
    selected_net2 <- selected_net1[!(selected_net1$MI==0),]
    selected_net3 <- as.matrix(selected_net2)
    dim(selected_net3)
    write.table(selected_net3,paste0("./Data/Results/Results_MI/BRCA_MI_",l,"N_LA.txt"),sep="\t",row.names = F,quote = F)

    ##Normalize Data
    max_MI <- max(as.numeric(selected_net3[,3]))
    min_MI <- min(as.numeric(selected_net3[,3]))

    Nor_MI1 <- apply(as.matrix(selected_net2[,3]),2,function(x) (as.numeric(x)-min_MI)/(max_MI-min_MI))

    Final_net_store <- cbind(selected_net3[,-3],Nor_MI1)
    colnames(Final_net_store) <- c("ID1","ID2","Normalized_MI")

    ### output results
    write.table(Final_net_store,paste0("./Data/Results/Results_MI/BRCA_MI_",l,"N_LA_Normalized.txt"),sep="\t",row.names = F,quote = F)
    ###
}
