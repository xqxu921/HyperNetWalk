library(limma)
library(dplyr)
library(tidyr)
library(stringr)

exp_data_file <- "../../data/processed_data/BRCA/exp_tpm_data.tsv"
phenotype_file <- "../../data/rawdata/TCGA-BRCA.clinical.tsv"
exp_data <- read.delim(exp_data_file, header = TRUE, as.is = TRUE, check.names = FALSE)
invalid_gene_idx <- which(is.na(rownames(exp_data)) | rownames(exp_data) == "" | rownames(exp_data) == ".")
if (length(invalid_gene_idx) > 0){
  exp_data <- exp_data[-invalid_gene_idx,]
}
colnames(exp_data) <- trimws(colnames(exp_data))
# exp_data <- as.matrix(exp_data[,grepl("-01A$|-11A$",colnames(exp_data))])
exp_data <- as.matrix(exp_data[which(rowSums(exp_data >= 1) >= ncol(exp_data) * 0.2), ])

phenotype_data <- read.delim(phenotype_file,header = TRUE, as.is = TRUE, check.names = FALSE)
phenotype_data <- phenotype_data %>%
    # filter(grepl("-01A$|-11A$",sample)) %>%
    mutate(
        ajcc_pathologic_stage.diagnoses = ifelse(grepl("-11A$|-11B$",sample),"Normal",ajcc_pathologic_stage.diagnoses),
        ajcc_pathologic_stage.diagnoses = ifelse(ajcc_pathologic_stage.diagnoses %in% c("Stage I", "Stage IA", "Stage IB"),"Stage I",ajcc_pathologic_stage.diagnoses),
        ajcc_pathologic_stage.diagnoses = ifelse(ajcc_pathologic_stage.diagnoses %in% c("Stage II", "Stage IIA", "Stage IIB"),"Stage II",ajcc_pathologic_stage.diagnoses),
        ajcc_pathologic_stage.diagnoses = ifelse(ajcc_pathologic_stage.diagnoses %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"),"Stage III",ajcc_pathologic_stage.diagnoses)) %>%
    filter(!ajcc_pathologic_stage.diagnoses %in% c("",NA,"Stage X"))
com_samples <- intersect(colnames(exp_data),phenotype_data$sample)
exp_data <- exp_data[,com_samples]
phenotype_data <- phenotype_data[match(com_samples,phenotype_data$sample),]
gr <- phenotype_data$ajcc_pathologic_stage.diagnoses
gr <- factor(gr)
design <- model.matrix(~0+gr)
colnames(design) <- gsub(" ","_",levels(gr))
fit <- lmFit(exp_data, design)
contrast.matrix <- makeContrasts(
    Stage_I_vs_Normal = Stage_I - Normal,
    Stage_II_vs_Normal = Stage_II - Normal,
    Stage_III_vs_Normal = Stage_III - Normal,
    Stage_IV_vs_Normal = Stage_IV - Normal,
    levels = design
)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2,0.01)
tT_All <- topTable(fit2, adjust = "fdr", sort.by = "B", number = Inf)
subsetList <- c("Stage_I_vs_Normal","Stage_II_vs_Normal","Stage_III_vs_Normal","Stage_IV_vs_Normal","adj.P.Val")
tT_All1s <- subset(tT_All,select=subsetList)

aml.UpandDown_All <- subset(tT_All1s,adj.P.Val < 0.05)

tT_All_sN <- aml.UpandDown_All
dim(tT_All_sN)

LFC_Layers <- data.frame(tT_All_sN[,c(1:4)])

fileGGI <- "../../data/RegNet_human_v2.txt"
GGIRN <- read.delim(
    fileGGI,
    header = T,
    sep = "\t"
)
colnames(GGIRN) <- c("SourceNode","TargetNode")

GGIRN1 <- distinct(GGIRN)
GGIRN2 <- GGIRN1[GGIRN1$SourceNode %in% rownames(tT_All_sN),]
GGIRN3 <- GGIRN2[GGIRN2$TargetNode %in% rownames(tT_All_sN),]

Network1 <- GGIRN3
write.table(Network1,"./Data/Results/DEG_Layers/BRCA_Network1_GGI_common.txt",quote = F,sep = "\t",row.names = F)
GNNodes1 <- distinct(data.frame(Network1$SourceNode))
GNNodes2 <- distinct(data.frame(Network1$TargetNode))
colnames(GNNodes1)<-"Genes"
colnames(GNNodes2)<-"Genes"
GNNodes <- rbind(GNNodes1,GNNodes2)
NetNodes<- distinct(GNNodes)
dim(NetNodes)

C_df1<-tT_All_sN[rownames(tT_All_sN) %in% (NetNodes$Genes),]
geneListID<- C_df1
cID <-c(rownames(geneListID))

####Export Data to compute MI for each Layer####
#######Layers#######
####################
##Normal
Normal_EXp <- exp_data[cID,which(gr=="Normal")]
########  L1  #######
L1_Exp <- exp_data[cID,which(gr=="Stage I")]
write.table(L1_Exp, "./Data/Results/DEG_Layers/BRCA_L1_Exp.txt", quote=F,sep="\t",row.names=TRUE)
ex_L1_N<-cbind(Genes=rownames(Normal_EXp),L1_Exp,Normal_EXp)
write.table(ex_L1_N, "./Data/Results/DEG_Layers/BRCA_ex_L1_N.txt", quote=F,sep="\t",row.names=FALSE)
###################
########  L2  #######
L2_Exp <- exp_data[cID,which(gr=="Stage II")]
write.table(L2_Exp, "./Data/Results/DEG_Layers/BRCA_L2_Exp.txt", quote=F,sep="\t",row.names=TRUE)
ex_L2_N<-cbind(Genes=rownames(Normal_EXp),L2_Exp,Normal_EXp)
write.table(ex_L2_N, "./Data/Results/DEG_Layers/BRCA_ex_L2_N.txt", quote=F,sep="\t",row.names=FALSE)
###################
########  L3  #######
L3_Exp <- exp_data[cID,which(gr=="Stage III")]
write.table(L3_Exp, "./Data/Results/DEG_Layers/BRCA_L3_Exp.txt", quote=F,sep="\t",row.names=TRUE)
ex_L3_N<-cbind(Genes=rownames(Normal_EXp),L3_Exp,Normal_EXp)
write.table(ex_L3_N, "./Data/Results/DEG_Layers/BRCA_ex_L3_N.txt", quote=F,sep="\t",row.names=FALSE)
###################
########  L4  #######
L4_Exp <- exp_data[cID,which(gr=="Stage IV")]
write.table(L4_Exp, "./Data/Results/DEG_Layers/BRCA_L4_Exp.txt", quote=F,sep="\t",row.names=TRUE)
ex_L4_N<-cbind(Genes=rownames(Normal_EXp),L4_Exp,Normal_EXp)
write.table(ex_L4_N, "./Data/Results/DEG_Layers/BRCA_ex_L4_N.txt", quote=F,sep="\t",row.names=FALSE)
###################
