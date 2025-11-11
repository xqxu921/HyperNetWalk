# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(data.table)
library(dplyr)
library(WGCNA)
library(Matrix)
source("src/model.R")
PPI_file <- "./data/NETWORK/9606.protein.links.v12.0.onlyAB.tsv"
protein_file <- "./data/metadata/9606.protein.info.v12.0.txt"

PPI_df <- fread(PPI_file)
protein_df <- fread(protein_file,
                    select = c("#string_protein_id", "preferred_name"))
setnames(protein_df, "#string_protein_id", "protein_id")
protein_map1 <- match(PPI_df$protein1, protein_df$protein_id)
protein_map2 <- match(PPI_df$protein2, protein_df$protein_id)
PPI_df$protein1 <- protein_df$preferred_name[protein_map1]
PPI_df$protein2 <- protein_df$preferred_name[protein_map2]
PPI_df$combined_score <- PPI_df$combined_score / 1000
fwrite(PPI_df,
       file = "./data/NETWORK/STRINGv12.txt",
       quote = FALSE,
       sep = "\t")

GRN_file <- "./data/NETWORK/human_TF_Target.txt"
GRN_df <- fread(GRN_file, header = F)
GRN_df <- GRN_df[, c(1, 3)]
if (any(duplicated(GRN_df))) {
  GRN_df <- unique(GRN_df)
}
colnames(GRN_df) <- c("TF", "Target")
fwrite(GRN_df,
       file = "./data/NETWORK/RegNet_human_V2.txt",
       quote = FALSE,
       sep = "\t")

# GC_file <- "./data/human_core_TF_Target.txt"
# GC_df <- fread(GC_file, header = F)
# GC_df <- GC_df[, c(1, 3)]
# if (any(duplicated(GC_df))) {
#   GC_df <- unique(GC_df)
# }
# colnames(GC_df) <- c("TF", "Target")
# fwrite(GC_df,
#        file = "./data/RegNet_human_core_V2.txt",
#        quote = FALSE,
#        sep = "\t")

cancers <- c(
  "PANCAN",
  "BRCA",
  "COAD",
  "HNSC",
  "KIRC",
  "KIRP",
  "LIHC",
  "LUAD",
  "LUSC",
  "PRAD",
  "STAD",
  "THCA",
  "UCEC"
)
for (cancer_type in cancers) {
  message("Start network construction for ",cancer_type)
  mut_data_file <- file.path("./data/processed", cancer_type, "mut_data.tsv")
  mut_data <- get_mut_data(mut_data_file)
  if (cancer_type == "PANCAN"){
    exp_data_file <- file.path("./data/processed", cancer_type, "exp_data.tsv")
  } else {
    exp_data_file <- file.path("./data/processed", cancer_type, "exp_tpm_data.tsv")
  }
  
  # exp_data <- get_exp_data(exp_data_file)
  # 
  # exp_genes <- rownames(exp_data)
  # exp_cor <- cor(t(exp_data), method = "pearson")
  # uptri_idx <- which(upper.tri(exp_cor), arr.ind = TRUE)
  # exp_cor_df <- data.frame(gene1 = exp_genes[uptri_idx[, 1]],
  #                          gene2 = exp_genes[uptri_idx[, 2]],
  #                          cor = exp_cor[uptri_idx])
  # CEN_gz_file <- gzfile(file.path("./data/processed_data", cancer_type, "CEN.txt.gz"),
  #                       "w")
  # write.table(
  #   exp_cor_df,
  #   file = CEN_gz_file,
  #   sep = "\t",
  #   quote = FALSE,
  #   row.names = FALSE
  # )
  # close(CEN_gz_file)
  # # adj_mat <- (abs(exp_cor) >= 0.7) + 0
  # # diag(adj_mat) <- 0
  # # save(adj_mat, file = "./data/processed_data/BRCA/co_exp_net.Rdata")
  # exp_data_scaled <- scale(t(exp_data), center = TRUE, scale = TRUE)
  # exp_data_scaled <- abs(exp_data_scaled)
  # exp_diff_z_cor <- cor(exp_data_scaled, method = "pearson")
  # exp_diff_z_cor_df <- data.frame(gene1 = exp_genes[uptri_idx[, 1]],
  #                                 gene2 = exp_genes[uptri_idx[, 2]],
  #                                 cor = exp_diff_z_cor[uptri_idx])
  # exp_diff_z_cor_df <- exp_diff_z_cor_df[exp_diff_z_cor_df$cor > 0, ]
  # DEN_z_gz_file <- gzfile(file.path("./data/processed_data", cancer_type, "DEN_z.txt.gz"),
  #                         "w")
  # write.table(
  #   exp_diff_z_cor_df,
  #   file = DEN_z_gz_file,
  #   sep = "\t",
  #   quote = FALSE,
  #   row.names = FALSE
  # )
  # close(DEN_z_gz_file)
  # adj_diff_z <- (exp_diff_z_cor >= 0.5) + 0
  # save(adj_diff_z, file = "./data/processed_data/BRCA/DEN_z.Rdata")
  exp_data <- read.delim(
    exp_data_file,
    header = TRUE,
    as.is = TRUE,
    check.names = FALSE
  )
  invalid_gene_index <- which(is.na(rownames(exp_data)) | rownames(exp_data) == "" | rownames(exp_data) == ".")
  if (length(invalid_gene_index) > 0){
    exp_data <- exp_data[-invalid_gene_index,]
  }
  if (cancer_type == "PANCAN"){
    exp_data <- exp_data[,grepl("-01$|-11$", colnames(exp_data))]
  } else {
    exp_data <- exp_data[,grepl("-01A$|-11A$", colnames(exp_data))]
    colnames(exp_data) <- substr(colnames(exp_data),1,15)
  }
  exp_data <- as.matrix(exp_data[which(rowSums(exp_data >= 1) >= ncol(exp_data) * 0.2 | rownames(exp_data) %in% rownames(mut_data)),])
  exp_data <- exp_data[rowSums(exp_data)!=0,]
  exp_data_tumor <- exp_data[, grepl("-01$", colnames(exp_data))]
  colnames(exp_data_tumor) <- substr(colnames(exp_data_tumor), 1, 12)
  exp_data_tumor <- as.matrix(exp_data_tumor)
  exp_cor_tumor <- WGCNA::cor(t(exp_data_tumor),method="pearson",nThreads = 32)
  exp_cor_tumor[is.na(exp_cor_tumor)] <- 0
  uptri_idx <- which(upper.tri(exp_cor_tumor), arr.ind = TRUE)
  
  if (sum(grepl("-11$", colnames(exp_data)))>1){
    exp_data_normal <- exp_data[, grepl("-11$", colnames(exp_data))]
    colnames(exp_data_normal) <- substr(colnames(exp_data_normal), 1, 12)
    exp_data_normal <- as.matrix(exp_data_normal)
    exp_cor_normal <- WGCNA::cor(t(exp_data_normal), method = "pearson",nThreads = 32)
    exp_cor_normal[is.na(exp_cor_normal)] <- 0
    # exp_cor_union_df <- data.frame(gene1 = rownames(exp_cor_tumor)[uptri_idx[, 1]],
    #                               gene2 = rownames(exp_cor_tumor)[uptri_idx[, 2]],
    #                               cor = pmax(abs(exp_cor_tumor[uptri_idx]),abs(exp_cor_normal[uptri_idx])))
    exp_cor_union <- pmax(abs(exp_cor_tumor),abs(exp_cor_normal))
  } else {
    # exp_cor_union_df <- data.frame(gene1 = rownames(exp_cor_tumor)[uptri_idx[, 1]],
    #                               gene2 = rownames(exp_cor_tumor)[uptri_idx[, 2]],
    #                               cor = abs(exp_cor_tumor[uptri_idx]))
    exp_cor_union <- abs(exp_cor_tumor) 
    next
  }
  # exp_cor_union_df <- exp_cor_union_df[exp_cor_union_df$cor != 0, ]
  # CEN_union_gz_file <- file.path("./data/processed_data",cancer_type,"CEN_union.txt.gz")
  # fwrite(
  #   exp_cor_union_df,
  #   file = CEN_union_gz_file,
  #   sep = "\t",
  #   quote = FALSE,
  #   row.names = FALSE,
  #   col.names = TRUE,
  #   compress = "gzip"
  # )
  exp_cor_upper <- matrix(0,nrow = nrow(exp_cor_tumor),ncol = ncol(exp_cor_tumor),dimnames = list(rownames(exp_cor_tumor),colnames(exp_cor_tumor)))
  exp_cor_upper[uptri_idx] <- exp_cor_union[uptri_idx]
  exp_cor_sparse <- Matrix(exp_cor_upper,sparse = TRUE)
  saveRDS(exp_cor_sparse, file = file.path("./data/processed",cancer_type,"CEN_union_matrix_uptri.rds"))
  
  # com_samples <- intersect(colnames(exp_data_normal), colnames(exp_data_tumor))
  # exp_data_tumor <- exp_data_tumor[, com_samples]
  # exp_data_normal <- exp_data_normal[, com_samples]
  # # exp_data_tumor <- exp_data_tumor[exp_genes, com_samples]
  # # exp_data_normal <- exp_data_normal[exp_genes, com_samples]
  # exp_diff <- exp_data_tumor - exp_data_normal
  # exp_diff_scaled <- scale(t(exp_diff), scale = TRUE, center = TRUE)
  # exp_diff_scaled <- abs(exp_diff_scaled)
  # exp_diff_scaled[is.na(exp_diff_scaled)] <- 0
  # exp_diff_cor <- WGCNA::cor(exp_diff_scaled, method = "pearson", nThreads = 32)
  # exp_diff_cor[is.na(exp_diff_cor)] <- 0
  # exp_diff_cor[exp_diff_cor < 0] <- 0
  # uptri_idx <- which(upper.tri(exp_diff_cor), arr.ind = TRUE)
  # exp_diff_cor_upper <- matrix(0,nrow = nrow(exp_diff_cor),ncol = ncol(exp_diff_cor),dimnames = list(rownames(exp_diff_cor),colnames(exp_diff_cor)))
  # exp_diff_cor_upper[uptri_idx] <- exp_diff_cor[uptri_idx]
  # exp_diff_cor_sparse <- Matrix(exp_diff_cor_upper, sparse = TRUE)
  # saveRDS(exp_diff_cor_sparse, file = file.path("./data/processed", cancer_type, "DEN_matrix_uptri.rds"))
  # exp_diff_cor <- cor(exp_diff_scaled, method = "pearson")
  # # uptri_idx <- which(upper.tri(exp_diff_cor), arr.ind = TRUE)
  # exp_diff_cor_df <- data.frame(gene1 = exp_genes[uptri_idx[, 1]],
  #                               gene2 = exp_genes[uptri_idx[, 2]],
  #                               cor = exp_diff_cor[uptri_idx])
  # exp_diff_cor_df <- exp_diff_cor_df[exp_diff_cor_df$cor > 0, ]
  # DEN_gz_file <- gzfile(file.path("./data/processed_data",cancer_type,"DEN.txt.gz"), "w")
  # write.table(
  #   exp_diff_cor_df,
  #   file = DEN_gz_file,
  #   sep = "\t",
  #   quote = FALSE,
  #   row.names = FALSE
  # )
  # close(DEN_gz_file)
  
}




# source('../IBTA-main/algorithms/CDN_fast.R')
# X <- matrix(list(), nrow = 1, ncol = 2)
# X[[1, 1]] <- scale(t(exp_data_tumor))
# X[[1, 1]][is.na(X[[1, 1]])] <- 0

# X[[1, 2]] <- scale(t(exp_data_normal))
# X[[1, 2]][is.na(X[[1, 2]])] <- 0
# delta <- CDN_fast(X, ifparallel = TRUE)