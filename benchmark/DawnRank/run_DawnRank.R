setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
### using a small subset of the TCGA dataset and a small KEGG gene
### interaction network, We will show how to get DawnRank Results
library(DawnRank)
library(Matrix)
library(igraph)
source("../../model_linux.R")

network_file <- "../../data/STRINGv12.txt"
edges <- read.delim(
  network_file,
  header = TRUE,
  as.is = TRUE,
  check.names = FALSE
)
edge_idx <- which(edges[, 3] >= 0.4)
edges <- edges[edge_idx, ]
network <- graph_from_data_frame(edges, directed = FALSE)
gene_metadata_file <- "../../data/gencode.v36.annotation.gtf.gene.probemap"
# network_file <- "../../data/RegNet_human_v2.txt"
network_file <- "../PersonaDrive/data/dawnrank_network.csv"
edges <- read.delim(
  network_file,
  header = TRUE,
  as.is = TRUE,
  check.names = FALSE
)
network <- graph_from_data_frame(edges, directed = TRUE)

adj_mat <- as.matrix(as_adjacency_matrix(network))

# cancer_type <- "PANCAN"
cancers <- c(
  # "PANCAN",
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
  mut_data_file <- paste0("../../data/processed_data/", cancer_type, "/mut_data.tsv")
  exp_data_file <- paste0("../../data/processed_data/",
                          cancer_type,
                          "/exp_tpm_data.tsv")
  
  outfile_dir <- paste0("../../results/DawnRank_hrw/", cancer_type, "/")
  if (!dir.exists(outfile_dir)) {
    dir.create(outfile_dir)
  }
  mut_data <- get_mut_data(mut_data_file)
  colnames(mut_data) <- substr(colnames(mut_data), 1, 12)
  mut_data <- as.matrix(mut_data)
  
  exp_data <- read.delim(
    exp_data_file,
    header = TRUE,
    as.is = TRUE,
    check.names = FALSE
  )
  exp_data_tumor <- exp_data[, grepl("-01A$|-01$", colnames(exp_data))]
  exp_data_normal <- exp_data[, grepl("-11A$|-11$", colnames(exp_data))]
  colnames(exp_data_tumor) <- substr(colnames(exp_data_tumor), 1, 12)
  colnames(exp_data_normal) <- substr(colnames(exp_data_normal), 1, 12)
  exp_data_tumor <- as.matrix(exp_data_tumor)
  exp_data_normal <- as.matrix(exp_data_normal)
  # normalize the tumor and normal data to get the differential expression
  normalizedDawn <- DawnNormalize(tumorMat = exp_data_tumor, normalMat = exp_data_normal)
  
  com_samples <- intersect(colnames(normalizedDawn), colnames(mut_data))
  com_genes <- intersect(rownames(normalizedDawn), rownames(mut_data))
  normalizedDawn <- normalizedDawn[com_genes, com_samples]
  mut_data <- mut_data[com_genes, com_samples]
  A <- matrix(
    0,
    nrow = length(com_genes),
    ncol = length(com_genes),
    dimnames = list(com_genes, com_genes)
  )
  A[intersect(rownames(adj_mat), com_genes), intersect(rownames(adj_mat), com_genes)] <- adj_mat[intersect(rownames(adj_mat), com_genes), intersect(rownames(adj_mat), com_genes)]
  # get the DawnRank Score Get some coffee, this might take a while!
  dawnRankScore <- DawnRank(
    adjMatrix = A,
    mutationMatrix = mut_data,
    expressionMatrix = normalizedDawn,
    mu = 3
  )
  dawnRankFrame <- dawnRankScore[[3]]
  samples <- unique(dawnRankFrame$Patient)
  genes_ranking <- logical(0)
  for (i in 1:length(samples)) {
    sample <- samples[i]
    sample_dawnRankFrame <- dawnRankFrame[dawnRankFrame$Patient == sample, ]
    sample_dawnRankFrame <- sample_dawnRankFrame[order(sample_dawnRankFrame$Rank, decreasing = TRUE), ]
    genes_ranking_i <- paste(sample_dawnRankFrame$Gene, collapse = ",")
    genes_ranking <- rbind(genes_ranking, c(sample, genes_ranking_i))
  }
  write.table(
    genes_ranking,
    file = paste0(outfile_dir, "genes_ranking.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  # get the aggregate DawnRank scores Get some coffee, this might take a
  # while!
  # aggregateDawnRankScore <- condorcetRanking(scoreMatrix = dawnRankScore[[2]], mutationMatrix = mut_data)
  # score_cohort <- as.matrix(aggregateDawnRankScore[[2]])
  # write.table(
  #   score_cohort,
  #   file = paste0(outfile_dir, "genes_ranking_cohort.txt"),
  #   sep = "\t",
  #   quote = FALSE,
  #   row.names = TRUE,
  #   col.names = FALSE
  # )

  res_pers <- dawnRankScore[[1]] * mut_data
  source("../pernalized2cohort.R")
  personalized2cohort_score_new(res_pers, output_dir = outfile_dir)
}

#######
library(doSNOW)
library(foreach)
num_cores <- min(100, parallel::detectCores())
cl <- makeSOCKcluster(num_cores)
registerDoSNOW(cl)
pb <- txtProgressBar(max = ncol(expressionMatrix), style = 3)
progress <- function(n)
  setTxtProgressBar(pb, n)
opts <- list(progress = progress)
results <- foreach(
  i = 1:ncol(expressionMatrix),
  .export = c("Dawn"),
  .options.snow = opts
) %dopar% {
  Dawn(
    adjMatrix,
    expressionMatrix[, i],
    mutationMatrix[, i],
    damping = damping,
    maxit = maxit,
    epsilon = epsilon,
    patientTag = colnames(expressionMatrix)[i],
    goldStandard = goldStandard
  )
}
close(pb)
stopCluster(cl)
for (i in 1:ncol(expressionMatrix)) {
  allPercentiles <- cbind(allPercentiles, results[[i]][[1]][, 2])
  allRanks <- cbind(allRanks, results[[i]][[1]][, 1])
  allMutated <- rbind(allMutated, results[[i]][[2]])
  iterations <- append(iterations, results[[i]][[3]])
}

cl <- makeSOCKcluster(num_cores)
registerDoSNOW(cl)
pb <- txtProgressBar(max = nrow(cmat), style = 3)
progress <- function(n)
  setTxtProgressBar(pb, n)
opts <- list(progress = progress)
results <- foreach(
  i = 1:nrow(cmat),
  .export = c("pwc"),
  .options.snow = opts
) %dopar% {
  result <- list()
  for (j in 1:ncol(cmat)) {
    if (j > i) {
      result[[j - i]] <- pwc(trunc[i], trunc[j], scoreMatrix, mutationMatrix)
    }
  }
  return(result)
}
for (i in 1:(length(results) - 1)) {
  for (k in 1:length(results[[i]])) {
    vres <- results[[i]][[k]]
    cmat[i, i + k] <- vres[1]
    cmat[i + k, i] <- vres[2]
  }
}
close(pb)
stopCluster(cl)
