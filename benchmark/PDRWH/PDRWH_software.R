#PDRWH
library(igraph)

get_mut_data <- function(data_file) {
  data <- read.delim(
    data_file,
    header = T,
    as.is = T,
    check.names = F
  )
  especial_gene_index <- which(is.na(rownames(data)) |
                                 rownames(data) == "")
  if (length(especial_gene_index) > 0) {
    data <- data[-especial_gene_index, ]
  } else {
    data <- data
  }
  data <- data[, grep("-01A$|-01$", colnames(data))]
  mut_data <- as.matrix(data[which(apply(data, 1, sum) != 0), which(apply(data, 2, sum) != 0)])
  
  return(mut_data)
}

get_exp_data <- function(data_file) {
  data <- read.delim(
    data_file,
    header = T,
    as.is = T,
    check.names = F
  )
  especial_gene_index <- which(is.na(rownames(data)) |
                                 rownames(data) == "" | rownames(data) == ".")
  if (length(especial_gene_index) > 0) {
    data <- data[-especial_gene_index, ]
  }
  data <- data[, grep("-01A$|-01$", colnames(data))]
  data <- data[which(rowSums(data > 0) >= 0.1 * ncol(data)), ]
  return(data)
}

get_normal_exp <- function(exp_data, com_samples) {
  sample_id <- colnames(exp_data)
  normal_idx <- which(substr(sample_id, 14, 14) == "1")
  normal_sample_id <- substr(sample_id[normal_idx], 1, 15)
  normal_exp <- exp_data[, normal_sample_id, drop = F]
  
  return(normal_exp)
}

get_ppi_network <- function(network_file, b_split) {
  Network <- read.table(network_file, header = T)
  edge_idx <- which(Network[, 3] >= b_split)
  Network <- Network[edge_idx, ]
  Graph <- graph.data.frame(Network)
  
  return(Graph)
}

combine_mut_exp <- function(data_mut_com, exp_data, cancer_type) {
  genes <- intersect(rownames(data_mut_com), rownames(exp_data))
  exp_data_scaled <- scale(t(exp_data[genes, ]), center = T, scale = T)
  exp_data_scaled[is.na(exp_data_scaled)] <- 0
  a <- as.vector(exp_data_scaled)
  z_thres <- min(abs(quantile(a, 0.95)), abs(quantile(a, 0.05)))
  outlying_idx <- (abs(exp_data_scaled) >= z_thres) + 0
  samples <- colnames(data_mut_com)
  
  data_mut_exp <- data_mut_com
  for (i in 1:length(samples)) {
    sample_i <- samples[i]
    outlying_gene <- genes[which(outlying_idx[sample_i, ] == 1)]
    data_mut_exp[outlying_gene, sample_i] <- 1
  }
  
  tmp <- which(rowSums(data_mut_exp) == 0)
  if (length(tmp) > 0) {
    data_mut_exp <- data_mut_exp[-tmp, ]
  }
  col_sum <- colSums(data_mut_exp)
  
  return(data_mut_exp)
}

get_edge_W <- function(sample_i, exp_data, H) {
  genes <- unique(rownames(exp_data))
  if (length(genes) <= 1) {
    edge_w <- matrix(1, dim(H)[2], 1)
  } else if (ncol(exp_data) == 1) {
    edge_w <- matrix(1, 1, 1)
  } else {
    exp_cor <- cor(exp_data[genes, ], method = "pearson")
    exp_cor_idx <- ((exp_cor) >= 0.4) + 0
    exp_cor <- abs(exp_cor) * exp_cor_idx
    edge_w <- exp_cor[, sample_i]
    edge_w_1 <- -(1 - edge_w)^2
    edge_w_2 <- edge_w_1 / (0.1^2 * 2)
    edge_w <- exp(edge_w_2)
  }
  edge_w[is.na(edge_w)] <- 0
  Edge_W <- diag(c(t(edge_w)))
  
  return(Edge_W)
}

get_vertex_W <- function(sample_i,
                         data_mut_exp,
                         com_mut_samples,
                         Graph,
                         sample_i_mut_genes) {
  mut_data_i <- data_mut_exp[, sample_i]
  vertex_W <- c()
  for (k in com_mut_samples) {
    vertex_w <- rep(0, length(sample_i_mut_genes))
    names(vertex_w) <- sample_i_mut_genes
    mut_data_k <- data_mut_exp[, k]
    mut_exp_gene <- names(mut_data_k[which(mut_data_k == 1)])
    
    v_i <- intersect(mut_exp_gene, V(Graph)$name) # questionable: 这里才考虑与PPI网络的点交集，前边sample_i_mut_genes还会包括不在PPI网络中的基因
    graph_i <- induced_subgraph(Graph, v_i)
    graph_i_degree <- degree(graph_i)[v_i]
    graph_i_degree[which(graph_i_degree == 0)] = 0.01
    
    co_mut_genes <- intersect(v_i, sample_i_mut_genes)
    vertex_w[co_mut_genes] <- graph_i_degree[co_mut_genes]
    vertex_w[setdiff(sample_i_mut_genes, co_mut_genes)] = 0.01 # questionable:在sample i中的突变基因，但是不在PPI或者当前neighbor sample中
    vertex_W <- cbind(vertex_W, vertex_w)
  }
  colnames(vertex_W) <- com_mut_samples
  
  return(vertex_W)
}

get_P <- function(H, Edge_W, Vertex_W) {
  H_W <- H %*% Edge_W
  Degree_v <- apply(H_W, 1, sum)
  D_v_inverse <- diag(1 / Degree_v)
  
  probability_hyperedge <- D_v_inverse %*% H_W
  rownames(probability_hyperedge) <- rownames(H)
  
  Degree_ve <- apply(Vertex_W, 2, sum)
  D_ve_inverse <- diag(1 / Degree_ve, nrow = length(Degree_ve))
  
  probability_vertex <- D_ve_inverse %*% t(Vertex_W)
  rownames(probability_vertex) <- colnames(H)
  
  P <- probability_hyperedge %*% probability_vertex
  
  return(as.matrix(P))
}

get_hyper_randomwalk <- function(P, H, theta = 0.85, mut_idx) {
  v0 <- rep(0, nrow(H))
  names(v0) <- rownames(H)
  v0[mut_idx] <- 1 / length(mut_idx)
  teleport <- v0
  
  Distance <- c()
  vi <- v0
  for (k in 1:20) {
    vj <- vi
    vi <- theta * t(P) %*% vi + (1 - theta) * teleport
    dis <- sum(abs(vj - vi))
    Distance <- append(Distance, dis)
    if (dis < 0.000001) {
      cat("Successfully converged after", k, "iterations!", "\n")
      break
    }
  }
  
  vi <- vi[mut_idx, , drop = F]
  Importance <- data.frame(vi[order(vi, decreasing = T), ])
  Importance <- Importance / sum(Importance)
  colnames(Importance) <- c("Importance")
  
  return(Importance)
}

hygraph_randomwalk <- function(sample_i,
                               data_mut_exp,
                               mut_idx,
                               co_mut,
                               exp_data,
                               Graph) {
  sample_i_mut_genes <- rownames(data_mut_exp)[which(data_mut_exp[, sample_i] == 1)]
  com_mut_sample_i <- which(co_mut[sample_i, ] >= 1) # neighbor samples
  com_mut_samples <- colnames(co_mut)[com_mut_sample_i]
  
  H <- data_mut_exp[sample_i_mut_genes, com_mut_samples, drop = F]
  E <- exp_data[, com_mut_samples, drop = F]
  Vertex_W <- get_vertex_W(sample_i,
                           data_mut_exp,
                           com_mut_samples,
                           Graph,
                           sample_i_mut_genes)
  Edge_W <- get_edge_W(sample_i, E, H)
  P <- get_P(H, Edge_W, Vertex_W)
  
  importance_score <- get_hyper_randomwalk(P, H, theta = 0.85, mut_idx)
  
  return(importance_score)
}

PDRWHscore <- function(cancer_type,
                       mut_data_file,
                       exp_data_file,
                       network_file,
                       outfile_dir) {
  mut_data <- get_mut_data(mut_data_file)
  mut_data_morethan_2 <- colnames(mut_data)[which(apply(mut_data, 2, sum) >=
                                                    3)] #突变基因数大于2的样本
  exp_data <- get_exp_data(exp_data_file)
  Graph <- get_ppi_network(network_file, b_split = 0.7)
  genes_degree <- degree(Graph)
  
  com_samples <- intersect(mut_data_morethan_2, colnames(exp_data))
  data_mut_com <- mut_data[, com_samples]
  
  exp_data <- exp_data[, com_samples]
  com_samples <- substr(com_samples, 1, 12) #结尾都是-01表示肿瘤样本
  colnames(exp_data) <- com_samples
  colnames(data_mut_com) <- com_samples
  
  data_mut_idx <- data_mut_com
  data_mut_exp <- combine_mut_exp(data_mut_com, exp_data, cancer_type)
  co_mut <- t(data_mut_com) %*% data_mut_com
  
  Importance_Score <- list()
  for (i in 1:length(com_samples)) {
    sample_i <- com_samples[i]
    print(paste0(i, " ", sample_i))
    mut_idx <- rownames(data_mut_idx)[which(data_mut_idx[, sample_i] == 1)]
    importance_score <- hygraph_randomwalk(sample_i, data_mut_exp, mut_idx, co_mut, exp_data, Graph)
    Importance_Score[[i]] <- importance_score
  }
  names(Importance_Score) <- com_samples
  
  save(Importance_Score,
       file = paste(outfile_dir, "Importance_Score.Rdata", sep = ""))
}

get_top_list <- function(PDRWH_score, top_num) {
  PDRWH_score <- as.matrix(PDRWH_score)
  top_list <- list()
  for (i in 1:ncol(PDRWH_score)) {
    temp <- PDRWH_score[order(PDRWH_score[, i], decreasing = T), i]
    target_genes <- names(temp)[which(temp != 0)]
    if (top_num < length(target_genes)) {
      top_list[[i]] <- target_genes[1:top_num]
    } else{
      top_list[[i]] <- target_genes
    }
  }
  
  return(top_list)
}

top_condorcet <- function(top_list) {
  genes <- unique(unlist(top_list))
  print(length(genes))
  condorcet_mat <- matrix(0, length(genes), length(genes), dimnames = list(genes, genes))
  
  for (K in 1:(length(top_list))) {
    gene_temp_list <- top_list[[K]]
    if (length(gene_temp_list) < 2) {
      next
    }
    for (i in 1:(length(gene_temp_list) - 1)) {
      gene_i <- gene_temp_list[i]
      for (j in (i + 1):length(gene_temp_list)) {
        gene_j <- gene_temp_list[j]
        condorcet_mat[gene_i, gene_j] <- condorcet_mat[gene_i, gene_j] + 1
      }
    }
  }
  
  return(condorcet_mat)
}

PDRWH_cohort_score <- function(input_dir) {
  PDRWH_score <- read.table(paste0(input_dir, "Importance_Score_mat.txt"))
  print("Get PDRWH-scores on the cohort-level... ")
  top_list <- get_top_list(PDRWH_score, 100)
  condorcet_mat <- top_condorcet(top_list)
  
  A <- (condorcet_mat > 0) + 0
  Degree_v <- apply(A, 2, sum)
  for (i in 1:length(Degree_v)) {
    if (Degree_v[i] == 0) {
      Degree_v[i] <- 1
    }
  }
  D_v_inverse <- diag(1 / Degree_v)
  X <- matrix(1 / length(rownames(A)), length(rownames(A)), 1)
  rownames(X) <- rownames(A)
  
  for (i in 1:30) {
    DX <- D_v_inverse %*% X
    X <- 0.9 * A %*% DX + 0.1 * X
  }
  gene_ordered <- X[order(X, decreasing = T), , drop = F]
  
  write.table(
    gene_ordered,
    paste0(input_dir, "Importance_Score_cohort.txt"),
    quote = F,
    sep = "\t",
    row.names = T,
    col.names = T
  )
  print(
    paste0(
      "Finished! The output file is saved in: ",
      input_dir,
      "Importance_Score_cohort.txt"
    )
  )
}

resultTransform <- function(cancer_type, input_dir, outfile_dir) {
  load(paste0(input_dir, "Importance_Score.Rdata"))
  all_genes <- c()
  for (i in (1:length(Importance_Score[]))) {
    genes <- rownames(Importance_Score[[i]])
    all_genes <- union(all_genes, genes)
  }
  
  res_matrix <- matrix(0, length(all_genes), length(Importance_Score[]))
  rownames(res_matrix) <- all_genes
  colnames(res_matrix) <- names(Importance_Score[])
  for (i in (1:length(Importance_Score[]))) {
    for (j in (1:length(rownames(Importance_Score[[i]])))) {
      res_matrix[rownames(Importance_Score[[i]])[j], i] <- Importance_Score[[i]][j, ]
    }
  }
  
  write.table(
    res_matrix,
    paste0(outfile_dir, "Importance_Score_mat.txt"),
    quote = T,
    sep = "\t"
  )
  print(paste0(
    "The output file is saved in: ",
    outfile_dir,
    "Importance_Score_mat.txt"
  ))
}
