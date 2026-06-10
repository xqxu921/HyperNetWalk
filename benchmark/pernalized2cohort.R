library(Matrix)
# source("../../model.R")
get_top_list <- function(genes_ranking, top_num) {
  top_list <- list()
  for (i in 1:length(genes_ranking)) {
    genes_ranking_i <- genes_ranking[[i]]
    if (top_num < length(genes_ranking_i)) {
      top_list[[i]] <- genes_ranking_i[1:top_num]
    } else {
      top_list[[i]] <- genes_ranking_i
    }
  }
  names(top_list) <- names(genes_ranking)
  return(top_list)
}

get_condorcet_mat <- function(top_list) {
  genes <- unique(unlist(top_list))
  condorcet_mat <- Matrix(
    0,
    length(genes),
    length(genes),
    dimnames = list(genes, genes),
    sparse = TRUE
  )
  
  for (k in 1:length(top_list)) {
    gene_temp_list <- top_list[[k]]
    if (length(gene_temp_list) < 2) {
      next
    }
    for (i in 1:(length(gene_temp_list) - 1)) {
      gene_i <- gene_temp_list[i]
      gene_j <- gene_temp_list[(i + 1):length(gene_temp_list)]
      condorcet_mat[gene_i, gene_j] <- condorcet_mat[gene_i, gene_j] + 1
    }
  }
  
  return(condorcet_mat)
}

personalized2cohort_score <- function(input_dir) {
  data <- read.table(
    file.path(input_dir, "genes_ranking.txt"),
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE
  )
  genes_ranking <- list()
  for (i in 1:nrow(data)) {
    sample_i <- data[i, 1]
    genes_ranking_i <- unlist(strsplit(data[i, 2], ","))
    genes_ranking[[sample_i]] <- genes_ranking_i
  }
  message("Aggregate patient-level results to cohort-level results...")
  top_list <- get_top_list(genes_ranking, 100)
  W <- get_condorcet_mat(top_list)
  
  A <- (W > 0) + 0
  Degree_v <- colSums(W)
  Degree_v[Degree_v == 0] <- 1
  D_v_inverse <- Diagonal(n = length(Degree_v), x = 1 / Degree_v)
  P <- W %*% D_v_inverse
  colnames(P) <- colnames(W)
  v0 <- matrix(1 / nrow(W), nrow(W), 1)
  rownames(v0) <- rownames(W)
  vt <- v0
  for (k in 1:50) {
    v_old <- vt
    vt <- 0.9 * P %*% vt + 0.1 * v0
    dis <- sum(abs(vt - v_old))
    if (dis < 1e-6) {
      break
    }
  }
  vt <- vt[order(vt, decreasing = T), , drop = F]
  vt <- as.matrix(vt)
  write.table(vt,file.path(input_dir,"genes_ranking_cohort.txt"),quote = F,sep = "\t",row.names = T,col.names = F)
}

get_hyper_P <- function(H, W_e, W_ve) {
  H_w <- H %*% W_e
  D_H_w <- normalize_rows(H_w)
  D_W_ve <- normalize_rows(t(W_ve))
  # D_W_ve <- t(W_ve)
  P <- D_H_w %*% D_W_ve
  return(P)
}

normalize_rows <- function(mat) {
  row_sums <- rowSums(mat)
  row_sums[row_sums == 0] <- 1
  return(mat / row_sums)
}

randomwalk_withrestart <-
  function(P_i,
           mut_genes_i,
           theta = 0.85,
           v0 = NULL) {
    # cat("Starting random walk on sample:", sample_i, "\n")
    if (is.null(v0)) {
      v0 <- rep(0, nrow(P_i))
      names(v0) <- rownames(P_i)
      v0[mut_genes_i] <- 1 / length(mut_genes_i)
    } else {
      if (!identical(names(v0), rownames(P_i))) {
        stop("The names of v0 and P_i are not consistent.")
      }
    }
    
    Distance_not_converged <- NA
    vt <- v0
    for (k in 1:10000) {
      v_old <- vt
      vt <- theta * t(P_i) %*% vt + (1 - theta) * v0
      dis <- sum(abs(vt - v_old))
      # Distance <- append(Distance,dis)
      
      if (dis < 0.000001) {
        break
      }
      if (k > 100){
        vt <- as.vector(vt)
        v_old <- as.vector(v_old)
        names(vt) <- names(v_old) <- names(v0)
        if (sum(abs(vt[mut_genes_i]-v_old[mut_genes_i])) < 0.000001){
          break
        }
      }
      if (k == 10000) {
        Distance_not_converged <- dis
      }
    }
    vt <- as.vector(vt)
    names(vt) <- names(v0)
    return(list(vt = vt, D_not_converged = Distance_not_converged))
  }

personalized2cohort_score_new <- function(res_pers,output_dir){
  message("Aggregate patient-level results to cohort-level results by hypergraph...")
  E <- colnames(res_pers)
  V <- rownames(res_pers)
  H <- Matrix(0, nrow(res_pers), ncol(res_pers), dimnames = list(V, E), sparse = TRUE)
  nz_idx <- which(res_pers > 0, arr.ind = TRUE)
  H[nz_idx] <- 1
  W_v <- res_pers
  W_e <- diag(length(E))
  P <- get_hyper_P(H, W_e, W_v)
  diag(P) <- 0
  P <- normalize_rows(P)
  v0 <- rep(1 / nrow(P), nrow(P))
  names(v0) <- rownames(P)
  rw <- randomwalk_withrestart(P, V,v0 = v0)
  scores <- rw$vt
  scores <- scores / sum(scores)
  scores <- sort(scores, decreasing = TRUE)
  scores <- as.matrix(scores)
  write.table(scores,file.path(output_dir,"genes_ranking_cohort.txt"),quote = F,sep = "\t",row.names = T,col.names = F)
}
