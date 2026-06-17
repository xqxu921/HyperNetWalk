# Load required libraries
library(igraph)
library(Matrix)
library(doSNOW)
library(foreach)
library(progress)
library(data.table)
library(dplyr)
library(readr)
library(purrr)

# Main function
run_hypernetwalk <- function(cancer_type,
                             mut_data_file,
                             exp_data_file,
                             output_dir,
                             max_degs = 500,
                             alpha = 0.5,
                             sigma = 0.1,
                             r = 0.15,
                             gamma = 0.15,
                             ifparallel = TRUE,
                             num_cores = 50) {
  # Start timing
  start_time <- proc.time()
  # Initialize logging
  logMessage("Initializing analysis for", cancer_type)
  
  # Load data
  logMessage("Loading data for", cancer_type)
  mut_data <- get_mut_data(mut_data_file)
  exp_data <- get_exp_data(exp_data_file)
  
  # Stage 1: Sample-specific RWR
  logMessage("\nStarting Stage 1: Sample-specific Random Walks")
  objs <- get_objs(cancer_type, mut_data, exp_data)
  mut_data <- objs$mut_mat
  # com_samples <- colnames(mut_data)
  exp_data <- objs$exp_data
  objs <- RW_on_layered_net(objs, ifparallel, num_cores, max_degs, alpha, r)
  
  # Stage 2: Hypergraph RWR
  logMessage("\nStarting Stage 2: Hypergraph Integration")
  logMessage("\nStarting Hypergraph Integration for Sample-specific Results")
  co_mut <- t(mut_data) %*% mut_data
  sim_mat <- cor(exp_data, method = "pearson")
  objs <- HRWR_sample_specific(objs, co_mut, sim_mat, cancer_type, sigma, gamma, ifparallel, num_cores)
  outfile <- file.path(output_dir, cancer_type)
  if (!dir.exists(outfile)) {
    dir.create(outfile, recursive = TRUE)
  }
  save(objs, file = file.path(outfile, "results.Rdata"))
  genes_ranking <- objs$mut_genes_ranks  
  ranking_df <- data.frame(
    sample = names(genes_ranking),
    ranking = sapply(genes_ranking, function(x)
      paste(x, collapse = ","))
  )
  write.table(
    ranking_df,
    file = file.path(outfile, "genes_ranking.txt"),
    sep = "\t",
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE
  )

  logMessage('\nStarting Hypergraph Integration for Cohort')
  scores <- HRWR_cohort(
    objs,
    gamma,
    num_cores
  )
  
  write.table(
    scores,
    file.path(outfile, "genes_ranking_cohort.txt"),
    quote = F,
    sep = "\t",
    row.names = T,
    col.names = F
  )
  logMessage("Analysis completed for", cancer_type)
  # End timing
  elapsed_time <- proc.time() - start_time
  logMessage("Total execution time for", cancer_type, ":", elapsed_time["elapsed"], "seconds")
}

# Data Loading Functions ----------------------------------------------------
get_mut_data <- function(data_file, min_mgs = 3) {
  logMessage("Loading mutation data from:", data_file)
  data <- tryCatch({
    read.delim(
      data_file,
      header = TRUE,
      as.is = TRUE,
      check.names = FALSE
    )
  }, error = function(e) {
    stop("Failed to load mutation data: ", e$message)
  })
  
  if (is.null(rownames(data))) {
    stop("The input data does not have row names (gene identifiers).")
  }
  # Filter invalid genes
  invalid_gene_index <- which(is.na(rownames(data)) |
                                rownames(data) == "" |
                                rownames(data) == ".")
  if (length(invalid_gene_index) > 0) {
    data <- data[-invalid_gene_index, ]
  }
  
  # Keep the samples with "-01A" suffix
  colnames(data) <- trimws(colnames(data))
  data <- data[, grepl("-01A$|-01$", colnames(data))]
  data <- Matrix(as.matrix(data), sparse = TRUE)
  data <- data[, colSums(data) >= min_mgs]
  # data <- data[, colSums(data) >= 3 & colSums(data) <= 500]
  data <- data[rowSums(data) != 0, ]
  logMessage(paste0(
    "Loaded mutation data with ",
    nrow(data),
    " genes and ",
    ncol(data),
    " samples"
  ))
  return(data)
}

# load exp data
get_exp_data <- function(data_file,filter=TRUE) {
  logMessage("Loading expression data from:", data_file)
  
  data <- tryCatch({
    read.delim(
      data_file,
      header = TRUE,
      as.is = TRUE,
      check.names = FALSE
    )
  }, error = function(e) {
    stop("Failed to load expression data: ", e$message)
  })
  if (is.null(rownames(data))) {
    stop("The input data does not have row names (gene identifiers).")
  }
  invalid_gene_index <- which(is.na(rownames(data)) |
                                rownames(data) == "" |
                                rownames(data) == ".")
  if (length(invalid_gene_index) > 0) {
    data <- data[-invalid_gene_index, ]
  }
  
  # only keep the samples with "-01A" suffix
  colnames(data) <- trimws(colnames(data))
  data <- as.matrix(data[, grepl("-01A$|-11A$|-01$|-11$", colnames(data))])
  # sample_id <- colnames(data)
  # colnames(exp_data) <- substr(sample_id, 1, 15)
  if (filter) {
    data <-
      as.matrix(data[which(rowSums(data >= 1) >= ncol(data) * 0.2), ])
  }
  logMessage(paste0(
    "Loaded expression data with ",
    nrow(data),
    " genes and ",
    ncol(data),
    " samples"
  ))
  return(data)
}

get_length_factor <- function(gene_length, l0, k) {
  if (gene_length <= l0) {
    return(1)
  } else {
    return(exp(-(gene_length - l0) / (k ^ 2)))
    # return(exp(-(gene_length-l0)^(1/k)))
  }
}

get_objs <- function(cancer_type, mut_data, exp_data) {
  mut_samples <- colnames(mut_data)
  exp_data_tumor <- exp_data[, grep("-01A$|-01$", colnames(exp_data)), drop =
                               F]
  exp_data_normal <- exp_data[, grep("-11A$|-11$", colnames(exp_data)), drop =
                                F]
  
  exp_samples <- colnames(exp_data_tumor)
  com_samples <- intersect(mut_samples, exp_samples)
  mut_data <- mut_data[, com_samples]
  logMessage(paste0(
    "Found ",
    length(com_samples),
    " common samples between mutation and expression data"
  ))
  
  deg_mat <- NULL
  if (ncol(exp_data_normal) > 1) {
    normal_mean <- rowMeans(exp_data_normal, na.rm = TRUE)
    normal_sd <- apply(exp_data_normal, 1, sd, na.rm = TRUE)
    z_normal <- (exp_data_tumor - normal_mean) / normal_sd
    z_normal[is.na(z_normal)] <- 0
    z_normal[is.infinite(z_normal)] <- 0
    deg_mat <- z_normal
  } else {
    z_tumor <- scale(t(exp_data_tumor), center = TRUE, scale = TRUE)
    z_tumor[is.na(z_tumor)] <- 0
    z_tumor[is.infinite(z_tumor)] <- 0
    z_tumor <- t(z_tumor)
    deg_mat <- z_tumor
    personalized <<- "lack of normal samples and use tumor mean and sd for z-score!"
  }
  exp_data_tumor <- exp_data_tumor[, com_samples]
  com_samples <- substr(com_samples, 1, 12)
  colnames(mut_data) <- com_samples
  colnames(exp_data_tumor) <- com_samples
  if (!is.null(deg_mat)) {
    colnames(deg_mat) <- substr(colnames(deg_mat), 1, 12)
    deg_mat <- deg_mat[, com_samples]
  }
  
  mut_genes <- list()
  for (i in 1:length(com_samples)) {
    sample_i <- com_samples[i]
    mut_genes_i <- rownames(mut_data)[which(mut_data[, sample_i] == 1)]
    mut_genes[[sample_i]] <- mut_genes_i
  }
  
  invalid_idx_m <- which(rowSums(mut_data) == 0)
  if (length(invalid_idx_m) > 0) {
    mut_data <- mut_data[-invalid_idx_m, ]
  }
  return(
    list(
      mut_mat = mut_data,
      exp_data = exp_data_tumor,
      exp_data_normal = exp_data_normal,
      deg_mat = deg_mat,
      mut_genes = mut_genes
    )
  )
}

RW_on_layered_net <- function(objs, ifparallel, num_cores, max_degs = 500, alpha = 0.5, theta = 0.15) {
  comb_adj_mat <- combine_STRING_and_omnipath()
  nodes_ppi <- rownames(comb_adj_mat)

  comb_ppi <- graph_from_adjacency_matrix(comb_adj_mat, weight = TRUE)
  comb_adj_mat0 <- comb_adj_mat
  comb_adj_mat0@x[comb_adj_mat0@x <= 1] <- 0
  comb_adj_mat0 <- drop0(comb_adj_mat0)
  comb_adj_lists <- build_adj_lists(comb_adj_mat0)
  # comb_adj_lists <- build_adj_lists(omni_dir_adj_mat)
  # comb_adj_mat@x[] <- 1
  deg_mat <- objs$deg_mat
  cor_mat <- WGCNA::cor(t(deg_mat))
  cor_mat <- abs(cor_mat)
  cor_mat[cor_mat < 0.1] <- 0.1
  idx <- summary(comb_adj_mat)
  g1 <- nodes_ppi[idx$i]
  g2 <- nodes_ppi[idx$j]
  w <- rep(0.1, length(g1))
  hit <- g1 %in% rownames(cor_mat) & g2 %in% colnames(cor_mat)
  w[hit] <- cor_mat[cbind(g1[hit], g2[hit])]
  comb_adj_mat@x <- w
  
  GRN_file <- "./data/RegNet_human_V2.txt"
  GRN_edges <- fread(GRN_file,
                     header = T,
                     sep = "\t",
                     select = 1:2)
  # 对GRN_edges添加一列全为1的权重
  GRN_edges$weight <- 1
  
  if (any(grep("/", GRN_edges$TF, fixed = TRUE))) {
    needs_split <- GRN_edges[grep("/", GRN_edges$TF, fixed = TRUE)]
    GRN_edges <- GRN_edges[!grepl("/", GRN_edges$TF, fixed = TRUE)]
    
    split_rows <- needs_split[, {
      tfs <- strsplit(TF, "/", fixed = TRUE)[[1]]
      .(TF = tfs,
        Target = Target,
        weight = weight)
    }, by = 1:nrow(needs_split)]
    split_rows[, nrow := NULL]
    GRN_edges <- rbind(GRN_edges, split_rows)
    GRN_edges <- unique(GRN_edges)
  }
  
  GRN_net <- graph_from_data_frame(GRN_edges, directed = TRUE)
  # GRN_net <- grn
  TF <- V(GRN_net)$name[which(igraph::degree(GRN_net, mode = "out") > 0)]
  
  comb_TF <- TF[grepl("::", TF, fixed = TRUE)]
  total_TF <- setdiff(TF, comb_TF)
  comb_TF_list <- strsplit(comb_TF, "::", fixed = TRUE)
  comb_TF_split <- unique(unlist(comb_TF_list))
  TF2comb <- lapply(comb_TF_split, function(tf) {
    comb_TF[sapply(comb_TF_list, function(x)
      tf %in% x)]
  })
  names(TF2comb) <- comb_TF_split
  total_TF <- union(total_TF, comb_TF_split)
  target <- setdiff(V(GRN_net)$name, TF)
  target <- target[which(igraph::degree(GRN_net, mode = "in")[target] >
                           0)]
  GRN_adj_mat <- as_adjacency_matrix(GRN_net, sparse = TRUE)
  GRN_adj_mat <- GRN_adj_mat[c(TF, target), c(TF, target)]
  
  GRN <- graph_from_adjacency_matrix(GRN_adj_mat)
  GRN_adj_lists <- build_adj_lists(GRN_adj_mat)
  
  
  mut_mat <- objs$mut_mat
  objs$exp_data <- objs$exp_data[which(rownames(objs$exp_data) %in% union(V(GRN)$name, rownames(mut_mat))), ]
  deg_mat <- objs$deg_mat
  # degs_coh <- objs$degs_coh
  mut_genes <- objs$mut_genes
  
  samples <- colnames(mut_mat)
  P0 <- matrix(
    0,
    nrow = nrow(mut_mat),
    ncol = ncol(mut_mat),
    dimnames = list(rownames(mut_mat), colnames(mut_mat))
  )
  P_tf <- matrix(
    0,
    nrow = length(TF),
    ncol = length(samples),
    dimnames = list(TF, samples)
  )
  
  if (ifparallel) {
    num_cores <- min(num_cores, parallel::detectCores() - 1)
    cl <- makeSOCKcluster(num_cores)
    registerDoSNOW(cl)
    on.exit(stopCluster(cl))
    
    pb <- txtProgressBar(max = length(samples), style = 3)
    progress <- function(n) {
      setTxtProgressBar(pb, n)
    }
    opts <- list(progress = progress)
    results <- foreach(
      i = 1:length(samples),
      .packages = c("igraph", "Matrix", "purrr"),
      .export = c(
        "randomwalk_withrestart",
        "normalize_rows",
        "multi_source_bfs_indices"
      ),
      .options.snow = opts
    ) %dopar% {
      sample_i <- samples[i]
      mut_genes_i <- mut_genes[[sample_i]]
      com_genes <- intersect(mut_genes_i, V(comb_ppi)$name)
      com_TF <- intersect(total_TF, V(comb_ppi)$name)
      d_out <- multi_source_bfs_indices(com_genes,
                                        comb_adj_lists$out,
                                        comb_adj_lists$idx_map,
                                        2)
      d_in <- multi_source_bfs_indices(com_TF, comb_adj_lists$'in', comb_adj_lists$idx_map, 2)
      sumd <- d_out + d_in
      candidate_idx <- which(!is.infinite(sumd) & (sumd <= 2))
      candidate_nodes <- names(d_out)[candidate_idx]
      
      if (length(intersect(com_genes, com_TF)) > 0) {
        # sub_d1 <- d1[,candidate_nodes]
        # sub_d2 <- d2[candidate_nodes,]
        keep_node <- sapply(1:length(candidate_nodes), function(k) {
          current_sum <- sumd[candidate_nodes[k]]
          if (current_sum != 2) {
            return(TRUE)
          }
          in_neighs_idx <- comb_adj_lists$'in'[[candidate_nodes[k]]]
          out_neighs_idx <- comb_adj_lists$out[[candidate_nodes[k]]]
          in_neighs <- names(comb_adj_lists$idx_map)[in_neighs_idx]
          out_neighs <- names(comb_adj_lists$idx_map)[out_neighs_idx]
          mut_in_neighs <- intersect(mut_genes_i, in_neighs)
          TF_out_neighs <- intersect(total_TF, out_neighs)
          if (length(mut_in_neighs) == 0 || length(TF_out_neighs) == 0) {
            #本身是TF或突变基因
            return(TRUE)
          }
          if (length(TF_out_neighs) > 1 || length(mut_in_neighs) > 1) {
            return(TRUE)
          }
          if (TF_out_neighs[1] == mut_in_neighs[1]) {
            return(FALSE)
          } else {
            return(TRUE)
          }
        })
        med_nodes <- candidate_nodes[keep_node]
      } else {
        med_nodes <- candidate_nodes
      }
      
      nodes_ppi_i <- union(mut_genes_i, med_nodes)
      # TF_i <- union(intersect(nodes_ppi_i,TF_act_i),intersect(mut_genes_i,total_TF))
      TF_i <- intersect(nodes_ppi_i, com_TF)
      if (length(TF_i) == 0) {
        P_i <- rep(1 / length(mut_genes_i), length(mut_genes_i))
        names(P_i) <- mut_genes_i
        return(list(P_i = P_i))
      }
      if (length(intersect(TF_i, comb_TF_split)) > 0) {
        comb_TF_split_i <- intersect(comb_TF_split, TF_i)
        comb_TF_i <- unique(unlist(TF2comb[comb_TF_split_i]))
        TF_nodes_i <- union(TF_i, comb_TF_i)
      } else {
        comb_TF_i <- NULL
        TF_nodes_i <- TF_i
      }
      
      d_out_grn <- multi_source_bfs_indices(TF_nodes_i, GRN_adj_lists$out, GRN_adj_lists$idx_map, 3)
      nodes_grn_i <- names(d_out_grn)[which(d_out_grn <= 3)]
      target_i <- intersect(nodes_grn_i, target)
      z_normal <- deg_mat[intersect(rownames(deg_mat), target_i), sample_i]
      thre <- max(sort(abs(z_normal), decreasing = TRUE)[min(max_degs, length(z_normal))], 1)
      deg_i <- names(z_normal)[which(abs(z_normal) >= thre)]
      deg_i <- intersect(target_i, deg_i)
      if (length(deg_i) == 0) {
        P_i <- rep(1 / length(mut_genes_i), length(mut_genes_i))
        names(P_i) <- mut_genes_i
        return(list(P_i = P_i))
      }
      
      d_in_grn <- multi_source_bfs_indices(deg_i, GRN_adj_lists$'in', GRN_adj_lists$idx_map, 3)
      d_grn <- d_out_grn + d_in_grn
      nodes_grn_i <- names(d_grn)[which(d_grn <= 3)]
      nodes_grn_i <- union(nodes_grn_i, TF_nodes_i)
      
      match_l1_TF <- match(TF_i, nodes_ppi_i)
      match_l2_TF <- match(TF_i, intersect(nodes_grn_i, TF))
      TF_all_i <- paste0("TF_", intersect(nodes_grn_i, TF))
      deg_i <- paste0("TAR_", deg_i)
      # deg_i[deg_i %in% target] <- paste0("TAR_",deg_i[deg_i %in% target])
      # deg_i[deg_i %in% TF] <- paste0("TF_",deg_i[deg_i %in% TF])
      target_i <- paste0("TAR_", intersect(nodes_grn_i, target))
      nodes_i <- c(nodes_ppi_i, TF_all_i, target_i)
      
      sub_GRN_adj <- GRN_adj_mat[c(substring(TF_all_i, 4), substring(target_i, 5)), c(substring(TF_all_i, 4), substring(target_i, 5))]
      sub_PPI_adj <- comb_adj_mat[intersect(nodes_ppi_i, nodes_ppi), intersect(nodes_ppi_i, nodes_ppi)]
      P <- Matrix(
        0,
        nrow = length(nodes_i),
        ncol = length(nodes_i),
        dimnames = list(nodes_i, nodes_i)
      )
      P[intersect(nodes_ppi_i, nodes_ppi), intersect(nodes_ppi_i, nodes_ppi)] <- sub_PPI_adj
      sub_PPI_adj <- P[nodes_ppi_i, nodes_ppi_i]
      if (alpha * (1-alpha) == 0) {
        message("Warning: alpha must be between 0 and 1, exclusive. Setting alpha to 0.5.") 
        alpha <- 0.5
      }
      P[cbind(match_l1_TF, length(nodes_ppi_i) + match_l2_TF)] <- alpha / (1-alpha) * pmax(1, colSums(sub_GRN_adj)[match_l2_TF])
      if (!is.null(intersect(comb_TF_i, nodes_grn_i))) {
        for (ctf in intersect(comb_TF_i, nodes_grn_i)) {
          ctf_l2_idx <- which(intersect(nodes_grn_i, TF) == ctf)
          tfs <- unique(strsplit(ctf, "::", fixed = TRUE)[[1]])
          ctf_l1_idx <- match(tfs, nodes_ppi_i)
          ctf_l1_idx <- ctf_l1_idx[!is.na(ctf_l1_idx)]
          if (length(ctf_l1_idx) > 0) {
            P[cbind(ctf_l1_idx, length(nodes_ppi_i) + ctf_l2_idx)] <- alpha / (1-alpha) * pmax(1, colSums(sub_GRN_adj)[ctf_l2_idx]) / length(ctf_l1_idx)
          }
        }
      }
      P[c(TF_all_i, target_i), c(TF_all_i, target_i)] <- sub_GRN_adj
      P[, mut_genes_i] <- 0
      mut_idx <- match(mut_genes_i, rownames(P))
      P[cbind(mut_idx, mut_idx)] <- 1
      P <- normalize_rows(t(P))
      
      zero_idx <- which(rowSums(P) == 0)
      P[cbind(zero_idx, zero_idx)] <- 1
      v0 <- numeric(length(nodes_i))
      names(v0) <- nodes_i
      v0[deg_i] <- 1
      # names(z_normal) <- paste0("TAR_",names(z_normal))
      # v0[deg_i] <- abs(z_normal[deg_i])
      v0 <- v0 / sum(v0)
      rw <- randomwalk_withrestart(P, mut_genes_i,theta = theta, v0 = v0)
      P_i <- rw$vt[intersect(mut_genes_i, nodes_i)]
      P_i <- P_i / sum(P_i)
      return(list(P_i = P_i, vt = rw$vt))
    }
    close(pb)
    for (i in 1:length(samples)) {
      result <- results[[i]]
      P0[names(result$P_i), samples[i]] <- result$P_i
      vt <- result$vt
      vt_TF <- numeric(length(TF))
      names(vt_TF) <- paste0("TF_", TF)
      vt_TF[intersect(names(vt), names(vt_TF))] <- vt[intersect(names(vt), names(vt_TF))]
      P_tf[, samples[i]] <- vt_TF
    }
  } else {
    pb <- progress_bar$new(
      format = "  Processing sample [:bar] :percent :elapsed",
      total = length(samples),
      clear = FALSE,
      width = 60
    )
    for (i in 1:length(samples)) {
      sample_i <- samples[i]
      mut_genes_i <- mut_genes[[sample_i]]
      com_genes <- intersect(mut_genes_i, V(comb_ppi)$name)
      com_TF <- intersect(total_TF, V(comb_ppi)$name)
      d_out <- multi_source_bfs_indices(com_genes, comb_adj_lists$out, comb_adj_lists$idx_map, 2)
      d_in <- multi_source_bfs_indices(com_TF, comb_adj_lists$'in', comb_adj_lists$idx_map, 2)
      sumd <- d_out + d_in
      candidate_idx <- which(!is.infinite(sumd) & (sumd <= 2))
      candidate_nodes <- names(d_out)[candidate_idx]
      # 进一步过滤候选节点
      if (length(intersect(com_genes, com_TF)) > 0) {
        keep_node <- sapply(1:length(candidate_nodes), function(k) {
          current_sum <- sumd[candidate_nodes[k]]
          if (current_sum != 2) {
            return(TRUE)
          }
          in_neighs_idx <- comb_adj_lists$'in'[[candidate_nodes[k]]]
          out_neighs_idx <- comb_adj_lists$out[[candidate_nodes[k]]]
          in_neighs <- names(comb_adj_lists$idx_map)[in_neighs_idx]
          out_neighs <- names(comb_adj_lists$idx_map)[out_neighs_idx]
          mut_in_neighs <- intersect(mut_genes_i, in_neighs)
          TF_out_neighs <- intersect(total_TF, out_neighs)
          if (length(mut_in_neighs) == 0 ||
              length(TF_out_neighs) == 0) {
            return(TRUE)
          }
          if (length(TF_out_neighs) > 1 ||
              length(mut_in_neighs) > 1) {
            return(TRUE)
          }
          if (TF_out_neighs[1] == mut_in_neighs[1]) {
            return(FALSE)
          } else {
            return(TRUE)
          }
        })
        med_nodes <- candidate_nodes[keep_node]
      } else {
        med_nodes <- candidate_nodes
      }
      nodes_ppi_i <- union(mut_genes_i, med_nodes)
      TF_i <- intersect(nodes_ppi_i, com_TF)
      if (length(TF_i) == 0) {
        P_i <- rep(1 / length(mut_genes_i), length(mut_genes_i))
        names(P_i) <- mut_genes_i
        P0[names(P_i), samples[i]] <- P_i
        next
      }
      if (length(intersect(TF_i, comb_TF_split)) > 0) {
        comb_TF_split_i <- intersect(comb_TF_split, TF_i)
        comb_TF_i <- unique(unlist(TF2comb[comb_TF_split_i]))
        TF_nodes_i <- union(TF_i, comb_TF_i)
      } else {
        comb_TF_i <- NULL
        TF_nodes_i <- TF_i
      }
      d_out_grn <- multi_source_bfs_indices(TF_nodes_i, GRN_adj_lists$out, GRN_adj_lists$idx_map, 3)
      nodes_grn_i <- names(d_out_grn)[which(d_out_grn <= 3)]
      target_i <- intersect(nodes_grn_i, target)
      z_normal <- deg_mat[intersect(rownames(deg_mat), target_i), sample_i]
      thre <- max(sort(abs(z_normal), decreasing = TRUE)[min(max_degs, length(z_normal))], 1)
      deg_i <- names(z_normal)[which(abs(z_normal) >= thre)]
      deg_i <- intersect(target_i, deg_i)
      if (length(deg_i) == 0) {
        P_i <- rep(1 / length(mut_genes_i), length(mut_genes_i))
        names(P_i) <- mut_genes_i
        P0[names(P_i), samples[i]] <- P_i
        next
      }
      d_in_grn <- multi_source_bfs_indices(deg_i, GRN_adj_lists$'in', GRN_adj_lists$idx_map, 3)
      d_grn <- d_out_grn + d_in_grn
      nodes_grn_i <- names(d_grn)[which(d_grn <= 3)]
      nodes_grn_i <- union(nodes_grn_i, TF_nodes_i)
      match_l1_TF <- match(TF_i, nodes_ppi_i)
      match_l2_TF <- match(TF_i, intersect(nodes_grn_i, TF))
      TF_all_i <- paste0("TF_", intersect(nodes_grn_i, TF))
      deg_i <- paste0("TAR_", deg_i)
      target_i <- paste0("TAR_", intersect(nodes_grn_i, target))
      nodes_i <- c(nodes_ppi_i, TF_all_i, target_i)
      sub_GRN_adj <- GRN_adj_mat[c(substring(TF_all_i, 4), substring(target_i, 5)), c(substring(TF_all_i, 4), substring(target_i, 5))]
      sub_PPI_adj <- comb_adj_mat[intersect(nodes_ppi_i, nodes_ppi), intersect(nodes_ppi_i, nodes_ppi)]
      P <- Matrix(
        0,
        nrow = length(nodes_i),
        ncol = length(nodes_i),
        dimnames = list(nodes_i, nodes_i)
      )
      P[intersect(nodes_ppi_i, nodes_ppi), intersect(nodes_ppi_i, nodes_ppi)] <- sub_PPI_adj
      sub_PPI_adj <- P[nodes_ppi_i, nodes_ppi_i]
      P[cbind(match_l1_TF, length(nodes_ppi_i) + match_l2_TF)] <- alpha / (1-alpha) * pmax(1, colSums(sub_GRN_adj)[match_l2_TF])
      if (!is.null(intersect(comb_TF_i, nodes_grn_i))) {
        for (ctf in intersect(comb_TF_i, nodes_grn_i)) {
          ctf_l2_idx <- which(intersect(nodes_grn_i, TF) == ctf)
          tfs <- unique(strsplit(ctf, "::", fixed = TRUE)[[1]])
          ctf_l1_idx <- match(tfs, nodes_ppi_i)
          ctf_l1_idx <- ctf_l1_idx[!is.na(ctf_l1_idx)]
          if (length(ctf_l1_idx) > 0) {
            P[cbind(ctf_l1_idx, length(nodes_ppi_i) + ctf_l2_idx)] <- alpha / (1-alpha) * pmax(1, colSums(sub_GRN_adj)[ctf_l2_idx]) / length(ctf_l1_idx)
          }
        }
      }
      P[c(TF_all_i, target_i), c(TF_all_i, target_i)] <- sub_GRN_adj
      P[, mut_genes_i] <- 0
      mut_idx <- match(mut_genes_i, rownames(P))
      P[cbind(mut_idx, mut_idx)] <- 1
      P <- normalize_rows(t(P))
      zero_idx <- which(rowSums(P) == 0)
      P[cbind(zero_idx, zero_idx)] <- 1
      v0 <- numeric(length(nodes_i))
      names(v0) <- nodes_i
      v0[deg_i] <- 1
      v0 <- v0 / sum(v0)
      rw <- randomwalk_withrestart(P, mut_genes_i, theta = theta, v0 = v0)
      P_i <- rw$vt[intersect(mut_genes_i, nodes_i)]
      P_i <- P_i / sum(P_i)
      P0[names(P_i), samples[i]] <- P_i
      vt <- rw$vt
      vt_TF <- numeric(length(TF))
      names(vt_TF) <- paste0("TF_", TF)
      vt_TF[intersect(names(vt), names(vt_TF))] <- vt[intersect(names(vt), names(vt_TF))]
      P_tf[, samples[i]] <- vt_TF
      pb$tick()
    }
  }
  
  objs$P0_mat <- P0
  objs$P_tf <- P_tf
  return(objs)
}

combine_STRING_and_omnipath <- function(){
  STRING_file <- "./data/STRINGv12.txt"
  STRING_edges <- fread(
      STRING_file,
      header = TRUE,
      sep = "\t",
      colClasses = c("character", "character", "numeric")
    )
  colnames(STRING_edges) <- c("protein1", "protein2", "score")
  STRING_edges <- STRING_edges %>% 
    filter(score >= 0.4) %>%
    filter(!grepl("_", protein1, fixed = TRUE) & !grepl("_", protein2, fixed = TRUE)) %>%
    filter(protein1 != "", protein2 != "", !is.na(protein1), !is.na(protein2)) %>%
    filter(protein1 != protein2) %>%
    mutate(weight = score ^ 2) %>%
    dplyr::select(protein1,protein2,weight) %>%
    distinct(protein1,protein2,.keep_all = TRUE)
  STRING_net <- graph_from_data_frame(STRING_edges, directed = FALSE)
  STRING_adj_mat <- as_adjacency_matrix(STRING_net, attr = "weight", sparse=TRUE)

  omnipath_url <- paste0(
    "https://omnipathdb.org/interactions?",
    "datasets=omnipath,pathwayextra,kinaseextra,ligrecextra",
    "&genesymbols=1",
    "&format=tsv"
  )
  data_dir <- "./data"
  omni_file <- file.path(data_dir,"omnipath_interactions.tsv")

  options(timeout = 1200)
  if (!file.exists(omni_file)){
    message("正在从 OmniPath 下载信号网络数据...")
    download.file(omnipath_url, omni_file, mode = "wb")
  } else {
    message("Using cached OmniPath file: ", omni_file)
  }

  message("Reading OmniPath data...")
  omni_raw <- fread(omni_file)
  required_cols <- c("source_genesymbol", "target_genesymbol")

  if (!all(required_cols %in% colnames(omni_raw))) {
    stop("OmniPath format changed — required columns missing.")
  }
  omni_raw <- omni_raw %>%
    # separate_rows(source_genesymbol, sep = "_") %>%
    # separate_rows(target_genesymbol, sep = "_") %>%
    filter(
      source_genesymbol != "", target_genesymbol != "",
      !is.na(source_genesymbol), !is.na(target_genesymbol),
      source_genesymbol != target_genesymbol
    )

  omni_edges <- omni_raw %>%
    dplyr::select(source = source_genesymbol, target = target_genesymbol) %>%
    distinct()
  message("Building directed signaling network...")
  omni_net <- graph_from_data_frame(omni_edges,directed = FALSE)
  omni_adj_mat <- as_adjacency_matrix(omni_net,sparse = TRUE)
  nodes_ppi <- union(V(STRING_net)$name,V(omni_net)$name)
  
  comb_adj_mat <- Matrix(0,nrow=length(nodes_ppi),ncol=length(nodes_ppi),dimnames=list(nodes_ppi,nodes_ppi))
  nz_STRING <- summary(STRING_adj_mat)
  row_idx_STRING <- match(rownames(STRING_adj_mat)[nz_STRING$i],nodes_ppi)
  col_idx_STRING <- match(colnames(STRING_adj_mat)[nz_STRING$j],nodes_ppi)
  comb_adj_mat[cbind(row_idx_STRING,col_idx_STRING)] <- nz_STRING$x
  # comb_adj_mat[cbind(row_idx_STRING,col_idx_STRING)] <- 0.1
  nz_omni <- summary(omni_adj_mat)
  row_idx_omni <- match(rownames(omni_adj_mat)[nz_omni$i],nodes_ppi)
  col_idx_omni <- match(colnames(omni_adj_mat)[nz_omni$j],nodes_ppi)
  comb_adj_mat[cbind(row_idx_omni,col_idx_omni)] <- 1

  omni_directed <- omni_raw %>%
    filter(consensus_direction == 1) %>%
    transmute(
      source = source_genesymbol,
      target = target_genesymbol
    ) %>%
    distinct(source, target, .keep_all = TRUE) %>%
    mutate(weight = 1.5)

  omni_directed_net <- graph_from_data_frame(omni_directed,directed = TRUE)
  omni_dir_adj_mat <- as_adjacency_matrix(omni_directed_net,sparse = TRUE, attr = "weight")
  nz_dir <- summary(omni_dir_adj_mat)
  row_idx_dir <- match(rownames(omni_dir_adj_mat)[nz_dir$i],nodes_ppi)
  col_idx_dir <- match(colnames(omni_dir_adj_mat)[nz_dir$j],nodes_ppi)
  comb_adj_mat[cbind(row_idx_dir,col_idx_dir)] <- 1.5
  reverse_exists <- omni_dir_adj_mat[cbind(nz_dir$j, nz_dir$i)] != 0
  one_way <- !reverse_exists
  comb_adj_mat[cbind(col_idx_dir[one_way],row_idx_dir[one_way])] <- 0
  comb_adj_mat <- drop0(comb_adj_mat)
  return(comb_adj_mat)
}

HRWR_sample_specific <- function(objs,
                                 co_mut,
                                 sim_mat,
                                 cancer_type,
                                 sigma,
                                 theta,
                                 ifparallel,
                                 num_cores) {
  cat("Starting patient-specific hypergraph random walks for samples in",
      cancer_type,
      "\n")
  samples <- colnames(co_mut)
  mut_mat <- objs$mut_mat
  mut_genes <- objs$mut_genes
  P0_mat <- objs$P0_mat[, 1:ncol(mut_mat), drop = F]
  P0_mat <- as(P0_mat, "sparseMatrix")
  N_mg <- colSums(mut_mat)
  IDF_mg <- log(ncol(mut_mat) / (rowSums(mut_mat) + 1))
  rank0_mat <- apply(P0_mat, 2, function(x)
    rank(-x, ties.method = "average"))
  
  rank_weights <- 1 / log2(rank0_mat + 1)
  rank_weights <- rank_weights * mut_mat
  
  prediction_all <-
    matrix(0,
           nrow(mut_mat),
           ncol(mut_mat),
           dimnames = list(rownames(mut_mat), samples))
  mut_genes_scores <- list()
  mut_genes_ranks <- list()
  
  hrw_not_converged <- c()
  if (ifparallel) {
    num_cores <- min(num_cores, parallel::detectCores() - 1)
    cl <- makeSOCKcluster(num_cores)
    registerDoSNOW(cl)
    on.exit(stopCluster(cl))
    
    pb <- txtProgressBar(max = length(samples), style = 3)
    progress <- function(n) {
      setTxtProgressBar(pb, n)
    }
    
    opts <- list(progress = progress)
    results <- foreach(
      i = 1:length(samples),
      .export = c(
        'get_hyperedge_weight',
        'get_hyper_P',
        'randomwalk_withrestart',
        'normalize_rows'
      ),
      .packages = c('igraph'),
      .options.snow = opts
    ) %dopar% {
      sample_i <- samples[i]
      mut_genes_i <- mut_genes[[sample_i]]
      V_i <- mut_genes_i
      candidate_neig <- colnames(co_mut)[which(co_mut[sample_i, ] >= 2)]
      E_i <- candidate_neig
      sim_mat_i <- as.matrix(sim_mat[sample_i, E_i, drop = F])
      H_i <- as.matrix(mut_mat[V_i, E_i, drop = F])
      W_ve_i <- as.matrix(P0_mat[V_i, E_i, drop = F])
      top_rsum <- apply(rank_weights[V_i, E_i, drop = F], 2, function(x)
        sum(head(sort(
          x, decreasing = TRUE
        ), 3)))
      W_e_i <- get_hyperedge_weight(sim_mat_i, sigma, ref_scores = top_rsum)
      P_i <- get_hyper_P(H_i, W_e_i, W_ve_i)
      diag(P_i) <- 0
      P_i <- normalize_rows(P_i)
      
      v0 <- numeric(length(V_i))
      names(v0) <- V_i
      v0[mut_genes_i] <- 1 / length(mut_genes_i)
      v0 <- v0 / sum(v0)
      hrw_i <- randomwalk_withrestart(P_i, mut_genes_i, theta = theta, v0 = v0)
      return(hrw_i)
    }
    close(pb)
    for (i in 1:length(samples)) {
      mut_genes_i <- mut_genes[[i]]
      prediction_i <- results[[i]]$vt
      prediction_all[names(prediction_i), i] <- prediction_i
      mut_genes_scores_i <- prediction_i[mut_genes_i]
      mut_genes_scores_i <-
        mut_genes_scores_i / sum(mut_genes_scores_i)
      mut_genes_scores_i <-
        sort(mut_genes_scores_i, decreasing = TRUE)
      mut_genes_scores[[samples[i]]] <- mut_genes_scores_i
      mut_genes_ranks[[samples[i]]] <- names(mut_genes_scores_i)
      
      if (!is.na(results[[i]]$D_not_converged)) {
        hrw_not_converged[samples[i]] <- results[[i]]$D_not_converged
      }
    }
  } else{
    pb <- progress_bar$new(
      format = "  Processing sample [:bar] :percent :elapsed",
      total = length(samples),
      clear = FALSE,
      width = 60
    )
    for (i in 1:length(samples)) {
      sample_i <- samples[i]
      mut_genes_i <- mut_genes[[sample_i]]
      V_i <- mut_genes_i
      candidate_neig <- colnames(co_mut)[which(co_mut[sample_i, ] >= 2)]
      E_i <- candidate_neig
      sim_mat_i <- as.matrix(sim_mat[sample_i, E_i, drop = F])
      H_i <- as.matrix(mut_mat[V_i, E_i, drop = F])
      W_ve_i <- as.matrix(P0_mat[V_i, E_i, drop = F])
      top_rsum <- apply(rank_weights[V_i, E_i, drop = F], 2, function(x)
        sum(head(sort(
          x, decreasing = TRUE
        ), 3)))
      W_e_i <- get_hyperedge_weight(sim_mat_i, sigma, ref_scores = top_rsum)
      P_i <- get_hyper_P(H_i, W_e_i, W_ve_i)
      diag(P_i) <- 0
      P_i <- normalize_rows(P_i)
      v0 <- numeric(length(V_i))
      names(v0) <- V_i
      v0[mut_genes_i] <- 1 / length(mut_genes_i)
      v0 <- v0 / sum(v0)
      hrw_i <- randomwalk_withrestart(P_i, mut_genes_i, theta = theta, v0 = v0)
      prediction_i <- hrw_i$vt
      
      prediction_all[names(prediction_i), i] <- prediction_i
      mut_genes_scores_i <- prediction_i[mut_genes_i]
      mut_genes_scores_i <-
        mut_genes_scores_i / sum(mut_genes_scores_i)
      mut_genes_scores_i <-
        sort(mut_genes_scores_i, decreasing = TRUE)
      mut_genes_scores[[samples[i]]] <- mut_genes_scores_i
      mut_genes_ranks[[samples[i]]] <- names(mut_genes_scores_i)
      pb$tick()
    }
  }
  
  if (length(hrw_not_converged) == 0) {
    cat(
      "Finished hypergraph random walks for samples in",
      cancer_type,
      "with all samples successfully converged",
      "\n"
    )
  } else{
    cat(
      "Finished hypergraph random walks for samples in",
      cancer_type,
      "with ",
      length(rw_not_converged),
      "samples not successfully converged",
      "\n"
    )
  }
  
  objs$prediction_mat <- prediction_all
  objs$mut_genes_scores <- mut_genes_scores
  objs$mut_genes_ranks <- mut_genes_ranks
  return(objs)
}

HRWR_cohort <- function(objs,
                        theta,
                        num_cores) {
  H <- objs$mut_mat
  V <- rownames(H)
  E <- colnames(H)
  W_e <- diag(ncol(H))
  W_ve <- objs$P0_mat[,1:ncol(H), drop = F]
  # prior <- objs$P0_mat[,ncol(H)+1]
  P <- get_hyper_P(H, W_e, W_ve)
  diag(P) <- 0
  P <- normalize_rows(P)
  
  v0 <- numeric(nrow(P))
  names(v0) <- rownames(P)
  v0[V] <- 1 / length(V)
  v0 <- v0 / sum(v0)
  rw <- randomwalk_withrestart(P, V, theta = theta, v0 = v0)
  scores <- rw$vt
  scores <- scores / sum(scores)
  scores <- sort(scores, decreasing = TRUE)
  scores <- as.matrix(scores)
  return(scores)
}

# compute hyperedge weights
get_hyperedge_weight <- function(sim_mat_i,
                                 sigma = 0.1,
                                 ref_scores = NULL) {
  W_e_i <- exp(-(1 - sim_mat_i) ^ 2 / (2 * sigma ^ 2))
  if (!is.null(ref_scores)) {
    # W_e_i <- matrix(ref_scores,nrow=1)
    W_e_i <- W_e_i * ref_scores
  }
  if (ncol(W_e_i) > 1) {
    W_e_i <- diag(c(W_e_i))
  }
  colnames(W_e_i) <- rownames(W_e_i) <- colnames(sim_mat_i)
  return(W_e_i)
}

# compute hypergraph transition matrix
get_hyper_P <- function(H, W_e, W_ve) {
  H_w <- H %*% W_e
  D_H_w <- normalize_rows(H_w)
  D_W_ve <- normalize_rows(t(W_ve))
  # D_W_ve <- t(W_ve)
  P <- D_H_w %*% D_W_ve
  return(P)
}

build_adj_lists <- function(adj_mat) {
  # adj_mat: sparse matrix (dgCMatrix/dgTMatrix) with rownames/colnames
  s <- summary(adj_mat)  # i, j, x
  n <- nrow(adj_mat)
  # create empty lists
  out_list <- vector("list", n)
  in_list  <- vector("list", n)
  for (k in seq_len(nrow(s))) {
    i <- s$i[k]
    j <- s$j[k]
    out_list[[i]] <- c(out_list[[i]], j)
    in_list[[j]]  <- c(in_list[[j]], i)
  }
  names(out_list) <- rownames(adj_mat)
  names(in_list)  <- rownames(adj_mat)
  idx_map <- setNames(seq_len(n), rownames(adj_mat)) # name -> index
  return(list(
    out = out_list,
    `in` = in_list,
    idx_map = idx_map
  ))
}

multi_source_bfs_indices <- function(source_names,
                                     adj_list,
                                     name_to_idx,
                                     max_dist) {
  # source_names: character vector of node names (may contain names not in name_to_idx)
  # adj_list: list of integer neighbor indices, names = node names
  # name_to_idx: named integer vector: name -> index
  n <- length(adj_list)
  dist <- rep(Inf, n)
  names(dist) <- names(adj_list)
  # map sources to indices (drop NA)
  src_idx <- name_to_idx[source_names]
  src_idx <- src_idx[!is.na(src_idx)]
  if (length(src_idx) == 0)
    return(dist)
  dist[src_idx] <- 0
  visited <- rep(FALSE, n)
  visited[src_idx] <- TRUE
  frontier <- unique(src_idx)
  for (d in 1:max_dist) {
    if (length(frontier) == 0)
      break
    next_frontier <- integer(0)
    # iterate frontier (usually small for small max_dist)
    for (u in frontier) {
      neigh <- adj_list[[u]]
      if (length(neigh) == 0)
        next
      new <- neigh[!visited[neigh]]
      if (length(new)) {
        dist[new] <- pmin(dist[new], d)
        visited[new] <- TRUE
        next_frontier <- c(next_frontier, new)
      }
    }
    frontier <- unique(next_frontier)
  }
  return(dist)
}

normalize_rows <- function(mat) {
  row_sums <- rowSums(mat)
  row_sums[row_sums == 0] <- 1
  return(mat / row_sums)
}

randomwalk_withrestart <-
  function(P_i,
           mut_genes_i,
           theta = 0.15,
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
      vt <- (1-theta) * t(P_i) %*% vt + theta * v0
      dis <- sum(abs(vt - v_old))
      # Distance <- append(Distance,dis)
      
      if (dis < 0.000001) {
        break
      }
      if (k > 100) {
        vt <- as.vector(vt)
        v_old <- as.vector(v_old)
        names(vt) <- names(v_old) <- names(v0)
        if (sum(abs(vt[mut_genes_i] - v_old[mut_genes_i])) < 0.000001) {
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

logMessage <- function(...,
                       level = "INFO",
                       timestamp = TRUE,
                       log_file = NULL,
                       append = TRUE) {
  # 参数说明：
  # ...: 可变参数，接收需要记录的日志内容
  # level: 日志级别（如 INFO, WARNING, ERROR）
  # timestamp: 是否添加时间戳（默认添加）
  # log_file: 日志文件路径（默认不写入文件）
  # append: 是否追加到文件（默认追加）
  
  # 组合日志内容
  msg <- paste(...)
  
  # 添加时间戳
  if (timestamp) {
    time_str <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    msg <- paste0("[", time_str, "] ", msg)
  }
  
  # 添加日志级别
  if (!is.null(level)) {
    msg <- paste0("[", level, "] ", msg)
  }
  
  # 输出到控制台
  message(msg)
  
  # 写入文件（如果指定了log_file）
  if (!is.null(log_file)) {
    if (!file.exists(log_file)) {
      dir.create(dirname(log_file),
                 showWarnings = FALSE,
                 recursive = TRUE)
    }
    write(msg, file = log_file, append = append)
  }
}
