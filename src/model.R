# Load required libraries
library(igraph)
library(progress)
library(limma)
library(Matrix)
library(doSNOW)
library(foreach)
library(data.table)
library(readr)
library(dplyr)
library(purrr)
library(WGCNA)

# Main function -----------------------------------------------------------
DGScore <- function(cancer_type,
                    outfile_dir,
                    mut_data_file,
                    exp_data_file,
                    PPI_file,
                    GRN_file,
                    DEN_file = NULL,
                    ME_file = NULL,
                    PPI_thre = 0.4,
                    DEG = TRUE,
                    AEG = FALSE,
                    isolated_restart = TRUE,
                    beta = 0.8,
                    # Weight for betweenness vs harmonic centrality
                    lambda1 = 0.8,
                    lambda2 = 0.8,
                    # Weight for PPI vs GRN networks
                    delta = 0.1,
                    # Hyperedge weight parameter
                    theta = 0.85,
                    ME_as_prior = TRUE,
                    # Whether to use ME network as prior
                    DEG_as_prior = TRUE,
                    reverse = TRUE,
                    max_dist = 3,
                    ifweighted_dist = FALSE,
                    # Whether to reverse GRN direction
                    ifparallel = TRUE,
                    num_cores = 20,
                    scale_genelength_stage1 = TRUE,
                    scale_genelength_stage2 = FALSE,
                    k = 100,
                    l0 = 100000,
                    no_selftrans = TRUE,
                    max_degs = 500,
                    use_cpdegs1 = FALSE,
                    use_cpdegs2 = FALSE,
                    exp_genes_in_GRN = TRUE,
                    weighted = FALSE,
                    absorbed = TRUE,
                    use_aegs = TRUE,
                    use_degs = TRUE,
                    prior_weighted = TRUE,
                    level = "all") {
  # 记录运行时间
  start_time <- Sys.time()
  # Initialize logging
  logMessage("Initializing analysis for", cancer_type)
  
  # Load and preprocess data -------------------------------------------------
  logMessage("Loading input data...")
  mut_data <- get_mut_data(mut_data_file)
  # mut_samples <- colnames(mut_data)
  exp_data <- get_exp_data(exp_data_file)
  # exp_data <- exp_data[which(rowSums(exp_data >= 1) >= ncol(exp_data) * 0.2 | rownames(exp_data) %in% rownames(mut_data)), ]
  exp_data <- exp_data[which(rowSums(exp_data >= 1) >= ncol(exp_data) * 0.2),]
  # exp_samples <- colnames(exp_data)[grep("-01A$", colnames(exp_data))]
  PPI_network <- get_network(PPI_file, cutoff = PPI_thre)
  rGRN <- get_reversed_GRN(GRN_file, reverse)
  # rGRN <- get_reversed_GRN_DoRothEA(reverse,levels = c("A","B","C"))
  # rGRN_dist <- get_reversed_GRN_dist(GRN_file, reverse, max_dist)
  if (exp_genes_in_GRN){
    exp_data <- exp_data[intersect(rownames(exp_data),union(V(rGRN)$name,rownames(mut_data))),]
  }
  # DEN_mat <- readRDS(DEN_file)
  # DEN_mat@x[DEN_mat@x < 0.5] <- 0
  # DEN_mat <- drop0(DEN_mat)

  gene_length_df <- read_tsv("./data/metadata/gencode.v36.annotation.gtf.gene.probemap",
                             show_col_types = FALSE) %>%
    mutate(length = chromEnd - chromStart + 1) %>%
    group_by(gene) %>%
    summarise(length = max(length), .groups = "drop")
  gene_lengths <- gene_length_df$length
  names(gene_lengths) <- gene_length_df$gene
  gene_length_factors <- map_dbl(gene_lengths, get_length_factor, l0, k)
  mut_data <- mut_data[rownames(mut_data) %in% gene_length_df$gene, ]
  exp_data <- exp_data[rownames(exp_data) %in% gene_length_df$gene, ]

  ME_centrality <- get_ME_centrality(ME_file, beta)
  
  # Stage 1: Initial Random Walk ---------------------------------------------
  logMessage("\nStarting Stage 1: Sample-specific Random Walks")
  nodes <- get_nodes(cancer_type,
                     mut_data,
                     exp_data,
                     rGRN,
                     AEG,
                     DEG,
                     max_degs,
                     use_cpdegs1)
  # genes <- union(rownames(mut_data),unique(unlist(nodes$degs)))
  # genes <- Reduce(union,list(rownames(mut_data),unique(unlist(nodes$degs)),unique(unlist(nodes$aegs))))
  DEN_mat <- nodes$cen_mat
  nodes$cen_mat <- NULL
  DEN_mat@x[DEN_mat@x < 0.5] <- 0
  DEN_mat <- drop0(DEN_mat)
  genes <- Reduce(union,list(rownames(mut_data),unique(unlist(nodes$degs)),nodes$degs_coh,unique(unlist(nodes$aegs))))
  DEN_mat <- DEN_mat[intersect(rownames(DEN_mat),genes), intersect(rownames(DEN_mat),genes)]
  nodes <- randomwalk_persample(
      nodes,
      PPI_network,
      rGRN,
      DEN_mat,
      ifweighted_dist,
      max_dist,
      ME_centrality,
      gene_length_factors,
      ME_as_prior,
      DEG_as_prior,
      scale_genelength_stage1,
      cancer_type,
      isolated_restart,
      weighted,
      absorbed,
      lambda1,
      lambda2,
      max_degs,
      ifparallel,
      num_cores,
      use_cpdegs2
  )
  # Stage 2: Hypergraph Random Walk ------------------------------------------
  logMessage("\nStarting Stage 2: Hypergraph Integration")
  # Save results -----------------------------------------------------------
  outfile <- file.path(outfile_dir, cancer_type)
  if (!dir.exists(outfile)) {
    dir.create(outfile, recursive = TRUE)
  }
  if (level == "all" | level == "individual"){
    mut_data <- nodes$mut_mat
    com_samples <- colnames(mut_data)
    # colnames(mut_data) <- substr(colnames(mut_data), 1, 12)
    # mut_data <- mut_data[, com_samples]
    exp_data <- exp_data[, grep("-01A$", colnames(exp_data))]
    colnames(exp_data) <- substr(colnames(exp_data), 1, 12)
    exp_data <- exp_data[, com_samples]
    co_mut <- t(mut_data) %*% mut_data
    sim_mat <- cor(exp_data, method = "pearson")
    nodes <- prediction_sample_specific(
      nodes,
      co_mut,
      sim_mat,
      gene_length_factors,
      scale_genelength_stage2,
      no_selftrans,
      cancer_type,
      delta,
      theta,
      ifparallel,
      num_cores
    )
    # save(nodes, file = file.path(outfile, "results.Rdata"))
    genes_ranking <- nodes$mut_genes_ranks
    ranking_df <- data.frame(
      sample = names(genes_ranking),
      ranking = sapply(genes_ranking, function(x)
        paste(x, collapse = ","))
    )
    write.table(
      ranking_df,
      file = file.path(outfile, "individual_results.txt"),
      sep = "\t",
      col.names = FALSE,
      row.names = FALSE,
      quote = FALSE
    )
  }
  
  if (level == "all" | level == "cohort"){
    logMessage('\nStarting Hypergraph Integration for Cohort')
    scores <- hypergraph_rw_cohort_singlecancer(
      nodes,
      scale_genelength_stage2,
      gene_length_factors,
      no_selftrans,
      theta,
      num_cores
    )
    write.table(
      scores,
      file.path(outfile, "cohort_results.txt"),
      quote = F,
      sep = "\t",
      row.names = T,
      col.names = F
    )  
  }
  logMessage("\nAnalysis completed successfully for", cancer_type)
  end_time <- Sys.time()
  logMessage("Total runtime:", round(difftime(end_time, start_time, units = "mins"),2), "minutes")
}

# Data Loading Functions ----------------------------------------------------
get_mut_data <- function(data_file) {
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
  data <- data[, colSums(data) >= 3]
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
get_exp_data <- function(data_file) {
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
  if (nchar(colnames(data)[1]) > 15){
    data <- as.matrix(data[, grepl("-01A$|-11A$", colnames(data))])
  } else {
    data <- as.matrix(data[, grepl("-01$|-11$", colnames(data))])
  }
  
  # sample_id <- colnames(data)
  # colnames(exp_data) <- substr(sample_id, 1, 15)
  # data <-
  #   as.matrix(data[which(rowSums(data >= 1) >= ncol(data) * 0.2), ])
  
  logMessage(paste0(
    "Loaded expression data with ",
    nrow(data),
    " genes and ",
    ncol(data),
    " samples"
  ))
  return(data)
}

# load PPI network
# get_PPI_network <- function(network_file, cutoff) {
#   PPI_edges <- read.delim(network_file, header = T)
#   edge_idx <- which(PPI_edges[, 3] >= cutoff)
#   PPI_edges <- PPI_edges[edge_idx, ]
#   PPI_network <- graph_from_data_frame(PPI_edges, directed = FALSE)
#   logMessage(paste0(
#     "Loaded PPI network with ",
#     vcount(PPI_network),
#     " nodes and ",
#     ecount(PPI_network),
#     " edges"
#   ))
#   return(PPI_network)
# }

# load undirected network
get_network <- function(network_file, cutoff=0) {
  edges <- tryCatch({
    fread(
      network_file,
      header = TRUE,
      sep = "\t",
      colClasses = c("character", "character", "numeric")
    )
  }, error = function(e) {
    stop("Failed to read the network file: ", e$message)
  })
  if(cutoff > 0){
    edge_idx <- which(abs(edges[, 3]) >= cutoff)
    edges <- edges[edge_idx, ]
    if (nrow(edges) == 0) {
      warning("No edges meet the cutoff criteria. Returning an empty graph.")
      return(make_empty_graph(directed = FALSE))
    }
  }
  colnames(edges)[1:3] <- c("node1", "node2", "weight")
  network <- graph_from_data_frame(edges, directed = FALSE)
  logMessage(paste0(
    "Loaded network with ",
    vcount(network),
    " nodes and ",
    ecount(network),
    " edges"
  ))
  return(network)
}

# load GRN network
get_reversed_GRN <- function(GRN_file, reverse = TRUE) {
  edges <- fread(
    GRN_file,
    header = T,
    sep = "\t",
    select = 1:2
  )
  net <- graph_from_data_frame(edges, directed = TRUE)
  logMessage(paste0(
    "Loaded GRN network with ",
    vcount(net),
    " nodes and ",
    ecount(net),
    " edges"
  ))
  if (reverse) {
    net <- reverse_edges(net)
  }
  return(net)
}


get_ME_centrality <- function(ME_file, beta) {
  if (is.null(ME_file)) {
    logMessage("No mutual exclusive network provided, skipping centrality calculation")
    ME_as_prior <<- FALSE
    scale_genelength_stage1 <<- FALSE
    return(NULL)
  }
  ME_edges <- fread(ME_file, header = T, sep = "\t")
  if (nrow(ME_edges) == 0) {
    logMessage("No mutual exclusive network provided, skipping centrality calculation")
    ME_as_prior <<- FALSE
    scale_genelength_stage1 <<- FALSE
    return(NULL)
  }
  ME_net <- graph_from_data_frame(ME_edges, directed = FALSE)
  if (vcount(ME_net) < 200){
    logMessage("ME network too small, skipping centrality calculation")
    ME_as_prior <<- FALSE
    scale_genelength_stage1 <<- FALSE
    return(NULL)
  }

  logMessage(paste0(
    "Loaded ME network with ",
    vcount(ME_net),
    " nodes and ",
    ecount(ME_net),
    " edges"
  ))
  bc <- betweenness(
    ME_net,
    v = V(ME_net),
    directed = FALSE,
    weights = NULL,
    normalized = TRUE,
    cutoff = -1
  )
  bc[bc == 0] <- min(bc[bc > 0]) / 2
  bc <- bc / sum(bc)
  # cc <- closeness(ME_net, weights = NULL, normalized = TRUE)
  # cc <- cc / sum(cc)
  hc <- harmonic_centrality(ME_net, weights = NULL, normalized = TRUE)
  hc <- hc / sum(hc)
  combine <- beta * bc + (1 - beta) * hc
  return(combine)
}

# combine mut genes and DEGs for each sample
get_nodes <- function(cancer_type,
                      mut_data,
                      exp_data,
                      rGRN,
                      AEG,
                      DEG,
                      max_degs,
                      use_cpdegs) {
  mut_samples <- colnames(mut_data)
  exp_data_tumor <- exp_data[, grep("-01A$|-01$", colnames(exp_data)),drop=F]
  exp_data_normal <- exp_data[, grep("-11A$|-11$", colnames(exp_data)),drop=F]
  
  exp_samples <- colnames(exp_data_tumor)
  com_samples <- intersect(mut_samples, exp_samples)
  mut_data <- mut_data[, com_samples]
  logMessage(paste0(
    "Found ",
    length(com_samples),
    " common samples between mutation and expression data"
  ))
  deg_mat <- NULL
  deg_genes <- NULL
  if (DEG) {
    logMessage("Conducting differential expression analysis...")
    group <- factor(c(rep("tumor", ncol(
      exp_data_tumor
    )), rep("normal", ncol(
      exp_data_normal
    ))))
    design <- model.matrix(~ 0 + group)
    colnames(design) <- levels(group)
    exp_combined <- cbind(exp_data_tumor, exp_data_normal)
    fit <- lmFit(exp_combined, design)
    contrast.matrix <- makeContrasts(tumor - normal, levels = design)
    fit <- contrasts.fit(fit, contrast.matrix)
    fit <- eBayes(fit)

    deg <- topTable(fit, number = Inf, adjust.method = "BH")
    deg_sig <- deg[deg$adj.P.Val < 0.05 & abs(deg$logFC) > 1, ]
    degs_coh <- rownames(deg_sig)

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

    # degs_coh <- intersect(degs_coh,rownames(deg_mat))
    if (use_cpdegs) {
      # thre <- apply(deg_mat, 2, function(x)
      #   max(quantile(abs(x), probs = 0.95), 2))
      thre <- apply(deg_mat,2,function(x) {
        x_abs <- abs(x)
        x_sorted <- sort(x_abs, decreasing = TRUE)
        cutoff <- if(length(x_sorted)>=max_degs) x_sorted[max_degs] else min(x_sorted)
        max(cutoff,2)
      })
      diff_idx <- sweep(abs(deg_mat), 2, thre, FUN = ">=") + 0
      diff_idx <- Matrix(diff_idx, sparse = TRUE)
    } else {
      deg_mat1 <- deg_mat[degs_coh,]
      thre <- apply(deg_mat1, 2, function(x) {
        x_abs <- abs(x)
        x_sorted <- sort(x_abs, decreasing = TRUE)
        cutoff <- if(length(x_sorted)>=max_degs) x_sorted[max_degs] else min(x_sorted)
        max(cutoff,2)
      })
      diff_idx <- sweep(abs(deg_mat1), 2, thre, FUN = ">=") + 0
      diff_idx <- Matrix(diff_idx, sparse = TRUE)
    }
    degs <- list()
    for (i in 1:length(com_samples)) {
      sample_i <- com_samples[i]
      degs[[substr(sample_i, 1, 12)]] <- rownames(diff_idx)[which(diff_idx[, sample_i] ==
                                                                        1)]
    }
  }
  logMessage("Calculating gene coexpression network...")
  exp_cor_tumor <- WGCNA::cor(t(exp_data_tumor), method = "pearson", nThreads = 32)
  exp_cor_tumor[is.na(exp_cor_tumor)] <- 0
  uptri_idx <- which(upper.tri(exp_cor_tumor), arr.ind = TRUE)

  if (ncol(exp_data_normal) > 1) {
    exp_cor_normal <- WGCNA::cor(t(exp_data_normal), method = "pearson",nThreads = 32)
    exp_cor_normal[is.na(exp_cor_normal)] <- 0
    exp_cor_union <- pmax(abs(exp_cor_tumor),abs(exp_cor_normal))
  } else {
    exp_cor_union <- abs(exp_cor_tumor)
  }
  
  exp_cor_upper <- matrix(0,nrow = nrow(exp_cor_tumor),ncol = ncol(exp_cor_tumor),dimnames = list(rownames(exp_cor_tumor),colnames(exp_cor_tumor)))
  exp_cor_upper[uptri_idx] <- exp_cor_union[uptri_idx]
  exp_cor_sparse <- Matrix(exp_cor_upper,sparse = TRUE)
  rm(exp_cor_union)
  gc()

  exp_data_tumor <- exp_data_tumor[, com_samples]
  com_samples <- substr(com_samples, 1, 12)
  colnames(mut_data) <- com_samples
  colnames(exp_data_tumor) <- com_samples
  if (!is.null(deg_mat)) {
    colnames(deg_mat) <- substr(colnames(deg_mat), 1, 12)
    deg_mat <- deg_mat[, com_samples]
  }
  logMessage("Calculating aberrant expression matrix...")
  exp_data_tumor_scaled <-
    scale(t(exp_data_tumor), center = T, scale = T)
  exp_data_tumor_scaled[is.na(exp_data_tumor_scaled)] <- 0
  z_thres <- 2
  aeg_idx <- (abs(exp_data_tumor_scaled) > z_thres) + 0
  aeg_idx <- Matrix(t(aeg_idx), sparse = TRUE)
  aeg_mat <- t(exp_data_tumor_scaled)
  
  mut_genes <- list()
  aegs <- list()
  for (i in 1:length(com_samples)) {
    sample_i <- com_samples[i]
    mut_genes_i <- rownames(mut_data)[which(mut_data[, sample_i] == 1)]
    mut_genes[[sample_i]] <- mut_genes_i
    if (AEG) {
      aegs_i <- rownames(aeg_idx)[which(aeg_idx[, sample_i] == 1)]
      aegs[[sample_i]] <- aegs_i
    }
  }
  
  invalid_idx_m <- which(rowSums(mut_data) == 0)
  if (length(invalid_idx_m) > 0) {
    mut_data <- mut_data[-invalid_idx_m, ]
  }
  return(
    list(
      mut_mat = mut_data,
      cen_mat = exp_cor_sparse,
      deg_mat = deg_mat,
      aeg_mat  = aeg_mat,
      mut_genes = mut_genes,
      degs_coh = degs_coh,
      degs = degs,
      aegs = aegs
    )
  )
}

# initial random walk for each sample
randomwalk_persample <- function(nodes,
                                 PPI_network,
                                 rGRN,
                                 DEN_mat,
                                 ifweighted_dist,
                                 max_dist,
                                 ME_centrality,
                                 gene_length_factors,
                                 ME_as_prior,
                                 DEG_as_prior,
                                 scale_genelength,
                                 cancer_type,
                                 isolated_restart,
                                 weighted,
                                 absorbed,
                                 lambda1,
                                 lambda2,
                                 max_degs,
                                 ifparallel,
                                 num_cores,
                                 use_cpdegs) {
  cat("Starting initial random walks for samples in",
      cancer_type,
      "\n")
  mut_mat <- nodes$mut_mat
  deg_mat <- nodes$deg_mat
  aeg_mat <- nodes$aeg_mat
  # nodes_mat <- nodes$nodes_mat
  mut_genes <- nodes$mut_genes
  degs_coh <- nodes$degs_coh
  degs <- nodes$degs
  aegs <- nodes$aegs
  # samples <- colnames(nodes_mat)
  samples <- colnames(mut_mat)
  P0 <- P1 <- P2 <- matrix(0, nrow(mut_mat), ncol(mut_mat))
  colnames(P0) <- colnames(P1) <- colnames(P2) <- samples
  rownames(P0) <- rownames(P1) <- rownames(P2) <- rownames(mut_mat)

  # rw_not_converged <- c()
  # start_time <- Sys.time()
  # exp_genes <- intersect(rownames(deg_mat),V(rGRN)$name)
  # rGRN_exp <- induced_subgraph(rGRN,exp_genes)
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
    results <- foreach(i = 1:length(samples), .packages = c("igraph", "Matrix"), .export = c("randomwalk_withrestart","get_sample_P","RW_on_PPI_and_rGRN","RW_on_DEN","normalize_rows","quantile_normalize"),.options.snow = opts) %dopar% {
      sample_i <- samples[i]
      mut_genes_i <- mut_genes[[sample_i]]
      if (identical(names(degs),samples)) {
        degs_i <- degs[[sample_i]]
      } else {
        degs_i <- degs
      }
      aegs_i <- aegs[[sample_i]]
      genes_i <- Reduce(union, list(mut_genes_i, degs_i, degs_coh))
      DEN_mat_i <- DEN_mat[intersect(rownames(DEN_mat),genes_i),intersect(rownames(DEN_mat),genes_i)]
      P_i1 <- RW_on_PPI_and_rGRN(mut_genes_i,aegs_i,PPI_network,rGRN,ME_as_prior,ME_centrality,aeg_mat[,sample_i],scale_genelength,gene_length_factors,isolated_restart,lambda1)
      # res_PPI <- RW_on_PPI(mut_genes_i, aegs_i, PPI_network, ME_as_prior, ME_centrality, scale_genelength, gene_length_factors, isolated_restart, alpha)
      # P_i2 <- RW_on_GRN_dist(mut_genes_i, degs_i, rGRN_dist, ifweighted_dist, DEG_as_prior, ME_as_prior=TRUE, ME_centrality, deg_mat[,sample_i], gene_length_factors, isolated_restart, alpha)
      # P_i2 <- RW_on_GRN_whole(mut_genes_i, degs_i, rGRN, DEG_as_prior, ME_as_prior = FALSE, ME_centrality, deg_mat[,sample_i], gene_length_factors, isolated_restart, alpha)
      # res_GRN <- RW_on_GRN_local(mut_genes_i, degs_i, rGRN, DEG_as_prior, ME_as_prior = FALSE, ME_centrality, deg_mat[, sample_i], gene_length_factors, max_dist, isolated_restart, alpha, max_degs)
      # P_i1 <- res_PPI$P_i
      # P_i2 <- res_GRN$P_i
      # P_i <- lambda * P_i1 + (1-lambda) * P_i2
      # vg_num <- length(res_GRN$vg_idx)
      P_i2 <- RW_on_DEN(mut_genes_i,degs_i,degs_coh,DEN_mat_i,DEG_as_prior,deg_mat[,sample_i],weighted,absorbed,use_cpdegs)

      qn <- quantile_normalize(P_i1,P_i2)
      P_i <- lambda2 * qn[[1]] + (1-lambda2) * qn[[2]]
      P_i <- P_i / sum(P_i)
      return(list(P_i1 = P_i1, P_i2 = P_i2, P_i = P_i))
    }

    close(pb)
    # stopCluster(cl)
    # P0 <- do.call(cbind, results)
    # colnames(P0) <- samples
    # rownames(P0) <- rownames(mut_mat)
    for (i in 1:length(samples)) {
      result <- results[[i]]
      # if (!is.na(result$D_not_converged)) {
      #   rw_not_converged[samples[i]] <- result$D_not_converged
      # }
      # P0[names(result$P0_i), samples[i]] <- result$P0_i
      # P0[names(result), samples[i]] <- result
      P0[names(result$P_i), samples[i]] <- result$P_i
      P1[names(result$P_i1), samples[i]] <- result$P_i1
      P2[names(result$P_i2), samples[i]] <- result$P_i2
      # vg_nums[samples[i]] <- result$vg_num
    }
    # end_time <- Sys.time()
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
      if (identical(names(degs),samples)) {
        degs_i <- degs[[sample_i]]
      } else {
        degs_i <- degs
      }
      aegs_i <- aegs[[sample_i]]
      genes_i <- Reduce(union, list(mut_genes_i, degs_i, degs_coh))
      DEN_mat_i <- DEN_mat[intersect(rownames(DEN_mat),genes_i),intersect(rownames(DEN_mat),genes_i)]
      P_i1 <- RW_on_PPI_and_rGRN(mut_genes_i,aegs_i,PPI_network,rGRN,ME_as_prior,ME_centrality,aeg_mat[,sample_i],scale_genelength,gene_length_factors,isolated_restart,lambda1)
      P_i2 <- RW_on_DEN(mut_genes_i,degs_i,degs_coh,DEN_mat_i,DEG_as_prior,deg_mat[,sample_i],weighted,absorbed,use_cpdegs)

      qn <- quantile_normalize(P_i1,P_i2)
      P_i <- lambda2 * qn[[1]] + (1-lambda2) * qn[[2]]
      P_i <- P_i / sum(P_i)
      P0[names(P_i),sample_i] <- P_i
      P1[names(P_i1),sample_i] <- P_i1
      P2[names(P_i2),sample_i] <- P_i2
      pb$tick()
    }
  }
  nodes$P1 <- P1
  nodes$P2 <- P2
  nodes$P0_mat <- P0
  return(nodes)
}


# initial transition matrix for each sample
get_sample_P <- function(genes, network) {
  subnetwork <-
    induced_subgraph(network, intersect(genes, V(network)$name))
  A <- matrix(
    0,
    nrow = length(genes),
    ncol = length(genes),
    dimnames = list(genes, genes)
  )
  adjacency_mat <- as.matrix(as_adjacency_matrix(subnetwork))
  A[rownames(adjacency_mat), colnames(adjacency_mat)] <-
    adjacency_mat
  A <- as(A, "sparseMatrix")
  
  P <- normalize_rows(A)
  return(P)
}

# Perform random walk with restart on the PPI subnetwork for a given sample.
# The subnetwork is induced by the union of mutated genes and abnormally expressed genes.
# The restart (prior) vector can be based on the centrality of mutated genes in the ME network,
# and can be further adjusted by gene length factors.
# For isolated nodes (nodes with no edges), the function supports either forced restart
# (assigning the restart vector to these nodes) or random transition to other nodes.
# The function returns the stationary distribution (steady-state probability) for the mutated genes.
RW_on_PPI <- function(mut_genes_i,
                      aegs_i,
                      PPI_network,
                      ME_as_prior,
                      ME_centrality,
                      scale_genelength,
                      gene_length_factors,
                      isolated_restart,
                      alpha) {
  # Induce subnetwork from PPI using mutated and abnormally expressed genes
  genes <- union(mut_genes_i, aegs_i)
  subnetwork <- induced_subgraph(PPI_network, intersect(genes, V(PPI_network)$name))
  A <- Matrix(
    0,
    nrow = length(genes),
    ncol = length(genes),
    dimnames = list(genes, genes),
    sparse = TRUE
  )
  adj_mat <- as_adjacency_matrix(subnetwork, sparse = TRUE)
  nz_idx <- summary(adj_mat)
  row_idx <- match(rownames(adj_mat)[nz_idx$i], rownames(A))
  col_idx <- match(colnames(adj_mat)[nz_idx$j], colnames(A))
  A[cbind(row_idx, col_idx)] <- nz_idx$x
  # A[rownames(adj_mat), colnames(adj_mat)] <- adj_mat
  P <- normalize_rows(A)
  
  # Construct the restart (prior) vector v0
  v0 <- numeric(length(genes))
  names(v0) <- genes
  if (ME_as_prior &&
      length(intersect(names(ME_centrality), mut_genes_i)) >= 1) {
    idx <- intersect(names(ME_centrality), mut_genes_i)
    v0[idx] <- ME_centrality[idx]
    # v0[setdiff(mut_genes_i, idx)] <- median(ME_centrality)
    # v0[setdiff(mut_genes_i, idx)] <- min(ME_centrality)
  } else {
    v0[mut_genes_i] <- 1 / length(mut_genes_i)
  }
  # Optionally scale the prior by gene length factors
  if (scale_genelength) {
    v0[mut_genes_i] <- v0[mut_genes_i] * gene_length_factors[mut_genes_i]
  }
  v0 <- v0 / sum(v0)
  # Handle isolated nodes: either force restart or use random transition
  if (isolated_restart) {
    zero_rows <- which(rowSums(P) == 0)
    if (length(zero_rows) > 0) {
      row_idx <- rep(zero_rows, each = ncol(P))
      col_idx <- rep(1:ncol(P), length(zero_rows))
      idx_matrix <- cbind(row_idx, col_idx)
      P[idx_matrix] <- rep(v0, length(zero_rows))
    }
    # P[rowSums(P) == 0, ] <- matrix(rep(v0, sum(rowSums(P) == 0)), ncol = ncol(P), byrow =
    #                                  TRUE)
  } else {
    P <- (1 - alpha) * P + alpha / ncol(P)
    P <- normalize_rows(P)
  }
  # Perform random walk with restart
  rw <- randomwalk_withrestart(P, mut_genes_i, v0 = v0)
  P_i <- rw$vt[mut_genes_i]
  P_i <- P_i / sum(P_i)
  # return(P_i)
  # vg_num <- sum(degree(subnetwork)[intersect(mut_genes_i,V(subnetwork)$name)]>0)
  degree_mut <- numeric(length(mut_genes_i))
  names(degree_mut) <- mut_genes_i
  degree_mut[V(subnetwork)$name] <- degree(subnetwork)
  # if (sum(degree_mut)==0){
  #   degree_mut <- degree_mut+1
  # }
  # kl <- kl_divergence(P_i,degree_mut/sum(degree_mut))
  vg_idx <- names(degree_mut)[which(degree_mut > 0)]
  comps <- components(subnetwork, mode = "weak")
  p <- tapply(P_i, comps$membership[names(P_i)],sum)
  entropy <- -sum(p * log(p + 1e-10))
  top_degree <- mean(sort(degree(subnetwork),decreasing = TRUE)[1:min(10,vcount(subnetwork))])
  # big_comp_ids <- which(comps$csize >= 3)
  # total_nodes_big <- sum(comps$csize[big_comp_ids])
  # vg_prop <- total_nodes_big / max(1, vcount(subnetwork))
  # vg_idx <- names(degree(subnetwork))[which(degree(subnetwork) > 0)]
  return(list(prior = v0, P_i = P_i, vg_idx = vg_idx, entropy = entropy, top_degree = top_degree, degree = degree_mut))
}

RW_on_PPI_and_rGRN <- function(mut_genes_i,
                               aegs_i,
                               PPI_network,
                               rGRN,
                               ME_as_prior,
                               ME_centrality,
                               z_tumor,
                               scale_genelength,
                               gene_length_factors,
                               isolated_restart,
                               lambda) {
  genes <- union(mut_genes_i, aegs_i)
  P_ppi <- get_sample_P(genes, PPI_network)
  P_rgrn <- get_sample_P(genes, rGRN)
  P <- lambda * P_ppi + (1 - lambda) * P_rgrn
  P <- normalize_rows(P)
  v0 <- numeric(length(genes))
  names(v0) <- genes
  if (ME_as_prior &&
      length(intersect(names(ME_centrality), mut_genes_i)) >= 1) {
    idx <- intersect(names(ME_centrality), mut_genes_i)
    v0[idx] <- ME_centrality[idx]
    v0[setdiff(mut_genes_i, idx)] <- median(ME_centrality)
  } else if (length(intersect(mut_genes_i,names(z_tumor))) >=1){
    # v0[mut_genes_i] <- 1 / length(mut_genes_i)
    v0[intersect(mut_genes_i,names(z_tumor))] <- abs(z_tumor[intersect(mut_genes_i,names(z_tumor))])
  } else {
    v0[mut_genes_i] <- 1 / length(mut_genes_i)
  }
  if (scale_genelength) {
    v0[mut_genes_i] <- v0[mut_genes_i] * gene_length_factors[mut_genes_i]
  }
  v0 <- v0 / sum(v0)
  if (isolated_restart) {
    zero_rows <- which(rowSums(P) == 0)
    if (length(zero_rows) > 0) {
      row_idx <- rep(zero_rows, each = ncol(P))
      col_idx <- rep(1:ncol(P), length(zero_rows))
      idx_matrix <- cbind(row_idx, col_idx)
      P[idx_matrix] <- rep(v0, length(zero_rows))
    }
  }
  rw <- randomwalk_withrestart(P, mut_genes_i, v0 = v0)
  P_i <- rw$vt[mut_genes_i]
  P_i <- P_i / sum(P_i)
  # subPPI <- induced_subgraph(PPI_network,intersect(V(PPI_network)$name,genes))
  # comps <- components(subPPI, mode = "weak")
  # p <- tapply(P_i, comps$membership[names(P_i)],sum)
  # entropy <- -sum(p * log(p + 1e-10))
  # return(list(P_i=P_i, entropy=entropy))
  return(P_i)
}

# Perform a random walk with restart on a local subgraph of the reversed GRN for a given sample.
# The subgraph is not simply induced by mutated and differentially expressed genes,
# but also includes all nodes on directed paths (in the reversed GRN) between these genes
# where the path length does not exceed a specified threshold (max_dist).
# The restart (prior) vector v0 is based on the differential expression scores of DEGs.
# For isolated nodes, the function supports either forced restart to mutated genes
# or random transition to any other gene.
# The function returns the stationary distribution (steady-state probability) for the mutated genes.
RW_on_GRN_local <- function(mut_genes_i,
                            degs_i,
                            rGRN,
                            DEG_as_prior,
                            ME_as_prior,
                            ME_centrality,
                            deg_score_i,
                            gene_length_factors,
                            max_dist,
                            isolated_restart,
                            alpha,
                            max_degs) {
  # Find all genes on short paths between mutated and DEG genes
  # d <- distances(rGRN, to = intersect(V(rGRN)$name, mut_genes_i), mode = "out")
  # d_min <- apply(d, 1, min)
  # ugs <- names(d_min)[which(d_min <= max_dist)]
  # ugs <- intersect(ugs,names(deg_score_i))
  # ugs_zscore <- deg_score_i[ugs]
  # degs_i <- names(sort(abs(ugs_zscore),decreasing=TRUE)[1:min(length(ugs_zscore),1000)])
  # degs_i <- names(ugs_zscore)[which(abs(ugs_zscore)>2)]
  genes <- union(mut_genes_i, degs_i)
  com_genes <- intersect(genes, V(rGRN)$name)
  # D1 <- distances(rGRN, v = intersect(V(rGRN)$name,degs_i), mode = "out")
  # D2 <- distances(rGRN, to = intersect(V(rGRN)$name,mut_genes_i), mode = "out")
  D1 <- distances(rGRN, v = com_genes, mode = "out")
  D2 <- distances(rGRN, to = com_genes, mode = "out")
  D1_min <- apply(D1, 2, min)
  D2_min <- apply(D2, 1, min)
  V <- names(D1_min)[which(D1_min + D2_min <= max_dist)]
  genes_total <- union(V, genes)
  # D <- distances(rGRN, v = setdiff(V(rGRN)$name,mut_genes_i), to = intersect(V(rGRN)$name,mut_genes_i), mode = "out")
  # D_min <- apply(D,1,min)
  # D_min_mut <- apply(D,2,min)
  # V <- names(D_min)[which(D_min <= max_dist)]
  # vg_idx <- names(D_min_mut)[which(D_min_mut < Inf)]
  # genes_total <- union(V,mut_genes_i)
  if (length(V) == 0){
    return(list(P_i = numeric(0), vg_idx = numeric(0)))
  }
  # Build the local subgraph and adjacency matrix
  subnet <- induced_subgraph(rGRN, V)
  d <- distances(subnet, v = intersect(V(subnet)$name,degs_i),to=intersect(V(subnet)$name,mut_genes_i),mode="out")
  d_min <- apply(d, 2, min)
  # d_min_deg <- apply(d, 1, min)
  # degs_i <- names(d_min_deg)[which(d_min_deg < Inf)]
  # vg_num <- sum(d_min < Inf)
  # vg_idx <- numeric(length(mut_genes_i))
  # names(vg_idx) <- mut_genes_i
  # vg_idx[names(d_min)] <- (d_min < Inf) + 0
  vg_idx <- names(d_min)[which(d_min < Inf)] # vg_idx <- colnames(d)

  adj_mat <- as_adjacency_matrix(subnet, sparse = TRUE) 
  A <- Matrix(0, nrow = length(genes_total), ncol = length(genes_total), dimnames = list(genes_total, genes_total), sparse = TRUE)
  A[rownames(adj_mat), colnames(adj_mat)] <- adj_mat
  mut_idx <- which(rownames(A) %in% mut_genes_i)
  A[mut_idx,]<-0
  A[cbind(mut_idx,mut_idx)] <- 1
  # zero_rows <- which(rowSums(A) == 0)
  # A[cbind(zero_rows, zero_rows)] <- 1
  P <- normalize_rows(A)
  v0 <- v1 <- numeric(length(genes_total))
  names(v0) <- names(v1) <- genes_total
  if (length(intersect(names(ME_centrality),mut_genes_i)) > 0) {
    idx <- intersect(names(ME_centrality),mut_genes_i)
    v1[idx] <- ME_centrality[idx]
    v1[setdiff(mut_genes_i,idx)] <- median(ME_centrality)
    v1[mut_genes_i] <- v1[mut_genes_i] * gene_length_factors[mut_genes_i]
  } else {
    # v1[mut_genes_i] <- gene_length_factors[mut_genes_i]
    v1[mut_genes_i] <- 1 / length(mut_genes_i)
  }
  if (DEG_as_prior && length(intersect(genes,names(deg_score_i))) >= 1) {
    # thre <- max(sort(abs(deg_score_i),decreasing = TRUE)[min(length(deg_score_i),max_degs)],2)
    # thre <- max(quantile(abs(deg_score_i),0.9),2)
    # idx <- abs(deg_score_i) >= thre
    # v0[degs_i] <- abs(deg_score_i[degs_i]) * idx
    # idx <- which(abs(deg_score_i) >= thre)
    # v0[degs_i[idx]] <- 1/length(idx)
    # v0[degs_i] <- abs(deg_score_i[degs_i])
    v0[intersect(genes,names(deg_score_i))] <- abs(deg_score_i[intersect(genes,names(deg_score_i))])
    # v0[intersect(genes_total,names(deg_score_i))] <- abs(deg_score_i[intersect(genes_total,names(deg_score_i))])
    # v0[setdiff(mut_genes_i,vg_idx)] <- 0
    # v0[setdiff(genes_total,names(deg_score_i))] <- median(abs(deg_score_i))
    # v0[degs_i] <- 1/length(degs_i)
    # v0[intersect(genes_total,names(deg_score_i))] <- abs(deg_score_i[intersect(genes_total,names(deg_score_i))])
    # v0[intersect(genes,names(deg_score_i))] <- abs(deg_score_i[intersect(genes,names(deg_score_i))])
  } else if (ME_as_prior) {
    v0 <- v1
  } else {
    v0[mut_genes_i] <- gene_length_factors[mut_genes_i]
  }
  v0 <- v0 / sum(v0)
  # Handle isolated nodes: forced restart or random transition
  # mut_idx <- which(rownames(P) %in% mut_genes_i)
  # P[cbind(mut_idx,mut_idx)] <- 1
  if (isolated_restart) {
    zero_rows <- which(rowSums(P) == 0)
    if (length(zero_rows) > 0) {
      row_idx <- rep(zero_rows, each = ncol(P))
      col_idx <- rep(1:ncol(P), length(zero_rows))
      idx_matrix <- cbind(row_idx, col_idx)
      P[idx_matrix] <- rep(v1, length(zero_rows))
    }
    # P[rowSums(P) == 0, mut_genes_i] <- 1 / length(mut_genes_i)
  } else {
    P <- (1 - alpha) * P + alpha / ncol(P)
  }
  P <- normalize_rows(P)
  # Perform random walk with restart
  rw <- randomwalk_withrestart(P, mut_genes_i, v0 = v0,theta = 1)
  # rw <- randomwalk_withrestart(P,mut_genes_i,v0=v0)
  P_i <- rw$vt[mut_genes_i]
  if (sum(P_i) == 0) {
    P_i <- v1[mut_genes_i]
  }
  # P_i[P_i==0] <- median(P_i[P_i!=0])
  P_i <- P_i / sum(P_i)
  # n <- length(P_i)
  # k <- 10
  # tail_k <- max(1, floor(n * 0.3))
  # sorted_P_i <- sort(P_i,decreasing = TRUE)
  # HTR <- mean(sorted_P_i[1:min(k,n)]) / mean(sorted_P_i[(n - tail_k + 1):n])
  degree_mut <- numeric(length(mut_genes_i))
  names(degree_mut) <- mut_genes_i
  degree_mut[intersect(mut_genes_i,V(subnet)$name)] <- degree(subnet,mode = "in")[intersect(mut_genes_i,V(subnet)$name)]
  # return(P_i)
  return(list(P_i = P_i, vg_idx = vg_idx, degree = degree_mut))
}

RW_on_DEN <- function(mut_genes_i,degs_i,degs_coh,DEN_mat_i,DEG_as_prior,deg_score_i,weighted,absorbed,use_cpdegs){
  if (use_cpdegs){
    genes <- Reduce(union, list(mut_genes_i, degs_i, degs_coh))
  } else {
    genes <- union(mut_genes_i,degs_i)
  }
  DEN_mat_i <- DEN_mat_i[intersect(genes,rownames(DEN_mat_i)),intersect(genes,rownames(DEN_mat_i))]
  v0 <- numeric(length(genes))
  names(v0) <- genes
  if (DEG_as_prior && length(intersect(union(mut_genes_i,degs_i), names(deg_score_i))) > 0){
    v0[intersect(union(mut_genes_i,degs_i),names(deg_score_i))] <- abs(deg_score_i[intersect(union(mut_genes_i,degs_i),names(deg_score_i))])
  } else {
    v0[union(mut_genes_i,degs_i)] <- 1
  }
  v0 <- v0/sum(v0)
  # thresholds <- c(0.5,0.4,0.3)
  # for (t in thresholds){
  #   mat <- DEN_mat_i
  #   mat@x[mat@x < t] <- 0
  #   mat <- drop0(mat)
  #   nz_idx <- summary(mat)
  #   nz_genes <- unique(union(rownames(mat)[nz_idx$i],colnames(mat)[nz_idx$j]))
  #   nz_mutgenes <- intersect(nz_genes,mut_genes_i)
  #   if (length(nz_mutgenes)/length(mut_genes_i) >= 0.8){
  #     break
  #   }
  # }
  # W <- as(mat, "generalMatrix") + t(as(mat, "generalMatrix"))
  W <- as(DEN_mat_i, "generalMatrix") + t(as(DEN_mat_i, "generalMatrix"))
  A <- (W>=0.5) * 1
  P <- Matrix(0,nrow = length(genes),ncol = length(genes),dimnames = list(genes,genes),sparse=TRUE)

  if (weighted){
    mat <- W
  } else {
    mat <- A
  }
  nzero_idx <- summary(mat)
  map <- match(rownames(mat),genes)
  new_i <- map[nzero_idx$i]
  new_j <- map[nzero_idx$j]
  # P <- sparseMatrix(i = new_i, j = new_j, x = nzero_idx$x, dims = c(length(genes),length(genes)), dimnames = list(genes,genes))
  P <- sparseMatrix(i = new_i, j = new_j, x = 1, dims = c(length(genes),length(genes)), dimnames = list(genes,genes))
  if (absorbed){
    mut_idx <- which(rownames(P) %in% mut_genes_i)
    P[mut_idx,] <- 0
    P[cbind(mut_idx,mut_idx)] <- 1
    theta <- 1
  } else {
    theta <- 0.85
  }
  if (any(rowSums(P)==0)){
    zero_idx <- which(rowSums(P)==0)
    # P[cbind(zero_idx,zero_idx)] <- 1
    row_idx <- rep(zero_idx, each = ncol(P))
    col_idx <- rep(1:ncol(P), length(zero_idx))
    idx_matrix <- cbind(row_idx, col_idx)
    P[idx_matrix] <- rep(v0, length(zero_idx))
  }
  P <- normalize_rows(P)
  rw <- randomwalk_withrestart(P, mut_genes_i, v0 = v0, theta=theta)
  P_i <- rw$vt[mut_genes_i]
  if (sum(P_i) == 0){
    if (sum(v0[mut_genes_i])!=0){
      P_i <- v0[mut_genes_i]
    } else {
      P_i <- rep(1,length(mut_genes_i))
      names(P_i) <- mut_genes_i
    }
  }
  P_i <- P_i / sum(P_i)
  return(P_i)
}

# compute stationary distribution of random walk with restart
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

# hypergraph random walks for each sample
prediction_sample_specific <- function(nodes,
                                       co_mut,
                                       sim_mat,
                                       gene_length_factors,
                                       scale_genelength,
                                       no_selftrans,
                                       cancer_type,
                                       delta,
                                       theta,
                                       ifparallel,
                                       num_cores) {
  cat("Starting patient-specific hypergraph random walks for samples in",
      cancer_type,
      "\n")
  samples <- colnames(co_mut)
  # nodes_mat <- as.matrix(nodes$nodes_mat)
  mut_mat <- nodes$mut_mat
  mut_genes <- nodes$mut_genes
  # degs <- nodes$degs
  # aegs <- nodes$aegs
  P0_mat <- nodes$P0_mat
  # P1 <- nodes$P1
  # P2 <- nodes$P2
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
      E_i <-
        colnames(co_mut)[which(co_mut[sample_i, ] >= 2)]# hyperedge
      # sim_mat_i <- sim_mat[sample_i, E_i, drop = F]
      # H_i <- mut_mat[V_i, E_i, drop = F]
      # W_ve_i <- P0_mat[V_i, E_i, drop = F]
      sim_mat_i <- as.matrix(sim_mat[sample_i, E_i, drop = FALSE])
      H_i <- as.matrix(mut_mat[V_i, E_i, drop = FALSE])
      W_ve_i <- as.matrix(P0_mat[V_i, E_i, drop = FALSE])
      W_e_i <- get_hyperedge_weight(sim_mat_i, delta)
      if (no_selftrans) {
        P_i <- get_hyper_P(H_i, W_e_i, W_ve_i)
        diag(P_i) <- 0
        P_i <- normalize_rows(P_i)
      } else{
        P_i <- get_hyper_P(H_i, W_e_i, W_ve_i)
      }
      v0 <- numeric(length(V_i))
      names(v0) <- V_i
      v0[mut_genes_i] <- 1 / length(mut_genes_i)
      if (scale_genelength) {
        v0[mut_genes_i] <- v0[mut_genes_i] * gene_length_factors[mut_genes_i]
      }
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
    # end_time <- Sys.time()
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
      E_i <-
        colnames(co_mut)[which(co_mut[sample_i, ] >= 2)]# hyperedge
      sim_mat_i <- sim_mat[sample_i, E_i, drop = F]
      V_i <- mut_genes_i
      H_i <- mut_mat[V_i, E_i, drop = F]
      W_e_i <- get_hyperedge_weight(sim_mat_i, delta)
      W_ve_i <- P0_mat[V_i, E_i, drop = F]
      if (no_selftrans) {
        P_i <- get_hyper_P(H_i, W_e_i, W_ve_i)
        diag(P_i) <- 0
        P_i <- normalize_rows(P_i)
      } else {
        P_i <- get_hyper_P(H_i, W_e_i, W_ve_i)
      }
      v0 <- rep(1, nrow(P_i))
      names(v0) <- rownames(P_i)
      if (scale_genelength) {
        v0[mut_genes_i] <- v0[mut_genes_i] * gene_length_factors[mut_genes_i]
      }
      v0 <- v0 / sum(v0)
      hrw_i <- randomwalk_withrestart(P_i, mut_genes_i, theta = theta, v0 = v0)
      prediction_i <- hrw_i$vt
      if (!is.na(hrw_i$D_not_converged)) {
        hrw_not_converged[sample_i] <- hrw_i$D_not_converged
      }
      prediction_all[names(prediction_i), sample_i] <- prediction_i
      mut_genes_scores_i <- prediction_i[mut_genes_i]
      mut_genes_scores_i <-
        mut_genes_scores_i / sum(mut_genes_scores_i)
      mut_genes_scores_i <-
        sort(mut_genes_scores_i, decreasing = TRUE)
      mut_genes_scores[[sample_i]] <- mut_genes_scores_i
      mut_genes_ranks[[sample_i]] <- names(mut_genes_scores_i)
      
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
  
  nodes$prediction_mat <- prediction_all
  nodes$mut_genes_scores <- mut_genes_scores
  nodes$mut_genes_ranks <- mut_genes_ranks
  return(nodes)
}

# compute hyperedge weights
get_hyperedge_weight <- function(sim_mat_i, delta = 0.1) {
  W_e_i <- exp(-(1 - sim_mat_i)^2 / (2 * delta^2))
  W_e_i <- diag(c(W_e_i))
  colnames(W_e_i) <- rownames(W_e_i) <- colnames(sim_mat_i)
  return(W_e_i)
}

# compute hypergraph transition matrix
get_hyper_P <- function(H, W_e, W_ve) {
  H_w <- H %*% W_e
  D_H_w <- normalize_rows(H_w)
  D_W_ve <- normalize_rows(t(W_ve))
  P <- D_H_w %*% D_W_ve
  return(P)
}

# 对行和非零的行进行归一化
normalize_rows <- function(mat) {
  row_sums <- rowSums(mat)
  row_sums[row_sums == 0] <- 1
  return(mat / row_sums)
}

hypergraph_rw_cohort_singlecancer <- function(nodes,
                                              scale_genelength,
                                              gene_length_factors,
                                              no_selftrans,
                                              theta,
                                              num_cores) {
  H <- nodes$mut_mat
  V <- rownames(H)
  E <- colnames(H)
  W_e <- diag(ncol(H))
  W_ve <- nodes$P0_mat
  if (no_selftrans) {
    P <- get_hyper_P(H, W_e, W_ve)
    diag(P) <- 0
    P <- normalize_rows(P)
  } else {
    P <- get_hyper_P(H, W_e, W_ve)
  }

  v0 <- numeric(nrow(P))
  names(v0) <- rownames(P)
  v0[V] <- 1 / length(V)
  if (scale_genelength) {
    v0[V] <- v0[V] * gene_length_factors[V]
  }
  v0 <- v0 / sum(v0)
  rw <- randomwalk_withrestart(P, V, theta = theta, v0 = v0)
  scores <- rw$vt
  scores <- scores / sum(scores)
  scores <- sort(scores, decreasing = TRUE)
  scores <- as.matrix(scores)
  return(scores)
}

get_length_factor <- function(gene_length, l0, k) {
  if (gene_length <= l0) {
    return(1)
  } else {
    return(exp(-(gene_length - l0) / (k^2)))
    # return(exp(-(gene_length-l0)^(1/k)))
  }
}

quantile_normalize <- function(x1,x2){
  # 对输入的向量进行分位数归一化
  if (length(x1) != length(x2)) {
    stop("The lengths of x1 and x2 must be equal.")
  }
  # 将两个向量合并为一个矩阵
  combined <- cbind(x1, x2)
  # 对每一列进行排序
  sorted_combined <- apply(combined, 2, sort)
  # 计算每一列的分位数
  quantiles <- apply(sorted_combined, 1, mean)
  # 将原始数据替换为分位数
  # normalized_x1 <- quantiles[rank(x1, ties.method = "first")]
  # normalized_x2 <- quantiles[rank(x2, ties.method = "first")]
  normalized_x1 <- approx(1:length(quantiles), quantiles, xout = rank(x1, ties.method = "average"))$y
  normalized_x2 <- approx(1:length(quantiles), quantiles, xout = rank(x2, ties.method = "average"))$y
  normalized_x1 <- normalized_x1 / sum(normalized_x1)
  normalized_x2 <- normalized_x2 / sum(normalized_x2)
  # 返回归一化后的向量
  names(normalized_x1) <- names(x1)
  names(normalized_x2) <- names(x2)
  return(list(normalized_x1 = normalized_x1, normalized_x2 = normalized_x2))
}

safe_normalize_rank <- function(deg) {
  r <- rank(-deg, ties.method = "max")
  if (max(r) == min(r)) {
    return(rep(0.5, length(r)))
  }
  r_norm <- (r - min(r)) / (max(r) - min(r))
  return(r_norm)
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