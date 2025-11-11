# This script processes mutation and expression data for multiple cancer types
# Author: xqxu
# Date: "2025-03-04"

# Set working directory to the script location
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(Matrix)
#' Process mutation data file
#'
#' @param data_file Path to the mutation data file
#' @return Processed mutation data as a data frame
#' @details Reads mutation data, removes invalid entries, and creates a gene-by-sample matrix
process_mut_data <- function(data_file) {
  message("Processing mutation data from: ", data_file)
  
  # Check if file exists
  if (!file.exists(data_file)) {
    stop("Mutation data file not found: ", data_file)
  }
  
  # Read mutation data
  data <- read.delim(
    data_file,
    header = TRUE,
    as.is = TRUE,
    check.names = FALSE
  )
  
  # Remove invalid entries
  invalid_gene_index <- which(is.na(data[, 2]) |
                                data[, 2] == "" |
                                data[, 2] == ".")
  if (length(invalid_gene_index) > 0) {
    data <- data[-invalid_gene_index, ]
  }
  
  # Keep only gene and sample columns
  data <- data[, c(1:2)]
  
  # Remove duplicated rows
  if (any(duplicated(data))) {
    data <- unique(data)
  }
  
  # Get unique genes and samples
  genes <- unique(data$gene)
  samples <- unique(data$sample)
  
  # Create mapping for sparse matrix construction
  row_map <- setNames(seq_along(genes), genes)
  col_map <- setNames(seq_along(samples), samples)
  row_idx <- row_map[data$gene]
  col_idx <- col_map[data$sample]
  
  # Create sparse matrix
  mut_data <- sparseMatrix(
    i = row_idx,
    j = col_idx,
    x = 1,
    dims = c(length(genes), length(samples)),
    dimnames = list(genes, samples)
  )
  
  # mut_data <- mut_data[, which(colSums(mut_data) >= 3)]
  #
  # # only samples with "-01A" suffix are kept
  # mut_data <- mut_data[, grepl("-01A$", colnames(mut_data))]
  
  mut_data <-
    mut_data[which(rowSums(mut_data) != 0), which(colSums(mut_data) != 0)]
  
  # Convert to data frame
  mut_data <- as.matrix(mut_data)
  mut_data <- as.data.frame(mut_data)
  return(mut_data)
}

#' Process expression data file
#'
#' @param data_file Path to the expression data file
#' @param metadata_file Path to the gene metadata file
#' @return Processed expression data as a data frame
#' @details Reads expression data, maps gene IDs to symbols, and handles duplicated genes
process_exp_data <- function(counts_file, tpm_file, metadata_file) {
  message("Processing expression counts data from: ", counts_file)
  message("Processing expression TPM data from: ", tpm_file)
  # Check if files exist
  if (!file.exists(counts_file)) {
    stop("Expression count data file not found: ", counts_file)
  }
  if (!file.exists(tpm_file)) {
    stop("Expression TPM data file not found: ", tpm_file)
  }
  if (!file.exists(metadata_file)) {
    stop("Gene metadata file not found: ", metadata_file)
  }
  
  # Read expression data and metadata
  exp_counts <- read.delim(
    counts_file,
    header = TRUE,
    as.is = TRUE,
    check.names = FALSE,
    row.names = 1
  )
  exp_tpm <- read.delim(
    tpm_file,
    header = TRUE,
    as.is = TRUE,
    check.names = FALSE,
    row.names = 1
  )
  metadata <- read.delim(
    metadata_file,
    header = TRUE,
    as.is = TRUE,
    check.names = FALSE
  )
  
  # Map gene IDs to gene symbols
  gene_id_map <- metadata[, c("id", "gene")]
  match_idx <- match(rownames(exp_counts), gene_id_map$id)
  exp_tpm <- exp_tpm[rownames(exp_counts), ]
  
  # Different ensembl ids may map to the same gene symbol
  if (any(duplicated(gene_id_map$gene))) {
    row_gene_id <- gene_id_map$gene[match_idx]
    table_results <- table(row_gene_id)
    duplicated_genes <- names(table_results)[table_results > 1]
    for (gene in duplicated_genes) {
      gene_idx <- which(row_gene_id == gene)
      exp_counts[gene_idx[1], ] <-
        colMeans(exp_counts[gene_idx, ], na.rm = T)
      exp_tpm[gene_idx[1], ] <-
        colMeans(exp_tpm[gene_idx, ], na.rm = T)
      exp_counts <- exp_counts[-gene_idx[-1], ]
      exp_tpm <- exp_tpm[-gene_idx[-1], ]
      row_gene_id <- row_gene_id[-gene_idx[-1]]
    }
    rownames(exp_counts) <- rownames(exp_tpm) <- row_gene_id
  } else{
    rownames(exp_counts) <- rownames(exp_tpm) <- gene_id_map$gene[match_idx]
  }
  
  # Remove invalid entries
  invalid_gene_index <- which(is.na(rownames(exp_counts)) |
                                rownames(exp_counts) == "" |
                                rownames(exp_counts) == ".")
  if (length(invalid_gene_index) > 0) {
    exp_counts <- exp_counts[-invalid_gene_index, ]
    exp_tpm <- exp_tpm[-invalid_gene_index, ]
  }
  # sample_id <- colnames(exp_data)
  # colnames(exp_data) <- substr(sample_id, 1, 15)
  # # only samples with "-01A" suffix are kept
  # exp_data <- exp_data[, grepl("-01A$", colnames(exp_data))]
  # exp_tpm <-
  #   as.matrix(exp_tpm[which(rowSums(exp_tpm >= 1) >= ncol(exp_tpm) * 0.2), ])
  # exp_counts <- as.matrix(exp_counts[rownames(exp_tpm), ])
  exp_tpm <- as.matrix(exp_tpm)
  exp_counts <- as.matrix(exp_counts)
  exp_counts <- 2^exp_counts - 1
  return(list(counts = exp_counts, tpm = exp_tpm))
}

get_DEG_data <- function(exp_data, theta = 1) {
  exp_data <- exp_data[, grepl("-01A$", colnames(exp_data))]
  exp_data_scaled <- scale(t(exp_data), center = T, scale = T)
  exp_data_scaled[is.na(exp_data_scaled)] <- 0
  
  # quantile_5 <- apply(exp_data_scaled, 2, function(x)
  #   quantile(x, probs = 0.05))
  # quantile_95 <- apply(exp_data_scaled, 2, function(x)
  #   quantile(x, probs = 0.95))
  # theta <- pmin(abs(quantile_5), abs(quantile_95))
  # deg_idx <- sweep(abs(exp_data_scaled), 2, theta, FUN = ">")
  deg_idx <- abs(exp_data_scaled) > theta
  return(deg_idx)
}

#' Main function to process all cancer types
#'
#' @param cancer_types Vector of cancer type codes to process
#' @param base_dir Base directory for data files
#' @return None, writes processed files to disk
main <- function(cancer_types, base_dir = "./data") {
  # Ensure base directories exist
  raw_data_dir <- file.path(base_dir, "raw")
  if (!dir.exists(raw_data_dir)) {
    stop("Raw data directory not found: ", raw_data_dir)
  }

  gene_metadata_file <- file.path(base_dir, "gene_mapping", "gencode.v36.annotation.gtf.gene.probemap")
  if (!file.exists(gene_metadata_file)) {
    stop("Gene metadata file not found: ", gene_metadata_file)
  }
  
  # Process each cancer type
  for (cancer_type in cancer_types) {
    message("\n=== Processing cancer type: ", cancer_type, " ===")
    
    # Define file paths
    mut_data_file <- file.path(raw_data_dir,
                               paste0("TCGA-", cancer_type, ".somaticmutation_wxs.tsv"))
    exp_counts_file <- file.path(raw_data_dir,
                                 paste0("TCGA-", cancer_type, ".star_counts.tsv"))
    exp_tpm_file <- file.path(raw_data_dir,
                              paste0("TCGA-", cancer_type, ".star_tpm.tsv"))
    output_dir <- file.path(base_dir, "processed", cancer_type)
    
    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
      message("Created output directory: ", output_dir)
    }
    
    # Process mutation data
    tryCatch({
      mut_data <- process_mut_data(mut_data_file)
      mut_output_file <- file.path(output_dir, "mut_data.tsv")
      write.table(
        mut_data,
        file = mut_output_file,
        sep = "\t",
        quote = F,
        row.names = T,
        col.names = T
      )
      message("Mutation data saved to: ", mut_output_file)
    }, error = function(e) {
      warning("Error processing mutation data for ",
              cancer_type,
              ": ",
              e$message)
    })
    
    # Process expression data
    tryCatch({
      exp_data <- process_exp_data(exp_counts_file, exp_tpm_file, gene_metadata_file)
      write.table(
        exp_data$counts,
        file = file.path(output_dir, "exp_counts_data.tsv"),
        sep = "\t",
        quote = F,
        row.names = T,
        col.names = T
      )
      write.table(
        exp_data$tpm,
        file = file.path(output_dir, "exp_tpm_data.tsv"),
        sep = "\t",
        quote = F,
        row.names = T,
        col.names = T
      )
      message(
        "Expression counts data saved to: ",
        file.path(output_dir, "exp_counts_data.tsv")
      )
      message("Expression TPM data saved to: ",
              file.path(output_dir, "exp_tpm_data.tsv"))
      # deg_idx <- get_DEG_data(exp_data$tpm)
      # write.table(
      #   deg_idx,
      #   file = file.path(output_dir, "DEGs.csv"),
      #   quote = TRUE,
      #   sep = ","
      # ) 
    }, error = function(e) {
      warning("Error processing expression data for ",
              cancer_type,
              ": ",
              e$message)
    })
    
    # Summarize results
    common_genes <- intersect(rownames(mut_data), rownames(exp_data$counts))
    common_samples <- intersect(colnames(mut_data), colnames(exp_data$counts))
    message(
      "Number of common samples: ",
      length(common_samples),
      " (",
      ncol(mut_data),
      ", ",
      ncol(exp_data$counts),
      ")"
    )
    message(
      "Number of common genes: ",
      length(common_genes),
      " (",
      nrow(mut_data),
      ", ",
      nrow(exp_data$counts),
      ")"
    )
    
    message("=== Completed processing for ", cancer_type, " ===")
  }
}

# List of cancer types to process
cancers <- c(
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

# Execute main function
main(cancers)

message("\n=== Processing PANCAN dataset ===")
mut_data_file <- "./data/raw/mc3.v0.2.8.PUBLIC.nonsilentGene.xena.gz"
exp_data_file <- "./data/raw/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz"
outfile <- "./data/processed/PANCAN"
if (!dir.exists(outfile)) {
  dir.create(outfile, recursive = TRUE)
  message("Created output directory: ", outfile)
}
mut_data <- read.delim(
  mut_data_file,
  header = TRUE,
  as.is = TRUE,
  check.names = FALSE
)
invalid_gene_index <- which(is.na(mut_data[, 1]) |
                                mut_data[, 1] == "" |
                                mut_data[, 1] == ".")
if (length(invalid_gene_index) > 0) {
  mut_data <- mut_data[-invalid_gene_index, ]
}
rownames(mut_data) <- mut_data[,1]
mut_data <- mut_data[,-1]
mut_data <-
  mut_data[which(rowSums(mut_data) != 0), which(colSums(mut_data) != 0)]
write.table(
  mut_data,
  file = file.path(outfile,"mut_data.tsv"),
  sep = "\t",
  quote = F,
  row.names = T,
  col.names = T
)
exp_data <- read.delim(
    exp_data_file,
    header = TRUE,
    as.is = TRUE,
    check.names = FALSE
  )
if (any(duplicated(exp_data[,1]))) {
  exp_genes_id <- exp_data[,1]
  exp_data <- exp_data[,-1]
  table_results <- table(exp_genes_id)
  duplicated_genes <- names(table_results)[table_results > 1]
  for (gene in duplicated_genes) {
    gene_idx <- which(exp_genes_id == gene)
    exp_data[gene_idx[1],] <- colMeans(exp_data[gene_idx,],na.rm = T)
    exp_data <- exp_data[-gene_idx[-1],]
    exp_genes_id <- exp_genes_id[-gene_idx[-1]]
  }
  rownames(exp_data) <- exp_genes_id
} else {
  rownames(exp_data) <- exp_data[,1]
  exp_data <- exp_data[,-1]
}
invalid_gene_index <- which(is.na(rownames(exp_data)) |
                                rownames(exp_data) == "" |
                                rownames(exp_data) == ".")
if (length(invalid_gene_index) > 0) {
  exp_data <- exp_data[-invalid_gene_index, ]
}
exp_data <- exp_data[rowSums(is.na(exp_data)) == 0, ]
exp_data <- as.matrix(exp_data)
write.table(
  exp_data,
  file = file.path(outfile, "exp_data.tsv"),
  sep = "\t",
  quote = F,
  row.names = T,
  col.names = T
)
