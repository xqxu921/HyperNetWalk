#!/usr/bin/env Rscript

library(argparse)
# 假设你有评估相关的函数在另一个文件中
# source("src/evaluation_functions.R")

# 创建解析器
parser <- ArgumentParser(description = "Cancer driver prediction evaluation")

# 添加参数
parser$add_argument("--mode", type="character", required=TRUE,
                    choices = c("pancancer", "single_cancer"),
                    help="Evaluation mode: pancancer or single_cancer")
parser$add_argument("--level", type="character", required=TRUE,
                    choices = c("cohort", "individual", "all"),
                    help="Prediction level: cohort, individual or both")
parser$add_argument("--cancer_type", type="character", default=NULL,
                    help="Cancer type, required for single_cancer mode")
parser$add_argument("--predicted", type="character", required=TRUE,
                    help="Predicted results file")
parser$add_argument("--benchmark", type="character", required=TRUE,
                    help="Benchmark gene set")
parser$add_argument("--output", type="character", required=TRUE,
                    help="Output evaluation results file")

# 解析命令行参数
args <- parser$parse_args()

# ----------------- 参数验证 -----------------
cat("Mode:", args$mode, "\n")
cat("Level:", args$level, "\n")
cat("Cancer type:", args$cancer_type, "\n")
cat("Predicted:", args$predicted, "\n")
cat("Benchmark:", args$benchmark, "\n")
cat("Output:", args$output, "\n")

# 验证必要参数
if (args$mode == "single_cancer" && is.null(args$cancer_type)) {
  stop("Cancer type must be specified for single_cancer mode")
}

# 验证文件存在
if (!file.exists(args$predicted)) {
  stop("Predicted file not found: ", args$predicted)
}

if (args$mode == "pancancer"){
  cancer_type <- "PANCAN"
} else {
  cancer_type <- args$cancer_type
}

# ----------------- 主评估逻辑 -----------------

if (!dir.exists(dirname(args$output))) {
  dir.create(dirname(args$output), recursive = TRUE)
}

# 读取基准基因集
if (args$benchmark == "CGC") {
    benchmark_genes <- read.delim("./data/DRIVER/CGC_Tier1.tsv",
                             header = T,
                             as.is = TRUE)[, 1]
} else {
    cgc_genes <- read.delim("./data/DRIVER/CGC_Tier1.tsv",
                             header = T,
                             as.is = TRUE)[, 1]
    intogen_genes <- read.delim("./data/DRIVER/Compendium_Cancer_Genes.tsv",
                                header = T,
                                as.is = TRUE)
    intogen_genes <- intogen_genes[intogen_genes$CANCER_TYPE == args$cancer_type,]$SYMBOL                            
    benchmark_genes <- union(cgc_genes, intogen_genes)
}

if (args$level == "cohort" | args$level == "all") {
    cohort_results <- read.table(
      file.path(args$predicted,"cohort_results.txt"),
      sep = "\t",
      header = FALSE,
      stringsAsFactors = FALSE
    )
    TP_top_50 <- length(intersect(cohort_results[1:50, 1], benchmark_genes))
    coh_prec_50 <- TP_top_50 / 50
    TP_top_100 <- length(intersect(cohort_results[1:100, 1], benchmark_genes))
    coh_prec_100 <- TP_top_100 / 100
    cat("For cancer type:", cancer_type, "\n")
    cat("Cohort Level - coh_precision@50:", coh_prec_50, "\n")
    cat("Cohort Level - coh_precision@100:", coh_prec_100, "\n")
    # 写入output文件 包括当前时间及结果，格式规整
    write(paste0(Sys.time(), "\nCohort Level - coh_precision@50: ", coh_prec_50,
                 "\nCohort Level - coh_precision@100: ", coh_prec_100, "\n"),
          file = args$output)
}

if (args$level == "individual" | args$level == "all") {
    individual_results <- read.table(
      file.path(args$predicted,"individual_results.txt"),
      sep = "\t",
      header = FALSE,
      stringsAsFactors = FALSE
    )
    ind_prec_10_list <- c()
    for (i in 1:nrow(individual_results)){
      sample_genes <- unlist(strsplit(individual_results[i,2],","))
      if (length(sample_genes) < 10){
        next
      }
      TP_num <- length(intersect(sample_genes, benchmark_genes))
      if (TP_num < 3){
        next
      }
      TP_top_10_ind <- length(intersect(sample_genes[1:10], benchmark_genes))
      ind_prec_10 <- TP_top_10_ind / 10
      ind_prec_10_list <- c(ind_prec_10_list, ind_prec_10)
    }
    avg_ind_prec_10 <- mean(ind_prec_10_list)
    cat("For cancer type:", cancer_type, "\n")
    cat("Individual Level - ind_precision@10:", avg_ind_prec_10, "\n")
    # 继续写入output文件 包括当前时间及结果，格式规整
    write(paste0(Sys.time(), "\nIndividual Level - ind_precision@10: ", avg_ind_prec_10, "\n"),
          file = args$output, append = TRUE)          
}
# -------------------------------------------
cat("Evaluation completed. Results saved to:", args$output, "\n")