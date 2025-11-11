#!/usr/bin/env Rscript

library(argparse)
source("src/model.R")

# 创建解析器
parser <- ArgumentParser(description = "Cancer driver prediction model")

# 添加参数
parser$add_argument("--mode", type="character", default="single_cancer",
                    help="Mode: pancancer or single_cancer")
parser$add_argument("--level", type="character", default="cohort",
                    help="Prediction level: cohort, individual or both")
parser$add_argument("--cancer_type", type="character", default=NULL,
                    help="Cancer type, required for single_cancer mode")
parser$add_argument("--input", type="character", required=TRUE,
                    help="Input data directory")
parser$add_argument("--ppi", type="character", required=TRUE,
                    help="PPI network file")
parser$add_argument("--grn", type="character", required=TRUE,
                    help="GRN network file")
parser$add_argument("--output", type="character", required=TRUE,
                    help="Output directory")
parser$add_argument("--cores", type="integer", default=1,
                    help="Number of CPU cores to use")

# 解析命令行参数
args <- parser$parse_args()

# ----------------- 使用参数 -----------------
cat("Mode:", args$mode, "\n")
cat("Level:", args$level, "\n")
cat("Cancer type:", args$cancer_type, "\n")
cat("Input dir:", args$input, "\n")
cat("PPI file:", args$ppi, "\n")
cat("GRN file:", args$grn, "\n")
cat("Output dir:", args$output, "\n")
cat("Cores:", args$cores, "\n")

# -------------------------------------------
if (args$mode == "pancancer") {
    mut_data_file <- file.path(args$input,"PANCAN","mut_data.tsv")
    exp_data_file <- file.path(args$input,"PANCAN","exp_data.tsv")
    cancer_type <- "PANCAN"
} else {
    mut_data_file <- file.path(args$input,args$cancer_type,"mut_data.tsv")
    exp_data_file <- file.path(args$input,args$cancer_type,"exp_tpm_data.tsv")
    cancer_type <- args$cancer_type
}
ME_file <- paste0("./data/NETWORK/", cancer_type,"_me_net.txt")
DEN_file <- file.path("./data/processed", cancer_type,"CEN_union_matrix_uptri.rds")
if (!file.exists(ME_file)){
  ME_file <- NULL
}
DGScore(cancer_type = cancer_type,
        outfile_dir = args$output,
        mut_data_file = mut_data_file,
        exp_data_file = exp_data_file,
        PPI_file = args$ppi,
        GRN_file = args$grn,
        DEN_file = DEN_file,
        ME_file = ME_file,
        num_cores = args$cores,
        level = args$level
        )