source("model.R")

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
start_time <- proc.time()
for (c in 1:length(cancers)){
    cancer_type <- cancers[c]
    mut_data_file <- file.path("./data/processed_data",cancer_type,"mut_data.tsv")
    exp_data_file <- file.path("./data/processed_data",cancer_type,"exp_tpm_data.tsv")
    output_dir <- "./results/HyperNetWalk"
    run_hypernetwalk(cancer_type, mut_data_file, exp_data_file, output_dir)
}
elapsed_time <- proc.time() - start_time
message("Total execution time: ", elapsed_time["elapsed"], " seconds")