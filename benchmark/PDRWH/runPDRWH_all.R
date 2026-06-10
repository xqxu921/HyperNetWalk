#run PDRWH in different cohorts
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("PDRWH_software.R")
source("../pernalized2cohort.R")

# cancers <- c(
#   "BRCA",
#   "BLCA",
#   "COAD",
#   "ESCA",
#   "GBM",
#   "HNSC",
#   "KIRC",
#   "LGG",
#   "LIHC",
#   "LUAD",
#   "LUSC",
#   "PAAD",
#   "SKCM",
#   "STAD",
#   "THCA",
#   "UCEC"
# )
# cancers <- c(
#   "KIRP",
#   "PRAD",
#   "READ"
# )
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

for (i in 1:length(cancers)) {
  cancer_type <- cancers[i]
  # print("----------------------------------------------------")
  # print(paste0("Run PDRWH for cancer type: ", cancer_type, "..."))
  # cat("\n\n")
  
  # snv_file <- paste0("../../data/processed_data/", cancer_type, "/mut_data.tsv")
  # exp_file <- paste0("../../data/processed_data/", cancer_type, "/exp_tpm_data.tsv")
  outfile_dir <- paste0("../../results/PDRWH/", cancer_type, "/")
  # if (!file.exists(outfile_dir)) {
  #   dir.create(outfile_dir,recursive = TRUE)
  # }
  # network_file <- "../../data/STRINGv12.txt"
  # PDRWHscore(cancer_type, snv_file, exp_file, network_file, outfile_dir)
  # resultTransform(cancer_type, outfile_dir, outfile_dir)
  # print("----------------------------------------------------")
  # print(paste0("PDRWH for ", cancer_type, ":  Achieved!"))
  # print("----------------------------------------------------")
  # cat("\n\n")
  
  # Reformatting outputs
  load(paste0(
    "../../results/PDRWH/",
    cancer_type,
    "/Importance_Score.Rdata"
  ))
  genes_ranking <- lapply(Importance_Score, rownames)
  # ranking_df <- data.frame(
  #   sample = names(Importance_Score),
  #   ranking = sapply(genes_ranking, function(x)
  #     paste(x, collapse = ","))
  # )
  # write.table(
  #   ranking_df,
  #   file = paste0(outfile_dir, "genes_ranking.txt"),
  #   sep = "\t",
  #   col.names = FALSE,
  #   row.names = FALSE,
  #   quote = FALSE
  # )
  # personalized2cohort_score(outfile_dir)

  genes_total <- unique(unlist(genes_ranking))
  res_pers <- Matrix(0, nrow = length(genes_total), ncol = length(genes_ranking),
    dimnames = list(genes_total, names(genes_ranking)),
    sparse = TRUE
  )

  for (j in 1:length(Importance_Score)) {
    sample_j <- names(Importance_Score)[j]
    genes_ranking_j <- genes_ranking[[j]]
    res_pers[genes_ranking_j, sample_j] <- Importance_Score[[j]][genes_ranking_j, 1]
  }

  personalized2cohort_score_new(res_pers, output_dir = outfile_dir)

  if (i == length(cancers)) {
    print("----------------------------------------------------")
    print(paste0("Congratrulations! All achieved!"))
    print("----------------------------------------------------")
  }
}
