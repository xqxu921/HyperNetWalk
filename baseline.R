library(Matrix)
source("model.R")
source("./benchmark/pernalized2cohort.R")

cancers <- c(
#   "PANCAN",
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

#### Random ####
for (cancer_type in cancers){
  message("Process ", cancer_type,"...\n")
  mut_file <- file.path("./data/processed_data",cancer_type,"mut_data.tsv")
  mut_data <- get_mut_data(mut_file)
  samples <- colnames(mut_data)
  outfile <- file.path("./results/Random",cancer_type)
  if (!dir.exists(outfile)){
    dir.create(outfile,recursive = TRUE)
  }
  
  # genes_ranking <- list()
  # for (sample in samples){
  #   mut_genes <- rownames(mut_data)[which(mut_data[,sample]!=0)]
  #   genes_ranking[[sample]] <- sample(mut_genes)
  # }
  # ranking_df <- data.frame(
  #   sample = names(genes_ranking),
  #   ranking = sapply(genes_ranking, function(x)
  #     paste(x, collapse = ","))
  # )
  # write.table(
  #   ranking_df,
  #   file = file.path(outfile, "genes_ranking.txt"),
  #   sep = "\t",
  #   col.names = FALSE,
  #   row.names = FALSE,
  #   quote = FALSE
  # )
  # personalized2cohort_score(outfile)
  # 生成一个nrow(mut_data)行1列的矩阵，值为0-1之间随机数
  scores <- matrix(runif(nrow(mut_data)), ncol = 1)
  rownames(scores) <- rownames(mut_data)
  scores <- sort(scores[,1], decreasing = TRUE)  
  scores <- scores/sum(scores)
  scores <- as.matrix(scores)
  write.table(scores,file.path(outfile,"genes_ranking_cohort.txt"),quote = F,sep = "\t",row.names = T,col.names = F)


  message("Finised ",cancer_type,"!\n")
}

#### Frequency-based ####
for (cancer_type in cancers){
  message("Process ", cancer_type,"...\n")
  mut_file <- file.path("./data/processed_data",cancer_type,"mut_data.tsv")
  mut_data <- get_mut_data(mut_file)
  samples <- colnames(mut_data)
  outfile <- file.path("./results/Frequency",cancer_type)
  if (!dir.exists(outfile)){
    dir.create(outfile,recursive = TRUE)
  }
  mut_fq <- rowSums(mut_data)/ncol(mut_data)
  
  genes_ranking <- list()
  for (sample in samples){
    mut_genes <- rownames(mut_data)[which(mut_data[,sample]!=0)]
    fq <- sort(mut_fq[mut_genes],decreasing = TRUE)
    genes_ranking[[sample]] <- names(fq)
  }
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

  mut_fq <- sort(mut_fq,decreasing = TRUE)
  mut_fq <- as.matrix(mut_fq,ncol=1)
  write.table(mut_fq,file.path(outfile,"genes_ranking_cohort.txt"),quote = F,sep = "\t",row.names = T,col.names = F)
}
