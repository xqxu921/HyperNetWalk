library(data.table)
library(dplyr)
library(WGCNA)
library(Matrix)
source("model.R")
PPI_file <- "./data/9606.protein.links.v12.0.onlyAB.tsv"
protein_file <- "./data/9606.protein.info.v12.0.txt"

PPI_df <- fread(PPI_file)
protein_df <- fread(protein_file,
                    select = c("#string_protein_id", "preferred_name"))
setnames(protein_df, "#string_protein_id", "protein_id")
protein_map1 <- match(PPI_df$protein1, protein_df$protein_id)
protein_map2 <- match(PPI_df$protein2, protein_df$protein_id)
PPI_df$protein1 <- protein_df$preferred_name[protein_map1]
PPI_df$protein2 <- protein_df$preferred_name[protein_map2]
PPI_df$combined_score <- PPI_df$combined_score / 1000
fwrite(PPI_df,
       file = "./data/STRINGv12.txt",
       quote = FALSE,
       sep = "\t")

GRN_file <- "./data/human_TF_Target.txt"
GRN_df <- fread(GRN_file, header = F)
#先只保留V6是TF或者Gene的行,然后只保留第1列和第3列,更改列名为TF和Target,最后去重

GRN_df <- GRN_df %>%
  # filter(V6 %in% c("TF", "Gene")) %>%
  select(V1, V3) %>%
  rename(TF = V1, Target = V3)

if (any(duplicated(GRN_df))) {
  GRN_df <- unique(GRN_df)
}
fwrite(GRN_df,
       file = "./data/RegNet_human_V2.txt",
       quote = FALSE,
       sep = "\t")

GC_file <- "./data/human_core_TF_Target.txt"
GC_df <- fread(GC_file, header = F)
GC_df <- GC_df %>%
  # filter(V6 %in% c("TF", "Gene")) %>%
  select(V1, V3) %>%
  rename(TF = V1, Target = V3)
if (any(duplicated(GC_df))) {
  GC_df <- unique(GC_df)
}
fwrite(GC_df,
       file = "./data/RegNet_human_core_V2.txt",
       quote = FALSE,
       sep = "\t")