library(dndscv)
library(dplyr)
cancers <- c("BRCA","COAD","HNSC","KIRC","KIRP","LIHC","LUAD","LUSC","PRAD","STAD","THCA","UCEC")
for (cancer in cancers) {
  mut_file <- paste0("../../data/rawdata/TCGA-", cancer, ".somaticmutation_wxs.tsv")
  data <- read.table(mut_file, header=TRUE, stringsAsFactors=FALSE, sep="\t")
  mutations <- data %>%
    mutate(
      sampleID = sample,
      chr = sub("^chr","", chrom),
      pos = start,
      ref = ref,
      mut = alt
    ) %>%
    select(sampleID, chr, pos, ref, mut) %>%
    distinct()
  dndsout <- dndscv(mutations,refdb = "hg38")
  sel_cv <- dndsout$sel_cv
  res_coh <- 1 / seq_along(sel_cv$gene_name)
  names(res_coh) <- sel_cv$gene_name
  res_coh <- as.matrix(res_coh)
  outfile <- paste0("../../results/dndscv/", cancer)
  dir.create(outfile, showWarnings = FALSE, recursive = TRUE)
  write.table(
    res_coh,
    file.path(outfile, "genes_ranking_cohort.txt"),
    quote = F,
    sep = "\t",
    row.names = T,
    col.names = F
  )
}

