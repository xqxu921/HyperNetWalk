# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
invisible(sapply(list.files("R", pattern = "\\.R$", full.names = TRUE), source))
library(parallel)
library(doSNOW)
library(foreach)
library(progress)
source('../pernalized2cohort.R')
source('../../model.R')
cancers <- c(
  # "PANCAN",
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

PPI_file <- "../../data/STRINGv12.txt"
for (i in 2:length(cancers)) {
  cancer_type <- cancers[i]
  message("Start Prodigy for ", cancer_type)
  mut_data_file <- paste0("../../data/processed_data/", cancer_type, "/mut_data.tsv")
  exp_data_file <- paste0("../../data/processed_data/", cancer_type, "/exp_counts_data.tsv")
  outfile_dir <- paste0("../../results/PRODIGY_hrw/", cancer_type, "/")
  if (!dir.exists(outfile_dir)) {
    dir.create(outfile_dir,recursive = TRUE)
  }
  
  mut_data <- get_mut_data(mut_data_file)
  colnames(mut_data) <- substr(colnames(mut_data), 1, 15)
  
  exp_data <- read.delim(
    exp_data_file,
    header = TRUE,
    as.is = TRUE,
    check.names = FALSE
  )
  if (sum(substr(colnames(exp_data),14,16)=="11A")==0){
    next()
  }
  exp_data <- exp_data[, grepl("-01A$|-11A$", colnames(exp_data))]
  colnames(exp_data) <- substr(colnames(exp_data), 1, 15)
  exp_mat <- as.matrix(exp_data)
  
  PPI_edges <- read.delim(PPI_file, header = T)
  edge_idx <- which(PPI_edges[, 3] >= 0.7)
  PPI_edges <- as.matrix(PPI_edges[edge_idx, ])
  
  samples <- intersect(colnames(exp_mat), colnames(mut_data))
  
  exp_mat <- exp_mat[which(rownames(exp_mat) %in% unique(c(PPI_edges[, 1], PPI_edges[, 2]))), ]
  # DEGs <- get_DEGs(
  #   exp_mat,
  #   samples,
  #   sample_origins = NULL,
  #   beta = 2,
  #   gamma = 0.05
  # )
  # # Identify sample origins (tumor or normal)
  # sample_origins <- rep("tumor", ncol(exp_mat))
  # sample_origins[substr(colnames(exp_mat), nchar(colnames(exp_mat)[1]) - 1, nchar(colnames(exp_mat)[1])) == "11"] <- "normal"
  # list_of_pathways <- get_pathway_list_from_graphite(
  #   source = "reactome",
  #   minimal_number_of_nodes = 10,
  #   num_of_cores = 100
  # )
  # save(DEGs,list_of_pathways, file = paste0(outfile_dir, "variants.Rdata"))
  # # Run PRODIGY
  # all_patients_scores <- PRODIGY_cohort(
  #   snv_matrix = mut_data,
  #   expression_matrix = exp_mat,
  #   network = PPI_edges,
  #   samples = samples,
  #   DEGs = DEGs,
  #   alpha = 0.05,
  #   pathway_list = list_of_pathways,
  #   num_of_cores = 100,
  #   sample_origins = sample_origins,
  #   write_results = F,
  #   results_folder = outfile_dir,
  #   beta = 2,
  #   gamma = 0.05,
  #   delta = 0.05
  # )
  # save(all_patients_scores,
  #      file = paste0(outfile_dir, "results.Rdata"))
  # Get driver gene rankings for all samples
  load(file.path("../../results/PRODIGY/", cancer_type, "results.Rdata"))
  # genes_ranking <- analyze_PRODIGY_results(all_patients_scores)
  # ranking_df <- data.frame(
  #   sample = substr(names(genes_ranking),1,12),
  #   ranking = sapply(genes_ranking, function(x)
  #     paste(x, collapse = ","))
  # )
  # ranking_df[ranking_df == ""] <- NA
  # ranking_df <- na.omit(ranking_df)
  # write.table(
  #   ranking_df,
  #   file = paste0(outfile_dir, "genes_ranking.txt"),
  #   sep = "\t",
  #   col.names = FALSE,
  #   row.names = FALSE,
  #   quote = FALSE
  # )
  # personalized2cohort_score(outfile_dir)
  # message("Analysis completed successfully for ", cancer_type)

  
  libraries = c("DESeq2","igraph","plyr","biomaRt","MASS","mixtools")
	for(j in 1:length(libraries)){
	try({library(libraries[j],character.only=T)})
	}
	if(class(all_patients_scores) == "matrix"){ all_patients_scores = list(all_patients_scores) }
	for(j in 1:length(all_patients_scores))
	{	
		all_patients_scores[[j]][which(is.na(all_patients_scores[[j]]))] = 0
		all_patients_scores[[j]][all_patients_scores[[j]] < 0] = 0
	}
	Prodigy_rankings = list()
	for(j in 1:length(all_patients_scores))
	{
		if(is.null(all_patients_scores[[j]])){Prodigy_rankings[[j]] = c(); next}
		pathways_to_take = c()
		#single pathway
		if(is.null(nrow(all_patients_scores[[j]])))
		{
			ranking = sort(all_patients_scores[[j]],decreasing=T)
		#less than 4 mutations, no pathway filtering
		} else if (nrow(all_patients_scores[[j]])==1) {
			sorted_idx = order(all_patients_scores[[j]],decreasing=T)
			ranking = all_patients_scores[[j]][,sorted_idx]

		} else if(ncol(all_patients_scores[[j]]) < 4) {
			ranking = sort(apply(all_patients_scores[[j]],2,function(x) sum(sort(x,decreasing=T))),decreasing=T)
		#filter pathways
		} else {
			pathways_to_take = names(which(apply(all_patients_scores[[j]],1,function(x) length(which(x > 0))) < ncol(all_patients_scores[[j]])/2))
			if(length(pathways_to_take) < 2)
			{
				ranking = sort(all_patients_scores[[j]][pathways_to_take,],decreasing=T)
				#if all pathways are filtered, abort pathway filtering
				if(all(ranking==0))
				{
					pathways_to_take = rownames(all_patients_scores[[j]])
					ranking = sort(apply(all_patients_scores[[j]][pathways_to_take,,drop=F],2,function(x) sum(sort(x,decreasing=T))),decreasing=T)
				}
			} else {
				ranking = sort(apply(all_patients_scores[[j]][pathways_to_take,],2,function(x) sum(sort(x,decreasing=T))),decreasing=T)
				if(all(ranking==0))
				{
					pathways_to_take = rownames(all_patients_scores[[j]])
					ranking = sort(apply(all_patients_scores[[j]][pathways_to_take,],2,function(x) sum(sort(x,decreasing=T))),decreasing=T)
				} else {
					pathways_to_take = pathways_to_take[which(apply(all_patients_scores[[j]][pathways_to_take,],1,sum)!=0)]
				}
			}
		}
		ranking = ranking[ranking > 0]
		if(length(ranking) == 0) {Prodigy_rankings[[j]] = c(); next}
		# single pathway was used, no bimodel distribution 
		if(length(pathways_to_take) < 2){Prodigy_rankings[[j]] = ranking;next}
		# check bimodel distribution
		possible_error = tryCatch(
		{
				bimodel_dist = normalmixEM(ranking,k=2)
				unimodel_dist = MASS::fitdistr(ranking,"normal")
		},
		error=function(cond) {
				cond
		})
		if(inherits(possible_error, "error")){ Prodigy_rankings[[j]] = ranking
		} else {
			# check which distribution is more likely
			if(bimodel_dist$loglik > unimodel_dist$loglik)
			{
					Prodigy_rankings[[j]] = ranking[bimodel_dist$posterior[,which(bimodel_dist$mu == max(bimodel_dist$mu))] > 0.5]
			} else {
					Prodigy_rankings[[j]] = ranking
			}
		}
	}

  genes <- unique(unlist(lapply(Prodigy_rankings, names)))
  res_pers <- Matrix(0, nrow = length(genes), ncol = length(Prodigy_rankings), dimnames = list(genes, names(all_patients_scores)), sparse = TRUE)
  for (j in 1:length(Prodigy_rankings)) {
    if (length(Prodigy_rankings[[j]]) == 0) {
      next
    }
    res_pers[names(Prodigy_rankings[[j]]), j] <- Prodigy_rankings[[j]]
  }
  zero_cix <- which(colSums(res_pers) == 0)
  if (length(zero_cix) > 0) {
    res_pers <- res_pers[, -zero_cix, drop = FALSE]
  }
  colnames(res_pers) <- substr(colnames(res_pers), 1, 12)
  personalized2cohort_score_new(res_pers, output_dir = outfile_dir)

  if (i==length(cancers)) {
    message("Congratrulations! All achieved!")
  }
}
