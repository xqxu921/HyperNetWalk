library(reshape2)
library(fgsea)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(forcats)
library(tidyr)
library(msigdbr)
library(cluster)
library(mclust)
library(aricode)
library(Rtsne)
library(stringr)

source("model.R")
source("plot_formal.R")
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
fig_dir <- "./figs/exp3_cancerspecific/"
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
}
# Cancer specific driver genes from IntOGen
tcga_map <- c(
  "BRCA" = "BRCA",
  "COAD" = "COAD",
  "HNSC" = "HNSC",
  "KIRC" = "CCRCC",
  "KIRP" = "PRCC",
  "LIHC" = "HCC",
  "LUAD" = "LUAD",
  "LUSC" = "LUSC",
  "PRAD" = "PRAD",
  "STAD" = "STAD",
  "THCA" = "WDTC",
  "UCEC" = "UCEC"
)
IntOGen_drivers <- read.delim(
"./reference_dg/IntOGen_Compendium_Cancer_Genes.tsv",
    header = T,
    as.is = TRUE
)

IntOGen_drivers$TCGA_CANCER <- names(tcga_map)[match(IntOGen_drivers$CANCER_TYPE, tcga_map)]
IntOGen_drivers <- IntOGen_drivers[!is.na(IntOGen_drivers$TCGA_CANCER), c("SYMBOL", "TCGA_CANCER")]
IntOGen_drivers <- IntOGen_drivers[!duplicated(IntOGen_drivers), ]
gene_cancer_counts <- IntOGen_drivers %>%
  group_by(SYMBOL) %>%
  summarise(n_cancers = n_distinct(TCGA_CANCER))
max_cancer_allowed <- 3
genes_to_keep <- gene_cancer_counts %>%
  filter(n_cancers <= max_cancer_allowed) %>%
  pull(SYMBOL)
IntOGen_drivers <- IntOGen_drivers[IntOGen_drivers$SYMBOL %in% genes_to_keep, ]
table(IntOGen_drivers$TCGA_CANCER)
all_pathways_list <- split(IntOGen_drivers$SYMBOL, IntOGen_drivers$TCGA_CANCER)

methods <- list(
    HyperNetWalk = "./results/HyperNetWalk/",
    DriverRWH = "./results/DriverRWH/",
    DriverMP = "./results/DriverMP/",
    DawnRank = "./results/DawnRank/"
)

res_coh <- list()
hit_df <- data.frame()
gsea_df <- data.frame()
ora_df <- data.frame()
for (cancer in cancers){
    cat("Processing", cancer, "...\n")
    mut_file <- file.path("./data/processed_data",cancer,"mut_data.tsv")
    mut_data <- get_mut_data(mut_file,1)
    colnames(mut_data) <- substr(colnames(mut_data), 1, 12)
    mut_genes <- rownames(mut_data)
    drivers_c <- IntOGen_drivers[IntOGen_drivers$TCGA_CANCER %in% cancer, ]$SYMBOL
    if (length(intersect(mut_genes, drivers_c)) == 0) {
        next
    }
    for (method in names(methods)){
        cat("  Method:", method, "...\n")
        res_file <- file.path(methods[[method]], cancer, "genes_ranking_cohort.txt")
        res <- read.table(
            res_file,
            row.names = 1
        )
        mut_genes_sup <- setdiff(mut_genes,rownames(res))
        if (length(mut_genes_sup) > 0) {
            padding <- matrix(0,
                                nrow = length(mut_genes_sup),
                                ncol = ncol(res))
            rownames(padding) <- mut_genes_sup
            colnames(padding) <- colnames(res)
            res <- rbind(res, padding)
        } 
        res_coh[[method]][[cancer]] <- res
        hit_res <- run_overlap_cross_cancer(res, IntOGen_drivers, top_n = 100) %>%
            mutate(
                Method = method,
                Predicted_Cancer = cancer
            )
        hit_df <- rbind(hit_df, hit_res)
        fgsea_res <- run_fgsea(res, all_pathways_list)
        if (nrow(fgsea_res) > 0) {
            # 把 fgsea 输出的 dataframe 直接拼接到大表里
            temp_df <- data.frame(
                Method = method,
                Predicted_Cancer = cancer,
                Reference_Cancer = fgsea_res$pathway, # fgsea 自动把 list 的名字放在 pathway 列
                NES = fgsea_res$NES,
                pval = fgsea_res$pval,
                padj = fgsea_res$padj
            )
            gsea_df <- rbind(gsea_df, temp_df)
        }
        #对cancers中所有cancer类型进行ORA分析，统计前100个基因中有多少是IntOGen驱动基因，并计算富集的p值，拼成一个大表
        ora_res <- lapply(cancers, function(ref_cancer){
            run_ora_topN(res, IntOGen_drivers, ref_cancer, top_n = 100) %>%
                mutate(
                    Method = method,
                    Predicted_Cancer = cancer,
                    Reference_Cancer = ref_cancer
                )
        }) %>% bind_rows()
        ora_df <- rbind(ora_df, ora_res)
    }
}

# Panel A: Cross-cancer overlap heatmap
res_cross_cancer <- lapply(names(methods),function(method){
    df <- hit_df %>% filter(Method == method) %>%
        group_by(Predicted_Cancer) %>%
        mutate(
            is_best = (Overlap_Proportion == max(Overlap_Proportion, na.rm = TRUE) & Overlap_Proportion > 0),
            label_text = if_else(Overlap_Proportion == 0, "0", sprintf("%.2f", Overlap_Proportion))
        ) %>%
        ungroup()
    p <- ggplot(df, aes(x = Reference_Cancer, y = Predicted_Cancer, fill = Overlap_Proportion)) +
        geom_tile(color = "white", linewidth = 0.8) + 
        geom_text(
            aes(
                label = label_text,
                color = is_best 
            ),
            size = 4, 
            fontface = "bold",
            show.legend = FALSE
        ) +
        scale_fill_gradient(
            low = "#F2F2F2", high = "#D94841", 
            limits = c(0, max(df$Overlap_Proportion, na.rm = TRUE)),
            name = "Overlap Proportion"
        ) +
        scale_color_manual(values = c("FALSE" = "#222222", "TRUE" = "white")) + 
        coord_fixed(ratio = 1) + 
        scale_x_discrete(limits = sort(unique(df$Reference_Cancer))) +
        scale_y_discrete(limits = rev(sort(unique(df$Predicted_Cancer)))) +
        theme_minimal(base_size = 14) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", face = "bold"),
            axis.text.y = element_text(color = "black", face = "bold"),
            axis.title = element_text(size = 14, face = "bold"),
            panel.grid = element_blank(), 
            legend.position = "right",
            legend.key.height = unit(1.0, "cm"),
            legend.key.width = unit(0.3, "cm"),
            legend.title = element_text(face = "bold", size = 11, angle = 90, hjust = 0.5),
            legend.title.position = "right",
            plot.title = element_blank() 
        ) +
        labs(
            x = "Reference Cancer Type (IntOGen Specific Drivers)",
            y = "Predicted Cancer Type (Model Scoring Cohort)"
        )
    ggsave(
        filename = file.path(fig_dir, paste0(method, "_cross_cancer_overlap.pdf")),
        plot = p,
        width = 6.5,     # 宽度 6.5 英寸 (约 16.5 厘米，几乎占满双栏宽度)
        height = 6.0,    # 高度 6.0 英寸 (留出 0.5 英寸的高度差给底部的长标签和右侧的图例)
        units = "in",    # 尺寸单位：英寸 (inches)
        dpi = 300,       # 矢量图的 dpi 其实无所谓，但如果以后要导 PNG，设为 300 是底线
        bg = "white" 
    )
    return(p)
})

# Panel B: comparison of different methods
#只保留predicted_cancer和reference_cancer相同的行,每种方法画boxplot+jitter图
library(tidyr)
library(ggsignif)
hit_df_filter <- hit_df %>% filter(Predicted_Cancer == Reference_Cancer)
median_scores <- hit_df_filter %>%
    filter(Method != "HyperNetWalk") %>%  # 注意：你的数据列名是首字母大写的 Method
    group_by(Method) %>%
    summarize(med_prec = median(Overlap_Proportion, na.rm = TRUE)) %>%
    arrange(desc(med_prec))

method_levels <- c("HyperNetWalk", median_scores$Method)
hit_df_filter$Method <- factor(hit_df_filter$Method, levels = method_levels)
hit_df_filter <- hit_df_filter %>%
    mutate(is_proposed = ifelse(Method == "HyperNetWalk", "Proposed", "Baseline"))

best_baseline <- as.character(median_scores$Method[1])
test_df <- hit_df_filter %>%
    filter(Method %in% c("HyperNetWalk", best_baseline)) %>%
    transmute(
        cancer = Reference_Cancer,
        Method = as.character(Method),
        value = Overlap_Proportion
    ) %>%
    pivot_wider(names_from = Method, values_from = value) %>%
    filter(
        !is.na(.data[["HyperNetWalk"]]),
        !is.na(.data[[best_baseline]])
    )

p_val <- wilcox.test(
    test_df[["HyperNetWalk"]],
    test_df[[best_baseline]],
    paired = TRUE,
    exact = FALSE
)$p.value

delta_median <- median(
    test_df[["HyperNetWalk"]] - test_df[[best_baseline]],
    na.rm = TRUE
)

p_label <- if (p_val < 0.001) {
    paste0("italic(P) < 0.001")
} else {
    paste0("italic(P) == ", sprintf("%.3f", p_val))
}

y_max <- max(hit_df_filter$Overlap_Proportion, na.rm = TRUE)
y_min <- min(hit_df_filter$Overlap_Proportion, na.rm = TRUE)
y_range <- y_max - y_min
if (y_range == 0) y_range <- 1

y_pos <- y_max + 0.06 * y_range

p_hit <- ggplot(hit_df_filter, aes(x = Method, y = Overlap_Proportion)) +
    # 箱线图层：半透明、去离群点、固定宽度
    geom_boxplot(
        aes(fill = is_proposed),
        outlier.shape = NA,
        alpha = 0.6,
        width = 0.62,
        color = "black",
        linewidth = 0.55
    ) +
    # 散点层：带白色描边的深灰弹珠
    geom_jitter(
        width = 0.16,
        size = 2.1,
        alpha = 0.8,
        shape = 21,
        fill = "#4D4D4D",
        color = "white",
        stroke = 0.5
    ) +
    
    # 手动指定颜色：Proposed(砖红) vs Baseline(高级灰)
    scale_fill_manual(values = c("Proposed" = "#D94841", "Baseline" = "#D9D9D9")) +
    
    # 留出顶部和底部的呼吸空间，防止散点贴边
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    
    # 替换为经典的坐标轴风格 (去背景网格，保留黑实线)
    theme_classic(base_size = 14) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", face = "bold"),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title = element_text(face = "bold", size = 14),
        axis.line = element_line(linewidth = 0.8, color = "black"),
        axis.ticks = element_line(linewidth = 0.8, color = "black"),
        legend.position = "none",
        plot.title = element_blank()
    ) +
    labs(
        x = "Methods",
        y = "Matched-reference overlap proportion"
    )
ggsave(
    filename = file.path(fig_dir, "method_comparison_overlap.pdf"), 
    plot = p_hit, 
    width = 3.8,   # 宽度 3.8 英寸 (约 9.6 厘米)
    height = 4.8,  # 高度 4.8 英寸 (留出底部空间给 45 度倾斜的长名字)
    units = "in", 
    dpi = 300,
    bg = "white"
)


#个体结果特异性分析
methods_pers <- list(
    HyperNetWalk = "./results/HyperNetWalk/",
    DawnRank = "./results/DawnRank/"
)

CGC <- read.delim("./reference_dg/CGC_Tier1.tsv",
                             header = T,
                             as.is = TRUE)[, 1]
res_pers_all <- list()
ref_list <- list()
top_pers_methods <- list()
hit_nums <- list()
for (cancer in cancers){
  mut_data <- get_mut_data(file.path("./data/processed_data",cancer,"mut_data.tsv"),1)
  ref <- get_filter_ref(mut_data, CGC, N_coh = 500)
  ref_list[[cancer]] <- ref
  driver_c <- IntOGen_drivers[IntOGen_drivers$TCGA_CANCER %in% cancer, ]$SYMBOL
  res_pers_c <- list()
  for (i in 1:length(methods_pers)){
    method <- names(methods_pers)[i]
    res_file <- file.path(methods_pers[[i]],cancer,"genes_ranking.txt")
    if(!file.exists(res_file)){
      next
    }
    res <- read.table(
      res_file,
      sep = "\t",
      header = FALSE,
      stringsAsFactors = FALSE
    )
    prediction <- list()
    hit_num <- c() 
    for (j in 1:nrow(res)) {
      sample_j <- substr(res[j, 1], 1, 12)
      genes_ranking_j <- unlist(strsplit(res[j, 2], ","))
      prediction[[sample_j]] <- genes_ranking_j
      top_df <- data.frame(
        cancer = cancer,
        sample = sample_j,
        gene = genes_ranking_j[1:min(ref$N_pers, length(genes_ranking_j))],
        rank = 1:min(ref$N_pers, length(genes_ranking_j))
      )
      if (!is.null(top_pers_methods[[method]])) {
        top_pers_methods[[method]] <- rbind(top_pers_methods[[method]], top_df)
      } else {
        top_pers_methods[[method]] <- top_df
      }
      hit_num <- c(hit_num, sum(genes_ranking_j[1:min(ref$N_pers, length(genes_ranking_j))] %in% driver_c))
    }
    res_pers_c[[method]] <- prediction
    hit_nums[[cancer]][[method]] <- hit_num
  }
  res_pers_all[[cancer]] <- res_pers_c
}

hit_stats_list <- list()

for (cancer in names(hit_nums)) {
  for (method in names(hit_nums[[cancer]])) {
    hits_vector <- hit_nums[[cancer]][[method]]
    
    # 防止某些方法在特定癌种运行失败导致 vector 为空
    if(length(hits_vector) == 0) next
    
    # 统计均值
    mean_hit <- mean(hits_vector, na.rm = TRUE)
    # 统计命中率 (hit_num > 0 的样本比例)
    prop_hit <- sum(hits_vector > 0, na.rm = TRUE) / length(hits_vector)
    
    hit_stats_list[[length(hit_stats_list) + 1]] <- data.frame(
      Cancer = cancer,
      Method = method,
      Mean_Hit = mean_hit,
      Prop_Hit = prop_hit
    )
  }
}


hit_summary_df <- bind_rows(hit_stats_list)
print("=== 样本命中率 (Hit Rate > 0) 统计 ===")
print(hit_summary_df %>% select(Cancer, Method, Prop_Hit) %>% pivot_wider(names_from = Method, values_from = Prop_Hit))
method_levels <- c("HyperNetWalk", setdiff(unique(hit_summary_df$Method), "HyperNetWalk"))
hit_summary_df$Method <- factor(hit_summary_df$Method, levels = method_levels)

custom_colors <- c("HyperNetWalk" = "#D94841", 
                   "DawnRank"     = "#D9D9D9") 
# 绘制样本比例分组柱状图，并且在柱子上显示具体比例数值
p_hit_prop <- ggplot(hit_summary_df, aes(x = Cancer, y = Prop_Hit, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.6) +
  geom_text(
    aes(label = sprintf("%.1f%%", Prop_Hit * 100)),
    position = position_dodge(width = 0.8),
    vjust = -0.5,
    size = 3.5,
    fontface = "bold",
    color = "black"
  ) +
  scale_fill_manual(values = custom_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", face = "bold"),
    axis.text.y = element_text(color = "black", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    axis.line = element_line(linewidth = 0.8, color = "black"),
    axis.ticks = element_line(linewidth = 0.8, color = "black"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(face = "bold", size = 11),
    legend.key.size = unit(0.5, "cm")
  ) +
  labs(
    x = "Cancer Type",
    y = "Proportion of Samples with ≥1 Hit"
  )
ggsave(file.path(fig_dir, "Fig_Personalized_HitProportion_Barplot.pdf"), plot = p_hit_prop, width = 8, height = 5, units = "in", dpi = 300, bg = "white")

# 3. 绘制均值分组柱状图
p_hit_mean <- ggplot(hit_summary_df, aes(x = Cancer, y = Mean_Hit, fill = Method)) +
  geom_bar(
    stat = "identity", 
    position = position_dodge(width = 0.8), 
    width = 0.7, 
    color = "black", 
    linewidth = 0.6
  ) +
  scale_fill_manual(values = custom_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", face = "bold"),
    axis.text.y = element_text(color = "black", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    axis.line = element_line(linewidth = 0.8, color = "black"),
    axis.ticks = element_line(linewidth = 0.8, color = "black"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(face = "bold", size = 11),
    legend.key.size = unit(0.5, "cm")
  ) +
  labs(
    x = "Cancer Type",
    y = "Mean Hit Number per Sample\n(Top Personalized Predictions)"
  )

print(p_hit_mean)
ggsave(file.path(fig_dir, "Fig_Personalized_MeanHit_Barplot.pdf"), plot = p_hit_mean, width = 8, height = 5, units = "in", dpi = 300, bg = "white")

# reference-independent analysis
jaccard_matrices_coh <- lapply(names(res_coh), function(method) {
  cancer_types <- names(res_coh[[method]])
  jaccard_matrix <- matrix(0, nrow = length(cancer_types), ncol = length(cancer_types))
  rownames(jaccard_matrix) <- cancer_types
  colnames(jaccard_matrix) <- cancer_types
  for (i in 1:length(cancer_types)) {
    for (j in 1:length(cancer_types)) {
      genes_i <- rownames(res_coh[[method]][[cancer_types[i]]])[1:100]
      genes_j <- rownames(res_coh[[method]][[cancer_types[j]]])[1:100]
      jaccard_matrix[i, j] <- calculate_jaccard(genes_i, genes_j)
    }
  }
  return(jaccard_matrix)
})  
names(jaccard_matrices_coh) <- names(res_coh)

jaccard_matrices_pers <- lapply(names(top_pers_methods),function(method){
    cancer_types <- unique(top_pers_methods[[method]]$cancer)
    jaccard_matrix <- matrix(0, nrow = length(cancer_types), ncol = length(cancer_types))
    top_union_genes <- lapply(cancer_types, function(cancer) {
        unique(unlist(top_pers_methods[[method]]$gene[top_pers_methods[[method]]$cancer == cancer]))
    })
    rownames(jaccard_matrix) <- cancer_types
    colnames(jaccard_matrix) <- cancer_types
    for (i in 1:length(cancer_types)) {
        for (j in 1:length(cancer_types)) {
            genes_i <- top_union_genes[[i]]
            genes_j <- top_union_genes[[j]]
            jaccard_matrix[i, j] <- calculate_jaccard(genes_i, genes_j)
        }
    }
    return(jaccard_matrix)
})
names(jaccard_matrices_pers) <- names(top_pers_methods)

for (i in 1:length(names(res_coh))) {
    method <- names(res_coh)[i]
    
    jaccard_matrix_coh <- jaccard_matrices_coh[[method]]
    jaccard_matrix_pers <- jaccard_matrices_pers[[method]]
    
    if (!is.null(jaccard_matrix_pers)) {
        p_combined <- plot_combined_jaccard(jaccard_matrix_coh, jaccard_matrix_pers, method)
        ggsave(filename = paste0(fig_dir, method, "_jaccard_heatmap_combined.pdf"), 
               plot = p_combined, width = 7.0, height = 7.0, bg = "white")
    } else {
        p_half <- plot_half_jaccard(jaccard_matrix_coh, method)
        ggsave(filename = paste0(fig_dir, method, "_jaccard_heatmap.pdf"), 
               plot = p_half, width = 7.0, height = 7.0, bg = "white")
    }
}

# Pathway enrichment analysis
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% select(gs_name, gene_symbol)
# m_t2g <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME") %>% 
#     select(gs_name, gene_symbol)
# 把数据框转成 fgsea 喜欢的 list 格式
pathways_hallmark <- split(x = m_t2g$gene_symbol, f = m_t2g$gs_name)

names(pathways_hallmark) <- gsub("HALLMARK_", "", names(pathways_hallmark))

# pathways_reactome <- split(x = m_t2g$gene_symbol, f = m_t2g$gs_name)
# names(pathways_reactome) <- gsub("REACTOME_", "", names(pathways_reactome))

top_union_genes <- lapply(cancers, function(cancer){
    unique(unlist(top_pers_methods[["HyperNetWalk"]]$gene[top_pers_methods[["HyperNetWalk"]]$cancer == cancer]))
})
names(top_union_genes) <- cancers
pathway_ora_df <- data.frame()
pathway_gsea_df <- data.frame()
pathway_ora_pers_df <- data.frame()
for (cancer in cancers) {
    res <- res_coh[["HyperNetWalk"]][[cancer]]
    if(is.null(res)) next
    ora_res <- enricher(
        gene = rownames(res)[1:100],
        universe = rownames(res),
        TERM2GENE = m_t2g,
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
    )
    if (!is.null(ora_res) && nrow(ora_res@result) > 0) {
        df <- ora_res@result %>% 
            filter(p.adjust < 0.05) %>% 
            mutate(Cancer = cancer)
        pathway_ora_df <- rbind(pathway_ora_df, df)
    }

    ora_pers_res <- enricher(
        gene = top_union_genes[[cancer]],
        universe = rownames(res),
        TERM2GENE = m_t2g,
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
    )
    if (!is.null(ora_pers_res) && nrow(ora_pers_res@result) > 0) {
        df_pers <- ora_pers_res@result %>% 
            filter(p.adjust < 0.05) %>% 
            mutate(Cancer = cancer)
        pathway_ora_pers_df <- rbind(pathway_ora_pers_df, df_pers)
    }
    
    gene_scores <- res[, 1]
    names(gene_scores) <- rownames(res)
    set.seed(42)
    gene_scores <- gene_scores + runif(length(gene_scores), 0, 1e-10)
    gene_scores <- sort(gene_scores, decreasing = TRUE)
    fgsea_res <- fgseaMultilevel(
        pathways = pathways_hallmark,
        stats = gene_scores,
        minSize = 15,
        maxSize = 500,
        eps = 0,
        scoreType = "pos"
    )
    #只保留显著的通路结果
    fgsea_res <- fgsea_res %>% filter(padj < 0.05 & pval < 0.05)
    
    if(nrow(fgsea_res) > 0) {
        fgsea_res$Cancer <- cancer
        pathway_gsea_df <- rbind(pathway_gsea_df, as.data.frame(fgsea_res))
    }
}

pathway_ora_df <- pathway_ora_df %>%
    mutate(ID = gsub("HALLMARK_", "", ID)) %>%
    # mutate(ID = gsub("REACTOME_", "", ID)) %>%
    mutate(ID = gsub("_", " ", ID))

pathway_ora_pers_df <- pathway_ora_pers_df %>%
    mutate(ID = gsub("HALLMARK_", "", ID)) %>%
    # mutate(ID = gsub("REACTOME_", "", ID)) %>%
    mutate(ID = gsub("_", " ", ID))

pathway_gsea_df <- pathway_gsea_df %>%
    mutate(pathway = gsub("HALLMARK_", "", pathway)) %>%
    # mutate(pathway = gsub("REACTOME_", "", pathway)) %>%
    mutate(pathway = gsub("_", " ", pathway))

#对每个癌症类型只保留最显著的5条通路结果
# pathway_ora_df_filter <- pathway_ora_df %>%
#     group_by(Cancer) %>%
#     slice_min(order_by = p.adjust, n = 5, with_ties = FALSE) %>%
#     ungroup()

# pathway_ora_pers_df_filter <- pathway_ora_pers_df %>%
#     group_by(Cancer) %>%
#     slice_min(order_by = p.adjust, n = 5, with_ties = FALSE) %>%
#     ungroup()

# pathway_gsea_df_filter <- pathway_gsea_df %>%
#     group_by(Cancer) %>%
#     slice_min(order_by = padj, n = 5, with_ties = FALSE) %>%
#     ungroup()

pathway_freq_ora <- pathway_ora_df %>%
    group_by(ID) %>%
    summarise(n_cancers = n()) 

pathway_freq_ora_pers <- pathway_ora_pers_df %>%
    group_by(ID) %>%
    summarise(n_cancers = n())

pathway_freq_gsea <- pathway_gsea_df %>%
    group_by(pathway) %>%
    summarise(n_cancers = n())

# cohort top 100 ORA bubble plot
plot_ora_df <- pathway_ora_df%>%
    mutate(
        GeneRatio_num = as.numeric(sub("/.*", "", GeneRatio)) / 
                        as.numeric(sub(".*/", "", GeneRatio)),
        log10_padj = -log10(p.adjust)
    )
ora_pathway_order <- pathway_freq_ora %>% 
    arrange(n_cancers, ID) %>% 
    pull(ID)

plot_ora_df$ID <- factor(plot_ora_df$ID, levels = ora_pathway_order)
plot_ora_df$Cancer <- factor(plot_ora_df$Cancer, levels = sort(unique(plot_ora_df$Cancer)))
p_ora_cohort <- ggplot(plot_ora_df, aes(x = Cancer, y = ID)) +
    geom_point(aes(size = GeneRatio_num, color = log10_padj)) +
    scale_color_gradient(low = "#E6E8EA", high = "#D94841", name = expression(bold("-log"[10]*"(FDR)"))) +
    scale_size_continuous(range = c(2.5, 7), name = "Gene Ratio") +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) +
    
    theme_bw(base_size = 14) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
        axis.text.y = element_text(face = "bold", color = "black", size = 10),
        axis.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype = "dashed", color = "grey85", linewidth = 0.5),
        panel.border = element_rect(color = "black", linewidth = 1),
        plot.title = element_text(face = "bold", hjust = 0.5)
    ) +
    labs(title = "Cohort ORA (Top 100)")

print(p_ora_cohort)
ggsave(file.path(fig_dir,"Fig_ORA_Cohort_Bubble.pdf"), p_ora_cohort, width = 7.5, height = 6.5)

# cohort GSEA bubble plot
plot_gsea_df <- pathway_gsea_df %>%
    mutate(log10_padj = -log10(padj))

gsea_pathway_order <- pathway_freq_gsea %>% 
    arrange(n_cancers, pathway) %>% 
    pull(pathway)

plot_gsea_df$pathway <- factor(plot_gsea_df$pathway, levels = gsea_pathway_order)
plot_gsea_df$Cancer <- factor(plot_gsea_df$Cancer, levels = sort(unique(plot_gsea_df$Cancer)))

p_gsea <- ggplot(plot_gsea_df, aes(x = Cancer, y = pathway)) +
    geom_point(aes(size = NES, color = log10_padj)) +
    scale_color_gradient(low = "#E6E8EA", high = "#D94841", name = expression(bold("-log"[10]*"(FDR)"))) +
    scale_size_continuous(range = c(2.5, 7), name = "NES") +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) +  
    
    
    theme_bw(base_size = 14) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
        axis.text.y = element_text(face = "bold", color = "black", size = 10),
        axis.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype = "dashed", color = "grey85", linewidth = 0.5),
        panel.border = element_rect(color = "black", linewidth = 1),
        plot.title = element_text(face = "bold", hjust = 0.5)
    ) +
    labs(title = "Cohort GSEA (Full Ranking)")

print(p_gsea)
ggsave(file.path(fig_dir,"Fig_GSEA_Cohort_Bubble.pdf"), p_gsea, width = 7.5, height = 7.5)

# personalized ORA bubble plot
plot_ora_pers_df <- pathway_ora_pers_df %>%
    mutate(
        GeneRatio_num = as.numeric(sub("/.*", "", GeneRatio)) / 
                        as.numeric(sub(".*/", "", GeneRatio)),
        log10_padj = -log10(p.adjust)
    )

ora_pers_pathway_order <- pathway_freq_ora_pers %>% 
    arrange(n_cancers, ID) %>% 
    pull(ID)

plot_ora_pers_df$ID <- factor(plot_ora_pers_df$ID, levels = ora_pers_pathway_order)
plot_ora_pers_df$Cancer <- factor(plot_ora_pers_df$Cancer, levels = sort(unique(plot_ora_pers_df$Cancer)))
p_ora_pers <- ggplot(plot_ora_pers_df, aes(x = Cancer, y = ID)) +
    geom_point(aes(size = GeneRatio_num, color = log10_padj)) +
    
    # 💡 由于这个图非常密集，我把低值颜色稍微调深了一点，防止看不清
    scale_color_gradient(low = "#B3C3D1", high = "#D94841", name = expression(bold("-log"[10]*"(FDR)"))) +
    scale_size_continuous(range = c(2, 6.5), name = "Gene Ratio") +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) + 
    theme_bw(base_size = 14) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
        axis.text.y = element_text(face = "bold", color = "black", size = 9), # 通路多，字号稍微缩小一点
        axis.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype = "dashed", color = "grey85", linewidth = 0.5),
        panel.border = element_rect(color = "black", linewidth = 1),
        plot.title = element_text(face = "bold", hjust = 0.5)
    ) +
    labs(title = "Personalized Union ORA")

print(p_ora_pers)
ggsave(file.path(fig_dir,"Fig_ORA_Pers_Bubble.pdf"), p_ora_pers, width = 8, height = 8)

# TF clusteing analysis
TF_mats <- list()
exp_mats <- list()
TFs <- c()
exp_genes <- c()
cancer_vec <- c()
for (cancer in cancers){
    load(file.path ("./results/HyperNetWalk",cancer,"results.Rdata"))
    TF_mat_c <- objs$P_tf
    TF_mats[[cancer]] <- TF_mat_c
    TFs <- union(TFs, rownames(TF_mat_c))

    exp_file <- file.path("./data/processed_data",cancer,"exp_tpm_data.tsv")
    exp_data <- get_exp_data(exp_file,filter = FALSE)
    exp_data <- exp_data[,grep("-01A$|-01$",colnames(exp_data))]
    colnames(exp_data) <- substr(colnames(exp_data),1,12)
    exp_data <- exp_data[,which(colnames(exp_data) %in% colnames(TF_mat_c))]
    exp_mats[[cancer]] <- exp_data
    exp_genes <- union(exp_genes, rownames(exp_data))

    cancer_vec <- c(cancer_vec, rep(cancer, ncol(exp_data)))
}

TF_matrix <- do.call(cbind, TF_mats)
exp_matrix <- do.call(cbind, exp_mats)

valid_samples1 <- which(colSums(TF_matrix) > 0)
valid_TFs <- which(rowSums(TF_matrix) > 0)
TF_matrix_valid <- TF_matrix[valid_TFs, valid_samples1]
cancer_vec1 <- cancer_vec[valid_samples1]

TF_mat_norm <- apply(TF_matrix_valid, 2, function(x) x / sum(x))
TF_mat_norm <- as.matrix(TF_mat_norm)
# TF_mat_log <- log10(TF_mat_norm+1e-10)
tf_vars <- apply(TF_mat_norm, 1, var,na.rm = TRUE)
thre <- quantile(tf_vars, 0.05)
keep_tfs <- names(tf_vars)[tf_vars > thre]
TF_mat_sel <- TF_mat_norm[keep_tfs, ]
TF_mat_z <- t(scale(t(TF_mat_sel)))
TF_mat_z[is.na(TF_mat_z)] <- 0

keep_expr <- rowMeans(exp_matrix >= 1) >= 0.2
exp_matrix_valid <- exp_matrix[keep_expr, ]
cancer_vec2 <- cancer_vec

vars <- apply(exp_matrix_valid, 1, var,na.rm = TRUE)
keep_exp_genes <- names(sort(vars, decreasing = TRUE))[1:2000]
exp_matrix_sel <- exp_matrix_valid[keep_exp_genes, ]
exp_matrix_z <- t(scale(t(exp_matrix_sel)))
exp_matrix_z[is.na(exp_matrix_z)] <- 0

library(umap)
TF_dr <- run_dimred(TF_mat_z,cancer_vec1)
exp_dr <- run_dimred(exp_matrix_z,cancer_vec2)

my_colors_palette <- c(
  "#8DD3C7", # 薄荷绿
  "#FFFFB3", # 奶油黄
  "#BEBADA", # 浅紫灰
  "#FB8072", # 珊瑚粉
  "#80B1D3", # 天空蓝
  "#FDB462", # 杏色
  "#B3DE69", # 草嫩绿
  "#FCCDE5", # 淡樱花
  "#D9D9D9", # 高级灰
  "#BC80BD", # 香芋紫
  "#CCEBC5", # 浅茶绿
  "#FFED6F"  # 柠檬黄
)
if(length(unique(cancer_vec)) > length(my_colors_palette)){
  get_palette <- colorRampPalette(my_colors_palette)
  my_colors <- get_palette(length(unique(cancer_vec)))
} else {
  my_colors <- my_colors_palette[1:length(unique(cancer_vec))]
}
names(my_colors) <- sort(unique(cancer_vec))

p_tf <- plot_umap_by_group(TF_dr$umap, cancer_vec1, my_colors, "UMAP of TF embeddings",metrics = TF_dr$metrics_pca)
ggsave(file.path(fig_dir,"Fig_TF_UMAP.pdf"), p_tf, width = 6.5, height = 6.5, bg = "white")
p_exp <- plot_umap_by_group(exp_dr$umap, cancer_vec2, my_colors, "UMAP of tumor expression profiles",metrics = exp_dr$metrics_pca)
ggsave(file.path(fig_dir,"Fig_Exp_UMAP.pdf"), p_exp, width = 6.5, height = 6.5, bg = "white")

tsne_tf <- run_tsne_embedding(TF_mat_z, group_vec = cancer_vec1, seed = 921)
p_tsne_tf <- plot_tsne_embedding(tsne_tf$tsne_df, my_colors,title = paste0("t-SNE of TF embeddings (", tsne_tf$n_pcs, " PCs)"), metrics = tsne_tf$metrics_pca
)
ggsave(file.path(fig_dir,"Fig_TF_tSNE.pdf"), p_tsne_tf, width = 6.5, height = 6.5, bg = "white")

tsne_exp <- run_tsne_embedding(
    mat_z = exp_matrix_z,
    group_vec = cancer_vec2,
    seed = 921
)

p_tsne_exp <- plot_tsne_embedding(
    emb_df = tsne_exp$tsne_df,
    color_map = my_colors,
    title = paste0("t-SNE of expression profiles (", tsne_exp$n_pcs, " PCs)"),
    metrics = tsne_exp$metrics_pca
)
ggsave(file.path(fig_dir,"Fig_Exp_tSNE.pdf"), p_tsne_exp, width = 6.5, height = 6.5, bg = "white")

metric_summary <- data.frame(
    Data = c("TF_embedding", "Expression", "TF_embedding", "Expression"),
    Method = c("UMAP", "UMAP", "tSNE", "tSNE"),
    Silhouette_2D = c(
        TF_dr$metrics_umap$silhouette_mean,
        exp_dr$metrics_umap$silhouette_mean,
        tsne_tf$metrics_tsne$silhouette_mean,
        tsne_exp$metrics_tsne$silhouette_mean
    ),
    ARI_2D = c(
        TF_dr$metrics_umap$ARI,
        exp_dr$metrics_umap$ARI,
        tsne_tf$metrics_tsne$ARI,
        tsne_exp$metrics_tsne$ARI
    ),
    NMI_2D = c(
        TF_dr$metrics_umap$NMI,
        exp_dr$metrics_umap$NMI,
        tsne_tf$metrics_tsne$NMI,
        tsne_exp$metrics_tsne$NMI
    ),
    Silhouette_PCA = c(
        TF_dr$metrics_pca$silhouette_mean,
        exp_dr$metrics_pca$silhouette_mean,
        tsne_tf$metrics_pca$silhouette_mean,
        tsne_exp$metrics_pca$silhouette_mean
    ),
    ARI_PCA = c(
        TF_dr$metrics_pca$ARI,
        exp_dr$metrics_pca$ARI,
        tsne_tf$metrics_pca$ARI,
        tsne_exp$metrics_pca$ARI
    ),
    NMI_PCA = c(
        TF_dr$metrics_pca$NMI,
        exp_dr$metrics_pca$NMI,
        tsne_tf$metrics_pca$NMI,
        tsne_exp$metrics_pca$NMI
    )
)

write.csv(metric_summary, file.path(fig_dir, "Embedding_metrics_summary.csv"), row.names = FALSE)

# cancer-specific top TFs
# TF_matrix汇总的是全部癌症类型全部样本的TF embedding，现在对同一癌症类型的样本取平均，得到每个癌症类型的平均TF embedding，可以用TF_matrix_valid
avg_TF_mat <- sapply(unique(cancer_vec1), function(cancer) {
    col_idx <- which(cancer_vec1 == cancer)
    rowMeans(TF_matrix_valid[, col_idx, drop = FALSE])
})
colnames(avg_TF_mat) <- unique(cancer_vec1)
#对每个癌症类型分别取top 20的TFs
top_TFs_by_cancer <- apply(avg_TF_mat, 2, function(x) {
    names(sort(x, decreasing = TRUE))[1:50]
})
top_TFs <- unique(as.vector(top_TFs_by_cancer))
#用TF_matrix_valid计算one-vs-rest differential TF scores，取top 5的TFs
diff_TF_scores <- sapply(unique(cancer_vec1), function(cancer) {
    col_idx <- which(cancer_vec1 == cancer)
    rest_idx <- which(cancer_vec1 != cancer)
    (rowMeans(TF_matrix_valid[, col_idx, drop = FALSE]) - rowMeans(TF_matrix_valid[, rest_idx, drop = FALSE])) / (apply(TF_matrix_valid, 1, sd) + 1e-10)
})
top_TFs_diff <- sapply(unique(cancer_vec1), function(cancer) {
    names(sort(diff_TF_scores[, cancer], decreasing = TRUE))[1:5]
})
GRN_file <- "./data/RegNet_human_V2.txt"
GRN_edges <- fread(
    GRN_file,
    header = T,
    sep = "\t",
    select = 1:2
)
hallmark_m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% select(gs_name, gene_symbol)
raw_universe <- unique(GRN_edges$Target)
valid_universe <- raw_universe[!str_detect(raw_universe, "^hsa-miR-|^MIR|^LINC|^LOC")]
enr_res <- list()
for (cancer in cancers){
    top_TFs_c <- rownames(diff_TF_scores)[diff_TF_scores[, cancer] > 1]
    GRN_edges_c <- GRN_edges %>%
        filter(TF %in% top_TFs_c) %>%
        distinct()
    # table(GRN_edges_c$TF)
    # table(GRN_edges_c$Target)
    targets <- unique(GRN_edges_c$Target)
    targets <- targets[!str_detect(targets, "^hsa-miR-|^MIR|^LINC|^LOC")]
    ora_targets <- enricher(
        gene = targets,
        universe = valid_universe, # 🚨 强制使用网络背景
        TERM2GENE = hallmark_m_t2g,
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
    )
    if (!is.null(ora_targets) && nrow(ora_targets@result) > 0) {
        ora_targets_filter <- ora_targets@result %>%
            filter(p.adjust < 0.05)
        
        enr_res[[cancer]] <- ora_targets_filter$ID
        
        if (nrow(ora_targets_filter) > 0) {
            plot_hallmark_bar(
                df = ora_targets_filter,
                cancer = cancer,
                outdir = fig_dir
            )
        }
    } else {
        enr_res[[cancer]] <- character(0)
    }
}

