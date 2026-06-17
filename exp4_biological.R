source("model.R")
source("plot_formal.R")
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggvenn)
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
GRN_file <- "./data/RegNet_human_V2.txt"
GRN_edges <- fread(
    GRN_file,
    header = T,
    sep = "\t",
    select = 1:2
)
GRN_edges$weight <- 1
if(any(grep("/",GRN_edges$TF,fixed = TRUE))){
    needs_split <- GRN_edges[grep("/",GRN_edges$TF,fixed = TRUE)]
    GRN_edges <- GRN_edges[!grepl("/", GRN_edges$TF, fixed = TRUE)]

    split_rows <- needs_split[, {
        tfs <- strsplit(TF, "/", fixed = TRUE)[[1]]
        .(TF = tfs, Target = Target, weight = weight)
    }, by = 1:nrow(needs_split)]
    split_rows[,nrow:=NULL]
    GRN_edges <- rbind(GRN_edges,split_rows)
    GRN_edges <- unique(GRN_edges)
}
GRN_net <- graph_from_data_frame(GRN_edges,directed = TRUE)
dg <- igraph::degree(GRN_net, mode = "out")
fig_dir <- "./figs/exp4_biological/"
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
}
CGC_df <- read.delim("./reference_dg/CGC_Tier1.tsv",
                             header = T,
                             as.is = TRUE)

CGC <- read.delim("./reference_dg/CGC_Tier1.tsv",
                             header = T,
                             as.is = TRUE)[, 1]
methods <- list(
    HyperNetWalk = "./results/HyperNetWalk/",
    DriverRWH = "./results/DriverRWH/",
    DriverMP = "./results/DriverMP/",
    DawnRank = "./results/DawnRank/",
    PDRWH = "./results/PDRWH/",
    Frequency = "./results/Frequency/"
)
res_coh <- list()
res_pers <- list()
top_pers <- list()
mfs <- list()
for (cancer in cancers){
    mut_file <- mut_file <- file.path("./data/processed_data",cancer,"mut_data.tsv")
    mut_data <- get_mut_data(mut_file,1)
    colnames(mut_data) <- substr(colnames(mut_data), 1, 12)
    mut_genes <- rownames(mut_data)
    mfs[[cancer]] <- rowSums(mut_data)/ncol(mut_data)
    ref <- get_filter_ref(mut_data, CGC, N_coh = 500)
    for (method in names(methods)){
        res_pers_file <- file.path(methods[[method]],cancer,"genes_ranking.txt")
        if (file.exists(res_pers_file)) {
            res <- read.table(
                res_pers_file,
                sep = "\t",
                header = FALSE,
                stringsAsFactors = FALSE
            )
            prediction <- list()
            for (j in 1:nrow(res)) {
                sample_j <- substr(res[j, 1], 1, 12)
                genes_ranking_j <- unlist(strsplit(res[j, 2], ","))
                prediction[[sample_j]] <- genes_ranking_j
                top_df <- data.frame(
                    cancer = cancer,
                    sample = sample_j,
                    gene = genes_ranking_j[1:min(ref$N_pers,length(genes_ranking_j))],
                    rank = 1:min(ref$N_pers,length(genes_ranking_j))
                )
                if (!is.null(top_pers[[method]])) {
                    top_pers[[method]] <- rbind(top_pers[[method]], top_df)
                } else {
                    top_pers[[method]] <- top_df
                }
            }
            res_pers[[method]][[cancer]] <- prediction
        }
        res_coh_file <- file.path(methods[[method]],cancer,"genes_ranking_cohort.txt")
        resc <- read.table(
            res_coh_file,
            row.names = 1
        )
        mut_genes_sup <- setdiff(mut_genes,rownames(resc))
        if (length(mut_genes_sup) > 0) {
            padding <- matrix(0,
                                nrow = length(mut_genes_sup),
                                ncol = ncol(resc))
            rownames(padding) <- mut_genes_sup
            colnames(padding) <- colnames(resc)
            resc <- rbind(resc, padding)
        }
        res_coh[[method]][[cancer]] <- resc
    }
}

# 对每个方法的cohort—level结果评估绘制柱状图，每个柱子代表一个癌症类型的top 100基因，柱子按个数分成四份，分别是已知高频驱动基因、已知低频驱动基因、其它高频突变基因和其它低频突变基因
thre <- 0.02
df_mf_coh <- data.frame()
for (method in names(res_coh)){
    for (cancer in cancers){
        top_c <- rownames(res_coh[[method]][[cancer]])[1:100]
        num_known_high <- sum(top_c %in% CGC & mfs[[cancer]][top_c] >= thre)
        num_known_low <- sum(top_c %in% CGC & mfs[[cancer]][top_c] < thre)
        num_other_high <- sum(!(top_c %in% CGC) & mfs[[cancer]][top_c] >= thre)
        num_other_low <- sum(!(top_c %in% CGC) & mfs[[cancer]][top_c] < thre)
        df_mf_coh <- rbind(df_mf_coh, data.frame(
            method = method,
            cancer = cancer,
            category = c("Known high-frequency", "Known low-frequency", "Other high-frequency", "Other low-frequency"),
            count = c(num_known_high, num_known_low, num_other_high, num_other_low)
        ))
    }
}
pastel_colors <- c(
  "Known high-frequency" = "#FB8072",  # 珊瑚粉
  "Known low-frequency"  = "#FDB462",  # 杏色
  "Other high-frequency" = "#80B1D3",  # 天空蓝
  "Other low-frequency"  = "#8DD3C7"   # 薄荷绿
)
cat_levels <- c("Known high-frequency", "Known low-frequency", 
                "Other high-frequency", "Other low-frequency")
all_methods <- unique(df_mf_coh$method)
other_methods <- sort(setdiff(all_methods, "HyperNetWalk"))
method_order <- c("HyperNetWalk", other_methods)
df_mf_coh_clean <- df_mf_coh %>%
  mutate(
    method = factor(method, levels = method_order),
    category = factor(category, levels = rev(cat_levels)) 
  )
bar_coh <- list()
for (cancer in cancers) {
  df_cancer <- df_mf_coh_clean %>% filter(cancer == !!cancer)
  if(nrow(df_cancer) == 0) next
  p <- ggplot(df_cancer, aes(x = method, y = count, fill = category)) +
    geom_col(color = "white", linewidth = 0.2, width = 0.65) +
    
    scale_fill_manual(
      values = pastel_colors,
      breaks = cat_levels, 
      name = "Gene Category"
    ) +
    
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    
    labs(
      title = cancer,
      x = NULL, 
      y = "Number of Genes in Top 100"
    ) +
    
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black", size = 11),
      axis.text.y = element_text(face = "bold", color = "black", size = 11),
      axis.title.y = element_text(face = "bold", size = 12, margin = margin(r = 5)),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black", linewidth = 0.6)
      # 这里不需要管 legend.position，交由 patchwork 全局处理
    )
  ggsave(
    filename = file.path(fig_dir,"coh",paste0(cancer, "_mutation_frequency_top100coh.pdf")), 
    plot = p, 
    width = 7.5, 
    height = 5, 
    bg = "white"
  )
  bar_coh[[cancer]] <- p
}

combined_plot <- wrap_plots(bar_coh, ncol = 3, nrow = 4) +
  
  # plot_layout(guides = "collect") 会自动提取所有小图的图例，并合并为一个公共图例
  plot_layout(guides = "collect") +
  
  # 为整张大图添加全局主标题和图例位置控制
  plot_annotation(
    title = "Mutation Frequency Distribution in Top 100 Predicted Drivers",
    theme = theme(
      plot.title = element_text(face = "bold", size = 20, hjust = 0.5, margin = margin(b = 20)),
      # 将公共图例放在大图的底部（或者 "right" 放在右侧），底部通常更节省空间
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(size = 13)
    )
  )

ggsave(
  filename = file.path(fig_dir, "All_Cancers_Mutation_Frequency_Combined_top100coh.pdf"), 
  plot = combined_plot, 
  width = 12,    # 总宽度
  height = 16,   # 总高度
  bg = "white",
  device = cairo_pdf
)

# 对每个方法的personalized结果评估绘制柱状图，每一个柱子代表所有样本中第K名预测结果中已知高频驱动基因、已知低频驱动基因、其它高频突变基因和其它低频突变基因的数量
df_mf_pers <- data.frame()
for (method in names(top_pers)){
    top_pers_m <- top_pers[[method]]
    #对于top_pers_m中的每一行对应于一个样本中的一个基因，添加一列category，表示该基因是已知高频驱动基因、已知低频驱动基因、其它高频突变基因还是其它低频突变基因
    top_pers_m <- top_pers_m %>%
    mutate(
        mf = purrr::map2_dbl(cancer, gene, ~{
        v <- mfs[[.x]]
        if (is.null(v)) return(NA_real_)
        out <- v[.y]
        if (length(out) == 0 || is.na(out)) NA_real_ else as.numeric(out)
        }),
        category = case_when(
        gene %in% CGC & mf >= thre ~ "Known high-frequency",
        gene %in% CGC & mf <  thre ~ "Known low-frequency",
        !(gene %in% CGC) & mf >= thre ~ "Other high-frequency",
        !(gene %in% CGC) & mf <  thre ~ "Other low-frequency",
        TRUE ~ NA_character_
        )
    )
    # 对于每个癌症类型和每个排名位置，统计不同category的基因百分比
    df_mf_pers <- rbind(df_mf_pers, top_pers_m %>%
    group_by(cancer, rank) %>%
    summarise(
        method = method,
        known_high_freq = sum(category == "Known high-frequency", na.rm = TRUE) / n(),
        known_low_freq = sum(category == "Known low-frequency", na.rm = TRUE) / n(),
        other_high_freq = sum(category == "Other high-frequency", na.rm = TRUE) / n(),
        other_low_freq = sum(category == "Other low-frequency", na.rm = TRUE) / n()
    )) 
}

df_plot_data <- df_mf_pers %>%
  pivot_longer(
    cols = ends_with("_freq"),
    names_to = "category_raw",
    values_to = "prop"
  ) %>%
  mutate(
    category = case_when(
      category_raw == "known_high_freq" ~ "Known high-frequency",
      category_raw == "known_low_freq"  ~ "Known low-frequency",
      category_raw == "other_high_freq" ~ "Other high-frequency",
      category_raw == "other_low_freq"  ~ "Other low-frequency"
    ),
    category = factor(category, levels = rev(cat_levels)),
    plot_prop = if_else(method == "DawnRank", -prop, prop)
  )

mirrored_plots <- list()
for (canc in cancers) {
  df_cancer <- df_plot_data %>% filter(cancer == !!canc)
  if (nrow(df_cancer) == 0) next
  p <- ggplot(df_cancer, aes(x = rank, y = plot_prop, fill = category)) +
    geom_col(color = "white", linewidth = 0.3, width = 0.8) +
    geom_hline(yintercept = 0, color = "black", linewidth = 1.2) +
    scale_fill_manual(
      values = pastel_colors,
      breaks = cat_levels,
      name = "Gene Category"
    ) +
    scale_y_continuous(
      labels = function(x) sprintf("%.2f", abs(x)),
      limits = c(-1, 1), # 强制 Y 轴范围为 -1 到 1，保证绝对对称
      breaks = seq(-1, 1, by = 0.25) # 每 0.25 划一个刻度
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    labs(
      title = canc,
      x = "Rank",
      y = "Proportion"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
      axis.text = element_text(face = "bold", color = "black"),
      axis.title = element_text(face = "bold"),
      axis.line.x = element_blank(), 
      axis.ticks.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black", linewidth = 0.6),
      axis.ticks.y = element_line(color = "black", linewidth = 0.6)
    ) +
    annotate("text", x = max(df_cancer$rank) * 0.95, y = 0.9, 
             label = "HyperNetWalk", fontface = "bold", color = "#4D4D4D", hjust = 1, size = 5) +
    annotate("text", x = max(df_cancer$rank) * 0.95, y = -0.9, 
             label = "DawnRank", fontface = "bold", color = "#4D4D4D", hjust = 1, size = 5)
  
  mirrored_plots[[canc]] <- p
}
combined_mirrored <- wrap_plots(mirrored_plots, ncol = 3) + 
plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "Mirrored_Pers_TopK_Combined.pdf"), combined_mirrored, width = 14, height = 18)

df_hnw_plot <- df_mf_pers %>%
  filter(method == "HyperNetWalk") %>%
  pivot_longer(
    cols = ends_with("_freq"),
    names_to = "category_raw",
    values_to = "prop"
  ) %>%
  mutate(
    category = case_when(
      category_raw == "known_high_freq" ~ "Known high-frequency",
      category_raw == "known_low_freq"  ~ "Known low-frequency",
      category_raw == "other_high_freq" ~ "Other high-frequency",
      category_raw == "other_low_freq"  ~ "Other low-frequency"
    ),
    category = factor(category, levels = rev(cat_levels))
  )

for (canc in unique(df_hnw_plot$cancer)) {
  
  df_cancer <- df_hnw_plot %>% filter(cancer == !!canc)
  
  p <- ggplot(df_cancer, aes(x = rank, y = prop, fill = category)) +
    geom_col(color = "white", linewidth = 0.3, width = 0.8) +
    
    scale_fill_manual(
      values = pastel_colors,
      breaks = cat_levels,
      name = "Gene Category"
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.05)),
      labels = scales::percent_format(),
      name = "Proportion of Samples"
    ) +
    scale_x_continuous(
      breaks = seq(1, max(df_cancer$rank), by = 2),
      expand = c(0.02, 0)
    ) +
    
    labs(
      title = canc,
      x = "Rank"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16, margin = margin(b = 15)),
      axis.text = element_text(face = "bold", color = "black"),
      axis.title = element_text(face = "bold"),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.8),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      panel.grid.major.y = element_line(color = "grey95")
    )
  ggsave(
    filename = file.path(fig_dir, "pers", paste0("HyperNetWalk_Pers_Profile_", canc, ".pdf")),
    plot = p,
    width = 7,
    height = 4,
    bg = "white"
  )
}

top_hnw_coh_df <- data.frame()
top_hnw_pers_df <- data.frame()
venn_coh_plots <- list()
venn_pers_plots <- list()
novel_res_coh <- list()
novel_res_pers <- list()
for (cancer in cancers){
    top_hnw_coh <- rownames(res_coh[["HyperNetWalk"]][[cancer]])[1:100]
    top_dawnrank_coh <- rownames(res_coh[["DawnRank"]][[cancer]])[1:100]
    top_drivermp_coh <- rownames(res_coh[["DriverMP"]][[cancer]])[1:100]
    top_driverrwh_coh <- rownames(res_coh[["DriverRWH"]][[cancer]])[1:100]
    list_coh <- list(
        HyperNetWalk = top_hnw_coh,
        DawnRank = top_dawnrank_coh,
        DriverMP = top_drivermp_coh,
        DriverRWH = top_driverrwh_coh
    )
    p_coh <- plot_venn_custom(list_coh, title = paste0(cancer, " Cohort Top 100"),colors = c("#FB8072", "#80B1D3", "#B3DE69", "#FDB462"))
    ggsave(
        filename = file.path(fig_dir, "coh", paste0(cancer, "_Venn_top100coh.pdf")),
        plot = p_coh,
        width = 6.5,
        height = 6,
        bg = "white"
    )
    venn_coh_plots[[cancer]] <- p_coh
    df_hnw_coh <- data.frame(
        gene = top_hnw_coh,
        cancer = cancer,
        rank = 1:length(top_hnw_coh),
        mf = mfs[[cancer]][top_hnw_coh],
        known = top_hnw_coh %in% CGC,
        dawnrank = match(top_hnw_coh, rownames(res_coh[["DawnRank"]][[cancer]])),
        drivermp = match(top_hnw_coh, rownames(res_coh[["DriverMP"]][[cancer]])),
        driverrwh = match(top_hnw_coh, rownames(res_coh[["DriverRWH"]][[cancer]]))
    )
    top_hnw_coh_df <- rbind(top_hnw_coh_df, df_hnw_coh)
    # 对于每个癌症类型，找出HyperNetWalk的top 100中在其它方法中排名在100以外的基因，认为这些基因是HyperNetWalk独有的预测结果，只保留gene名
    novel_res_coh[[cancer]] <- df_hnw_coh %>%
        filter(dawnrank > 100 & drivermp > 100 & driverrwh > 100) %>%
        select(gene) %>%
        pull()
    
    # 对于每个癌症类型，对于top_pers[["HyperNetWalk"]]中cancer对应的行，合并gene列的元素相同的行并计数，保留cancer列、gene列和count列，按照count列降序排列
    top_hnw_pers <- top_pers[["HyperNetWalk"]] %>%
        filter(cancer == !!cancer) %>%
        filter(rank <= 5) %>%
        group_by(gene) %>%
        summarise(count = n()) %>%
        arrange(desc(count)) %>%
        mutate(mf = mfs[[cancer]][gene],
               known = gene %in% CGC)
    top_dawnrank_pers <- top_pers[["DawnRank"]] %>%
        filter(cancer == !!cancer) %>%
        filter(rank <= 5) %>%
        select(gene) %>%
        unique() %>%
        pull()
    top_pdrwh_pers <- top_pers[["PDRWH"]] %>%
        filter(cancer == !!cancer) %>%
        filter(rank <= 5) %>%
        select(gene) %>%
        unique() %>%
        pull()
    list_pers <- list(
        HyperNetWalk = top_hnw_pers$gene,
        DawnRank = top_dawnrank_pers,
        PDRWH = top_pdrwh_pers
    )
    p_pers <- plot_venn_custom(list_pers, title = paste0(cancer, " Personalized Top 5 Union"), colors = c("#FB8072", "#80B1D3", "#CCEBC5"))
    ggsave(
        filename = file.path(fig_dir, "pers", paste0(cancer, "_Venn_top5pers.pdf")),
        plot = p_pers,
        width = 6,
        height = 6,
        bg = "white"
    )
    venn_pers_plots[[cancer]] <- p_pers
    top_hnw_pers <- top_hnw_pers %>%
        mutate(cancer = cancer,
               dawnrank_top5 = gene %in% top_dawnrank_pers,
               pdrwh_top5 = gene %in% top_pdrwh_pers)
              
    top_hnw_pers_df <- rbind(top_hnw_pers_df, top_hnw_pers)
    # 对于每个癌症类型，找出HyperNetWalk的top 5中在DawnRank和PDRWH的top 5中都没有出现的基因，认为这些基因是HyperNetWalk独有的预测结果，只保留gene名
    novel_res_pers[[cancer]] <- top_hnw_pers %>%
        filter(!dawnrank_top5 & !pdrwh_top5) %>%
        select(gene) %>%
        pull()
}

for (cancer in cancers){
    message("For cancer ", cancer, "...")
    # 打印novel基因的数量，打印novel基因中已知驱动基因的数量
    message("Novel genes in cohort-level: ", length(novel_res_coh[[cancer]]))
    message("Known drivers in novel cohort-level genes: ", sum(novel_res_coh[[cancer]] %in% CGC))
    message("Novel genes in personalized-level: ", length(novel_res_pers[[cancer]]))
    message("Known drivers in novel personalized-level genes: ", sum(novel_res_pers[[cancer]] %in% CGC))
}

combined_venn_coh <- wrap_plots(venn_coh_plots, ncol = 3, nrow = 4) +
    plot_annotation(
        title = "Cohort-level Top 100 Predicted Drivers Overlap",
        theme = theme(
            plot.title = element_text(face = "bold", size = 24, hjust = 0.5, margin = margin(b = 20))
        )
    )

ggsave(
    filename = file.path(fig_dir, "All_Cancers_Venn_Cohort_Combined.pdf"),
    plot = combined_venn_coh,
    width = 15,
    height = 20,
    bg = "white",
    device = cairo_pdf # 🌟 强烈建议使用 cairo_pdf，处理 Venn 图的半透明重叠部分效果远超默认 PDF 引擎
)

combined_venn_pers <- wrap_plots(venn_pers_plots, ncol = 3, nrow = 4) +
    plot_annotation(
        title = "Personalized-level Top 5 Union Drivers Overlap",
        theme = theme(
            plot.title = element_text(face = "bold", size = 22, hjust = 0.5, margin = margin(b = 20))
        )
    )

ggsave(
    filename = file.path(fig_dir, "All_Cancers_Venn_Pers_Combined.pdf"),
    plot = combined_venn_pers,
    width = 12,
    height = 16,
    bg = "white",
    device = cairo_pdf
)

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
CGC2_df <- read.delim("./reference_dg/CGC_Tier2.tsv",
                             header = T,
                             as.is = TRUE)
NCG_df <- read.delim("./reference_dg/NCG.tsv",
                             header = T,
                             as.is = TRUE)
CM_df <- read.delim("./reference_dg/cancermine_collated.tsv",
                             header = T,
                             as.is = TRUE)

# BRCA
top_hnw_coh_brca <- top_hnw_coh_df %>% filter(cancer == "BRCA") %>% arrange(desc(mf))
novel_coh_brca <- top_hnw_coh_brca %>% 
  filter(gene %in% novel_res_coh[["BRCA"]])
top_hnw_pers_brca <- top_hnw_pers_df %>% filter(cancer == "BRCA") %>% arrange(desc(count))
novel_pers_brca <- top_hnw_pers_brca %>% 
  filter(gene %in% novel_res_pers[["BRCA"]])

cgc1_brca <- CGC_df[grepl("breast", CGC_df$Tumour.Types.Somatic., ignore.case = TRUE), ]$Gene.Symbol
levelA_brca <- intersect(novel_coh_brca$gene, cgc1_brca)
levelB_brca <- setdiff(intersect(novel_coh_brca$gene, CGC),levelA_brca)
levelA_brca_pers <- intersect(novel_pers_brca$gene, cgc1_brca)
levelB_brca_pers <- setdiff(intersect(novel_pers_brca$gene, CGC),levelA_brca_pers)

cgc2_brca <- CGC2_df[grepl("breast", CGC2_df$Tumour.Types.Somatic., ignore.case = TRUE), ]$Gene.Symbol
levelB_brca <- union(levelB_brca, intersect(novel_coh_brca$gene, cgc2_brca))
levelC_brca <- intersect(novel_coh_brca$gene, CGC2_df$Gene.Symbol)
levelB_brca_pers <- union(levelB_brca_pers, intersect(novel_pers_brca$gene, cgc2_brca))
levelC_brca_pers <- intersect(novel_pers_brca$gene, CGC2_df$Gene.Symbol)

ncg_brca <- unique(NCG_df[grepl("breast", NCG_df$cancer_type, ignore.case = TRUE), ]$symbol)
levelB_brca <- union(levelB_brca, intersect(novel_coh_brca$gene, ncg_brca))
levelC_brca <- union(levelC_brca, intersect(novel_coh_brca$gene, NCG_df$symbol))
levelB_brca_pers <- union(levelB_brca_pers, intersect(novel_pers_brca$gene, ncg_brca))
levelC_brca_pers <- union(levelC_brca_pers, intersect(novel_pers_brca$gene, NCG_df$symbol))

brca_cm_terms <- c(
  "breast cancer", "breast carcinoma", "breast adenocarcinoma", 
  "sporadic breast cancer", "male breast cancer", "female breast cancer",
  "breast ductal carcinoma", "invasive ductal carcinoma", 
  "breast lobular carcinoma", "invasive lobular carcinoma",
  "luminal breast carcinoma A", "luminal breast carcinoma B", 
  "estrogen-receptor negative breast cancer",
  "inflammatory breast carcinoma", "breast secretory carcinoma", 
  "mammary Paget's disease", "breast metaplastic carcinoma", 
  "breast angiosarcoma", "breast malignant phyllodes tumor", 
  "breast adenomyoepithelioma"
)
cm_brca_genes <- CM_df %>%
  filter(cancer_normalized %in% brca_cm_terms) %>%
  # 最好加上这一步，过滤掉只有 1 篇文献支持的孤证，保证数据的硬核程度
  # filter(citation_count >= 2) %>% 
  pull(gene_normalized) %>%
  unique()
levelB_brca <- union(levelB_brca, intersect(novel_coh_brca$gene, cm_brca_genes))
levelC_brca <- union(levelC_brca,intersect(novel_coh_brca$gene, CM_df$gene_normalized))
levelB_brca_pers <- union(levelB_brca_pers, intersect(novel_pers_brca$gene, cm_brca_genes))
levelC_brca_pers <- union(levelC_brca_pers, intersect(novel_pers_brca$gene, CM_df$gene_normalized))

levelB_brca <- setdiff(levelB_brca, levelA_brca)
levelC_brca <- setdiff(levelC_brca, union(levelA_brca, levelB_brca))
levelB_brca_pers <- setdiff(levelB_brca_pers, levelA_brca_pers)
levelC_brca_pers <- setdiff(levelC_brca_pers, union(levelA_brca_pers, levelB_brca_pers))
length(levelA_brca)
length(levelB_brca)
length(levelC_brca)
setdiff(novel_coh_brca$gene, c(levelA_brca, levelB_brca, levelC_brca))
length(levelA_brca_pers)
length(levelB_brca_pers)
length(levelC_brca_pers)
setdiff(novel_pers_brca$gene, c(levelA_brca_pers, levelB_brca_pers, levelC_brca_pers))
#在novel_coh_brca中添加一列"evidence_level"，如果gene在levelA_brca中则为"A"，如果gene在levelB_brca中则为"B"，如果gene在levelC_brca中则为"C"，否则为"D";再加一列out_degree
novel_coh_brca <- novel_coh_brca %>%
  mutate(out_degree = dg[gene]) %>%
  mutate(evidence_level = case_when(
    gene %in% levelA_brca ~ "A",
    gene %in% levelB_brca ~ "B",
    gene %in% levelC_brca ~ "C",
    TRUE ~ "D"
  )) 

write.csv(novel_coh_brca, file.path(fig_dir, "coh", "BRCA_novel_coh_genes.csv"), row.names = FALSE)

novel_coh_df <- data.frame(
  cancer = "BRCA",
  num_novel = length(novel_coh_brca$gene),
  num_levelA = length(levelA_brca),
  num_levelB = length(levelB_brca),
  num_levelC = length(levelC_brca),
  num_levelD = length(setdiff(novel_coh_brca$gene, c(levelA_brca, levelB_brca, levelC_brca)))
)

novel_pers_df <- data.frame(
  cancer = "BRCA",
  num_novel = length(novel_pers_brca$gene),
  num_levelA = length(levelA_brca_pers),
  num_levelB = length(levelB_brca_pers),
  num_levelC = length(levelC_brca_pers),
  num_levelD = length(setdiff(novel_pers_brca$gene, c(levelA_brca_pers, levelB_brca_pers, levelC_brca_pers)))
)

#COAD
top_hnw_coh_coad <- top_hnw_coh_df %>% filter(cancer == "COAD") %>% arrange(desc(mf))
novel_coh_coad <- top_hnw_coh_coad %>% 
  filter(gene %in% novel_res_coh[["COAD"]])
top_hnw_pers_coad <- top_hnw_pers_df %>% filter(cancer == "COAD") %>% arrange(desc(count))
novel_pers_coad <- top_hnw_pers_coad %>% 
  filter(gene %in% novel_res_pers[["COAD"]])

cgc1_coad <- CGC_df %>%
  filter(grepl("colorectal|colon|large intestine|\\bCRC\\b", Tumour.Types.Somatic., ignore.case = TRUE)) %>%
  pull(Gene.Symbol) %>%
  unique()
levelA_coad <- intersect(novel_coh_coad$gene, cgc1_coad)
levelB_coad <- setdiff(intersect(novel_coh_coad$gene, CGC),levelA_coad)
levelA_coad_pers <- intersect(novel_pers_coad$gene, cgc1_coad)
levelB_coad_pers <- setdiff(intersect(novel_pers_coad$gene, CGC),levelA_coad_pers)

cgc2_coad <- CGC2_df %>%
  filter(grepl("colorectal|colon|large intestine|\\bCRC\\b", Tumour.Types.Somatic., ignore.case = TRUE)) %>%
  pull(Gene.Symbol) %>%
  unique()
levelB_coad <- union(levelB_coad, intersect(novel_coh_coad$gene, cgc2_coad))
levelC_coad <- intersect(novel_coh_coad$gene, CGC2_df$Gene.Symbol)
levelB_coad_pers <- union(levelB_coad_pers, intersect(novel_pers_coad$gene, cgc2_coad))
levelC_coad_pers <- intersect(novel_pers_coad$gene, CGC2_df$Gene.Symbol)

ncg_coad <- NCG_df %>%
  filter(cancer_type == "colorectal_adenocarcinoma") %>%
  pull(symbol) %>%
  unique()
levelB_coad <- union(levelB_coad, intersect(novel_coh_coad$gene, ncg_coad))
levelC_coad <- union(levelC_coad, intersect(novel_coh_coad$gene, NCG_df$symbol))
levelB_coad_pers <- union(levelB_coad_pers, intersect(novel_pers_coad$gene, ncg_coad))
levelC_coad_pers <- union(levelC_coad_pers, intersect(novel_pers_coad$gene, NCG_df$symbol))

coad_cm_terms <- c("colon cancer", "colon carcinoma", "colon adenocarcinoma", "colorectal cancer", "colorectal carcinoma", "colorectal adenocarcinoma", "rectum cancer", "rectum adenocarcinoma")
cm_coad_genes <- CM_df %>%
  filter(cancer_normalized %in% coad_cm_terms) %>%
  # filter(citation_count >= 2) %>% 
  pull(gene_normalized) %>%
  unique()
levelB_coad <- union(levelB_coad, intersect(novel_coh_coad$gene, cm_coad_genes))
levelC_coad <- union(levelC_coad, intersect(novel_coh_coad$gene, CM_df$gene_normalized))
levelB_coad_pers <- union(levelB_coad_pers, intersect(novel_pers_coad$gene, cm_coad_genes))
levelC_coad_pers <- union(levelC_coad_pers, intersect(novel_pers_coad$gene, CM_df$gene_normalized))
levelB_coad <- setdiff(levelB_coad, levelA_coad)
levelC_coad <- setdiff(levelC_coad, union(levelA_coad, levelB_coad))
levelB_coad_pers <- setdiff(levelB_coad_pers, levelA_coad_pers)
levelC_coad_pers <- setdiff(levelC_coad_pers, union(levelA_coad_pers, levelB_coad_pers))
length(levelA_coad)
length(levelB_coad)
length(levelC_coad)
setdiff(novel_coh_coad$gene, c(levelA_coad, levelB_coad, levelC_coad))
length(levelA_coad_pers)
length(levelB_coad_pers)
length(levelC_coad_pers)
setdiff(novel_pers_coad$gene, c(levelA_coad_pers, levelB_coad_pers, levelC_coad_pers))
novel_coh_coad <- novel_coh_coad %>%
  mutate(out_degree = dg[gene]) %>%
  mutate(evidence_level = case_when(
    gene %in% levelA_coad ~ "A",
    gene %in% levelB_coad ~ "B",
    gene %in% levelC_coad ~ "C",
    TRUE ~ "D"
  ))
write.csv(novel_coh_coad, file.path(fig_dir, "coh", "COAD_novel_coh_genes.csv"), row.names = FALSE)
novel_coh_df <- rbind(novel_coh_df, data.frame(
  cancer = "COAD",
  num_novel = length(novel_coh_coad$gene),
  num_levelA = length(levelA_coad),
  num_levelB = length(levelB_coad),
  num_levelC = length(levelC_coad),
  num_levelD = length(setdiff(novel_coh_coad$gene, c(levelA_coad, levelB_coad, levelC_coad)))
))

novel_pers_df <- rbind(novel_pers_df, data.frame(
  cancer = "COAD",
  num_novel = length(novel_pers_coad$gene),
  num_levelA = length(levelA_coad_pers),
  num_levelB = length(levelB_coad_pers),
  num_levelC = length(levelC_coad_pers),
  num_levelD = length(setdiff(novel_pers_coad$gene, c(levelA_coad_pers, levelB_coad_pers, levelC_coad_pers)))
))  

#HNSC
top_hnw_coh_hnsc <- top_hnw_coh_df %>% filter(cancer == "HNSC") %>% arrange(desc(mf))
novel_coh_hnsc <- top_hnw_coh_hnsc %>% 
  filter(gene %in% novel_res_coh[["HNSC"]])
top_hnw_pers_hnsc <- top_hnw_pers_df %>% filter(cancer == "HNSC") %>% arrange(desc(count))
novel_pers_hnsc <- top_hnw_pers_hnsc %>% 
  filter(gene %in% novel_res_pers[["HNSC"]])

cgc1_hnsc <- CGC_df %>%
  filter(grepl("HNSCC|head and neck|head-neck|oral|tongue|larynx|pharynx", Tumour.Types.Somatic., ignore.case = TRUE)) %>%
  pull(Gene.Symbol) %>%
  unique()
levelA_hnsc <- intersect(novel_coh_hnsc$gene, cgc1_hnsc)
levelB_hnsc <- setdiff(intersect(novel_coh_hnsc$gene, CGC),levelA_hnsc)
levelA_hnsc_pers <- intersect(novel_pers_hnsc$gene, cgc1_hnsc)
levelB_hnsc_pers <- setdiff(intersect(novel_pers_hnsc$gene, CGC),levelA_hnsc_pers)

cgc2_hnsc <- CGC2_df %>%
  filter(grepl("HNSCC|head and neck|head-neck|oral|tongue|larynx|pharynx", Tumour.Types.Somatic., ignore.case = TRUE)) %>%
  pull(Gene.Symbol) %>%
  unique()
levelB_hnsc <- union(levelB_hnsc, intersect(novel_coh_hnsc$gene, cgc2_hnsc))
levelC_hnsc <- intersect(novel_coh_hnsc$gene, CGC2_df$Gene.Symbol)
levelB_hnsc_pers <- union(levelB_hnsc_pers, intersect(novel_pers_hnsc$gene, cgc2_hnsc))
levelC_hnsc_pers <- intersect(novel_pers_hnsc$gene, CGC2_df$Gene.Symbol)

ncg_hnsc <- NCG_df %>%
  filter(cancer_type == "squamous_head_and_neck_cancer") %>%
  pull(symbol) %>%
  unique()
levelB_hnsc <- union(levelB_hnsc, intersect(novel_coh_hnsc$gene, ncg_hnsc))
levelC_hnsc <- union(levelC_hnsc, intersect(novel_coh_hnsc$gene, NCG_df$symbol))
levelB_hnsc_pers <- union(levelB_hnsc_pers, intersect(novel_pers_hnsc$gene, ncg_hnsc))
levelC_hnsc_pers <- union(levelC_hnsc_pers, intersect(novel_pers_hnsc$gene, NCG_df$symbol))

hnsc_cm_terms <- c("head and neck cancer", "head and neck carcinoma", "head and neck squamous cell carcinoma", "oral squamous cell carcinoma", "tongue cancer", "tongue carcinoma", "tongue squamous cell carcinoma", "larynx cancer", "laryngeal carcinoma", "laryngeal squamous cell carcinoma", "oropharynx cancer", "hypopharynx cancer", "tonsil cancer", "lip cancer")
cm_hnsc_genes <- CM_df %>%
  filter(cancer_normalized %in% hnsc_cm_terms) %>%
  # filter(citation_count >= 2) %>% 
  pull(gene_normalized) %>%
  unique()
levelB_hnsc <- union(levelB_hnsc, intersect(novel_coh_hnsc$gene, cm_hnsc_genes))
levelC_hnsc <- union(levelC_hnsc, intersect(novel_coh_hnsc$gene, CM_df$gene_normalized))
levelB_hnsc_pers <- union(levelB_hnsc_pers, intersect(novel_pers_hnsc$gene, cm_hnsc_genes))
levelC_hnsc_pers <- union(levelC_hnsc_pers, intersect(novel_pers_hnsc$gene, CM_df$gene_normalized))
levelB_hnsc <- setdiff(levelB_hnsc, levelA_hnsc)
levelC_hnsc <- setdiff(levelC_hnsc, union(levelA_hnsc, levelB_hnsc))
levelB_hnsc_pers <- setdiff(levelB_hnsc_pers, levelA_hnsc_pers)
levelC_hnsc_pers <- setdiff(levelC_hnsc_pers, union(levelA_hnsc_pers, levelB_hnsc_pers))
length(levelA_hnsc)
length(levelB_hnsc)
length(levelC_hnsc)
setdiff(novel_coh_hnsc$gene, c(levelA_hnsc, levelB_hnsc, levelC_hnsc))
length(levelA_hnsc_pers)
length(levelB_hnsc_pers)
length(levelC_hnsc_pers)
setdiff(novel_pers_hnsc$gene, c(levelA_hnsc_pers, levelB_hnsc_pers, levelC_hnsc_pers))
novel_coh_hnsc <- novel_coh_hnsc %>%
  mutate(out_degree = dg[gene]) %>%
  mutate(evidence_level = case_when(
    gene %in% levelA_hnsc ~ "A",
    gene %in% levelB_hnsc ~ "B",
    gene %in% levelC_hnsc ~ "C",
    TRUE ~ "D"
  ))
write.csv(novel_coh_hnsc, file.path(fig_dir, "coh", "HNSC_novel_coh_genes.csv"), row.names = FALSE)
novel_coh_df <- rbind(novel_coh_df, data.frame(
  cancer = "HNSC",
  num_novel = length(novel_coh_hnsc$gene),
  num_levelA = length(levelA_hnsc),
  num_levelB = length(levelB_hnsc),
  num_levelC = length(levelC_hnsc),
  num_levelD = length(setdiff(novel_coh_hnsc$gene, c(levelA_hnsc, levelB_hnsc, levelC_hnsc)))
))  
novel_pers_df <- rbind(novel_pers_df, data.frame(
  cancer = "HNSC",
  num_novel = length(novel_pers_hnsc$gene),
  num_levelA = length(levelA_hnsc_pers),
  num_levelB = length(levelB_hnsc_pers),
  num_levelC = length(levelC_hnsc_pers),
  num_levelD = length(setdiff(novel_pers_hnsc$gene, c(levelA_hnsc_pers, levelB_hnsc_pers, levelC_hnsc_pers)))
))

#KIRC
top_hnw_coh_kirc <- top_hnw_coh_df %>% filter(cancer == "KIRC") %>% arrange(desc(mf))
novel_coh_kirc <- top_hnw_coh_kirc %>% 
  filter(gene %in% novel_res_coh[["KIRC"]])
top_hnw_pers_kirc <- top_hnw_pers_df %>% filter(cancer == "KIRC") %>% arrange(desc(count))
novel_pers_kirc <- top_hnw_pers_kirc %>% 
  filter(gene %in% novel_res_pers[["KIRC"]])

cgc1_kirc <- CGC_df %>%
  filter(grepl("CCRCC|clear cell renal|kidney|\\brenal\\b", Tumour.Types.Somatic., ignore.case = TRUE)) %>%
  pull(Gene.Symbol) %>%
  unique()
levelA_kirc <- intersect(novel_coh_kirc$gene, cgc1_kirc)
levelB_kirc <- setdiff(intersect(novel_coh_kirc$gene, CGC),levelA_kirc)
levelA_kirc_pers <- intersect(novel_pers_kirc$gene, cgc1_kirc)
levelB_kirc_pers <- setdiff(intersect(novel_pers_kirc$gene, CGC),levelA_kirc_pers) 

cgc2_kirc <- CGC2_df %>%
  filter(grepl("CCRCC|clear cell renal|kidney|\\brenal\\b", Tumour.Types.Somatic., ignore.case = TRUE)) %>%
  pull(Gene.Symbol) %>%
  unique()
levelB_kirc <- union(levelB_kirc, intersect(novel_coh_kirc$gene, cgc2_kirc))
levelC_kirc <- intersect(novel_coh_kirc$gene, CGC2_df$Gene.Symbol)
levelB_kirc_pers <- union(levelB_kirc_pers, intersect(novel_pers_kirc$gene, cgc2_kirc))
levelC_kirc_pers <- intersect(novel_pers_kirc$gene, CGC2_df$Gene.Symbol)

ncg_kirc <- NCG_df %>%
  filter(cancer_type %in% c("clear_cell_renal_cancer", "renal_cancer_(all_histologies)")) %>%
  pull(symbol) %>%
  unique()
levelB_kirc <- union(levelB_kirc, intersect(novel_coh_kirc$gene, ncg_kirc))
levelC_kirc <- union(levelC_kirc, intersect(novel_coh_kirc$gene, NCG_df$symbol))
levelB_kirc_pers <- union(levelB_kirc_pers, intersect(novel_pers_kirc$gene, ncg_kirc))
levelC_kirc_pers <- union(levelC_kirc_pers, intersect(novel_pers_kirc$gene, NCG_df$symbol))

kirc_cm_terms <- c("clear cell renal cell carcinoma", "kidney cancer", "renal carcinoma", "renal cell carcinoma")
cm_kirc_genes <- CM_df %>%
  filter(cancer_normalized %in% kirc_cm_terms) %>%
  # filter(citation_count >= 2) %>% 
  pull(gene_normalized) %>%
  unique()
levelB_kirc <- union(levelB_kirc, intersect(novel_coh_kirc$gene, cm_kirc_genes))
levelC_kirc <- union(levelC_kirc, intersect(novel_coh_kirc$gene, CM_df$gene_normalized))
levelB_kirc_pers <- union(levelB_kirc_pers, intersect(novel_pers_kirc$gene, cm_kirc_genes))
levelC_kirc_pers <- union(levelC_kirc_pers, intersect(novel_pers_kirc$gene, CM_df$gene_normalized))
levelB_kirc <- setdiff(levelB_kirc, levelA_kirc)
levelC_kirc <- setdiff(levelC_kirc, union(levelA_kirc, levelB_kirc))
levelB_kirc_pers <- setdiff(levelB_kirc_pers, levelA_kirc_pers)
levelC_kirc_pers <- setdiff(levelC_kirc_pers, union(levelA_kirc_pers, levelB_kirc_pers))
length(levelA_kirc)
length(levelB_kirc)
length(levelC_kirc)
setdiff(novel_coh_kirc$gene, c(levelA_kirc, levelB_kirc, levelC_kirc))
length(levelA_kirc_pers)
length(levelB_kirc_pers)
length(levelC_kirc_pers)
setdiff(novel_pers_kirc$gene, c(levelA_kirc_pers, levelB_kirc_pers, levelC_kirc_pers))
novel_coh_kirc <- novel_coh_kirc %>%
  mutate(out_degree = dg[gene]) %>%
  mutate(evidence_level = case_when(
    gene %in% levelA_kirc ~ "A",
    gene %in% levelB_kirc ~ "B",
    gene %in% levelC_kirc ~ "C",
    TRUE ~ "D"
  ))
write.csv(novel_coh_kirc, file.path(fig_dir, "coh", "KIRC_novel_coh_genes.csv"), row.names = FALSE) 
novel_coh_df <- rbind(novel_coh_df, data.frame(
  cancer = "KIRC",
  num_novel = length(novel_coh_kirc$gene),
  num_levelA = length(levelA_kirc),
  num_levelB = length(levelB_kirc), 
  num_levelC = length(levelC_kirc),
  num_levelD = length(setdiff(novel_coh_kirc$gene, c(levelA_kirc, levelB_kirc, levelC_kirc)))
))
novel_pers_df <- rbind(novel_pers_df, data.frame(
  cancer = "KIRC",
  num_novel = length(novel_pers_kirc$gene),
  num_levelA = length(levelA_kirc_pers),
  num_levelB = length(levelB_kirc_pers),
  num_levelC = length(levelC_kirc_pers),  
  num_levelD = length(setdiff(novel_pers_kirc$gene, c(levelA_kirc_pers, levelB_kirc_pers, levelC_kirc_pers)))
))

#KIRP
top_hnw_coh_kirp <- top_hnw_coh_df %>% filter(cancer == "KIRP") %>% arrange(desc(mf))
novel_coh_kirp <- top_hnw_coh_kirp %>% 
  filter(gene %in% novel_res_coh[["KIRP"]])
top_hnw_pers_kirp <- top_hnw_pers_df %>% filter(cancer == "KIRP") %>% arrange(desc(count))
novel_pers_kirp <- top_hnw_pers_kirp %>% 
  filter(gene %in% novel_res_pers[["KIRP"]])

cgc1_kirp <- CGC_df %>%
  filter(grepl("papillary renal", Tumour.Types.Somatic., ignore.case = TRUE)) %>%
  pull(Gene.Symbol) %>%
  unique()
levelA_kirp <- intersect(novel_coh_kirp$gene, cgc1_kirp)
levelB_kirp <- setdiff(intersect(novel_coh_kirp$gene, CGC),levelA_kirp)
levelA_kirp_pers <- intersect(novel_pers_kirp$gene, cgc1_kirp)
levelB_kirp_pers <- setdiff(intersect(novel_pers_kirp$gene, CGC),levelA_kirp_pers)

cgc2_kirp <- CGC2_df %>%
  filter(grepl("papillary renal", Tumour.Types.Somatic., ignore.case = TRUE)) %>%
  pull(Gene.Symbol) %>%
  unique()
levelB_kirp <- union(levelB_kirp, intersect(novel_coh_kirp$gene, cgc2_kirp))
levelC_kirp <- intersect(novel_coh_kirp$gene, CGC2_df$Gene.Symbol)
levelB_kirp_pers <- union(levelB_kirp_pers, intersect(novel_pers_kirp$gene, cgc2_kirp))
levelC_kirp_pers <- intersect(novel_pers_kirp$gene, CGC2_df$Gene.Symbol)

ncg_kirp <- NCG_df %>%
  filter(cancer_type %in% c("papillary_renal_cell_carcinoma", "renal_cancer_(all_histologies)")) %>%
  pull(symbol) %>%
  unique()
levelB_kirp <- union(levelB_kirp, intersect(novel_coh_kirp$gene, ncg_kirp))
levelC_kirp <- union(levelC_kirp, intersect(novel_coh_kirp$gene, NCG_df$symbol))
levelB_kirp_pers <- union(levelB_kirp_pers, intersect(novel_pers_kirp$gene, ncg_kirp))
levelC_kirp_pers <- union(levelC_kirp_pers, intersect(novel_pers_kirp$gene, NCG_df$symbol))

kirp_cm_terms <- c("papillary renal cell carcinoma", "familial renal papillary carcinoma")
cm_kirp_genes <- CM_df %>%
  filter(cancer_normalized %in% kirp_cm_terms) %>%
  # filter(citation_count >= 2) %>% 
  pull(gene_normalized) %>%
  unique()
levelB_kirp <- union(levelB_kirp, intersect(novel_coh_kirp$gene, cm_kirp_genes))
levelC_kirp <- union(levelC_kirp, intersect(novel_coh_kirp$gene, CM_df$gene_normalized))
levelB_kirp_pers <- union(levelB_kirp_pers, intersect(novel_pers_kirp$gene, cm_kirp_genes))
levelC_kirp_pers <- union(levelC_kirp_pers, intersect(novel_pers_kirp$gene, CM_df$gene_normalized))
levelB_kirp <- setdiff(levelB_kirp, levelA_kirp)
levelC_kirp <- setdiff(levelC_kirp, union(levelA_kirp, levelB_kirp))
levelB_kirp_pers <- setdiff(levelB_kirp_pers, levelA_kirp_pers)
levelC_kirp_pers <- setdiff(levelC_kirp_pers, union(levelA_kirp_pers, levelB_kirp_pers))
length(levelA_kirp)
length(levelB_kirp)
length(levelC_kirp)
setdiff(novel_coh_kirp$gene, c(levelA_kirp, levelB_kirp, levelC_kirp))
length(levelA_kirp_pers)
length(levelB_kirp_pers)
length(levelC_kirp_pers)
setdiff(novel_pers_kirp$gene, c(levelA_kirp_pers, levelB_kirp_pers, levelC_kirp_pers))
novel_coh_kirp <- novel_coh_kirp %>%
  mutate(out_degree = dg[gene]) %>%
  mutate(evidence_level = case_when(
    gene %in% levelA_kirp ~ "A",
    gene %in% levelB_kirp ~ "B",
    gene %in% levelC_kirp ~ "C",
    TRUE ~ "D"
  ))
write.csv(novel_coh_kirp, file.path(fig_dir, "coh", "KIRP_novel_coh_genes.csv"), row.names = FALSE) 
novel_coh_df <- rbind(novel_coh_df, data.frame(
  cancer = "KIRP",
  num_novel = length(novel_coh_kirp$gene),
  num_levelA = length(levelA_kirp),
  num_levelB = length(levelB_kirp),
  num_levelC = length(levelC_kirp),
  num_levelD = length(setdiff(novel_coh_kirp$gene, c(levelA_kirp, levelB_kirp, levelC_kirp)))
))
novel_pers_df <- rbind(novel_pers_df, data.frame(
  cancer = "KIRP",
  num_novel = length(novel_pers_kirp$gene),
  num_levelA = length(levelA_kirp_pers),
  num_levelB = length(levelB_kirp_pers),
  num_levelC = length(levelC_kirp_pers),
  num_levelD = length(setdiff(novel_pers_kirp$gene, c(levelA_kirp_pers, levelB_kirp_pers, levelC_kirp_pers)))
))

#LIHC
top_hnw_coh_lihc <- top_hnw_coh_df %>% filter(cancer == "LIHC") %>% arrange(desc(mf))
novel_coh_lihc <- top_hnw_coh_lihc %>% 
  filter(gene %in% novel_res_coh[["LIHC"]])
top_hnw_pers_lihc <- top_hnw_pers_df %>% filter(cancer == "LIHC") %>% arrange(desc(count))
novel_pers_lihc <- top_hnw_pers_lihc %>% 
  filter(gene %in% novel_res_pers[["LIHC"]])

cgc1_lihc <- CGC_df %>%
  filter(grepl("liver|hepatocellular|hepatoblastoma", Tumour.Types.Somatic., ignore.case = TRUE)) %>%
  pull(Gene.Symbol) %>%
  unique()
levelA_lihc <- intersect(novel_coh_lihc$gene, cgc1_lihc)
levelB_lihc <- setdiff(intersect(novel_coh_lihc$gene, CGC),levelA_lihc)
levelA_lihc_pers <- intersect(novel_pers_lihc$gene, cgc1_lihc)
levelB_lihc_pers <- setdiff(intersect(novel_pers_lihc$gene, CGC),levelA_lihc_pers)

cgc2_lihc <- CGC2_df %>%
  filter(grepl("liver|hepatocellular|hepatoblastoma", Tumour.Types.Somatic., ignore.case = TRUE)) %>%
  pull(Gene.Symbol) %>%
  unique()
levelB_lihc <- union(levelB_lihc, intersect(novel_coh_lihc$gene, cgc2_lihc))
levelC_lihc <- intersect(novel_coh_lihc$gene, CGC2_df$Gene.Symbol)
levelB_lihc_pers <- union(levelB_lihc_pers, intersect(novel_pers_lihc$gene, cgc2_lihc))
levelC_lihc_pers <- intersect(novel_pers_lihc$gene, CGC2_df$Gene.Symbol)

ncg_lihc <- NCG_df %>%
  filter(cancer_type %in% c("hepatocellular_carcinoma", "hepatocellular_carcinoma, cholangiocarcinoma")) %>%
  pull(symbol) %>%
  unique()
levelB_lihc <- union(levelB_lihc, intersect(novel_coh_lihc$gene, ncg_lihc))
levelC_lihc <- union(levelC_lihc, intersect(novel_coh_lihc$gene, NCG_df$symbol))
levelB_lihc_pers <- union(levelB_lihc_pers, intersect(novel_pers_lihc$gene, ncg_lihc))
levelC_lihc_pers <- union(levelC_lihc_pers, intersect(novel_pers_lihc$gene, NCG_df$symbol))

lihc_cm_terms <- c("liver cancer", "liver carcinoma", "hepatocellular carcinoma", "fibrolamellar carcinoma")
cm_lihc_genes <- CM_df %>%
  filter(cancer_normalized %in% lihc_cm_terms) %>%
  # filter(citation_count >= 2) %>% 
  pull(gene_normalized) %>%
  unique()
levelB_lihc <- union(levelB_lihc, intersect(novel_coh_lihc$gene, cm_lihc_genes))
levelC_lihc <- union(levelC_lihc, intersect(novel_coh_lihc$gene, CM_df$gene_normalized))
levelB_lihc_pers <- union(levelB_lihc_pers, intersect(novel_pers_lihc$gene, cm_lihc_genes))
levelC_lihc_pers <- union(levelC_lihc_pers, intersect(novel_pers_lihc$gene, CM_df$gene_normalized))
levelB_lihc <- setdiff(levelB_lihc, levelA_lihc)
levelC_lihc <- setdiff(levelC_lihc, union(levelA_lihc, levelB_lihc))
levelB_lihc_pers <- setdiff(levelB_lihc_pers, levelA_lihc_pers)
levelC_lihc_pers <- setdiff(levelC_lihc_pers, union(levelA_lihc_pers, levelB_lihc_pers))
length(levelA_lihc)
length(levelB_lihc)
length(levelC_lihc)
setdiff(novel_coh_lihc$gene, c(levelA_lihc, levelB_lihc, levelC_lihc))
length(levelA_lihc_pers)  
length(levelB_lihc_pers)
length(levelC_lihc_pers)
setdiff(novel_pers_lihc$gene, c(levelA_lihc_pers, levelB_lihc_pers, levelC_lihc_pers))
novel_coh_lihc <- novel_coh_lihc %>%
  mutate(out_degree = dg[gene]) %>%
  mutate(evidence_level = case_when(
    gene %in% levelA_lihc ~ "A",  
    gene %in% levelB_lihc ~ "B",
    gene %in% levelC_lihc ~ "C",
    TRUE ~ "D"  
  ))
write.csv(novel_coh_lihc, file.path(fig_dir, "coh", "LIHC_novel_coh_genes.csv"), row.names = FALSE) 
novel_coh_df <- rbind(novel_coh_df, data.frame(
  cancer = "LIHC",
  num_novel = length(novel_coh_lihc$gene),
  num_levelA = length(levelA_lihc),
  num_levelB = length(levelB_lihc),
  num_levelC = length(levelC_lihc),
  num_levelD = length(setdiff(novel_coh_lihc$gene, c(levelA_lihc, levelB_lihc, levelC_lihc)))
))
novel_pers_df <- rbind(novel_pers_df, data.frame(
  cancer = "LIHC",
  num_novel = length(novel_pers_lihc$gene), 
  num_levelA = length(levelA_lihc_pers),
  num_levelB = length(levelB_lihc_pers),
  num_levelC = length(levelC_lihc_pers),
  num_levelD = length(setdiff(novel_pers_lihc$gene, c(levelA_lihc_pers, levelB_lihc_pers, levelC_lihc_pers))) 
))

#LUAD
top_hnw_coh_luad <- top_hnw_coh_df %>% filter(cancer == "LUAD") %>% arrange(desc(mf))
novel_coh_luad <- top_hnw_coh_luad %>% 
  filter(gene %in% novel_res_coh[["LUAD"]])
top_hnw_pers_luad <- top_hnw_pers_df %>% filter(cancer == "LUAD") %>% arrange(desc(count))
novel_pers_luad <- top_hnw_pers_luad %>% 
  filter(gene %in% novel_res_pers[["LUAD"]])

cgc1_luad <- CGC_df %>%
  filter(!grepl("small cell lung|lung SCC", Tumour.Types.Somatic., ignore.case = TRUE)) %>%
  filter(grepl("lung adenocarcinoma|NSCLC|lung cancer|lung carcinoma|\\blung\\b", Tumour.Types.Somatic., ignore.case = TRUE)) %>%
  pull(Gene.Symbol) %>%
  unique()
levelA_luad <- intersect(novel_coh_luad$gene, cgc1_luad)
levelB_luad <- setdiff(intersect(novel_coh_luad$gene, CGC),levelA_luad)
levelA_luad_pers <- intersect(novel_pers_luad$gene, cgc1_luad)
levelB_luad_pers <- setdiff(intersect(novel_pers_luad$gene, CGC),levelA_luad_pers)

cgc2_luad <- CGC2_df %>%
  # 第一步：排除小细胞肺癌(SCLC)和肺鳞癌(Lung SCC)
  filter(!grepl("SCLC|small cell lung|lung SCC", Tumour.Types.Somatic., ignore.case = TRUE)) %>%
  # 第二步：匹配 lung adenocarcinoma, NSCLC, lung cancer, 以及单独的 lung
  filter(grepl("lung adenocarcinoma|NSCLC|lung cancer|lung carcinoma|\\blung\\b", Tumour.Types.Somatic., ignore.case = TRUE)) %>%
  pull(Gene.Symbol) %>%
  unique()
levelB_luad <- union(levelB_luad, intersect(novel_coh_luad$gene, cgc2_luad))
levelC_luad <- intersect(novel_coh_luad$gene, CGC2_df$Gene.Symbol)
levelB_luad_pers <- union(levelB_luad_pers, intersect(novel_pers_luad$gene, cgc2_luad))
levelC_luad_pers <- intersect(novel_pers_luad$gene, CGC2_df$Gene.Symbol)

luad_ncg_terms <- c(
  "lung_adenocarcinoma", 
  "non-small_cell_lung_cancer", 
  "lung_cancer_(all_histologies)"
)
ncg_luad <- NCG_df %>%
  filter(cancer_type %in% luad_ncg_terms) %>%
  pull(symbol) %>% 
  unique()
levelB_luad <- union(levelB_luad, intersect(novel_coh_luad$gene, ncg_luad))
levelC_luad <- union(levelC_luad, intersect(novel_coh_luad$gene, NCG_df$symbol))
levelB_luad_pers <- union(levelB_luad_pers, intersect(novel_pers_luad$gene, ncg_luad))
levelC_luad_pers <- union(levelC_luad_pers, intersect(novel_pers_luad$gene, NCG_df$symbol))

luad_cm_terms <- c(
  "lung adenocarcinoma", 
  "lung papillary adenocarcinoma", 
  "lung non-small cell carcinoma", 
  "lung cancer", 
  "lung carcinoma"
)
cm_luad_genes <- CM_df %>%
  filter(cancer_normalized %in% luad_cm_terms) %>%
  # filter(citation_count >= 2) %>% 
  pull(gene_normalized) %>%
  unique()

levelB_luad <- union(levelB_luad, intersect(novel_coh_luad$gene, cm_luad_genes))
levelC_luad <- union(levelC_luad, intersect(novel_coh_luad$gene, CM_df$gene_normalized))
levelB_luad_pers <- union(levelB_luad_pers, intersect(novel_pers_luad$gene, cm_luad_genes))
levelC_luad_pers <- union(levelC_luad_pers, intersect(novel_pers_luad$gene, CM_df$gene_normalized))
levelB_luad <- setdiff(levelB_luad, levelA_luad)
levelC_luad <- setdiff(levelC_luad, union(levelA_luad, levelB_luad))
length(levelA_luad)
length(levelB_luad)
length(levelC_luad)
setdiff(novel_coh_luad$gene, c(levelA_luad, levelB_luad, levelC_luad))
length(levelA_luad_pers)
length(levelB_luad_pers)
length(levelC_luad_pers)
setdiff(novel_pers_luad$gene, c(levelA_luad_pers, levelB_luad_pers, levelC_luad_pers))
novel_coh_luad <- novel_coh_luad %>%
  mutate(out_degree = dg[gene]) %>%
  mutate(evidence_level = case_when(
    gene %in% levelA_luad ~ "A",
    gene %in% levelB_luad ~ "B",
    gene %in% levelC_luad ~ "C",
    TRUE ~ "D"
  ))
write.csv(novel_coh_luad, file.path(fig_dir, "coh", "LUAD_novel_coh_genes.csv"), row.names = FALSE) 
novel_coh_df <- rbind(novel_coh_df, data.frame(
  cancer = "LUAD",
  num_novel = length(novel_coh_luad$gene),
  num_levelA = length(levelA_luad),
  num_levelB = length(levelB_luad),
  num_levelC = length(levelC_luad),
  num_levelD = length(setdiff(novel_coh_luad$gene, c(levelA_luad, levelB_luad, levelC_luad)))
))
novel_pers_df <- rbind(novel_pers_df, data.frame(
  cancer = "LUAD",
  num_novel = length(novel_pers_luad$gene),
  num_levelA = length(levelA_luad_pers),
  num_levelB = length(levelB_luad_pers),
  num_levelC = length(levelC_luad_pers),
  num_levelD = length(setdiff(novel_pers_luad$gene, c(levelA_luad_pers, levelB_luad_pers, levelC_luad_pers)))
))

#LUSC
top_hnw_coh_lusc <- top_hnw_coh_df %>% filter(cancer == "LUSC") %>% arrange(desc(mf))
novel_coh_lusc <- top_hnw_coh_lusc %>% 
  filter(gene %in% novel_res_coh[["LUSC"]])
top_hnw_pers_lusc <- top_hnw_pers_df %>% filter(cancer == "LUSC") %>% arrange(desc(count))
novel_pers_lusc <- top_hnw_pers_lusc %>% 
  filter(gene %in% novel_res_pers[["LUSC"]])

cgc1_lusc <- CGC_df %>%
  filter(grepl("lung SCC", Tumour.Types.Somatic., ignore.case = TRUE)) %>%
  pull(Gene.Symbol) %>%
  unique()
levelA_lusc <- intersect(novel_coh_lusc$gene, cgc1_lusc)
levelB_lusc <- setdiff(intersect(novel_coh_lusc$gene, CGC),levelA_lusc)
levelA_lusc_pers <- intersect(novel_pers_lusc$gene, cgc1_lusc)
levelB_lusc_pers <- setdiff(intersect(novel_pers_lusc$gene, CGC),levelA_lusc_pers)

cgc2_lusc <- CGC2_df %>%
  filter(grepl("lung SCC", Tumour.Types.Somatic., ignore.case = TRUE)) %>%
  pull(Gene.Symbol) %>%
  unique()
levelB_lusc <- union(levelB_lusc, intersect(novel_coh_lusc$gene, cgc2_lusc))
levelC_lusc <- intersect(novel_coh_lusc$gene, CGC2_df$Gene.Symbol)
levelB_lusc_pers <- union(levelB_lusc_pers, intersect(novel_pers_lusc$gene, cgc2_lusc))
levelC_lusc_pers <- intersect(novel_pers_lusc$gene, CGC2_df$Gene.Symbol)

ncg_lusc <- NCG_df %>%
  filter(cancer_type %in% c("lung_squamous_cell_carcinoma", "non-small_cell_lung_cancer", "lung_cancer_(all_histologies)")) %>%
  pull(symbol) %>%
  unique()
levelB_lusc <- union(levelB_lusc, intersect(novel_coh_lusc$gene, ncg_lusc))
levelC_lusc <- union(levelC_lusc, intersect(novel_coh_lusc$gene, NCG_df$symbol))
levelB_lusc_pers <- union(levelB_lusc_pers, intersect(novel_pers_lusc$gene, ncg_lusc))
levelC_lusc_pers <- union(levelC_lusc_pers, intersect(novel_pers_lusc$gene, NCG_df$symbol))

lusc_cm_terms <- c("lung squamous cell carcinoma")
cm_lusc_genes <- CM_df %>%
  filter(cancer_normalized %in% lusc_cm_terms) %>%
  # filter(citation_count >= 2) %>% 
  pull(gene_normalized) %>%
  unique()
levelB_lusc <- union(levelB_lusc, intersect(novel_coh_lusc$gene, cm_lusc_genes))
levelC_lusc <- union(levelC_lusc, intersect(novel_coh_lusc$gene, CM_df$gene_normalized))
levelB_lusc_pers <- union(levelB_lusc_pers, intersect(novel_pers_lusc$gene, cm_lusc_genes))
levelC_lusc_pers <- union(levelC_lusc_pers, intersect(novel_pers_lusc$gene, CM_df$gene_normalized))
levelB_lusc <- setdiff(levelB_lusc, levelA_lusc)
levelC_lusc <- setdiff(levelC_lusc, union(levelA_lusc, levelB_lusc))
levelB_lusc_pers <- setdiff(levelB_lusc_pers, levelA_lusc_pers)
levelC_lusc_pers <- setdiff(levelC_lusc_pers, union(levelA_lusc_pers, levelB_lusc_pers))
length(levelA_lusc)
length(levelB_lusc)
length(levelC_lusc)
setdiff(novel_coh_lusc$gene, c(levelA_lusc, levelB_lusc, levelC_lusc))
length(levelA_lusc_pers)
length(levelB_lusc_pers)
length(levelC_lusc_pers)
setdiff(novel_pers_lusc$gene, c(levelA_lusc_pers, levelB_lusc_pers, levelC_lusc_pers))
novel_coh_lusc <- novel_coh_lusc %>%
  mutate(out_degree = dg[gene]) %>%
  mutate(evidence_level = case_when(
    gene %in% levelA_lusc ~ "A",
    gene %in% levelB_lusc ~ "B",
    gene %in% levelC_lusc ~ "C",
    TRUE ~ "D"
  ))
write.csv(novel_coh_lusc, file.path(fig_dir, "coh", "LUSC_novel_coh_genes.csv"), row.names = FALSE) 
novel_coh_df <- rbind(novel_coh_df, data.frame(
  cancer = "LUSC",
  num_novel = length(novel_coh_lusc$gene),
  num_levelA = length(levelA_lusc),
  num_levelB = length(levelB_lusc),
  num_levelC = length(levelC_lusc),
  num_levelD = length(setdiff(novel_coh_lusc$gene, c(levelA_lusc, levelB_lusc, levelC_lusc)))
))
novel_pers_df <- rbind(novel_pers_df, data.frame(
  cancer = "LUSC",
  num_novel = length(novel_pers_lusc$gene),
  num_levelA = length(levelA_lusc_pers),
  num_levelB = length(levelB_lusc_pers),
  num_levelC = length(levelC_lusc_pers),
  num_levelD = length(setdiff(novel_pers_lusc$gene, c(levelA_lusc_pers, levelB_lusc_pers, levelC_lusc_pers)))
))

#PRAD
top_hnw_coh_prad <- top_hnw_coh_df %>% filter(cancer == "PRAD") %>% arrange(desc(mf))
novel_coh_prad <- top_hnw_coh_prad %>% 
  filter(gene %in% novel_res_coh[["PRAD"]])
top_hnw_pers_prad <- top_hnw_pers_df %>% filter(cancer == "PRAD") %>% arrange(desc(count))
novel_pers_prad <- top_hnw_pers_prad %>% 
  filter(gene %in% novel_res_pers[["PRAD"]])

cgc1_prad <- CGC_df %>%
  filter(grepl("prostate", Tumour.Types.Somatic., ignore.case = TRUE)) %>%
  pull(Gene.Symbol) %>%
  unique()
levelA_prad <- intersect(novel_coh_prad$gene, cgc1_prad)
levelB_prad <- setdiff(intersect(novel_coh_prad$gene, CGC),levelA_prad)
levelA_prad_pers <- intersect(novel_pers_prad$gene, cgc1_prad)
levelB_prad_pers <- setdiff(intersect(novel_pers_prad$gene, CGC),levelA_prad_pers)

cgc2_prad <- CGC2_df %>%
  filter(grepl("prostate", Tumour.Types.Somatic., ignore.case = TRUE)) %>%
  pull(Gene.Symbol) %>%
  unique()
levelB_prad <- union(levelB_prad, intersect(novel_coh_prad$gene, cgc2_prad))
levelC_prad <- intersect(novel_coh_prad$gene, CGC2_df$Gene.Symbol)
levelB_prad_pers <- union(levelB_prad_pers, intersect(novel_pers_prad$gene, cgc2_prad))
levelC_prad_pers <- intersect(novel_pers_prad$gene, CGC2_df$Gene.Symbol)

ncg_prad <- NCG_df %>%
  filter(cancer_type == "prostate_cancer") %>%
  pull(symbol) %>%
  unique()
levelB_prad <- union(levelB_prad, intersect(novel_coh_prad$gene, ncg_prad))
levelC_prad <- union(levelC_prad, intersect(novel_coh_prad$gene, NCG_df$symbol))
levelB_prad_pers <- union(levelB_prad_pers, intersect(novel_pers_prad$gene, ncg_prad))
levelC_prad_pers <- union(levelC_prad_pers, intersect(novel_pers_prad$gene, NCG_df$symbol))

prad_cm_terms <- c("prostate cancer", "prostate carcinoma", "prostate adenocarcinoma")
cm_prad_genes <- CM_df %>%
  filter(cancer_normalized %in% prad_cm_terms) %>%
  # filter(citation_count >= 2) %>% 
  pull(gene_normalized) %>%
  unique()
levelB_prad <- union(levelB_prad, intersect(novel_coh_prad$gene, cm_prad_genes))
levelC_prad <- union(levelC_prad, intersect(novel_coh_prad$gene, CM_df$gene_normalized))
levelB_prad_pers <- union(levelB_prad_pers, intersect(novel_pers_prad$gene, cm_prad_genes))
levelC_prad_pers <- union(levelC_prad_pers, intersect(novel_pers_prad$gene, CM_df$gene_normalized))
levelB_prad <- setdiff(levelB_prad, levelA_prad)
levelC_prad <- setdiff(levelC_prad, union(levelA_prad, levelB_prad))
levelB_prad_pers <- setdiff(levelB_prad_pers, levelA_prad_pers)
levelC_prad_pers <- setdiff(levelC_prad_pers, union(levelA_prad_pers, levelB_prad_pers))
length(levelA_prad)
length(levelB_prad)
length(levelC_prad)
setdiff(novel_coh_prad$gene, c(levelA_prad, levelB_prad, levelC_prad))
length(levelA_prad_pers)
length(levelB_prad_pers)
length(levelC_prad_pers)
setdiff(novel_pers_prad$gene, c(levelA_prad_pers, levelB_prad_pers, levelC_prad_pers))
novel_coh_prad <- novel_coh_prad %>%
  mutate(out_degree = dg[gene]) %>%
  mutate(evidence_level = case_when(
    gene %in% levelA_prad ~ "A",
    gene %in% levelB_prad ~ "B",
    gene %in% levelC_prad ~ "C",
    TRUE ~ "D"
  ))
write.csv(novel_coh_prad, file.path(fig_dir, "coh", "PRAD_novel_coh_genes.csv"), row.names = FALSE) 
novel_coh_df <- rbind(novel_coh_df, data.frame(
  cancer = "PRAD",
  num_novel = length(novel_coh_prad$gene),
  num_levelA = length(levelA_prad),
  num_levelB = length(levelB_prad),
  num_levelC = length(levelC_prad),
  num_levelD = length(setdiff(novel_coh_prad$gene, c(levelA_prad, levelB_prad, levelC_prad)))
))
novel_pers_df <- rbind(novel_pers_df, data.frame(
  cancer = "PRAD",
  num_novel = length(novel_pers_prad$gene),
  num_levelA = length(levelA_prad_pers),
  num_levelB = length(levelB_prad_pers),
  num_levelC = length(levelC_prad_pers),
  num_levelD = length(setdiff(novel_pers_prad$gene, c(levelA_prad_pers, levelB_prad_pers, levelC_prad_pers))) 
))

#STAD
top_hnw_coh_stad <- top_hnw_coh_df %>% filter(cancer == "STAD") %>% arrange(desc(mf))
novel_coh_stad <- top_hnw_coh_stad %>% 
  filter(gene %in% novel_res_coh[["STAD"]])
top_hnw_pers_stad <- top_hnw_pers_df %>% filter(cancer == "STAD") %>% arrange(desc(count))
novel_pers_stad <- top_hnw_pers_stad %>% 
  filter(gene %in% novel_res_pers[["STAD"]])

cgc1_stad <- CGC_df %>%
  filter(grepl("stomach|gastric", Tumour.Types.Somatic., ignore.case = TRUE)) %>%
  pull(Gene.Symbol) %>%
  unique()
levelA_stad <- intersect(novel_coh_stad$gene, cgc1_stad)
levelB_stad <- setdiff(intersect(novel_coh_stad$gene, CGC),levelA_stad)
levelA_stad_pers <- intersect(novel_pers_stad$gene, cgc1_stad)
levelB_stad_pers <- setdiff(intersect(novel_pers_stad$gene, CGC),levelA_stad_pers)

cgc2_stad <- CGC2_df %>%
  filter(grepl("stomach|gastric", Tumour.Types.Somatic., ignore.case = TRUE)) %>%
  pull(Gene.Symbol) %>%
  unique()
levelB_stad <- union(levelB_stad, intersect(novel_coh_stad$gene, cgc2_stad))
levelC_stad <- intersect(novel_coh_stad$gene, CGC2_df$Gene.Symbol)
levelB_stad_pers <- union(levelB_stad_pers, intersect(novel_pers_stad$gene, cgc2_stad))
levelC_stad_pers <- intersect(novel_pers_stad$gene, CGC2_df$Gene.Symbol)

ncg_stad <- NCG_df %>%
  filter(cancer_type %in% c("gastric_cancer", "gastric_adenocarcinoma", "diffuse_gastric_adenocarcinoma", "mucinous_gastric_cancer", "pan-gastric")) %>%
  pull(symbol) %>%
  unique()
levelB_stad <- union(levelB_stad, intersect(novel_coh_stad$gene, ncg_stad))
levelC_stad <- union(levelC_stad, intersect(novel_coh_stad$gene, NCG_df$symbol))
levelB_stad_pers <- union(levelB_stad_pers, intersect(novel_pers_stad$gene, ncg_stad))
levelC_stad_pers <- union(levelC_stad_pers, intersect(novel_pers_stad$gene, NCG_df$symbol))

stad_cm_terms <- c("stomach cancer", "stomach carcinoma", "gastric cancer", "gastric adenocarcinoma", "diffuse gastric cancer", "hereditary diffuse gastric cancer", "gastric cardia adenocarcinoma", "microinvasive gastric cancer")
cm_stad_genes <- CM_df %>%
  filter(cancer_normalized %in% stad_cm_terms) %>%
  # filter(citation_count >= 2) %>% 
  pull(gene_normalized) %>%
  unique()
levelB_stad <- union(levelB_stad, intersect(novel_coh_stad$gene, cm_stad_genes))
levelC_stad <- union(levelC_stad, intersect(novel_coh_stad$gene, CM_df$gene_normalized))
levelB_stad_pers <- union(levelB_stad_pers, intersect(novel_pers_stad$gene, cm_stad_genes))
levelC_stad_pers <- union(levelC_stad_pers, intersect(novel_pers_stad$gene, CM_df$gene_normalized))
levelB_stad <- setdiff(levelB_stad, levelA_stad)
levelC_stad <- setdiff(levelC_stad, union(levelA_stad, levelB_stad))
levelB_stad_pers <- setdiff(levelB_stad_pers, levelA_stad_pers)
levelC_stad_pers <- setdiff(levelC_stad_pers, union(levelA_stad_pers, levelB_stad_pers))
length(levelA_stad)
length(levelB_stad)
length(levelC_stad)
setdiff(novel_coh_stad$gene, c(levelA_stad, levelB_stad, levelC_stad))
length(levelA_stad_pers)
length(levelB_stad_pers)
length(levelC_stad_pers)
setdiff(novel_pers_stad$gene, c(levelA_stad_pers, levelB_stad_pers, levelC_stad_pers))
novel_coh_stad <- novel_coh_stad %>%
  mutate(out_degree = dg[gene]) %>%
  mutate(evidence_level = case_when(
    gene %in% levelA_stad ~ "A",
    gene %in% levelB_stad ~ "B",
    gene %in% levelC_stad ~ "C",
    TRUE ~ "D"
  ))
write.csv(novel_coh_stad, file.path(fig_dir, "coh", "STAD_novel_coh_genes.csv"), row.names = FALSE) 
novel_coh_df <- rbind(novel_coh_df, data.frame(
  cancer = "STAD",
  num_novel = length(novel_coh_stad$gene),
  num_levelA = length(levelA_stad),
  num_levelB = length(levelB_stad),
  num_levelC = length(levelC_stad),
  num_levelD = length(setdiff(novel_coh_stad$gene, c(levelA_stad, levelB_stad, levelC_stad)))
))
novel_pers_df <- rbind(novel_pers_df, data.frame(
  cancer = "STAD",
  num_novel = length(novel_pers_stad$gene), 
  num_levelA = length(levelA_stad_pers),
  num_levelB = length(levelB_stad_pers),  
  num_levelC = length(levelC_stad_pers),
  num_levelD = length(setdiff(novel_pers_stad$gene, c(levelA_stad_pers, levelB_stad_pers, levelC_stad_pers)))
))

#THCA
top_hnw_coh_thca <- top_hnw_coh_df %>% filter(cancer == "THCA") %>% arrange(desc(mf))
novel_coh_thca <- top_hnw_coh_thca %>% 
  filter(gene %in% novel_res_coh[["THCA"]])
top_hnw_pers_thca <- top_hnw_pers_df %>% filter(cancer == "THCA") %>% arrange(desc(count))
novel_pers_thca <- top_hnw_pers_thca %>% 
  filter(gene %in% novel_res_pers[["THCA"]])

cgc1_thca <- CGC_df %>%
  filter(grepl("thyroid", Tumour.Types.Somatic., ignore.case = TRUE)) %>%
  pull(Gene.Symbol) %>%
  unique()
levelA_thca <- intersect(novel_coh_thca$gene, cgc1_thca)
levelB_thca <- setdiff(intersect(novel_coh_thca$gene, CGC),levelA_thca)
levelA_thca_pers <- intersect(novel_pers_thca$gene, cgc1_thca)
levelB_thca_pers <- setdiff(intersect(novel_pers_thca$gene, CGC),levelA_thca_pers)

cgc2_thca <- CGC2_df %>%
  filter(grepl("thyroid", Tumour.Types.Somatic., ignore.case = TRUE)) %>%
  pull(Gene.Symbol) %>%
  unique()
levelB_thca <- union(levelB_thca, intersect(novel_coh_thca$gene, cgc2_thca))
levelC_thca <- intersect(novel_coh_thca$gene, CGC2_df$Gene.Symbol)
levelB_thca_pers <- union(levelB_thca_pers, intersect(novel_pers_thca$gene, cgc2_thca))
levelC_thca_pers <- intersect(novel_pers_thca$gene, CGC2_df$Gene.Symbol)

ncg_thca <- NCG_df %>%
  filter(cancer_type %in% c("papillary_thyroid_cancer", "anaplastic_thyroid_carcinoma")) %>%
  pull(symbol) %>%
  unique()
levelB_thca <- union(levelB_thca, intersect(novel_coh_thca$gene, ncg_thca))
levelC_thca <- union(levelC_thca, intersect(novel_coh_thca$gene, NCG_df$symbol))
levelB_thca_pers <- union(levelB_thca_pers, intersect(novel_pers_thca$gene, ncg_thca))
levelC_thca_pers <- union(levelC_thca_pers, intersect(novel_pers_thca$gene, NCG_df$symbol))

thca_cm_terms <- c("thyroid gland cancer", "thyroid gland carcinoma", "thyroid gland papillary carcinoma", "thyroid gland follicular carcinoma", "thyroid gland medullary carcinoma", "thyroid gland anaplastic carcinoma", "differentiated thyroid gland carcinoma")
cm_thca_genes <- CM_df %>%
  filter(cancer_normalized %in% thca_cm_terms) %>%
  # filter(citation_count >= 2) %>% 
  pull(gene_normalized) %>%
  unique()
levelB_thca <- union(levelB_thca, intersect(novel_coh_thca$gene, cm_thca_genes))
levelC_thca <- union(levelC_thca, intersect(novel_coh_thca$gene, CM_df$gene_normalized))
levelB_thca_pers <- union(levelB_thca_pers, intersect(novel_pers_thca$gene, cm_thca_genes))
levelC_thca_pers <- union(levelC_thca_pers, intersect(novel_pers_thca$gene, CM_df$gene_normalized))
levelB_thca <- setdiff(levelB_thca, levelA_thca)
levelC_thca <- setdiff(levelC_thca, union(levelA_thca, levelB_thca))
levelB_thca_pers <- setdiff(levelB_thca_pers, levelA_thca_pers)
levelC_thca_pers <- setdiff(levelC_thca_pers, union(levelA_thca_pers, levelB_thca_pers))
length(levelA_thca)
length(levelB_thca)
length(levelC_thca)
setdiff(novel_coh_thca$gene, c(levelA_thca, levelB_thca, levelC_thca))
length(levelA_thca_pers)
length(levelB_thca_pers)
length(levelC_thca_pers)
setdiff(novel_pers_thca$gene, c(levelA_thca_pers, levelB_thca_pers, levelC_thca_pers))
novel_coh_thca <- novel_coh_thca %>%
  mutate(out_degree = dg[gene]) %>%
  mutate(evidence_level = case_when(
    gene %in% levelA_thca ~ "A",
    gene %in% levelB_thca ~ "B",
    gene %in% levelC_thca ~ "C",
    TRUE ~ "D"
  ))
write.csv(novel_coh_thca, file.path(fig_dir, "coh", "THCA_novel_coh_genes.csv"), row.names = FALSE) 
novel_coh_df <- rbind(novel_coh_df, data.frame(
  cancer = "THCA",  
  num_novel = length(novel_coh_thca$gene),
  num_levelA = length(levelA_thca),
  num_levelB = length(levelB_thca),
  num_levelC = length(levelC_thca),
  num_levelD = length(setdiff(novel_coh_thca$gene, c(levelA_thca, levelB_thca, levelC_thca)))
))
novel_pers_df <- rbind(novel_pers_df, data.frame(
  cancer = "THCA",
  num_novel = length(novel_pers_thca$gene),
  num_levelA = length(levelA_thca_pers),  
  num_levelB = length(levelB_thca_pers),
  num_levelC = length(levelC_thca_pers),
  num_levelD = length(setdiff(novel_pers_thca$gene, c(levelA_thca_pers, levelB_thca_pers, levelC_thca_pers)))
))

#UCEC
top_hnw_coh_ucec <- top_hnw_coh_df %>% filter(cancer == "UCEC") %>% arrange(desc(mf))
novel_coh_ucec <- top_hnw_coh_ucec %>% 
  filter(gene %in% novel_res_coh[["UCEC"]])
top_hnw_pers_ucec <- top_hnw_pers_df %>% filter(cancer == "UCEC") %>% arrange(desc(count))
novel_pers_ucec <- top_hnw_pers_ucec %>% 
  filter(gene %in% novel_res_pers[["UCEC"]])

cgc1_ucec <- CGC_df %>%
  filter(grepl("endometrial|endometrium|endometrioid", Tumour.Types.Somatic., ignore.case = TRUE)) %>%
  pull(Gene.Symbol) %>%
  unique()
levelA_ucec <- intersect(novel_coh_ucec$gene, cgc1_ucec)
levelB_ucec <- setdiff(intersect(novel_coh_ucec$gene, CGC),levelA_ucec)
levelA_ucec_pers <- intersect(novel_pers_ucec$gene, cgc1_ucec)
levelB_ucec_pers <- setdiff(intersect(novel_pers_ucec$gene, CGC),levelA_ucec_pers)

cgc2_ucec <- CGC2_df %>%
  filter(grepl("endometrial|endometrium|endometrioid", Tumour.Types.Somatic., ignore.case = TRUE)) %>%
  pull(Gene.Symbol) %>%
  unique()
levelB_ucec <- union(levelB_ucec, intersect(novel_coh_ucec$gene, cgc2_ucec))
levelC_ucec <- intersect(novel_coh_ucec$gene, CGC2_df$Gene.Symbol)
levelB_ucec_pers <- union(levelB_ucec_pers, intersect(novel_pers_ucec$gene, cgc2_ucec))
levelC_ucec_pers <- intersect(novel_pers_ucec$gene, CGC2_df$Gene.Symbol)
ncg_ucec <- NCG_df %>%
  filter(cancer_type %in% c("endometrial_cancer", "serous_endometrial_cancer", "clear_cell_endometrial_cancer")) %>%
  pull(symbol) %>%
  unique()
levelB_ucec <- union(levelB_ucec, intersect(novel_coh_ucec$gene, ncg_ucec))
levelC_ucec <- union(levelC_ucec, intersect(novel_coh_ucec$gene, NCG_df$symbol))
levelB_ucec_pers <- union(levelB_ucec_pers, intersect(novel_pers_ucec$gene, ncg_ucec))
levelC_ucec_pers <- union(levelC_ucec_pers, intersect(novel_pers_ucec$gene, NCG_df$symbol))
ucec_cm_terms <- c("endometrial cancer", "endometrial carcinoma", "endometrial adenocarcinoma", "endometrial serous adenocarcinoma", "endometrial clear cell adenocarcinoma", "endometrial squamous cell carcinoma", "uterine cancer", "uterine corpus cancer", "uterine corpus endometrial carcinoma")
cm_ucec_genes <- CM_df %>%
  filter(cancer_normalized %in% ucec_cm_terms) %>%
  # filter(citation_count >= 2) %>% 
  pull(gene_normalized) %>%
  unique()
levelB_ucec <- union(levelB_ucec, intersect(novel_coh_ucec$gene, cm_ucec_genes))
levelC_ucec <- union(levelC_ucec, intersect(novel_coh_ucec$gene, CM_df$gene_normalized))
levelB_ucec_pers <- union(levelB_ucec_pers, intersect(novel_pers_ucec$gene, cm_ucec_genes))
levelC_ucec_pers <- union(levelC_ucec_pers, intersect(novel_pers_ucec$gene, CM_df$gene_normalized))
levelB_ucec <- setdiff(levelB_ucec, levelA_ucec)
levelC_ucec <- setdiff(levelC_ucec, union(levelA_ucec, levelB_ucec))
levelB_ucec_pers <- setdiff(levelB_ucec_pers, levelA_ucec_pers)
levelC_ucec_pers <- setdiff(levelC_ucec_pers, union(levelA_ucec_pers, levelB_ucec_pers))
length(levelA_ucec)
length(levelB_ucec)
length(levelC_ucec)
setdiff(novel_coh_ucec$gene, c(levelA_ucec, levelB_ucec, levelC_ucec))
length(levelA_ucec_pers)
length(levelB_ucec_pers)
length(levelC_ucec_pers)
setdiff(novel_pers_ucec$gene, c(levelA_ucec_pers, levelB_ucec_pers, levelC_ucec_pers))
novel_coh_ucec <- novel_coh_ucec %>%
  mutate(out_degree = dg[gene]) %>%
  mutate(evidence_level = case_when(
    gene %in% levelA_ucec ~ "A",
    gene %in% levelB_ucec ~ "B",
    gene %in% levelC_ucec ~ "C",
    TRUE ~ "D"
  ))
write.csv(novel_coh_ucec, file.path(fig_dir, "coh", "UCEC_novel_coh_genes.csv"), row.names = FALSE) 
novel_coh_df <- rbind(novel_coh_df, data.frame(
  cancer = "UCEC",
  num_novel = length(novel_coh_ucec$gene),
  num_levelA = length(levelA_ucec),
  num_levelB = length(levelB_ucec),
  num_levelC = length(levelC_ucec),
  num_levelD = length(setdiff(novel_coh_ucec$gene, c(levelA_ucec, levelB_ucec, levelC_ucec)))
))
novel_pers_df <- rbind(novel_pers_df, data.frame(
  cancer = "UCEC",
  num_novel = length(novel_pers_ucec$gene),
  num_levelA = length(levelA_ucec_pers),
  num_levelB = length(levelB_ucec_pers),
  num_levelC = length(levelC_ucec_pers),
  num_levelD = length(setdiff(novel_pers_ucec$gene, c(levelA_ucec_pers, levelB_ucec_pers, levelC_ucec_pers)))
))


# subnet
comb_adj_mat <- combine_STRING_and_omnipath()
comb_ppi <- graph_from_adjacency_matrix(comb_adj_mat, weight = TRUE)
nodes_ppi <- rownames(comb_adj_mat)

comb_adj_mat0 <- comb_adj_mat
comb_adj_mat0@x[comb_adj_mat0@x <= 1] <- 0
comb_adj_mat0 <- drop0(comb_adj_mat0)
comb_adj_lists <- build_adj_lists(comb_adj_mat0)

TF <- V(GRN_net)$name[which(igraph::degree(GRN_net,mode="out") > 0)]

comb_TF <- TF[grepl("::",TF,fixed=TRUE)]
total_TF <- setdiff(TF,comb_TF)
comb_TF_list <- strsplit(comb_TF, "::", fixed = TRUE)
comb_TF_split <- unique(unlist(comb_TF_list))
TF2comb <- lapply(comb_TF_split, function(tf) {
  comb_TF[sapply(comb_TF_list, function(x) tf %in% x)]
})
names(TF2comb) <- comb_TF_split
total_TF <- union(total_TF,comb_TF_split)
target <- setdiff(V(GRN_net)$name,TF)
target <- target[which(igraph::degree(GRN_net,mode="in")[target]>0)]
GRN_adj_mat <- as_adjacency_matrix(GRN_net, sparse = TRUE)
GRN_adj_mat <- GRN_adj_mat[c(TF,target),c(TF,target)]
GRN <- graph_from_adjacency_matrix(GRN_adj_mat)
GRN_adj_lists <- build_adj_lists(GRN_adj_mat)
GRN_list <- list(adj_mat = GRN_adj_mat, adj_lists = GRN_adj_lists, TF = TF, target = target, total_TF = total_TF, comb_TF = comb_TF, TF2comb = TF2comb, comb_TF_split = comb_TF_split)

# BRCA
cancer <- "BRCA"
res_file <- file.path("./results/HyperNetWalk", cancer, "results.Rdata")
load(res_file)
deg_mat <- objs$deg_mat
cor_mat <- WGCNA::cor(t(deg_mat))
cor_mat <- abs(cor_mat)
cor_mat[cor_mat < 0.1] <- 0.1
idx <- summary(comb_adj_mat)
g1 <- nodes_ppi[idx$i]
g2 <- nodes_ppi[idx$j]
w <- rep(0.1, length(g1))
hit <- g1 %in% rownames(cor_mat) & g2 %in% colnames(cor_mat)
w[hit] <- cor_mat[cbind(g1[hit], g2[hit])]
comb_adj_mat_brca <- comb_adj_mat
comb_adj_mat_brca@x <- w

gene1 <- "TBX3"
downstream1 <- get_downstream_TFs(gene1,comb_adj_mat0,total_TF)
TF1 <- intersect(downstream1, TF)
dg_TF1 <- degree(GRN_net, mode="out")[TF1]
sub_comb_adj_mat1 <- comb_adj_mat_brca[union(gene1, downstream1), union(gene1, downstream1)]
sub_net1 <- get_subgraph(sub_comb_adj_mat1,comb_adj_mat0,dg_TF1,gene1) 

gene2 <- "ZNF263"
downstream2 <- get_downstream_TFs(gene2,comb_adj_mat0,total_TF)
TF2 <- intersect(downstream2, TF)
# dg_TF2 <- degree(GRN_net, mode="out")[TF2]

d_out2 <- multi_source_bfs_indices(gene2,GRN_adj_lists$out,GRN_adj_lists$idx_map,3)
downstream_TF2 <- intersect(names(d_out2)[which(d_out2 <= 3)], TF)
TF_brca <- nodes$P2
TF_brca <- apply(TF_brca,2,function(x) x / sum(x))
TF_mean_brca <- rowMeans(TF_brca)
top_downstream_TF2 <- names(sort(TF_mean_brca[downstream_TF2], decreasing = TRUE))[1:10]
dg_TF2 <- degree(GRN_net, mode="out")[top_downstream_TF2]
sub_adj_mat2 <- GRN_adj_mat[union(gene2, top_downstream_TF2), union(gene2, top_downstream_TF2)]
sub_net2 <- get_subgraph(sub_adj_mat2,GRN_adj_mat,TF_mean_brca[top_downstream_TF2],gene2)

neighbor_2 <- colnames(comb_adj_mat_brca)[which(comb_adj_mat_brca[gene2,] > 0)]

gene3 <- "ERG"
downstream3 <- get_downstream_TFs(gene3,comb_adj_mat0,total_TF)
TF3 <- intersect(downstream3, TF)
dg_TF3 <- degree(GRN_net, mode="out")[TF3]
sub_comb_adj_mat3 <- comb_adj_mat_brca[union(gene3, downstream3), union(gene3, downstream3)]
sub_net3 <- get_subgraph(sub_comb_adj_mat3,comb_adj_mat0,dg_TF3,gene3)

gene4 <- "MAZ"
downstream4 <- get_downstream_TFs(gene4,comb_adj_mat0,total_TF)
TF4 <- intersect(downstream4, TF)
dg_TF4 <- degree(GRN_net, mode="out")[TF4]
sub_comb_adj_mat4 <- comb_adj_mat_brca[union(gene4, downstream4), union(gene4, downstream4)]
sub_net4 <- get_subgraph(sub_comb_adj_mat4,comb_adj_mat0,dg_TF4,gene4)

gene5 <- "KDM4C"
downstream5 <- get_downstream_TFs(gene5,comb_adj_mat0,total_TF)
TF5 <- intersect(downstream5, TF)
dg_TF5 <- degree(GRN_net, mode="out")[TF5]
sub_comb_adj_mat5 <- comb_adj_mat_brca[union(gene5, downstream5), union(gene5, downstream5)]
sub_net5 <- get_subgraph(sub_comb_adj_mat5,comb_adj_mat0,dg_TF5,gene5)

gene6 <- "EBF1"
downstream6 <- get_downstream_TFs(gene6,comb_adj_mat0,total_TF)
TF6 <- intersect(downstream6, TF)
dg_TF6 <- degree(GRN_net, mode="out")[TF6]

gene9 <- "SPI1"
downstream9 <- get_downstream_TFs(gene9,comb_adj_mat0,total_TF)
TF9 <- intersect(downstream9, TF)
dg_TF9 <- degree(GRN_net, mode="out")[TF9]
sub_comb_adj_mat9 <- comb_adj_mat_brca[union(gene9, downstream9), union(gene9, downstream9)]
sub_net9 <- get_subgraph(sub_comb_adj_mat9,comb_adj_mat0,dg_TF9,gene9)

gene10 <- "TCF12"
downstream10 <- get_downstream_TFs(gene10,comb_adj_mat0,total_TF)
TF10 <- intersect(downstream10, TF)
dg_TF10 <- degree(GRN_net, mode="out")[TF10]

gene11 <- "SP1"
downstream11 <- get_downstream_TFs(gene11,comb_adj_mat0,total_TF)
TF11 <- intersect(downstream11, TF)
dg_TF11 <- degree(GRN_net, mode="out")[TF11]
sub_comb_adj_mat11 <- comb_adj_mat_brca[union(gene11, downstream11), union(gene11, downstream11)]
sub_net11 <- get_subgraph(sub_comb_adj_mat11,comb_adj_mat0,dg_TF11,gene11)

gene12 <- "PHF8"
downstream12 <- get_downstream_TFs(gene12,comb_adj_mat0,total_TF)
TF12 <- intersect(downstream12, TF)
dg_TF12 <- degree(GRN_net, mode="out")[TF12]
sub_comb_adj_mat12 <- comb_adj_mat_brca[union(gene12, downstream12), union(gene12, downstream12)]
sub_net12 <- get_subgraph(sub_comb_adj_mat12,comb_adj_mat0,dg_TF12,gene12)

gene13 <- "STAG1"
downstream13 <- get_downstream_TFs(gene13,comb_adj_mat0,total_TF)
TF13 <- intersect(downstream13, TF)
dg_TF13 <- degree(GRN_net, mode="out")[TF13]

gene14 <- "BCOR"
downstream14 <- get_downstream_TFs(gene14,comb_adj_mat0,total_TF)
TF14 <- intersect(downstream14, TF)
dg_TF14 <- degree(GRN_net, mode="out")[TF14]
sub_comb_adj_mat14 <- comb_adj_mat_brca[union(gene14, downstream14), union(gene14, downstream14)]
sub_net14 <- get_subgraph(sub_comb_adj_mat14,comb_adj_mat0,dg_TF14,gene14)

# LUAD
cancer <- "LUAD"
res_file <- file.path("./results/HyperNetWalk", cancer, "results.Rdata")
load(res_file)
deg_mat <- objs$deg_mat
cor_mat <- WGCNA::cor(t(deg_mat))
cor_mat <- abs(cor_mat)
cor_mat[cor_mat < 0.1] <- 0.1
idx <- summary(comb_adj_mat)
g1 <- nodes_ppi[idx$i]
g2 <- nodes_ppi[idx$j]
w <- rep(0.1, length(g1))
hit <- g1 %in% rownames(cor_mat) & g2 %in% colnames(cor_mat)
w[hit] <- cor_mat[cbind(g1[hit], g2[hit])]
comb_adj_mat_luad <- comb_adj_mat
comb_adj_mat_luad@x <- w

gene1 <- "FLI1"
downstream1 <- get_downstream_TFs(gene1,comb_adj_mat0,total_TF)
TF1 <- intersect(downstream1, TF)
dg_TF1 <- degree(GRN_net, mode="out")[TF1]
sub_comb_adj_mat1 <- comb_adj_mat_luad[union(gene1, downstream1), union(gene1, downstream1)]
sub_net1 <- get_subgraph(sub_comb_adj_mat1,comb_adj_mat0,dg_TF1,gene1)

gene2 <- "KDM4C"
downstream2 <- get_downstream_TFs(gene2,comb_adj_mat0,total_TF)
TF2 <- intersect(downstream2, TF)
dg_TF2 <- degree(GRN_net, mode="out")[TF2]
sub_comb_adj_mat2 <- comb_adj_mat_luad[union(gene2, downstream2), union(gene2, downstream2)]
sub_net2 <- get_subgraph(sub_comb_adj_mat2,comb_adj_mat0,dg_TF2,gene2)

gene3 <- "TFAP2A"
downstream3 <- get_downstream_TFs(gene3,comb_adj_mat0,total_TF)
TF3 <- intersect(downstream3, TF)
dg_TF3 <- degree(GRN_net, mode="out")[TF3]
sub_comb_adj_mat3 <- comb_adj_mat_luad[union(gene3, downstream3), union(gene3, downstream3)]
sub_net3 <- get_subgraph(sub_comb_adj_mat3,comb_adj_mat0,dg_TF3,gene3)

gene4 <- "GABPA"
downstream4 <- get_downstream_TFs(gene4,comb_adj_mat0,total_TF)
TF4 <- intersect(downstream4, TF)
dg_TF4 <- degree(GRN_net, mode="out")[TF4]
sub_comb_adj_mat4 <- comb_adj_mat_luad[union(gene4, downstream4), union(gene4, downstream4)]
sub_net4 <- get_subgraph(sub_comb_adj_mat4,comb_adj_mat0,dg_TF4,gene4)

gene5 <- "BRD2"
downstream5 <- get_downstream_TFs(gene5,comb_adj_mat0,total_TF)
TF5 <- intersect(downstream5, TF)
dg_TF5 <- degree(GRN_net, mode="out")[TF5]
sub_comb_adj_mat5 <- comb_adj_mat_luad[union(gene5, downstream5), union(gene5, downstream5)]
sub_net5 <- get_subgraph(sub_comb_adj_mat5,comb_adj_mat0,dg_TF5,gene5)

gene6 <- "MAX"
downstream6 <- get_downstream_TFs(gene6,comb_adj_mat0,total_TF)
TF6 <- intersect(downstream6, TF)
dg_TF6 <- degree(GRN_net, mode="out")[TF6]
sub_comb_adj_mat6 <- comb_adj_mat_luad[union(gene6, downstream6), union(gene6, downstream6)]
sub_net6 <- get_subgraph(sub_comb_adj_mat6,comb_adj_mat0,dg_TF6,gene6)

gene7 <- "EBF1"
downstream7 <- get_downstream_TFs(gene7,comb_adj_mat0,total_TF)
TF7 <- intersect(downstream7, TF)
dg_TF7 <- degree(GRN_net, mode="out")[TF7]

gene8 <- "KDM5B"
downstream8 <- get_downstream_TFs(gene8,comb_adj_mat0,total_TF)
TF8 <- intersect(downstream8, TF)
dg_TF8 <- degree(GRN_net, mode="out")[TF8]
sub_comb_adj_mat8 <- comb_adj_mat_luad[union(gene8, downstream8), union(gene8, downstream8)]
sub_net8 <- get_subgraph(sub_comb_adj_mat8,comb_adj_mat0,dg_TF8,gene8)

gene9 <- "BCOR"
downstream9 <- get_downstream_TFs(gene9,comb_adj_mat0,total_TF)
TF9 <- intersect(downstream9, TF)
dg_TF9 <- degree(GRN_net, mode="out")[TF9]
sub_comb_adj_mat9 <- comb_adj_mat_luad[union(gene9, downstream9), union(gene9, downstream9)]
sub_net9 <- get_subgraph(sub_comb_adj_mat9,comb_adj_mat0,dg_TF9,gene9)

gene10 <- "GATA1"
downstream10 <- get_downstream_TFs(gene10,comb_adj_mat0,total_TF)
TF10 <- intersect(downstream10, TF)
dg_TF10 <- degree(GRN_net, mode="out")[TF10]
sub_comb_adj_mat10 <- comb_adj_mat_luad[union(gene10, downstream10), union(gene10, downstream10)]
sub_net10 <- get_subgraph(sub_comb_adj_mat10,comb_adj_mat0,dg_TF10,gene10)

gene11 <- "TAF3"
downstream11 <- get_downstream_TFs(gene11,comb_adj_mat0,total_TF)
TF11 <- intersect(downstream11, TF)
dg_TF11 <- degree(GRN_net, mode="out")[TF11]
