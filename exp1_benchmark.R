source("model.R")
source("plot_formal.R")
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsci)

# Panel A: box plot
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
fig_dir <- "./figs/exp1_benchmark/"
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
}
CGC <- read.delim("./reference_dg/CGC_Tier1.tsv",
                             header = T,
                             as.is = TRUE)[, 1]
# Personalized results
# Personalized methods: HyperNetWalk, DawnRank, Prodigy, PersonaDrive, PDRWH, PITCH, Frequency, Random
# Metrics: Precision@k, Recall@k, F1@k (k=1,...,Nc)
methods_pers <- list(
  HyperNetWalk = "./results/HyperNetWalk/",
  DawnRank = "./results/DawnRank/",
  PRODIGY = "./results/PRODIGY/",
  PersonaDrive = "./results/PersonaDrive/",
  PDRWH = "./results/PDRWH/",
  PITCH = "./results/PITCH_ad/",
  Frequency = "./results/Frequency/",
  Random = "./results/Random/"
)
res_pers_all <- list()
df_pers <- list()
ref_list <- list()
Nc_min <- 10
top_pers_methods <- list()
for (cancer in cancers){
  mut_data <- get_mut_data(file.path("./data/processed_data",cancer,"mut_data.tsv"),1)
  ref <- get_filter_ref(mut_data, CGC, N_coh = 500)
  ref_list[[cancer]] <- ref
  Nc_min <- min(Nc_min, ref$N_pers)
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
    for (j in 1:nrow(res)) {
      sample_j <- substr(res[j, 1], 1, 12)
      genes_ranking_j <- unlist(strsplit(res[j, 2], ","))
      prediction[[sample_j]] <- genes_ranking_j
      top_df <- data.frame(
        sample = sample_j,
        gene = genes_ranking_j[1:min(ref$N_pers, length(genes_ranking_j))],
        rank = 1:min(ref$N_pers, length(genes_ranking_j))
      )
      if (!is.null(top_pers_methods[[method]])) {
        top_pers_methods[[method]] <- rbind(top_pers_methods[[method]], top_df)
      } else {
        top_pers_methods[[method]] <- top_df
      }
    }
    res_pers_c[[method]] <- prediction
  
  }
  res_pers_all[[cancer]] <- res_pers_c
  PRF_c <- lapply(names(res_pers_c),function(method){
    compute_PRF(res_pers_c[[method]],ref$selected_samples,ref$reference_genes_pers,ref$mut_genes_pers,ref$N_pers)
  })
  names(PRF_c) <- names(res_pers_c)
  df_pers_c <- do.call(rbind, lapply(seq_along(res_pers_c), function(i) {
    data.frame(
      method = names(res_pers_c)[i],
      n_gene = 1:ref$N_pers,
      precision = PRF_c[[i]][1, ],
      recall = PRF_c[[i]][2, ],
      f1 = PRF_c[[i]][3, ]
    )
  }))
  df_pers_c <- df_pers_c %>%
    # pivot_wider(id_cols = method, names_from = n_gene, values_from = c(precision, recall, f1),names_glue = "avg_{.value}@{n_gene}") %>%
    mutate(cancer_type = cancer) %>%
    select(cancer_type, everything())
  df_pers[[cancer]] <- df_pers_c
}

get_prec_nums <- function(df,Nc_min,CGC){
  #添加一列0，1表示gene是否在CGC中
  df <- df %>%
    mutate(is_CGC = ifelse(gene %in% CGC, 1, 0))
  #从1：Nc_min,计算rank=i的行中有多少is_CGC=1的gene
  prec_nums <- sapply(1:Nc_min, function(i) {
    sum(df$rank == i & df$is_CGC == 1)
  })
  # total_nums <- sapply(1:Nc_min, function(i) {
  #   sum(df$rank == i)
  # })
  return(prec_nums)
}
#对top_pers_methods中每个df,用get_prec_nums函数计算prec_nums,并添加一列n_gene表示rank,合并成一个df_comb_num
df_comb_num <- do.call(rbind, lapply(names(top_pers_methods), function(method) {
  df_method <- top_pers_methods[[method]]
  prec_nums <- get_prec_nums(df_method, Nc_min, CGC)
  df <- data.frame(
    method = method,
    n_gene = 1:Nc_min,
    prec_nums = prec_nums
  )
  return(df)
}))
#画折线图，x轴为n_gene，y轴为prec_nums，不同method用不同颜色表示
ggplot(df_comb_num, aes(x = n_gene, y = prec_nums, color = method)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  scale_color_manual(values = c("HyperNetWalk" = "#D94841", "DawnRank" = "#1F77B4", "PRODIGY" = "#2CA02C", "PersonaDrive" = "#FF7F0E", "PDRWH" = "#9467BD", "PITCH" = "#8C564B", "Frequency" = "#E377C2", "Random" = "#7F7F7F")) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    axis.line = element_line(linewidth = 0.8, color = "black"),
    axis.ticks = element_line(linewidth = 0.8, color = "black"),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  labs(
    x = paste0("Top N Genes (N ≤ ", Nc_min, ")"),
    y = "Number of CGC Genes"
  )


df_all_pers <- lapply(df_pers, function(df){
  df %>%
    filter(n_gene == Nc_min) %>%
    select(cancer_type, method, precision, recall, f1) %>%
    rename(avg_precision = precision, avg_recall = recall, avg_f1 = f1)
}) %>%
  bind_rows()

get_all_cancers_plot <- function(df_all,metric,y_label) {
  median_scores <- df_all %>%
    filter(method != "HyperNetWalk") %>%
    group_by(method) %>%
    summarize(med_prec = median(!!sym(metric), na.rm = TRUE)) %>%
    arrange(desc(med_prec)) # 从高到低排序

  method_levels <- c("HyperNetWalk", median_scores$method)
  df_all$method <- factor(df_all$method, levels = method_levels)

  df_all <- df_all %>%
    mutate(is_proposed = ifelse(method == "HyperNetWalk", "Proposed", "Baseline"))

  p <- ggplot(df_all, aes(x = method, y = !!sym(metric))) +
    geom_boxplot(
      aes(fill = is_proposed),
      outlier.shape = NA,
      alpha = 0.6,
      width = 0.62,
      color = "black",
      linewidth = 0.55
    ) +
    # 顶级质感：带白色描边的深灰弹珠
    geom_jitter(
      width = 0.16,
      size = 2.1,
      alpha = 0.8,
      shape = 21,
      fill = "#4D4D4D",
      color = "white",
      stroke = 0.5
    ) +
    scale_fill_manual(values = c("Proposed" = "#D94841", "Baseline" = "#D9D9D9")) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", face = "bold"),
      axis.text.y = element_text(color = "black", size = 12),
      axis.title = element_text(face = "bold", size = 14),
      axis.line = element_line(linewidth = 0.8, color = "black"),
      axis.ticks = element_line(linewidth = 0.8, color = "black"),
      legend.position = "none"
    ) +
    labs(
      x = "Methods",
      y = y_label
    )
  return(p)
}

panel_A1 <- get_all_cancers_plot(df_all_pers,"avg_precision", paste0("Precision@", Nc_min))

# Cohort results
# Cohort methods: HyperNetWalk, dndscv, DriverNet, Subdyquency, DriverRWH, DrvierMP
# Metrics: Precision@k, Recall@k, F1@k, AUROC, AUPRC, pAUROC, pAUPRC
methods_coh <- list(
  HyperNetWalk = "./results/HyperNetWalk/",
  dndscv = "./results/dndscv/",
  DriverNet = "./results/DriverNet/",
  Subdyquency = "./results/Subdyquency/",
  DriverRWH = "./results/DriverRWH/",
  DriverMP = "./results/DriverMP/",
  Frequency = "./results/Frequency/",
  Random = "./results/Random/"
)
res_coh_all <- list()
df_coh <- list()
for (cancer in cancers){
  mut_data <- get_mut_data(file.path("./data/processed_data",cancer,"mut_data.tsv"),1)
  ref <- get_filter_ref(mut_data, CGC, N_coh = 500)
  res_coh_c <- list()
  for (i in 1:length(methods_coh)){
    method <- names(methods_coh)[i]
    res_file <- file.path(methods_coh[[i]],cancer,"genes_ranking_cohort.txt")
    if(!file.exists(res_file)){
      next
    }
    res <- read.table(
      res_file,
      row.names = 1
    )
    mut_genes_sup <- setdiff(ref$mut_genes_coh,rownames(res))
    if (length(mut_genes_sup) > 0) {
      padding <- matrix(0,
                        nrow = length(mut_genes_sup),
                        ncol = ncol(res))
      rownames(padding) <- mut_genes_sup
      colnames(padding) <- colnames(res)
      res <- rbind(res, padding)
    } 
    res_coh_c[[method]] <- res
  }
  res_coh_all[[cancer]] <- res_coh_c
  PRF_c <- lapply(names(res_coh_c),function(method){
    compute_PRF(res_coh_c[[method]],ref$selected_samples,ref$reference_genes_coh,ref$mut_genes_coh,ref$N_coh)
  })
  names(PRF_c) <- names(res_coh_c)
  df_coh_c <- do.call(rbind, lapply(seq_along(res_coh_c), function(i) {
    data.frame(
      method = names(res_coh_c)[i],
      n_gene = 1:ref$N_coh,
      precision = PRF_c[[i]][1, ],
      recall = PRF_c[[i]][2, ],
      f1 = PRF_c[[i]][3, ]
    )
  }))
  df_coh_c <- df_coh_c %>%
    mutate(cancer_type = cancer) %>%
    select(cancer_type, everything())
  df_coh[[cancer]] <- df_coh_c
}
df_all_coh <- lapply(df_coh, function(df){
  df %>%
    filter(n_gene == 100) %>%
    select(cancer_type, method, precision, recall, f1)
}) %>% 
  bind_rows()

panel_A2 <- get_all_cancers_plot(df_all_coh,"precision", "Precision@100")

df_all_coh50 <- lapply(df_coh, function(df){
  df %>%
    filter(n_gene == 50) %>%
    select(cancer_type, method, precision, recall, f1)
}) %>% 
  bind_rows()
panel_sup1 <- get_all_cancers_plot(df_all_coh50,"precision", "Precision@50") 
ggsave(
  filename = file.path(fig_dir, "PanelSup1.pdf"),
  plot = panel_sup1,
  width = 5.5,      # 给底部倾斜的 X 轴标签留足横向空间
  height = 5.5,    # 压低高度，显得精干且大气
  device = "pdf",
  useDingbats = FALSE # 修复 AI 打开 PDF 时图例可能乱码的问题
)
df_all_coh500 <- lapply(df_coh, function(df){
  df %>%
    filter(n_gene == 500) %>%
    select(cancer_type, method, precision, recall, f1)
}) %>% 
  bind_rows()
panel_sup2 <- get_all_cancers_plot(df_all_coh500,"precision", "Precision@500")
ggsave(
  filename = file.path(fig_dir, "PanelSup2.pdf"),
  plot = panel_sup2,
  width = 5.5,      # 给底部倾斜的 X 轴标签留足横向空间
  height = 5.5,    # 压低高度，显得精干且大气
  device = "pdf",
  useDingbats = FALSE # 修复 AI 打开 PDF 时图例可能乱码的问题
)

delta_pers <- df_all_pers %>%
  group_by(cancer_type) %>%
  summarize(
    hyper_score = avg_precision[method == "HyperNetWalk"],
    best_base_score = max(avg_precision[method != "HyperNetWalk"]),
    delta = hyper_score - best_base_score
  ) %>%
  arrange(desc(delta))

delta_coh <- df_all_coh %>%
  group_by(cancer_type) %>%
  summarize(
    hyper_score = precision[method == "HyperNetWalk"],
    best_base_score = max(precision[method != "HyperNetWalk"]),
    delta = hyper_score - best_base_score
  ) %>%
  arrange(desc(delta))

print("Personalized Level Advantage (Delta):")
print(delta_pers)
print("Cohort Level Advantage (Delta):")
print(delta_coh)

# Panel B: heatmap of cohort 
# compute AUPRC, AUROC, pAUPRC, pAUROC for cohort methods
library(pROC)
ROC_data_all <- lapply(names(res_coh_all),function(cancer){
  res_coh_c <- res_coh_all[[cancer]]
  get_ROC(res_coh_c,CGC) %>%
    mutate(cancer_type = cancer)
}) %>% bind_rows()
PR_data_all <- lapply(names(res_coh_all),function(cancer){
  res_coh_c <- res_coh_all[[cancer]]
  get_PR(res_coh_c,CGC) %>%
    mutate(cancer_type = cancer)
}) %>% bind_rows()

auroc_summary <- ROC_data_all %>%
  distinct(cancer_type, method, AUC, pAUC)
auprc_summary <- PR_data_all %>%
  distinct(cancer_type, method, AUPR, pAUPR)

library(patchwork)
library(RColorBrewer) 

df_heat <- auroc_summary %>%
  left_join(auprc_summary, by = c("cancer_type", "method"))

method_ranking <- df_heat %>%
  filter(method != "HyperNetWalk") %>% # 先把你的方法拿出来
  group_by(method) %>%
  summarize(mean_score = mean(pAUPR, na.rm = TRUE)) %>%
  arrange(mean_score) %>% # 升序排列（因为 ggplot 的 Y 轴是从下往上画的）
  pull(method)

method_levels_heat <- c(as.character(method_ranking), "HyperNetWalk")

df_heat$method <- factor(df_heat$method, levels = method_levels_heat)
df_heat$cancer_type <- factor(df_heat$cancer_type, levels = cancers)

heat_colors <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100)

library(forcats)
plot_benchmark_heatmap <- function(data, metric_col, title_text) {
  
  # 数据预处理小魔法：按癌种分组，找出每行的最高分，打上 "is_best" 标签
  plot_data <- data %>%
    group_by(cancer_type) %>%
    mutate(
      is_best = (.data[[metric_col]] == max(.data[[metric_col]], na.rm = TRUE))
    ) %>%
    ungroup()
  
  ggplot(plot_data, aes(x = method, y = fct_rev(cancer_type), fill = .data[[metric_col]])) +
    
    # 矩形色块与白色切割线
    geom_tile(color = "white", linewidth = 0.8) + 
    
    # 🌟 冠军高光逻辑：是本行第一名 (is_best == TRUE) 变白，其余全黑
    geom_text(
      aes(
        label = sprintf("%.1f", .data[[metric_col]]),
        color = is_best 
      ),
      size = 3.5, 
      fontface = "bold",
      show.legend = FALSE
    ) +
    
    # 单色系渐变：浅灰过渡到红
    scale_fill_gradient(low = "#F2F2F2", high = "#D94841", name = title_text) +
    
    # 定义字体颜色：FALSE (非冠军) = 深碳黑，TRUE (冠军) = 纯白
    scale_color_manual(values = c("FALSE" = "#222222", "TRUE" = "white")) + 
    
    coord_fixed(ratio = 1) + 
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", face = "bold"),
      axis.text.y = element_text(color = "black", size = 12),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right",
      legend.key.height = unit(1.0, "cm"),
      legend.key.width = unit(0.3, "cm"),
      legend.title = element_text(face = "bold", size = 11, angle = 90, hjust = 0.5),
      legend.title.position = "right",
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
    ) +
    labs(title = title_text)
}

# ==========================================
# 绘制 4 张图
# ==========================================
p_auc   <- plot_benchmark_heatmap(df_heat, "AUC", "AUROC")
p_pauc  <- plot_benchmark_heatmap(df_heat, "pAUC", "pAUROC")
p_aupr  <- plot_benchmark_heatmap(df_heat, "AUPR", "AUPRC")
p_paupr <- plot_benchmark_heatmap(df_heat, "pAUPR", "pAUPRC")

ggsave(
  filename = file.path(fig_dir, "PanelB.pdf"),
  plot = p_paupr,
  width = 8.5, 
  height = 16,
  device = "pdf",
  useDingbats = FALSE # 修复 AI 打开 PDF 时图例可能乱码的问题
)

pauc_comb <- (p_auc + p_pauc) / (p_aupr + p_paupr) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 18))
ggsave(
  filename = file.path(fig_dir, "AUC_combined.pdf"),
  plot = pauc_comb,
  width = 16, 
  height = 16,
  device = "pdf",
  useDingbats = FALSE # 修复 AI 打开 PDF 时图例可能乱码的问题
)

roc_plots <- plot_roc(ROC_data_all)
pr_plots <- plot_pr(PR_data_all)

selected_cancers <- c("BRCA", "LUAD", "THCA")
# 把选定癌种的 ROC 图放在一起排成一行
roc_selected <- roc_plots[selected_cancers] %>%
  wrap_plots(nrow = 1) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 18))
ggsave(
  filename = file.path(fig_dir, "ROC_selected.pdf"),
  plot = roc_selected,
  width = 16, 
  height = 5.5,
  device = "pdf",
  useDingbats = FALSE # 修复 AI 打开 PDF 时图例可能乱码的问题
)

# 把选定癌种的 PR 图放在一起排成一行
pr_selected <- pr_plots[selected_cancers] %>%
  wrap_plots(nrow = 1) + plot_annotation(tag_levels = "A")
ggsave(
  filename = file.path(fig_dir, "PR_selected.pdf"),
  plot = pr_selected,
  width = 16, 
  height = 5.5,
  device = "pdf",
  useDingbats = FALSE # 修复 AI 打开 PDF 时图例可能乱码的问题
)

# 把所有癌种的 ROC 图放在一起排成4行3列
roc_all <- roc_plots[cancers] %>%
  wrap_plots(nrow = 4, ncol = 3) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 18))
ggsave(
  filename = file.path(fig_dir, "ROC_all.pdf"),
  plot = roc_all,
  width = 16, 
  height = 22,
  device = "pdf",
  useDingbats = FALSE # 修复 AI 打开 PDF 时图例可能乱码的问题
)

# 把所有癌种的 PR 图放在一起排成4行3列
pr_all <- pr_plots[cancers] %>%
  wrap_plots(nrow = 4, ncol = 3) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 18))
ggsave(
  filename = file.path(fig_dir, "PR_all.pdf"),
  plot = pr_all,
  width = 16, 
  height = 22,
  device = "pdf",
  useDingbats = FALSE # 修复 AI 打开 PDF 时图例可能乱码的问题
)

# Panel C: Precision@k curve
prf_pers <- plot_PRF_curve(df_pers)
prf_coh <- plot_PRF_curve(df_coh)
prec_pers_all <- prf_pers$prec_curves %>%
  wrap_plots(nrow = 4, ncol = 3) + 
  plot_layout(guides = "collect") +  # 🌟 关键：一键收集 12 张图的图例
  plot_annotation(tag_levels = "A") & 
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    legend.position = "bottom"       # 🌟 确保唯一图例在最下方
  )
ggsave(
  filename = file.path(fig_dir, "Precision_curve_all.pdf"),
  plot = prec_pers_all,
  width = 16, 
  height = 22,
  device = "pdf",
  useDingbats = FALSE # 修复 AI 打开 PDF 时图例可能乱码的问题
)
prec_coh_all <- prf_coh$prec_curves %>%
  wrap_plots(nrow = 4, ncol = 3) + 
  plot_layout(guides = "collect") +  # 🌟 关键：一键收集 12 张图的图例
  plot_annotation(tag_levels = "A") & 
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    legend.position = "bottom"       # 🌟 确保唯一图例在最下方
  )
ggsave(
  filename = file.path(fig_dir, "Precision_curve_coh_all.pdf"),
  plot = prec_coh_all,
  width = 16, 
  height = 22,
  device = "pdf",
  useDingbats = FALSE # 修复 AI 打开 PDF 时图例可能乱码的问题
)
recall_pers_all <- prf_pers$rec_curves %>%
  wrap_plots(nrow = 4, ncol = 3) + 
  plot_layout(guides = "collect") +  # 🌟 关键：一键收集 12 张图的图例
  plot_annotation(tag_levels = "A") & 
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    legend.position = "bottom"       # 🌟 确保唯一图例在最下方
  )
ggsave(
  filename = file.path(fig_dir, "Recall_curve_all.pdf"),
  plot = recall_pers_all,
  width = 16, 
  height = 22,
  device = "pdf",
  useDingbats = FALSE # 修复 AI 打开 PDF 时图例可能 乱码的问题
)
recall_coh_all <- prf_coh$rec_curves %>%
  wrap_plots(nrow = 4, ncol = 3) + 
  plot_layout(guides = "collect") +  # 🌟 关键：一键收集 12 张图的图例
  plot_annotation(tag_levels = "A") & 
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    legend.position = "bottom"       # 🌟 确保唯一图例在最下方
  )
ggsave(
  filename = file.path(fig_dir, "Recall_curve_coh_all.pdf"),
  plot = recall_coh_all,
  width = 16, 
  height = 22,
  device = "pdf",
  useDingbats = FALSE # 修复 AI 打开 PDF 时图例可能乱码的问题
)
f1_pers_all <- prf_pers$f1_curves %>%
  wrap_plots(nrow = 4, ncol = 3) + 
  plot_layout(guides = "collect") +  # 🌟 关键：一键收集 12 张图的图例
  plot_annotation(tag_levels = "A") & 
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    legend.position = "bottom"       # 🌟 确保唯一图例在最下方
  )
ggsave(
  filename = file.path(fig_dir, "F1_curve_all.pdf"), 
  plot = f1_pers_all,
  width = 16, 
  height = 22,
  device = "pdf",
  useDingbats = FALSE # 修复 AI 打开 PDF 时图例可能乱码的问题
)
f1_coh_all <- prf_coh$f1_curves %>%
  wrap_plots(nrow = 4, ncol = 3) + 
  plot_layout(guides = "collect") +  # 🌟 关键：一键收集 12 张图的图例
  plot_annotation(tag_levels = "A") & 
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    legend.position = "bottom"       # 🌟 确保唯一图例在最下方
  )
ggsave(
  filename = file.path(fig_dir, "F1_curve_coh_all.pdf"),
  plot = f1_coh_all,
  width = 16, 
  height = 22,
  device = "pdf",
  useDingbats = FALSE # 修复 AI 打开 PDF 时图例可能乱码的问题
)


# 把panel A,B,C,D放在一起，其中A，B在一行，C（prec_pers_selected）占一行，D(prec_coh_selected)占一行；注意A中包含两张子图
library(patchwork)
library(Cairo)
panel_A <- wrap_elements(panel_A1 | panel_A2)
ggsave(
  filename = file.path(fig_dir, "PanelA.pdf"),
  plot = panel_A,
  width = 11.5, 
  height = 7,
  units = "cm",
  device = cairo_pdf
)
panel_B <- wrap_elements(p_paupr)
ggsave(
  filename = file.path(fig_dir, "PanelB.pdf"),
  plot = panel_B,
  width = 6.5, 
  height = 7,
  units = "cm",
  device = cairo_pdf
)
panel_C <- wrap_elements(
  (
    prf_coh$prec_curves[["BRCA"]] | 
    prf_coh$prec_curves[["LUAD"]] | 
    prf_coh$prec_curves[["THCA"]]
  ) + 
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.box = "horizontal")
)
ggsave(
  filename = file.path(fig_dir, "PanelC.pdf"),
  plot = panel_C,
  width = 18, 
  height = 6.5,
  units = "cm",
  device = cairo_pdf
)
panel_D <- wrap_elements(
  (
    prf_pers$prec_curves[["BRCA"]] | 
    prf_pers$prec_curves[["LUAD"]] | 
    prf_pers$prec_curves[["THCA"]]
  ) + 
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",legend.box = "horizontal")
)
ggsave(
  filename = file.path(fig_dir, "PanelD.pdf"),
  plot = panel_D,
  width = 18, 
  height = 6.5,
  units = "cm",
  device = cairo_pdf
)

library(patchwork)

# 依然保留这个基础补丁，确保文字旋转的锚点正确，防止在 AI 里错位
base_ratio_theme <- theme(
  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
  plot.margin = margin(t = 10, r = 10, b = 20, l = 10, unit = "pt") # 给大图留足边距
)

# 【Panel A】 大尺寸，比例约 1.7 : 1
panel_A <- wrap_elements((panel_A1 + base_ratio_theme) | (panel_A2 + base_ratio_theme))
ggsave(
  filename = file.path(fig_dir, "PanelA.pdf"),
  plot = panel_A,
  width = 14, height = 8, units = "in", # 使用大英寸
  device = cairo_pdf 
)

# 【Panel B】 热图，比例约 0.8 : 1 (假设你的热图比较窄)
panel_B <- wrap_elements(p_paupr)
ggsave(
  filename = file.path(fig_dir, "PanelB.pdf"),
  plot = panel_B,
  width = 6.5, height = 8, units = "in",
  device = cairo_pdf 
)

# 【Panel C & D】 长横幅，比例约 3 : 1
panel_C <- wrap_elements(
  (
    prf_coh$prec_curves[["BRCA"]] | 
    prf_coh$prec_curves[["LUAD"]] | 
    prf_coh$prec_curves[["THCA"]]
  ) + 
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.box = "horizontal")
)
ggsave(
  filename = file.path(fig_dir, "PanelC.pdf"),
  plot = panel_C,
  width = 18, height = 6, units = "in",
  device = cairo_pdf 
)

# Panel D 同理...
panel_D <- wrap_elements(
  (
    prf_pers$prec_curves[["BRCA"]] | 
    prf_pers$prec_curves[["LUAD"]] | 
    prf_pers$prec_curves[["THCA"]]
  ) + 
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.box = "horizontal")
)
ggsave(
  filename = file.path(fig_dir, "PanelD.pdf"),
  plot = panel_D,
  width = 18, height = 6, units = "in",
  device = cairo_pdf 
)

# Personalized results combined as cohort results
res_coh0_all <- list()
df_coh0 <- list()
methods_pers <- list(
  HyperNetWalk = "./results/HyperNetWalk/",
  DawnRank = "./results/DawnRank/",
  `PRODIGY'(P)'` = "./results/PRODIGY/",
  `PersonaDrive'(H)'` = "./results/PersonaDrive_hrw/",
  PDRWH = "./results/PDRWH/",
  `PITCH'(H)'` = "./results/PITCH_ad_hrw/",
  DriverMP = "./results/DriverMP/"
)
for (cancer in cancers){
   mut_data <- get_mut_data(file.path("./data/processed_data",cancer,"mut_data.tsv"),1)
  ref <- get_filter_ref(mut_data, CGC, N_coh = 500)
  res_coh_c <- list()
  for (i in 1:length(methods_pers)){
    method <- names(methods_pers)[i]
    if (method %in% c("Frequency","Random")){
      next
    }
    res_file <- file.path(methods_pers[[i]],cancer,"genes_ranking_cohort.txt")
    if(!file.exists(res_file)){
      next
    }
    res <- read.table(
      res_file,
      row.names = 1
    )
    mut_genes_sup <- setdiff(ref$mut_genes_coh,rownames(res))
    if (length(mut_genes_sup) > 0) {
      padding <- matrix(0,
                        nrow = length(mut_genes_sup),
                        ncol = ncol(res))
      rownames(padding) <- mut_genes_sup
      colnames(padding) <- colnames(res)
      res <- rbind(res, padding)
    } 
    res_coh_c[[method]] <- res
  }
  res_coh0_all[[cancer]] <- res_coh_c
  PRF_c <- lapply(names(res_coh_c),function(method){
    compute_PRF(res_coh_c[[method]],ref$selected_samples,ref$reference_genes_coh,ref$mut_genes_coh,ref$N_coh)
  })
  names(PRF_c) <- names(res_coh_c)
  df_coh_c <- do.call(rbind, lapply(seq_along(res_coh_c), function(i) {
    data.frame(
      method = names(res_coh_c)[i],
      n_gene = 1:ref$N_coh,
      precision = PRF_c[[i]][1, ],
      recall = PRF_c[[i]][2, ],
      f1 = PRF_c[[i]][3, ]
    )
  }))
  df_coh_c <- df_coh_c %>%
    mutate(cancer_type = cancer) %>%
    select(cancer_type, everything())
  df_coh0[[cancer]] <- df_coh_c
}
df_all_coh0 <- lapply(df_coh0, function(df){
  df %>%
    filter(n_gene == 100) %>%
    select(cancer_type, method, precision, recall, f1)
}) %>% 
  bind_rows()
panel_sup3 <- get_all_cancers_plot(df_all_coh0,"precision", "Precision@100")
ggsave(
  filename = file.path(fig_dir, "PanelSup3.pdf"),
  plot = panel_sup3,
  width = 5.5,      # 给底部倾斜的 X 轴标签留足横向空间
  height = 5.5,    # 压低高度，显得精干且大气
  device = "pdf",
  useDingbats = FALSE # 修复 AI 打开 PDF 时图例可能乱码的问题
)
df_all_coh0_50 <- lapply(df_coh0, function(df){
  df %>%
    filter(n_gene == 50) %>%
    select(cancer_type, method, precision, recall, f1)
}) %>% 
  bind_rows()
panel_sup4 <- get_all_cancers_plot(df_all_coh0_50,"precision", "Precision@50")
ggsave(
  filename = file.path(fig_dir, "PanelSup4.pdf"),
  plot = panel_sup4,
  width = 5.5,      # 给底部倾斜的 X 轴标签留足横向空间
  height = 5.5,    # 压低高度，显得精干且大气
  device = "pdf",
  useDingbats = FALSE # 修复 AI 打开 PDF 时图例可能乱码的问题
)
df_all_coh0_500 <- lapply(df_coh0, function(df){
  df %>%
    filter(n_gene == 500) %>%
    select(cancer_type, method, precision, recall, f1)
}) %>% 
  bind_rows()
panel_sup5 <- get_all_cancers_plot(df_all_coh0_500,"precision", "Precision@500")
ggsave(
  filename = file.path(fig_dir, "PanelSup5.pdf"),
  plot = panel_sup5,
  width = 5.5,      # 给底部倾斜的 X 轴标签留足横向空间
  height = 5.5,    # 压低高度，显得精干且大气
  device = "pdf",
  useDingbats = FALSE # 修复 AI 打开 PDF 时图例可能乱码的问题
)

#把panel_sup3, panel_sup4, panel_sup5放在一行
panel_sup_comb <- panel_sup3 | panel_sup4 | panel_sup5 + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 18))
ggsave(
  filename = file.path(fig_dir, "PanelSup_comb.pdf"),
  plot = panel_sup_comb,
  width = 18,
  height = 5.5,
  device = "pdf",
  useDingbats = FALSE # 修复 AI 打开 PDF 时图例可能乱码的问题
)

ROC_data_all0 <- lapply(names(res_coh0_all),function(cancer){
  res_coh_c <- res_coh0_all[[cancer]]
  get_ROC(res_coh_c,CGC) %>%
    mutate(cancer_type = cancer)
}) %>% bind_rows()
PR_data_all0 <- lapply(names(res_coh0_all),function(cancer){
  res_coh_c <- res_coh0_all[[cancer]]
  get_PR(res_coh_c,CGC) %>%
    mutate(cancer_type = cancer)
}) %>% bind_rows()
auroc_summary0 <- ROC_data_all0 %>%
  distinct(cancer_type, method, AUC, pAUC)
auprc_summary0 <- PR_data_all0 %>%
  distinct(cancer_type, method, AUPR, pAUPR)
df_heat0 <- auroc_summary0 %>%
  left_join(auprc_summary0, by = c("cancer_type", "method"))
method_ranking0 <- df_heat0 %>%
  filter(method != "HyperNetWalk") %>% # 先把你的方法拿出来
  group_by(method) %>%
  summarize(mean_score = mean(pAUPR, na.rm = TRUE)) %>%
  arrange(mean_score) %>% # 升序排列（因为 ggplot 的 Y 轴是从下往上画的）
  pull(method)

method_levels_heat0 <- c(as.character(method_ranking0), "HyperNetWalk")

df_heat0$method <- factor(df_heat0$method, levels = method_levels_heat0)
df_heat0$cancer_type <- factor(df_heat0$cancer_type, levels = cancers)
p_auc0   <- plot_benchmark_heatmap(df_heat0, "AUC", "AUROC")
p_pauc0  <- plot_benchmark_heatmap(df_heat0, "pAUC", "pAUROC")
p_aupr0  <- plot_benchmark_heatmap(df_heat0, "AUPR", "AUPRC")
p_paupr0 <- plot_benchmark_heatmap(df_heat0, "pAUPR", "pAUPRC")

pauc_comb0 <- (p_auc0 + p_pauc0) / (p_aupr0 + p_paupr0) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 18))
ggsave(
  filename = file.path(fig_dir, "AUC_combined_per2coh.pdf"),
  plot = pauc_comb0,
  width = 16, 
  height = 16,
  device = "pdf",
  useDingbats = FALSE # 修复 AI 打开 PDF 时图例可能乱码的问题
)

roc_plots0 <- plot_roc(ROC_data_all0)
pr_plots0 <- plot_pr(PR_data_all0)

roc_all0 <- roc_plots0[cancers] %>%
  wrap_plots(nrow = 4, ncol = 3) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 18))
ggsave(
  filename = file.path(fig_dir, "ROC_all_pers2coh.pdf"),
  plot = roc_all0,
  width = 16, 
  height = 22,
  device = "pdf",
  useDingbats = FALSE # 修复 AI 打开 PDF 时图例可能乱码的问题
)

# 把所有癌种的 PR 图放在一起排成4行3列
pr_all0 <- pr_plots0[cancers] %>%
  wrap_plots(nrow = 4, ncol = 3) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 18))
ggsave(
  filename = file.path(fig_dir, "PR_all_pers2coh.pdf"),
  plot = pr_all0,
  width = 16, 
  height = 22,
  device = "pdf",
  useDingbats = FALSE # 修复 AI 打开 PDF 时图例可能乱码的问题
)

prf_coh0 <- plot_PRF_curve(df_coh0)
prec_coh_all0 <- prf_coh0$prec_curves %>%
  wrap_plots(nrow = 4, ncol = 3) + 
  plot_layout(guides = "collect") +  # 🌟 关键：一键收集 12 张图的图例
  plot_annotation(tag_levels = "A") & 
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    legend.position = "bottom"       # 🌟 确保唯一图例在最下方
  )
ggsave(
  filename = file.path(fig_dir, "Precision_curve_coh_all_pers2coh.pdf"),
  plot = prec_coh_all0,
  width = 16, 
  height = 22,
  device = "pdf",
  useDingbats = FALSE # 修复 AI 打开 PDF 时图例可能乱码的问题
)
recall_coh_all0 <- prf_coh0$rec_curves %>%
  wrap_plots(nrow = 4, ncol = 3) + 
  plot_layout(guides = "collect") +  # 🌟 关键：一键收集 12 张图的图例
  plot_annotation(tag_levels = "A") & 
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    legend.position = "bottom"       # 🌟 确保唯一图例在最下方
  )
ggsave(
  filename = file.path(fig_dir, "Recall_curve_coh_all_pers2coh.pdf"),
  plot = recall_coh_all0,
  width = 16, 
  height = 22,
  device = "pdf",
  useDingbats = FALSE # 修复 AI 打开 PDF 时图例可能乱码的问题
)
f1_coh_all0 <- prf_coh0$f1_curves %>%
  wrap_plots(nrow = 4, ncol = 3) + 
  plot_layout(guides = "collect") +  # 🌟 关键：一键收集 12 张图的图例
  plot_annotation(tag_levels = "A") & 
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    legend.position = "bottom"       # 🌟 确保唯一图例在最下方
  )
ggsave(
  filename = file.path(fig_dir, "F1_curve_coh_all_pers2coh.pdf"),
  plot = f1_coh_all0,
  width = 16, 
  height = 22,
  device = "pdf",
  useDingbats = FALSE # 修复 AI 打开 PDF 时图例可能乱码的问题
)