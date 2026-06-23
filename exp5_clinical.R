source("model.R")
source("plot_formal.R")
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggsci)
library(bindrcpp)

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
fig_dir <- "./figs/exp5_clinical"
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir,recursive = TRUE)
}
CGC <- read.delim("./reference_dg/CGC_Tier1.tsv",
                             header = T,
                             as.is = TRUE)[, 1]
methods <- list(
    HyperNetWalk = "./results/HyperNetWalk/",
    DriverRWH = "./results/DriverRWH/",
    DriverMP = "./results/DriverMP/",
    DawnRank = "./results/DawnRank/",
    PDRWH = "./results/PDRWH/"
)

res_coh <- list()
res_pers <- list()
top_pers <- list()
for (cancer in cancers){
    mut_file <- mut_file <- file.path("./data/processed_data",cancer,"mut_data.tsv")
    mut_data <- get_mut_data(mut_file,1)
    colnames(mut_data) <- substr(colnames(mut_data), 1, 12)
    mut_genes <- rownames(mut_data)
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

#读入dgidb文件
dgidb_df <- read.delim("./data/dgidb_interactions.tsv", header = TRUE, stringsAsFactors = FALSE)
#dag=druggable genes
broad_dag <- unique(dgidb_df$gene_name)
approved_dag <- dgidb_df %>% 
    filter(approved==TRUE) %>%
    pull(gene_name) %>%
    unique()
anticancer_dag <- dgidb_df %>% 
    filter(anti_neoplastic==TRUE) %>%
    pull(gene_name) %>%
    unique()

dag_df <- data.frame()
dag_hrw_vecs <- list()
for (method in names(methods)){
    for (cancer in cancers){
        top_c <- rownames(res_coh[[method]][[cancer]])[1:100]
        prop_broad_c <- sum(top_c %in% broad_dag)/length(top_c)
        prop_approved_c <- sum(top_c %in% approved_dag)/length(top_c)
        prop_anticancer_c <- sum(top_c %in% anticancer_dag)/length(top_c)
        if (method %in% names(top_pers)){
            top_p <- top_pers[[method]] %>%
                filter(cancer == !!cancer) %>%
                mutate(Is_Broad = gene %in% broad_dag,
                       Is_Approved = gene %in% approved_dag,
                       Is_Anticancer = gene %in% anticancer_dag)
            #把样本相同的行合并成一行,只保留cancer,sample,num_broad,num_approved,num_anticancer三列，num_broad是该样本中命中的broad基因数量，num_approved是该样本中命中的approved基因数量，num_anticancer是该样本中命中的anticancer基因数量
            if (method == "HyperNetWalk"){
                top_p_summary <- top_p %>%
                    group_by(sample) %>%
                    summarise(num_broad = sum(Is_Broad),
                            num_approved = sum(Is_Approved),
                            num_anticancer = sum(Is_Anticancer))
                dag_hrw_vecs[[cancer]] <- top_p_summary
            }
            #统计top_p中broad、approved、anticancer的数量，注意不同样本间基因可能重复
            prop_broad_p <- top_p %>%
                filter(Is_Broad) %>%
                distinct(gene) %>%
                nrow() / length(unique(top_p$gene))
            prop_approved_p <- top_p %>%
                filter(Is_Approved) %>%
                distinct(gene) %>%
                nrow() / length(unique(top_p$gene))
            prop_anticancer_p <- top_p %>%
                filter(Is_Anticancer) %>%
                distinct(gene) %>%
                nrow() / length(unique(top_p$gene))
            top_p_5union <- top_p %>%
                filter(rank <= 5)
            prop_broad_p_5union <- top_p_5union %>%
                filter(Is_Broad) %>%
                distinct(gene) %>%
                nrow() / length(unique(top_p_5union$gene))
            prop_approved_p_5union <- top_p_5union %>%
                filter(Is_Approved) %>%
                distinct(gene) %>%
                nrow() / length(unique(top_p_5union$gene))
            prop_anticancer_p_5union <- top_p_5union %>%
                filter(Is_Anticancer) %>%
                distinct(gene) %>%
                nrow() / length(unique(top_p_5union$gene))
            
            hit_broad_samples <- top_p %>%
                filter(Is_Broad) %>%
                distinct(sample) %>%
                nrow()
            hit_approved_samples <- top_p %>%
                filter(Is_Approved) %>%
                distinct(sample) %>%
                nrow()
            hit_anticancer_samples <- top_p %>%
                filter(Is_Anticancer) %>%
                distinct(sample) %>%
                nrow()
            hit_broad_samples_5union <- top_p_5union %>%
                filter(Is_Broad) %>%
                distinct(sample) %>%
                nrow()
            hit_approved_samples_5union <- top_p_5union %>%
                filter(Is_Approved) %>%
                distinct(sample) %>%
                nrow()
            hit_anticancer_samples_5union <- top_p_5union %>%
                filter(Is_Anticancer) %>%
                distinct(sample) %>%
                nrow()
            # 计算所有样本中命中的broad基因的平均数量，注意有的样本可能命中多个broad基因，有的样本可能没有命中broad基因，平均数量应该是所有样本中命中broad基因的总数量除以样本总数
            hit_broad_mean <- top_p %>%
                filter(Is_Broad) %>%
                group_by(sample) %>%
                summarise(num_broad = n()) %>%
                summarise(mean_broad = sum(num_broad)/length(unique(top_p$sample))) %>%
                pull(mean_broad)
            hit_approved_mean <- top_p %>%
                filter(Is_Approved) %>%
                group_by(sample) %>%
                summarise(num_approved = n()) %>%
                summarise(mean_approved = sum(num_approved)/length(unique(top_p$sample))) %>%
                pull(mean_approved)
            hit_anticancer_mean <- top_p %>%
                filter(Is_Anticancer) %>%
                group_by(sample) %>%
                summarise(num_anticancer = n()) %>%
                summarise(mean_anticancer = sum(num_anticancer)/length(unique(top_p$sample))) %>%
                pull(mean_anticancer)
            dag_df <- rbind(dag_df,data.frame(
                method = method,
                cancer = cancer,
                prop_broad_c = prop_broad_c,
                prop_approved_c = prop_approved_c,
                prop_anticancer_c = prop_anticancer_c,
                prop_broad_p = prop_broad_p,
                prop_approved_p = prop_approved_p,
                prop_anticancer_p = prop_anticancer_p,
                prop_broad_p_5union = prop_broad_p_5union,
                prop_approved_p_5union = prop_approved_p_5union,
                prop_anticancer_p_5union = prop_anticancer_p_5union,
                hit_broad_samples = hit_broad_samples,
                hit_approved_samples = hit_approved_samples,
                hit_anticancer_samples = hit_anticancer_samples,
                hit_broad_samples_5union = hit_broad_samples_5union,
                hit_approved_samples_5union = hit_approved_samples_5union,
                hit_anticancer_samples_5union = hit_anticancer_samples_5union,
                hit_broad_mean = hit_broad_mean,
                hit_approved_mean = hit_approved_mean,
                hit_anticancer_mean = hit_anticancer_mean
            ))
        } else {
           dag_df <- rbind(dag_df,data.frame(
                method = method,
                cancer = cancer,
                prop_broad_c = prop_broad_c,
                prop_approved_c = prop_approved_c,
                prop_anticancer_c = prop_anticancer_c,
                prop_broad_p = NA,
                prop_approved_p = NA,
                prop_anticancer_p = NA,
                prop_broad_p_5union = NA,
                prop_approved_p_5union = NA,
                prop_anticancer_p_5union = NA,
                hit_broad_samples = NA,
                hit_approved_samples = NA,
                hit_anticancer_samples = NA,
                hit_broad_samples_5union = NA,
                hit_approved_samples_5union = NA,
                hit_anticancer_samples_5union = NA,
                hit_broad_mean = NA,
                hit_approved_mean = NA,
                hit_anticancer_mean = NA
            ))
        }
    }
}
write.table(dag_df, file = file.path(fig_dir,"dag_enrichment_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

dag_hrw_df <- dag_df %>%
    filter(method == "HyperNetWalk")

plot_data <- dag_hrw_df %>%
  select(cancer, prop_broad_c, prop_approved_c, prop_anticancer_c) %>%
  rename(
    Broad = prop_broad_c,
    Approved = prop_approved_c,
    `Anti-neoplastic` = prop_anticancer_c  # 修正为 DGIdb 原本标注
  ) %>%
  pivot_longer(cols = c(Broad, Approved, `Anti-neoplastic`), 
               names_to = "Category", 
               values_to = "Proportion")
cancer_order <- sort(unique(plot_data$cancer))
plot_data <- plot_data %>%
  mutate(
    cancer = factor(cancer, levels = cancer_order),
    # 调整图例顺序，符合逻辑递进
    Category = factor(Category, levels = c("Broad", "Anti-neoplastic", "Approved"))
  )
p_bar <- ggplot(plot_data, aes(x = cancer, y = Proportion, fill = Category)) +
  # width = 0.7 让柱子变瘦，position_dodge(0.8) 拉开组内间距，极其提升高级感
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c(
    "Approved" = "#B3DE69", 
    "Anti-neoplastic" = "#80B1D3", 
    "Broad" = "#FB8072"
  )) +
  labs(
    x = "Cancer Types",
    y = "Proportion of Driver Genes (%)",
    fill = "Drug-Gene Interaction"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8),
    legend.position = "top",          # 顶刊通常把宽图例放在上方
    legend.title = element_blank(),   # 隐去多余的图例标题
    legend.background = element_blank()
  )
p_dot <- ggplot(plot_data, aes(x = cancer, y = Proportion)) +
  geom_linerange(
    data = plot_data %>% group_by(cancer) %>% summarise(min = min(Proportion), max = max(Proportion)),
    aes(x = cancer, ymin = min, ymax = max),
    inherit.aes = FALSE, color = "grey80", linewidth = 1.2
  ) +
  geom_point(aes(fill = Category, shape = Category), size = 4.5, color = "white", stroke = 0.6) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0.1, 0.1))
  ) +
  scale_shape_manual(values = c(21, 22, 24)) +
  scale_fill_manual(values = c(
    "Approved" = "#B3DE69", 
    "Anti-neoplastic" = "#80B1D3", 
    "Broad" = "#FB8072"
  )) +
  labs(
    x = "Cancer Types",
    y = "DGIdb-supported genes in top 100 (%)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold"),
    axis.text.y = element_text(color = "black"),
    panel.grid.major.y = element_line(color = "grey92", linetype = "dashed"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.background = element_blank(),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8)
  )
ggsave(file.path(fig_dir,"DGIdb_Proportion_DotPlot.pdf"), p_dot, width = 8.5, height = 5.5)

plot_data_p <- dag_hrw_df %>%
  select(cancer, prop_broad_p, prop_approved_p, prop_anticancer_p) %>%
  rename(
    Broad = prop_broad_p,
    Approved = prop_approved_p,
    `Anti-neoplastic` = prop_anticancer_p
  ) %>%
  pivot_longer(cols = c(Broad, Approved, `Anti-neoplastic`), 
               names_to = "Category", 
               values_to = "Proportion")
plot_data_p <- plot_data_p %>%
  mutate(
    cancer = factor(cancer, levels = cancer_order),
    Category = factor(Category, levels = c("Broad", "Anti-neoplastic", "Approved"))
  )
p_dot_p <- ggplot(plot_data_p, aes(x = cancer, y = Proportion)) +
  geom_linerange(
    data = plot_data_p %>% group_by(cancer) %>% summarise(min = min(Proportion, na.rm = TRUE), max = max(Proportion, na.rm = TRUE)),
    aes(x = cancer, ymin = min, ymax = max),
    inherit.aes = FALSE, color = "grey80", linewidth = 1.2
  ) +
  geom_point(aes(fill = Category, shape = Category), size = 4.5, color = "white", stroke = 0.6) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0.1, 0.1))
  ) +
  scale_shape_manual(values = c(21, 22, 24)) +
  scale_fill_manual(values = c(
    "Approved" = "#B3DE69", 
    "Anti-neoplastic" = "#80B1D3", 
    "Broad" = "#FB8072"
  )) +
  labs(
    x = "Cancer Types",
    y = "DGIdb-supported genes in top 100 (%)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold"),
    axis.text.y = element_text(color = "black"),
    panel.grid.major.y = element_line(color = "grey92", linetype = "dashed"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.background = element_blank(),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8)
  )
ggsave(file.path(fig_dir,"DGIdb_Proportion_PatientLevel_DotPlot.pdf"), p_dot_p, width = 8.5, height = 5.5)

patient_coverage_df <- bind_rows(dag_hrw_vecs, .id = "cancer") %>%
  mutate(
    Hit_Category = case_when(
      num_anticancer == 0 ~ "0",
      num_anticancer == 1 ~ "1",
      num_anticancer == 2 ~ "2",
      num_anticancer >= 3 ~ "\u2265 3"
    ),
    Hit_Category = factor(Hit_Category, levels = c("\u2265 3", "2", "1", "0"))
  )
coverage_summary <- patient_coverage_df %>%
  group_by(cancer, Hit_Category) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(cancer) %>%
  mutate(Proportion = Count / sum(Count))
coverage_summary$cancer <- factor(coverage_summary$cancer, levels = cancer_order)
p_stacked <- ggplot(coverage_summary, aes(x = cancer, y = Proportion, fill = Hit_Category)) +
  # position = "fill" 自动将高度拉伸至 100%
  geom_bar(stat = "identity", position = "fill", width = 0.75, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c(
    "0" = "#D9D9D9",        # 你的色卡：高级灰 (无靶点)
    "1" = "#CCEBC5",        # 你的色卡：浅茶绿
    "2" = "#8DD3C7",        # 你的色卡：薄荷绿
    "\u2265 3" = "#80B1D3"  # 你的色卡：天空蓝 (靶点极其丰富)
  )) +
  
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  labs(
    x = "Cancer Types",
    y = "Proportion of Patients (%)",
    fill = "Number of Actionable\nTargets per Patient"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold"),
    axis.text.y = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 11)
  )

plot_data_vln <- bind_rows(dag_hrw_vecs, .id = "cancer") %>%
  rename(
    Broad = num_broad,
    Approved = num_approved,
    `Anti-neoplastic` = num_anticancer # 保持之前严谨的 DGIdb 术语
  ) %>%
  pivot_longer(
    cols = c(Broad, Approved, `Anti-neoplastic`), 
    names_to = "Category", 
    values_to = "Hit_Count"
  )
plot_data_vln <- plot_data_vln %>%
  mutate(
    cancer = factor(cancer, levels = cancer_order),
    Category = factor(Category, levels = c("Broad", "Anti-neoplastic", "Approved"))
  )
p_boxplot_clinical <- ggplot(plot_data_vln, aes(x = cancer, y = Hit_Count, fill = Category)) +
  
  # 【核心修改】：直接使用箱线图，去除所有多余的花哨设计
  geom_boxplot(
    position = position_dodge(width = 0.75), 
    width = 0.6,             # 箱子不要太胖，保持高级的纤细感
    alpha = 0.85,            # 极轻微的透明度，让重叠的边缘更柔和
    color = "black",         # 边框为黑色
    linewidth = 0.5,         # 边框粗细
    outlier.size = 1,        # 缩小离群点的大小，防止太抢眼
    outlier.alpha = 0.5,     # 离群点半透明，减弱噪音感
    outlier.shape = 21       # 离群点为空心圆，看起来更透气
  ) +
  
  scale_fill_manual(values = c(
    "Approved" = "#B3DE69", 
    "Anti-neoplastic" = "#80B1D3", 
    "Broad" = "#FB8072"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "Cancer Types",
    y = "Number of Actionable Targets per Patient",
    fill = "Drug-Gene Interaction Level"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold"),
    axis.text.y = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8),
    # 保留虚线，这对于箱线图比较中位数极其有帮助
    panel.grid.major.y = element_line(color = "grey90", linetype = "dashed"),
    legend.position = "top",          
    legend.title = element_text(face = "bold"), 
    legend.background = element_blank()
  )
ggsave(file.path(fig_dir,"DGIdb_Hit_Count_BoxPlot.pdf"), p_boxplot_clinical, width = 8.5, height = 5.5)

#读入OncoKB文件
oncokb_df <- read.table("./data/oncokb_biomarker_drug_associations.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
pan_solid_tumors <- c(
    "All Solid Tumors", 
    "All Solid Tumors (excluding Bladder Cancer)", # 你的12个里没有BLCA，所以适用
    "All Solid Tumors (excluding Cholangiocarcinoma, Bladder Cancer)" # 同样适用
)
oncokb_dict <- list(
  BRCA = c("Breast Cancer", pan_solid_tumors, "All Solid Tumors (excluding Colorectal Cancer)"),
  
  COAD = c("Colorectal Cancer", 
           "Esophagogastric Cancer, Tubular Adenoma of the Colon, Gastrointestinal Neuroendocrine Tumors of the Esophagus/Stomach, Anal Cancer", 
           "All Solid Tumors"), # COAD必须排除 "excluding Colorectal Cancer" 的泛癌标签
  
  HNSC = c("Head and Neck Squamous Cell Carcinoma", pan_solid_tumors, "All Solid Tumors (excluding Colorectal Cancer)"),
  
  KIRC = c("Renal Cell Carcinoma", pan_solid_tumors, "All Solid Tumors (excluding Colorectal Cancer)"),
  
  KIRP = c("Renal Cell Carcinoma", pan_solid_tumors, "All Solid Tumors (excluding Colorectal Cancer)"),
  
  LIHC = c("Hepatobiliary Cancer", "Biliary Tract, Hepatobiliary Cancer", pan_solid_tumors, "All Solid Tumors (excluding Colorectal Cancer)"),
  
  LUAD = c("Non-Small Cell Lung Cancer", 
           "Non-Small Cell Lung Cancer (excluding Lung Squamous Cell Carcinoma)", 
           "Melanoma, Non-Small Cell Lung Cancer, Low-Grade Serous Ovarian Cancer", 
           pan_solid_tumors, "All Solid Tumors (excluding Colorectal Cancer)"),
  
  LUSC = c("Non-Small Cell Lung Cancer", 
           "Melanoma, Non-Small Cell Lung Cancer, Low-Grade Serous Ovarian Cancer", 
           pan_solid_tumors, "All Solid Tumors (excluding Colorectal Cancer)"),
  
  PRAD = c("Prostate Cancer", pan_solid_tumors, "All Solid Tumors (excluding Colorectal Cancer)"),
  
  STAD = c("Esophagogastric Cancer", "Esophagogastric Adenocarcinoma", 
           "Esophagogastric Cancer, Tubular Adenoma of the Colon, Gastrointestinal Neuroendocrine Tumors of the Esophagus/Stomach, Anal Cancer", 
           pan_solid_tumors, "All Solid Tumors (excluding Colorectal Cancer)"),
  
  THCA = c("Thyroid Cancer", "Medullary Thyroid Cancer", "Anaplastic Thyroid Cancer", 
           pan_solid_tumors, "All Solid Tumors (excluding Colorectal Cancer)"),
  
  UCEC = c("Endometrial Cancer", "Endometrial Cancer, Ovarian Cancer", 
           "Uterine Serous Carcinoma/Uterine Papillary Serous Carcinoma", 
           pan_solid_tumors, "All Solid Tumors (excluding Colorectal Cancer)")
)
abg_summary_df <- data.frame()
abg_res_df <- data.frame()
abg_hrw_vecs <- list()
for (cancer in cancers){
    abg_genes <- oncokb_df %>%
        filter(Cancer.Types %in% oncokb_dict[[cancer]]) %>%
        filter(!grepl("R", Level)) %>%
        group_by(Gene) %>%
        summarise(Highest_Evidence_Level = min(Level,na.rm = TRUE)) %>%
        ungroup()
    abg_summary_c <- abg_genes %>%
        group_by(Highest_Evidence_Level) %>%
        summarise(Num_Genes = n(), .groups = "drop") %>%
        complete(Highest_Evidence_Level = c("1", "2", "3", "4"), fill = list(Num_Genes = 0)) %>%
        pivot_wider(names_from = Highest_Evidence_Level, values_from = Num_Genes) %>%
        mutate(cancer = .env$cancer) %>% 
        select(cancer, `1`, `2`, `3`, `4`)
    abg_summary_df <- rbind(abg_summary_df, abg_summary_c)

    for (method in names(methods)){
        top_c <- rownames(res_coh[[method]][[cancer]])[1:100]
        num_level1_c <- sum(top_c %in% abg_genes$Gene[abg_genes$Highest_Evidence_Level == "1"])
        num_level2_c <- sum(top_c %in% abg_genes$Gene[abg_genes$Highest_Evidence_Level == "2"])
        num_level3_c <- sum(top_c %in% abg_genes$Gene[abg_genes$Highest_Evidence_Level == "3"])
        num_level4_c <- sum(top_c %in% abg_genes$Gene[abg_genes$Highest_Evidence_Level == "4"])
        if (method %in% names(top_pers)){
            top_p <- top_pers[[method]] %>%
                filter(cancer == !!cancer) %>%
                mutate(Highest_Evidence_Level = case_when(
                    gene %in% abg_genes$Gene[abg_genes$Highest_Evidence_Level == "1"] ~ "1",
                    gene %in% abg_genes$Gene[abg_genes$Highest_Evidence_Level == "2"] ~ "2",
                    gene %in% abg_genes$Gene[abg_genes$Highest_Evidence_Level == "3"] ~ "3",
                    gene %in% abg_genes$Gene[abg_genes$Highest_Evidence_Level == "4"] ~ "4",
                    TRUE ~ NA_character_
                ))
            if (method == "HyperNetWalk"){
                top_p_summary <- top_p %>%
                    group_by(sample) %>%
                    summarise(num_level1 = sum(Highest_Evidence_Level == "1", na.rm = TRUE),
                            num_level2 = sum(Highest_Evidence_Level == "2", na.rm = TRUE),
                            num_level3 = sum(Highest_Evidence_Level == "3", na.rm = TRUE),
                            num_level4 = sum(Highest_Evidence_Level == "4", na.rm = TRUE))
                abg_hrw_vecs[[cancer]] <- top_p_summary
            }
            num_level1_p <- top_p %>%
                filter(Highest_Evidence_Level == "1") %>%
                distinct(gene) %>%
                nrow()
            num_level2_p <- top_p %>%
                filter(Highest_Evidence_Level == "2") %>%
                distinct(gene) %>%
                nrow()
            num_level3_p <- top_p %>%   
                filter(Highest_Evidence_Level == "3") %>%
                distinct(gene) %>%
                nrow()
            num_level4_p <- top_p %>%
                filter(Highest_Evidence_Level == "4") %>%
                distinct(gene) %>%
                nrow()
            hit_level1_samples <- top_p %>%
                filter(Highest_Evidence_Level == "1") %>%
                distinct(sample) %>%
                nrow()
            hit_level2_samples <- top_p %>%
                filter(Highest_Evidence_Level == "2") %>%
                distinct(sample) %>%
                nrow()
            hit_level3_samples <- top_p %>%
                filter(Highest_Evidence_Level == "3") %>%
                distinct(sample) %>%
                nrow()
            hit_level4_samples <- top_p %>%
                filter(Highest_Evidence_Level == "4") %>%
                distinct(sample) %>%
                nrow()
            abg_res_df <- rbind(abg_res_df, data.frame(
                cancer = cancer,
                method = method,
                num_level1_c = num_level1_c,
                num_level2_c = num_level2_c,
                num_level3_c = num_level3_c,
                num_level4_c = num_level4_c,
                num_level1_p = num_level1_p,
                num_level2_p = num_level2_p,
                num_level3_p = num_level3_p,
                num_level4_p = num_level4_p,
                hit_level1_samples = hit_level1_samples,
                hit_level2_samples = hit_level2_samples,
                hit_level3_samples = hit_level3_samples,
                hit_level4_samples = hit_level4_samples
            ))
        } else {
            abg_res_df <- rbind(abg_res_df, data.frame(
                cancer = cancer,
                method = method,
                num_level1_c = num_level1_c,
                num_level2_c = num_level2_c,
                num_level3_c = num_level3_c,
                num_level4_c = num_level4_c,
                num_level1_p = NA,
                num_level2_p = NA,
                num_level3_p = NA,
                num_level4_p = NA,
                hit_level1_samples = NA,
                hit_level2_samples = NA,
                hit_level3_samples = NA,
                hit_level4_samples = NA
            ))
        }
    }
}
write.table(abg_summary_df, file = file.path(fig_dir,"abg_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
#把abg_res_df按照cancer和method排序
abg_res_df <- abg_res_df %>%
    arrange(method)

write.table(abg_res_df, file = file.path(fig_dir,"abg_res.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

abg_total <- abg_summary_df %>%
  rename(
    Total_L1 = `1`,
    Total_L2 = `2`,
    Total_L3 = `3`,
    Total_L4 = `4`
  )

plot_data_recall <- abg_res_df %>%
  filter(method == "HyperNetWalk") %>%
  left_join(abg_total, by = "cancer") %>%
  mutate(
    denom_tier1 = Total_L1 + Total_L2,
    # 分子：群体 Top 100 中实际命中的 L1 和 L2 基因数
    num_tier1_c = num_level1_c + num_level2_c,
    # 正确的召回率
    Recall_Tier1 = num_tier1_c / denom_tier1,
    
    # Tier 2: Investigational (Level 3 + 4)
    denom_tier2 = Total_L3 + Total_L4,
    num_tier2_c = num_level3_c + num_level4_c,
    Recall_Tier2 = num_tier2_c / denom_tier2
  ) %>%
  select(cancer, Recall_Tier1, Recall_Tier2) %>%
  # 宽表转长表，用于画图
  pivot_longer(
    cols = c(Recall_Tier1, Recall_Tier2), 
    names_to = "Tier", 
    values_to = "Recall"
  ) %>%
  mutate(
    Tier = case_when(
      Tier == "Recall_Tier1" ~ "Tier 1: Standard of Care (Level 1/2)",
      Tier == "Recall_Tier2" ~ "Tier 2: Investigational (Level 3/4)"
    ),
    Tier = factor(Tier, levels = c(
      "Tier 1: Standard of Care (Level 1/2)",
      "Tier 2: Investigational (Level 3/4)"
    ))
  )
cancer_order <- sort(unique(plot_data_recall$cancer))
plot_data_recall$cancer <- factor(plot_data_recall$cancer, levels = cancer_order)
p_recall_dot <- ggplot(plot_data_recall, aes(x = cancer, y = Recall)) +
  geom_linerange(
    data = plot_data_recall %>% group_by(cancer) %>% summarise(min = min(Recall), max = max(Recall)),
    aes(x = cancer, ymin = min, ymax = max),
    inherit.aes = FALSE, color = "grey80", linewidth = 1.2
  ) +
  geom_point(aes(fill = Tier, shape = Tier), size = 4.5, color = "white", stroke = 0.6) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0.1, 0.1)) 
  ) +
  scale_shape_manual(values = c(21, 22)) +
  scale_fill_manual(values = c(
    "Tier 1: Standard of Care (Level 1/2)" = "#FB8072",  # 珊瑚粉
    "Tier 2: Investigational (Level 3/4)"  = "#80B1D3"   # 天空蓝
  )) +
  labs(
    x = "Cancer Types",
    y = "Cohort Recall of OncoKB Targets (%)" # 明确标注为 Cohort Recall
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold"),
    axis.text.y = element_text(color = "black"),
    panel.grid.major.y = element_line(color = "grey92", linetype = "dashed"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.direction = "vertical",
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8)
  )
ggsave(file.path(fig_dir,"OncoKB_Recall_DotPlot.pdf"), p_recall_dot, width = 8.5, height = 5.5)

patient_oncokb_tier <- bind_rows(abg_hrw_vecs, .id = "cancer") %>%
  mutate(
    Clinical_Tier = case_when(
      (num_level1 >= 1 | num_level2 >= 1) ~ "Tier 1: Standard of Care (Level 1/2)",
      (num_level3 >= 1 | num_level4 >= 1) ~ "Tier 2: Investigational (Level 3/4)",
      TRUE ~ "Tier 3: No Actionable Target"
    ),
    Clinical_Tier = factor(Clinical_Tier, levels = c(
      "Tier 3: No Actionable Target",
      "Tier 2: Investigational (Level 3/4)",
      "Tier 1: Standard of Care (Level 1/2)"
    ))
  )
tier_summary_oncokb <- patient_oncokb_tier %>%
  group_by(cancer, Clinical_Tier) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(cancer) %>%
  mutate(Proportion = Count / sum(Count))

cancer_order <- sort(unique(tier_summary_oncokb$cancer))
tier_summary_oncokb$cancer <- factor(tier_summary_oncokb$cancer, levels = cancer_order)
p_oncokb_stacked <- ggplot(tier_summary_oncokb, aes(x = cancer, y = Proportion, fill = Clinical_Tier)) +
  geom_bar(stat = "identity", position = "fill", width = 0.75, color = "white", linewidth = 0.4) +
  scale_fill_manual(values = c(
    "Tier 3: No Actionable Target"         = "#D9D9D9",  # 高级灰
    "Tier 2: Investigational (Level 3/4)"  = "#80B1D3",  # 天空蓝
    "Tier 1: Standard of Care (Level 1/2)" = "#FB8072"   # 珊瑚粉
  )) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  labs(
    x = "Cancer Types",
    y = "Proportion of Patients (%)",
    fill = "Highest OncoKB Actionability Tier"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold"),
    axis.text.y = element_text(color = "black"),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.direction = "vertical" 
  )
ggsave(file.path(fig_dir,"OncoKB_Patient_Tier_StackedBarPlot.pdf"), p_oncokb_stacked, width = 8.5, height = 5.5)

plot_data_oncokb_bg <- abg_summary_df %>%
  pivot_longer(
    cols = c(`1`, `2`, `3`, `4`),
    names_to = "Level",
    values_to = "Count"
  ) %>%
  mutate(
    Level = paste0("Level ", Level),
    Level = factor(Level, levels = c("Level 1", "Level 2", "Level 3", "Level 4"))
  )

cancer_order <- sort(unique(plot_data_oncokb_bg$cancer))
plot_data_oncokb_bg$cancer <- factor(plot_data_oncokb_bg$cancer, levels = cancer_order)

p_oncokb_background <- ggplot(plot_data_oncokb_bg, aes(x = cancer, y = Count, fill = Level)) +
  # 堆叠柱状图，加上极其纤细的白边让分界更加清爽
  geom_bar(stat = "identity", width = 0.75, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c(
    "Level 1" = "#FB8072", # 珊瑚粉
    "Level 2" = "#FDB462", # 杏色 (作为粉色的柔和延伸)
    "Level 3" = "#80B1D3", # 天空蓝
    "Level 4" = "#BEBADA"  # 浅紫灰 (存在感较低的冷色)
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "Cancer Types",
    y = "Number of Genes in OncoKB",
    fill = "Evidence Level"
  ) +
  
  # 顶级简约学术主题
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold"),
    axis.text.y = element_text(color = "black"),
    # 加一条淡淡的网格线方便审稿人横向对比绝对数值
    panel.grid.major.y = element_line(color = "grey92", linetype = "dashed"),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.background = element_blank(),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8)
  )
ggsave(file.path(fig_dir,"OncoKB_Evidence_Level_BarPlot.pdf"), p_oncokb_background, width = 8.5, height = 5.5)

#survival analysis
p_brca <- plot_survival_curve_RBO("BRCA",CGC,fig_dir)
p_coad <- plot_survival_curve_RBO("COAD",CGC,fig_dir)
p_hnsc <- plot_survival_curve_RBO("HNSC",CGC,fig_dir)
p_kirc <- plot_survival_curve_RBO("KIRC",CGC,fig_dir)
p_kirp <- plot_survival_curve_RBO("KIRP",CGC,fig_dir)
p_lihc <- plot_survival_curve_RBO("LIHC",CGC,fig_dir)
p_luad <- plot_survival_curve_RBO("LUAD",CGC,fig_dir)
p_lusc <- plot_survival_curve_RBO("LUSC",CGC,fig_dir)
p_prad <- plot_survival_curve_RBO("PRAD",CGC,fig_dir)
p_stad <- plot_survival_curve_RBO("STAD",CGC,fig_dir)
p_thca <- plot_survival_curve_RBO("THCA",CGC,fig_dir)
p_ucec <- plot_survival_curve_RBO("UCEC",CGC,fig_dir)

p_brca_tf <- plot_survival_curve_TF("BRCA",fig_dir)
p_coad_tf <- plot_survival_curve_TF("COAD",fig_dir)
p_hnsc_tf <- plot_survival_curve_TF("HNSC",fig_dir)
p_kirc_tf <- plot_survival_curve_TF("KIRC",fig_dir)
p_kirp_tf <- plot_survival_curve_TF("KIRP",fig_dir)
p_lihc_tf <- plot_survival_curve_TF("LIHC",fig_dir)
p_luad_tf <- plot_survival_curve_TF("LUAD",fig_dir)
p_lusc_tf <- plot_survival_curve_TF("LUSC",fig_dir)
p_prad_tf <- plot_survival_curve_TF("PRAD",fig_dir)
p_stad_tf <- plot_survival_curve_TF("STAD",fig_dir)
p_thca_tf <- plot_survival_curve_TF("THCA",fig_dir)
p_ucec_tf <- plot_survival_curve_TF("UCEC",fig_dir)

pvals <- numeric()
pvals["BRCA"] <- p_brca_tf$pvals["K=3"]
pvals["COAD"] <- p_coad_tf$pvals["K=3"]
pvals["HNSC"] <- p_hnsc_tf$pvals["K=3"]
pvals["KIRC"] <- p_kirc_tf$pvals["K=3"]
pvals["KIRP"] <- p_kirp_tf$pvals["K=3"]
pvals["LIHC"] <- p_lihc_tf$pvals["K=3"]
pvals["LUAD"] <- p_luad_tf$pvals["K=3"]
pvals["LUSC"] <- p_lusc_tf$pvals["K=4"]
pvals["PRAD"] <- p_prad_tf$pvals["K=2"]
pvals["STAD"] <- p_stad_tf$pvals["K=4"]
pvals["THCA"] <- p_thca_tf$pvals["K=2"]
pvals["UCEC"] <- p_ucec_tf$pvals["K=4"]

pvals_adj <- p.adjust(pvals, method = "BH")

library(ggpubr)
p_lists <- list(
    BRCA = p_brca_tf$plots[["K=3"]],
    COAD = p_coad_tf$plots[["K=3"]],
    HNSC = p_hnsc_tf$plots[["K=3"]],
    KIRC = p_kirc_tf$plots[["K=3"]],
    KIRP = p_kirp_tf$plots[["K=3"]],
    LIHC = p_lihc_tf$plots[["K=3"]],
    LUAD = p_luad_tf$plots[["K=3"]],
    LUSC = p_lusc_tf$plots[["K=4"]],
    PRAD = p_prad_tf$plots[["K=2"]],
    STAD = p_stad_tf$plots[["K=4"]],
    THCA = p_thca_tf$plots[["K=2"]],
    UCEC = p_ucec_tf$plots[["K=4"]]
)
idx_row_by_row <- as.vector(matrix(1:12, ncol = 3, byrow = TRUE))
p_lists_reordered <- p_lists[idx_row_by_row]

arranged_grid <- arrange_ggsurvplots(
  p_lists_reordered, 
  print = FALSE, 
  ncol = 3, 
  nrow = 4,
  risk.table.height = 0.25
)
ggsave(filename = file.path(fig_dir,"survival_analysis_hypernetwalk_subtypes.pdf"), plot = arranged_grid, width = 18, height = 22,device="pdf")

p_kirp <- plot_visualize_subtype("KIRP",3,fig_dir)
p_lihc <- plot_visualize_subtype("LIHC",3,fig_dir)
p_ucec <- plot_visualize_subtype("UCEC",4,fig_dir)

library(patchwork)
#拼成一张图，各占一行
combined_plot <- (p_kirp / p_lihc / p_ucec) + plot_layout(ncol = 1, heights = c(1, 1, 1))
ggsave(filename = file.path(fig_dir,"subtype_visualization.pdf"), plot = combined_plot, width = 12.5, height = 18, device = "pdf")
