library(ConsensusClusterPlus)
library(survival)
library(survminer)
get_filter_ref <- function(mut_data, genelist, N_coh) {
  selected_samples <- c()
  
  reference_sizes <- c()
  reference_genes_pers <- list()
  mut_genes_pers <- list()
  for (i in 1:ncol(mut_data)) {
    #只保留sample_i的前12个字符
    sample_i <- substr(colnames(mut_data)[i], 1, 12)
    mut_genes_i <- rownames(mut_data)[mut_data[, i] != 0]
    reference_i <- intersect(mut_genes_i, genelist)
    if (length(reference_i) >= 3) {
      selected_samples <- c(selected_samples, sample_i)
      reference_sizes <- c(reference_sizes, length(reference_i))
      reference_genes_pers[[sample_i]] <- reference_i
      mut_genes_pers[[sample_i]] <- mut_genes_i
    }
  }
  N_pers <- 2 * median(reference_sizes)
  # N_pers <- 10
  
  mut_genes_coh <- rownames(mut_data)[rowSums(mut_data) != 0]
  reference_genes_coh <- intersect(mut_genes_coh, genelist)
  
  return(
    list(
      N_pers = N_pers,
      selected_samples = selected_samples,
      reference_genes_pers = reference_genes_pers,
      mut_genes_pers = mut_genes_pers,
      N_coh = N_coh,
      reference_genes_coh = reference_genes_coh,
      mut_genes_coh = mut_genes_coh
    )
  )
}

compute_PRF <- function(prediction,
                        selected_samples,
                        reference_genes,
                        mut_genes,
                        N) {
  if (length(names(prediction)) > 1) {
    samples <- names(prediction)
    P <- c()
    R <- c()
    F1 <- c()
    null_samples <- c()
    valid_samples <- c()
    for (i in 1:length(samples)) {
      sample_i <- substr(samples[i], 1, 12)
      if (!sample_i %in% selected_samples) {
        next
      }
      valid_samples <- c(valid_samples, sample_i)
      prediction_i <- prediction[[i]]
      reference_i <- reference_genes[[sample_i]]
      num_pos_i <- length(reference_i)
      if (num_pos_i != 0) {
        P_i <- c()
        R_i <- c()
        F1_i <- c()
        mut_genes_i <- mut_genes[[sample_i]]
        for (n in 1:N) {
          if (n > length(mut_genes_i)) {
            p <- c <- f1 <- -1
          } else{
            num_tp <- length(intersect(prediction_i[1:min(n, length(prediction_i))], reference_i))
            p <- num_tp / n
            r <- num_tp / num_pos_i
            if (p + r != 0) {
              f1 <- 2 * p * r / (p + r)
            } else{
              f1 <- 0
            }
          }
          P_i <- c(P_i, p)
          R_i <- c(R_i, r)
          F1_i <- c(F1_i, f1)
        }
        P <- rbind(P, P_i)
        R <- rbind(R, R_i)
        F1 <- rbind(F1, F1_i)
      } else {
        null_samples <- c(null_samples, sample_i)
      }
    }
    rownames(P) <- valid_samples
    P_idx <- (P > -1) + 0
    R_idx <- (R > -1) + 0
    F1_idx <- (F1 > -1) + 0
    average_P <- colSums(P * P_idx) / colSums(P_idx)
    average_R <- colSums(R * R_idx) / colSums(R_idx)
    average_F1 <- colSums(F1 * F1_idx) / colSums(F1_idx)
    PRF <- rbind(average_P, average_R, average_F1)
    if (length(null_samples) > 0) {
      cat(length(null_samples),
          "sample(s) didn't predict known cancer genes.",
          "\n")
    }
  } else {
    P <- c()
    R <- c()
    F1 <- c()
    num_pos <- length(reference_genes)
    ###
    # num_pos <- length(intersect(reference_genes, rownames(prediction)))
    
    for (n in 1:N) {
      if (n > length(mut_genes)) {
        p <- c <- f1 <- -1
      } else{
        num_tp <- length(intersect(rownames(prediction)[1:n], reference_genes))
        p <- num_tp / n
        r <- num_tp / num_pos
        if (p + r != 0) {
          f1 <- 2 * p * r / (p + r)
        } else{
          f1 <- 0
        }
      }
      P <- c(P, p)
      R <- c(R, r)
      F1 <- c(F1, f1)
    }
    PRF <- rbind(P, R, F1)
  }
  rownames(PRF) <- c("Precision", "Recall", "F1_score")
  return(PRF)
}

get_ROC <- function(res_coh,genelist){
    roc_data <- lapply(names(res_coh),function(method){
        res <- res_coh[[method]]
        labels <- ifelse(rownames(res) %in% genelist, 1, 0)
        roc_obj <- roc(labels,res[,1],percent=TRUE,quiet=TRUE)
        auc_full <- as.numeric(roc_obj$auc)
        auc_partial <- as.numeric(auc(roc_obj, partial.auc = c(100, 90), partial.auc.focus = "specificity"))
        data.frame(
            Sensitivity = c(roc_obj$sensitivities),
            Specificity = c(roc_obj$specificities),
            AUC = c(auc_full),
            pAUC = c(auc_partial),
            method = factor(method,levels = names(res_coh))
        )
    }) %>% bind_rows()
    roc_data$method <- factor(roc_data$method, levels = names(res_coh))
    return(roc_data)
}

library(dplyr)

get_PR <- function(res_coh, genelist){
    pr_data <- lapply(names(res_coh), function(method){
        res <- res_coh[[method]]
        labels <- ifelse(rownames(res) %in% genelist, 1, 0)
        scores <- res[,1]
        
        ord <- order(scores, decreasing = TRUE)
        labels <- labels[ord]
        scores <- scores[ord]
        
        total_positives <- sum(labels)
        
        tp <- cumsum(labels)
        fp <- seq_along(labels) - tp
        
        recall <- tp / total_positives
        precision <- tp / (tp + fp)
        precision[is.na(precision)] <- 0

        pr_points <- data.frame(
            Precision = precision,
            Recall = recall
        )

        pr_points <- rbind(data.frame(Precision = 1, Recall = 0), pr_points)
        aupr <- 0
        for (i in 2:nrow(pr_points)) {
            aupr <- aupr + (pr_points$Recall[i] - pr_points$Recall[i - 1]) * (pr_points$Precision[i] + pr_points$Precision[i - 1]) / 2
        }
        max_recall_threshold <- 0.1
        idx_over <- which(pr_points$Recall > max_recall_threshold)[1]
        
        pr_partial <- pr_points[pr_points$Recall <= max_recall_threshold, ]
        if (!is.na(idx_over)) {
            idx_under <- idx_over - 1
            r1 <- pr_points$Recall[idx_under]
            r2 <- pr_points$Recall[idx_over]
            p1 <- pr_points$Precision[idx_under]
            p2 <- pr_points$Precision[idx_over]
            if (r1 < max_recall_threshold) {
                p_cut <- p1 + (p2 - p1) * (max_recall_threshold - r1) / (r2 - r1)
                pr_partial <- rbind(pr_partial, data.frame(Precision = p_cut, Recall = max_recall_threshold))
            }
        }
        raw_paupr <- 0
        for (i in 2:nrow(pr_partial)) {
            raw_paupr <- raw_paupr + (pr_partial$Recall[i] - pr_partial$Recall[i - 1]) * (pr_partial$Precision[i] + pr_partial$Precision[i - 1]) / 2
        }
        
        pr_points$method <- method
        pr_points$AUPR <- aupr * 100
        pr_points$pAUPR <- raw_paupr * 100 # 存储归一化后的百分比
        
        return(pr_points)
    }) %>% bind_rows()
    
    # 强制因子排序，保证画图时图例顺序与输入一致
    pr_data$method <- factor(pr_data$method, levels = names(res_coh))
    
    return(pr_data)
}

plot_roc <- function(roc_df_all){
    methods <- levels(factor(roc_df_all$method))
    cancers <- unique(roc_df_all$cancer_type)
    
    # 使用你精心挑选的高级调色盘
    color_palette <- c(
        "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", 
        "#91D1C2", "#DC0000", "#7E6148", "#B09C85", "#5081C3", "#E07B91"
    )[1:length(methods)]
    names(color_palette) <- methods
    
    # 极简虚线逻辑：HyperNetWalk 实线，其余全用简单的长虚线
    linetype_vector <- setNames(
        c("solid", rep("longdash", length(methods) - 1)),
        methods
    )
    
    # 统一主题
    pub_theme <- theme_classic(base_size = 12) +
        theme(
            text = element_text(color = "black", family = "sans"),
            axis.text = element_text(color = "black", size = 11),
            axis.title = element_text(color = "black", size = 13, face = "bold"),
            axis.line = element_line(linewidth = 0.8, color = "black"),
            axis.ticks = element_line(linewidth = 0.8, color = "black"),
            plot.title = element_text(face = "bold", size = 15, hjust = 0.5, margin = margin(b = 10)),
            plot.margin = margin(t = 10, r = 20, b = 10, l = 10),
            legend.position = "none" # 彻底干掉图例
        )
        
    ROC_plots <- list()
    
    for (current_cancer in cancers){
        # 1. 安全过滤数据
        roc_df <- roc_df_all %>% filter(cancer_type == current_cancer)
        
        # 2. 🌟 Z-index 置顶魔法：把 HyperNetWalk 的数据移到最后，确保画在最上面
        roc_df <- roc_df %>%
            arrange(method == "HyperNetWalk")
        
        # 3. 🌟 智能文本标签排版：按 AUC 降序排列，给个递增的索引用于计算 Y 坐标
        auc_summary <- roc_df %>%
            distinct(method, AUC) %>%
            arrange(desc(AUC)) %>%
            mutate(
                rank = row_number(),
                # 文字的 Y 坐标：从 45% 向下递减，间隔 6%
                # 文字的 X 坐标：固定在 50% 的位置，向右对齐 (hjust=0)
                text_y = 45 - (rank - 1) * 6
            )
            
        p <- ggplot(roc_df, aes(x = 100 - Specificity, y = Sensitivity, color = method)) +
            
            # 对角线放在最底层
            geom_abline(
                slope = 1, intercept = 0, 
                linetype = "dashed", color = "#A9A9A9", linewidth = 0.8
            ) +
            
            # 画 ROC 曲线
            geom_line(aes(linetype = method), linewidth = 1.0, alpha = 0.85) +
            
            # 🌟 智能写入 AUC 文字
            geom_text(
                data = auc_summary,
                aes(
                    x = 55, # 放在图的右下角空白区
                    y = text_y, 
                    label = sprintf("%s: %.1f%%", method, AUC),
                    color = method
                ),
                hjust = 0, # 左对齐
                size = 4.0,
                fontface = "bold", # 文字加粗增加清晰度
                show.legend = FALSE
            ) +
            
            # 坐标轴设置
            scale_x_continuous(name = "False Positive Rate (%)", limits = c(0, 100), expand = c(0.01, 0.01)) +
            scale_y_continuous(name = "True Positive Rate (%)", limits = c(0, 100), expand = c(0.01, 0.01)) +
            
            scale_color_manual(values = color_palette) +
            scale_linetype_manual(values = linetype_vector) +
            
            ggtitle(current_cancer) +
            pub_theme
            
        ROC_plots[[current_cancer]] <- p
    }
    
    return(ROC_plots)
}

plot_pr <- function(pr_df_all){
    methods <- levels(factor(pr_df_all$method))
    cancers <- unique(pr_df_all$cancer_type)
    
    # 高级调色盘
    color_palette <- c(
        "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", 
        "#91D1C2", "#DC0000", "#7E6148", "#B09C85", "#5081C3", "#E07B91"
    )[1:length(methods)]
    names(color_palette) <- methods
    
    # 极简虚线逻辑：HyperNetWalk 实线，其余全用长虚线
    linetype_vector <- setNames(
        c("solid", rep("longdash", length(methods) - 1)),
        methods
    )
    
    # 统一主题
    pub_theme <- theme_classic(base_size = 12) +
        theme(
            text = element_text(color = "black", family = "sans"),
            axis.text = element_text(color = "black", size = 11),
            axis.title = element_text(color = "black", size = 13, face = "bold"),
            axis.line = element_line(linewidth = 0.8, color = "black"),
            axis.ticks = element_line(linewidth = 0.8, color = "black"),
            plot.title = element_text(face = "bold", size = 15, hjust = 0.5, margin = margin(b = 10)),
            plot.margin = margin(t = 10, r = 20, b = 10, l = 10),
            legend.position = "none" # 干掉图例
        )
        
    PR_plots <- list()
    
    for (current_cancer in cancers){
        # 1. 安全过滤
        pr_df <- pr_df_all %>% filter(cancer_type == current_cancer)
        
        # 2. 🌟 Z-index 置顶魔法：把 HyperNetWalk 移到最后画
        pr_df <- pr_df %>%
            arrange(method == "HyperNetWalk")
            
        # 3. 🌟 智能文本排版：按 AUPR 降序排列
        aupr_summary <- pr_df %>%
            distinct(method, AUPR) %>%
            arrange(desc(AUPR)) %>%
            mutate(
                rank = row_number(),
                # 文字的 Y 坐标：PR曲线右上角是空的，所以我们从 95% 开始往下排
                text_y = 95 - (rank - 1) * 6
            )
            
        p <- ggplot(pr_df, aes(x = Recall * 100, y = Precision * 100, color = method)) +
            
            # 画 PR 曲线
            geom_line(aes(linetype = method), linewidth = 1.0, alpha = 0.85) +
            
            # 🌟 智能写入 AUPR 文字
            geom_text(
                data = aupr_summary,
                aes(
                    x = 60, # 放在图的右上偏中区域
                    y = text_y, 
                    label = sprintf("%s: %.1f%%", method, AUPR),
                    color = method
                ),
                hjust = 0, # 左对齐
                size = 4.0,
                fontface = "bold",
                show.legend = FALSE
            ) +
            
            # 坐标轴设置：加上 0.01 的 expand 防止粗线条被切断
            scale_x_continuous(name = "Recall (%)", limits = c(0, 100), expand = c(0.01, 0.01)) +
            scale_y_continuous(name = "Precision (%)", limits = c(0, 100), expand = c(0.01, 0.01)) +
            
            scale_color_manual(values = color_palette) +
            scale_linetype_manual(values = linetype_vector) +
            
            ggtitle(current_cancer) +
            pub_theme
            
        PR_plots[[current_cancer]] <- p
    }
    
    return(PR_plots)
}

plot_PRF_curve <- function(df_list){
    cancers <- names(df_list)
    methods <- unique(unlist(lapply(df_list, function(df) df$method)))
    methods <- setdiff(methods, "Random") 
    # 采用高级基因组学调色盘，你的方法固定使用 "#E64B35" (顶刊红)
    color_palette <- c(
        "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", 
        "#91D1C2", "#DC0000", "#7E6148", "#B09C85", "#5081C3", "#E07B91",
        "#76B041", "#D3A03C", "#6C6EAF", "#A8584E"
    )[1:length(methods)]
    names(color_palette) <- methods
    
    # 极简虚线逻辑：HyperNetWalk 独享实线，其余全部使用干净的长虚线
    linetype_vector <- setNames(
        c("solid", rep("longdash", length(methods) - 1)),
        methods
    )
    
    # 2. 定义出版级主题
    pub_theme <- theme_classic(base_size = 12) +
        theme(
            text = element_text(color = "black", family = "sans"), 
            axis.text = element_text(color = "black", size = 11),
            axis.title = element_text(color = "black", size = 13, face = "bold"),
            axis.line = element_line(linewidth = 0.8, color = "black"),
            axis.ticks = element_line(linewidth = 0.8, color = "black"),
            plot.title = element_text(hjust = 0.5, face = "bold", size = 15, margin = margin(b = 12)),
            # 图例设置 (在单图中隐藏，后期可用 patchwork 提取公共图例)
            legend.title = element_blank(),
            legend.text = element_text(size = 11),
            legend.background = element_blank(),
            legend.key = element_blank(),
            legend.key.width = unit(1.8, "cm"), 
            legend.spacing.x = unit(0.2, "cm"),
            legend.position = "bottom"
        )
    
    # 3. 定义独立的动态 Y 轴生成器
    get_dynamic_y <- function(metric_values) {
        val_min <- min(metric_values, na.rm = TRUE)
        val_max <- max(metric_values, na.rm = TRUE)
        
        # 动态 padding：取数据极差的 15%，且至少保证有 0.05 的绝对缓冲空间
        padding <- max((val_max - val_min) * 0.15, 0.05) 
        
        scale_y_continuous(
            # 严格确保放大后的坐标系不突破 0 和 1 的物理极限
            limits = c(max(0, val_min - padding), min(1, val_max + padding)),
            expand = c(0, 0)
        )
    }
    
    # 初始化存储列表
    prec_curves <- list()
    rec_curves <- list()
    f1_curves <- list()
    
    # 4. 开始循环绘图
    for (i in seq_along(cancers)){
        cancer_name <- cancers[i]
        
        # 🌟 Z-index 置顶魔法：强制将 HyperNetWalk 排到最后画，确保红线永远覆盖在灰线上
        df <- df_list[[i]] %>%
            arrange(method == "HyperNetWalk") %>%
            filter(method != "Random") %>%
            group_by(method) %>%
            filter(any(precision > 0) || any(recall > 0) || any(f1 > 0)) %>%
            ungroup()
        df$method <- factor(df$method, levels = methods)

        K <- max(as.numeric(df$n_gene), na.rm = TRUE)
        if (K <= 30){
            x_breaks <- pretty(1:K, n = min(10, K))
            x_breaks <- x_breaks[x_breaks >= 1]
        } else {
            x_breaks <- seq(0,K,length.out=5)
        }
        
        # 提取公共绘图层，减少代码冗余
        base_layers <- list(
            geom_line(aes(linetype = method), linewidth = 1.0, alpha = 0.85),
            scale_x_continuous(breaks = x_breaks, expand = c(0.02, 0.02)),
            scale_color_manual(values = color_palette),
            scale_linetype_manual(values = linetype_vector),
            pub_theme,
            theme(legend.position = "bottom", legend.box = "horizontal"),
            ggtitle(cancer_name) 
        )
        
        # 绘制 Precision 曲线，并注入其专属动态 Y 轴
        prec_curve <- ggplot(df, aes(x = n_gene, y = precision, color = method, group = method)) +
            base_layers +
            get_dynamic_y(df$precision) + 
            labs(x = "Top K genes", y = "Precision")
            
        # 绘制 Recall 曲线，并注入其专属动态 Y 轴
        rec_curve <- ggplot(df, aes(x = n_gene, y = recall, color = method, group = method)) +
            base_layers +
            get_dynamic_y(df$recall) + 
            labs(x = "Top K genes", y = "Recall")

        # 绘制 F1 曲线，并注入其专属动态 Y 轴
        f1_curve <- ggplot(df, aes(x = n_gene, y = f1, color = method, group = method)) +
            base_layers +
            get_dynamic_y(df$f1) + 
            labs(x = "Top K genes", y = "F1 Score")
            
        prec_curves[[cancer_name]] <- prec_curve
        rec_curves[[cancer_name]] <- rec_curve
        f1_curves[[cancer_name]] <- f1_curve
    }
    
    return(list(prec_curves = prec_curves, rec_curves = rec_curves, f1_curves = f1_curves))
}

run_overlap_cross_cancer <- function(res, IntOGen_drivers, top_n = 100) {
    overlap_df <- data.frame()
    for (ref_cancer in cancers) {
        ref_genes <- IntOGen_drivers[IntOGen_drivers$TCGA_CANCER %in% ref_cancer, ]$SYMBOL
        ref_genes <- intersect(ref_genes, rownames(res))
        if (length(ref_genes) == 0) {
            next
        }
        top_genes <- rownames(res)[1:top_n]
        overlap_genes <- intersect(top_genes, ref_genes)
        overlap_count <- length(overlap_genes)
        overlap_prop <- overlap_count / length(ref_genes)
        overlap_df <- rbind(overlap_df, data.frame(
            Reference_Cancer = ref_cancer,
            Overlap_Count = overlap_count,
            Overlap_Proportion = overlap_prop
        ))
    }
    return(overlap_df)
}

run_fgsea <- function(res, all_pathways_list){
    gene_scores <- res[,1]
    names(gene_scores) <- rownames(res)
    set.seed(921)
    noise <- runif(length(gene_scores), min = 0, max = 1e-10)
    gene_scores <- gene_scores + noise
    gene_scores <- sort(gene_scores, decreasing = TRUE)
    fgsea_res <- fgseaMultilevel(
        pathways = all_pathways_list,
        stats = gene_scores,
        minSize = 10,
        maxSize = 500,
        eps = 0,
        scoreType = "pos",
        gseaParam = 0
    )
    return(fgsea_res)
}

analysis_gsea_cross_cancer <- function(df) {
  df <- df %>%
    group_by(Predicted_Cancer) %>%
    mutate(
      p_rank     = rank(pval, ties.method = "min"),
      Is_Best    = if_else(p_rank == 1, "Best", "Other"),
      neg_log_padj = pmin(-log10(padj), 12)  # 直接在这里 clamp，语义更清晰
    ) %>%
    ungroup()

  # 对角线标记（预测癌种 == 参考癌种）
  df <- df %>%
    mutate(Is_Diagonal = (Predicted_Cancer == Reference_Cancer))

  p <- ggplot(df, aes(x = Reference_Cancer, y = Predicted_Cancer)) +
    geom_point(aes(size = NES, color = neg_log_padj)) +
    # 黑框：最显著的那个
    geom_point(
      data = subset(df, Is_Best == "Best"),
      aes(size = NES),
      shape = 21, fill = NA, color = "black", stroke = 1.5
    ) +
    scale_color_gradientn(
      colors = c("#4575B4", "#E0F3F8", "#FDAE61", "#D73027"),
      values  = scales::rescale(c(1, 3, 6, 12)),
      breaks  = c(1, 3, 6, 12),
      limits  = c(1, 12),
      oob     = scales::squish,
      name    = "-log10(padj)"
    ) +
    scale_size_continuous(range = c(2, 9), name = "NES") +
    scale_x_discrete(limits = sort(unique(df$Reference_Cancer))) +
    scale_y_discrete(limits = rev(sort(unique(df$Predicted_Cancer)))) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x    = element_text(angle = 45, hjust = 1, color = "black",
                                    size = 11, face = "bold"),
      axis.text.y    = element_text(color = "black", size = 11, face = "bold"),
      axis.title     = element_text(size = 12, face = "bold"),
      legend.position = "right",
      panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      panel.border   = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    labs(
      x = "Reference Cancer Type (IntOGen Specific Drivers)",
      y = "Predicted Cancer Type (Model Scoring Cohort)"
    )
  return(p)
}

run_ora_topN <- function(res, IntOGen_drivers, cancer, top_n = 100) {
  
  universe_genes <- rownames(res)
  N_total <- length(universe_genes)
  
  drivers_c <- unique(IntOGen_drivers[IntOGen_drivers$TCGA_CANCER == cancer, ]$SYMBOL)
  drivers_in_universe <- intersect(drivers_c, universe_genes)
  M_pathway <- length(drivers_in_universe)
  
  res_sorted <- res[order(res[, 1], decreasing = TRUE), , drop = FALSE]
  top_genes <- rownames(res_sorted)[1:top_n]
  
  overlap_genes <- intersect(top_genes, drivers_in_universe)
  x_hit <- length(overlap_genes)
  
  # ---------------------------------------------------------
  # 方式 A：使用 Fisher 精确检验 (最推荐，清晰不易错)
  # ---------------------------------------------------------
  # 构建 2x2 列联表
  #               在真实驱动名单内  不在名单内
  # 在 Top N 内      x_hit            top_n - x_hit
  # 不在 Top N 内    M_pathway-x_hit  剩余的普通基因
  
  contingency_matrix <- matrix(c(
    x_hit,                            # Top 100 且是 Driver
    top_n - x_hit,                    # Top 100 但不是 Driver (假阳性)
    M_pathway - x_hit,                # 不是 Top 100 但是 Driver (漏报)
    N_total - top_n - (M_pathway - x_hit) # 不是 Top 100 也不是 Driver (真阴性)
  ), nrow = 2, byrow = TRUE)
  
  # 🚨 极度重要：必须使用 alternative = "greater"！
  # 因为我们只关心“过表达/富集 (Over-representation)”，不关心“欠表达”。
  fisher_res <- fisher.test(contingency_matrix, alternative = "greater")
  p_val_fisher <- fisher_res$p.value
  
  # ---------------------------------------------------------
  # 方式 B：使用超几何检验 phyper (与 Fisher 结果完全等价，可做验证)
  # phyper(命中数-1, 宇宙中的总红球, 宇宙中的总白球, 抽出的球数, lower.tail=FALSE)
  # p_val_hyper <- phyper(x_hit - 1, M_pathway, N_total - M_pathway, top_n, lower.tail = FALSE)
  # ---------------------------------------------------------
  
  # 计算命中率 (Precision) 和 召回率 (Recall)
  precision <- x_hit / top_n
  recall <- x_hit / M_pathway
  
  # 返回一个整理好的结果框
  return(data.frame(
    Hit_Count = x_hit,
    Pathway_Size = M_pathway,
    Precision = precision,
    Recall = recall,
    P_value = p_val_fisher
  ))
}

# Jaccard index
calculate_jaccard <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  if (union == 0) {
    return(0)
  } else {
    return(intersection / union)
  }
}
plot_combined_jaccard <- function(mat_coh, mat_pers, method) {
  
  cancers <- rownames(mat_coh)
  df <- expand.grid(Cancer_Y = cancers, Cancer_X = cancers, stringsAsFactors = FALSE)
  off_diag_vals <- c(mat_coh[lower.tri(mat_coh)], mat_pers[upper.tri(mat_pers)])
  max_val <- max(off_diag_vals, na.rm = TRUE)
  if (max_val <= 0 || is.na(max_val)) max_val <- 1
  
  df_plot <- df %>%
    mutate(
      RowIdx = match(Cancer_Y, cancers),
      ColIdx = match(Cancer_X, cancers),
      Type = case_when(
        RowIdx > ColIdx ~ "Cohort",         # 左下红
        RowIdx < ColIdx ~ "Personalized",   # 右上蓝
        RowIdx == ColIdx ~ "Diagonal"       # 对角灰
      ),
      Jaccard = case_when(
        Type == "Cohort" ~ mat_coh[cbind(Cancer_Y, Cancer_X)],
        Type == "Personalized" ~ mat_pers[cbind(Cancer_Y, Cancer_X)],
        Type == "Diagonal" ~ NA_real_ 
      ),
      label_text = case_when(
        is.na(Jaccard) ~ "", 
        Jaccard == 0 ~ "0",
        TRUE ~ sprintf("%.2f", Jaccard)
      )
    )
  
  df_plot$Cancer_X <- factor(df_plot$Cancer_X, levels = cancers)
  df_plot$Cancer_Y <- factor(df_plot$Cancer_Y, levels = rev(cancers))
  
  p <- ggplot(df_plot, aes(x = Cancer_X, y = Cancer_Y)) +
    geom_tile(aes(fill = Type, alpha = Jaccard), color = "white", linewidth = 0.8) +
    
    geom_text(aes(label = label_text), color = "#222222", size = 3.5, fontface = "bold") +
    scale_fill_manual(
      values = c("Cohort" = "#D94841", "Personalized" = "#4575B4", "Diagonal" = "#E6E8EA"),
      name = "Data Source"
    ) +
    scale_alpha_continuous(
      range = c(0.1, 1),
      limits = c(0, max_val),
      breaks = scales::pretty_breaks(n = 4),
      name = "Jaccard\nIndex",
      guide = guide_legend(override.aes = list(fill = "#4D4D4D")) 
    ) +
    
    coord_fixed(ratio = 1) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", face = "bold"),
      axis.text.y = element_text(color = "black", face = "bold"),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 11),
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
    ) +
    labs(title = method)
  
  return(p)
}
plot_half_jaccard <- function(mat_coh, method) {
  
  cancers <- rownames(mat_coh)
  mat_coh[upper.tri(mat_coh, diag = FALSE)] <- NA
  
  df <- reshape2::melt(mat_coh, na.rm = TRUE)
  colnames(df) <- c("Cancer_Y", "Cancer_X", "Jaccard")

  df <- df %>%
    mutate(
      is_diag = (as.character(Cancer_X) == as.character(Cancer_Y)),
      plot_val = if_else(is_diag, NA_real_, Jaccard)
    )
  max_val <- max(df$plot_val, na.rm = TRUE)
  if (max_val == 0 || is.na(max_val)) max_val <- 1
  
  df_plot <- df %>%
    mutate(
      Cancer_X = factor(Cancer_X, levels = cancers),
      Cancer_Y = factor(Cancer_Y, levels = rev(cancers)),
      label_text = case_when(
        is_diag ~ "",
        plot_val == 0 ~ "0",
        TRUE ~ sprintf("%.2f", plot_val)
      )
    )
  
  p <- ggplot(df_plot, aes(x = Cancer_X, y = Cancer_Y, fill = plot_val)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = label_text), color = "#222222", size = 3.5, fontface = "bold") +
    
    scale_fill_gradient(
      low = "#F2F2F2", high = "#D94841",
      limits = c(0, max_val),
      na.value = "#E6E8EA",
      name = "Jaccard Index"
    ) +
    
    coord_fixed(ratio = 1) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", face = "bold"),
      axis.text.y = element_text(color = "black", face = "bold"),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = c(0.8, 0.8),
      legend.title = element_text(face = "bold", size = 11, angle = 90, hjust = 0.5),
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
    ) +
    labs(title = method)
  
  return(p)
}

compute_embedding_metrics <- function(
    emb_mat,
    group_vec,
    k = NULL,
    seed = 921,
    do_kmeans = TRUE
) {
    # emb_mat: samples x dims
    emb_mat <- as.matrix(emb_mat)
    group_vec <- as.factor(group_vec)

    if (nrow(emb_mat) != length(group_vec)) {
        stop("nrow(emb_mat) must equal length(group_vec)")
    }

    if (is.null(k)) {
        k <- length(unique(group_vec))
    }

    out <- list()

    # 1) silhouette using true labels
    # silhouette 要求 cluster label 为 integer
    group_int <- as.integer(group_vec)

    # 若只有一个组或有异常情况，则返回 NA
    if (length(unique(group_int)) > 1) {
        d_emb <- dist(emb_mat)
        sil <- cluster::silhouette(group_int, d_emb)
        out$silhouette_mean <- mean(sil[, "sil_width"], na.rm = TRUE)
        sil_df <- as.data.frame(sil)
        out$silhouette_by_group <- data.frame(
            Group = levels(group_vec),
            MeanSilhouette = tapply(sil_df$sil_width, group_vec, mean, na.rm = TRUE)
        )
    } else {
        out$silhouette_mean <- NA_real_
        out$silhouette_by_group <- NULL
    }

    # 2) kmeans + ARI/NMI
    if (do_kmeans) {
        set.seed(seed)
        km <- stats::kmeans(emb_mat, centers = k, nstart = 50)

        pred_cluster <- km$cluster
        true_group <- as.integer(group_vec)

        out$kmeans_cluster <- pred_cluster
        out$ARI <- mclust::adjustedRandIndex(pred_cluster, true_group)
        out$NMI <- aricode::NMI(pred_cluster, true_group)

        # 再给一个聚类后的 silhouette（看 cluster 本身是否紧）
        if (length(unique(pred_cluster)) > 1) {
            sil_km <- cluster::silhouette(pred_cluster, dist(emb_mat))
            out$kmeans_silhouette_mean <- mean(sil_km[, "sil_width"], na.rm = TRUE)
        } else {
            out$kmeans_silhouette_mean <- NA_real_
        }
    } else {
        out$kmeans_cluster <- NULL
        out$ARI <- NA_real_
        out$NMI <- NA_real_
        out$kmeans_silhouette_mean <- NA_real_
    }

    out
}

run_dimred <- function(
    mat_z,
    group_vec = NULL,
    seed = 921,
    var_explained = 0.9,
    umap_neighbors = 30,
    umap_min_dist = 0.3,
    compute_metrics = TRUE
) {
    set.seed(seed)
    pca_res <- prcomp(t(mat_z), center = FALSE, scale. = FALSE)

    cum_var <- cumsum(pca_res$sdev^2) / sum(pca_res$sdev^2)
    n_pcs <- which(cum_var >= var_explained)[1]
    if (is.na(n_pcs)) n_pcs <- min(50, ncol(pca_res$x))

    pca_mat <- pca_res$x[, 1:n_pcs, drop = FALSE]

    umap_config <- umap::umap.defaults
    umap_config$n_neighbors <- umap_neighbors
    umap_config$min_dist <- umap_min_dist

    umap_fit <- umap::umap(pca_mat, config = umap_config)
    umap_mat <- umap_fit$layout

    metrics_pca <- NULL
    metrics_umap <- NULL

    if (!is.null(group_vec) && compute_metrics) {
        metrics_pca  <- compute_embedding_metrics(pca_mat, group_vec, seed = seed)
        metrics_umap <- compute_embedding_metrics(umap_mat, group_vec, seed = seed)
    }

    list(
        pca = pca_res,
        n_pcs = n_pcs,
        pca_mat = pca_mat,
        umap = umap_fit,
        umap_mat = umap_mat,
        metrics_pca = metrics_pca,
        metrics_umap = metrics_umap
    )
}

plot_umap_by_group <- function(
    umap_fit,
    group_vec,
    color_map,
    title,
    metrics = NULL,
    show_metrics = TRUE
) {
    df <- data.frame(
        UMAP1 = umap_fit$layout[, 1],
        UMAP2 = umap_fit$layout[, 2],
        Group = group_vec
    )

    subtitle_text <- NULL
    if (!is.null(metrics) && show_metrics) {
        subtitle_text <- paste0(
            "Sil(true labels) = ", sprintf("%.3f", metrics$silhouette_mean),
            " | ARI(kmeans) = ", sprintf("%.3f", metrics$ARI),
            " | NMI(kmeans) = ", sprintf("%.3f", metrics$NMI)
        )
    }

    ggplot(df, aes(UMAP1, UMAP2)) +
        geom_point(aes(fill = Group), shape = 21, color = "grey95",
                   size = 1.8, alpha = 0.85, stroke = 0.15) +
        scale_fill_manual(values = color_map) +
        labs(
            title = title,
            subtitle = subtitle_text,
            x = "UMAP 1",
            y = "UMAP 2",
            fill = NULL
        ) +
        theme_bw(base_size = 12) +
        theme(
            plot.title = element_text(face = "bold", hjust = 0.5),
            plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey35"),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_rect(color = "grey60", linewidth = 0.5),
            legend.position = "right",
            legend.text = element_text(size = 9),
            legend.key = element_blank(),
            aspect.ratio = 1
        ) +
        guides(fill = guide_legend(override.aes = list(size = 4, alpha = 1)))
}

run_tsne_embedding <- function(
    mat_z,
    group_vec = NULL,
    seed = 921,
    var_explained = 0.90,
    perplexity = NULL,
    dims = 2,
    compute_metrics = TRUE
) {
    stopifnot(is.matrix(mat_z) || is.data.frame(mat_z))
    mat_z <- as.matrix(mat_z)

    set.seed(seed)
    pca_res <- prcomp(t(mat_z), center = FALSE, scale. = FALSE)

    cum_var <- cumsum(pca_res$sdev^2) / sum(pca_res$sdev^2)
    n_pcs <- which(cum_var >= var_explained)[1]
    if (is.na(n_pcs)) n_pcs <- min(50, ncol(pca_res$x))

    pca_input <- pca_res$x[, 1:n_pcs, drop = FALSE]

    n_samples <- nrow(pca_input)
    if (is.null(perplexity)) {
        perplexity <- min(30, floor((n_samples - 1) / 3))
        perplexity <- max(perplexity, 5)
    }

    set.seed(seed)
    tsne_out <- Rtsne(
        pca_input,
        dims = dims,
        perplexity = perplexity,
        check_duplicates = FALSE,
        verbose = TRUE
    )

    tsne_df <- data.frame(
        Dim1 = tsne_out$Y[, 1],
        Dim2 = tsne_out$Y[, 2]
    )

    if (!is.null(group_vec)) {
        tsne_df$Group <- group_vec
    }

    metrics_pca <- NULL
    metrics_tsne <- NULL

    if (!is.null(group_vec) && compute_metrics) {
        metrics_pca  <- compute_embedding_metrics(pca_input, group_vec, seed = seed)
        metrics_tsne <- compute_embedding_metrics(tsne_out$Y, group_vec, seed = seed)
    }

    list(
        tsne_df = tsne_df,
        tsne_out = tsne_out,
        pca_res = pca_res,
        pca_input = pca_input,
        n_pcs = n_pcs,
        perplexity = perplexity,
        metrics_pca = metrics_pca,
        metrics_tsne = metrics_tsne
    )
}

plot_tsne_embedding <- function(
    emb_df,
    color_map,
    title = "t-SNE",
    metrics = NULL,
    show_metrics = TRUE,
    point_size = 1.8,
    alpha = 0.85,
    legend_position = "right",
    label_centroids = FALSE
) {
    subtitle_text <- NULL
    if (!is.null(metrics) && show_metrics) {
        subtitle_text <- paste0(
            "Sil(true labels) = ", sprintf("%.3f", metrics$silhouette_mean),
            " | ARI(kmeans) = ", sprintf("%.3f", metrics$ARI),
            " | NMI(kmeans) = ", sprintf("%.3f", metrics$NMI)
        )
    }

    p <- ggplot(emb_df, aes(x = Dim1, y = Dim2)) +
        geom_point(
            aes(fill = Group),
            shape = 21,
            color = "grey95",
            stroke = 0.15,
            size = point_size,
            alpha = alpha
        ) +
        scale_fill_manual(values = color_map) +
        labs(
            title = title,
            subtitle = subtitle_text,
            x = "t-SNE 1",
            y = "t-SNE 2",
            fill = NULL
        ) +
        theme_bw(base_size = 12) +
        theme(
            plot.title = element_text(face = "bold", hjust = 0.5),
            plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey35"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "grey60", linewidth = 0.5),
            axis.text = element_text(color = "black"),
            legend.position = legend_position,
            legend.key = element_blank(),
            legend.text = element_text(size = 9),
            aspect.ratio = 1
        ) +
        guides(fill = guide_legend(override.aes = list(size = 4, alpha = 1)))

    if (label_centroids) {
        centroid_df <- emb_df %>%
            group_by(Group) %>%
            summarise(
                Dim1 = median(Dim1, na.rm = TRUE),
                Dim2 = median(Dim2, na.rm = TRUE),
                .groups = "drop"
            )

        p <- p + ggrepel::geom_text_repel(
            data = centroid_df,
            aes(label = Group),
            color = "grey20",
            size = 3.2,
            fontface = "bold",
            inherit.aes = FALSE,
            max.overlaps = Inf
        )
    }

    p
}

clean_hallmark_name <- function(x) {
  x <- gsub("^HALLMARK_", "", x)
  x <- gsub("_", " ", x)
  x
}

plot_hallmark_bar <- function(df, cancer, outdir,
                              width = 7,
                              base_height = 2.2,
                              per_pathway_height = 0.32,
                              max_height = 10) {
  
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  plot_df <- df %>%
    dplyr::filter(p.adjust < 0.05) %>%
    dplyr::mutate(
      Pathway = clean_hallmark_name(ID),
      logFDR = -log10(p.adjust)
    ) %>%
    dplyr::arrange(logFDR) %>%
    dplyr::mutate(
      Pathway = factor(Pathway, levels = Pathway)
    )
  
  if (nrow(plot_df) == 0) return(NULL)
  
  p <- ggplot(plot_df, aes(x = logFDR, y = Pathway, fill=logFDR)) +
    geom_col(color = "white", linewidth=0.4, width = 0.7) +
    scale_fill_gradient(low = "#F3D6D2", high = "#D94841", guide = "none") +
    labs(
      title = paste0(cancer, " Hallmark enrichment"),
      x = expression(-log[10]~"(FDR)"),
      y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.y = element_text(face = "bold", color = "black", size = 10),
      axis.text.x = element_text(color = "black"),
      axis.title.x = element_text(face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black")
    )
  
  fig_height <- min(max_height, base_height + nrow(plot_df) * per_pathway_height)
  
  ggsave(
    filename = file.path(outdir, paste0("Hallmark_bar_", cancer, ".pdf")),
    plot = p,
    width = width,
    height = fig_height,
    bg = "white"
  )
  
  return(p)
}

plot_venn_custom <- function(
  gene_lists,
  title = NULL,
  colors = NULL,
  highlight_set = "HyperNetWalk",
  title_size = 14,
  set_name_size = 4.8,
  count_size = 4.0,
  fill_alpha = 0.50,
  stroke_size = 0.8
) {
  library(ggvenn)
  library(ggplot2)

  # 默认配色：优先让 HyperNetWalk 更醒目，其余颜色偏柔和
  if (is.null(colors)) {
    colors <- c(
      "HyperNetWalk" = "#E64B35FF",  # 更醒目的红
      "DawnRank"     = "#4DBBD5FF",  # 蓝
      "DriverMP"     = "#00A087FF",  # 绿
      "DriverRWH"    = "#F39B7FFF",  # 橙
      "PDRWH"        = "#7E6148FF",  # 棕
      "PersonaDrive" = "#B09C85FF"   # 浅棕
    )
  }

  set_names <- names(gene_lists)
  n_sets <- length(gene_lists)

  # 如果输入列表没有名字，自动补名字
  if (is.null(set_names)) {
    set_names <- paste0("Set", seq_len(n_sets))
    names(gene_lists) <- set_names
  }

  # 生成 fill 颜色：优先按名字匹配，否则按顺序取
  fill_colors <- colors[set_names]
  if (any(is.na(fill_colors))) {
    unnamed_cols <- unname(colors)
    fill_colors[is.na(fill_colors)] <- unnamed_cols[seq_len(sum(is.na(fill_colors)))]
  }

  # 如果指定高亮集合存在，则把它放到第一位，保证视觉中心更突出
  if (!is.null(highlight_set) && highlight_set %in% set_names) {
    ord <- c(which(set_names == highlight_set), which(set_names != highlight_set))
    gene_lists <- gene_lists[ord]
    fill_colors <- fill_colors[ord]
    set_names <- set_names[ord]
    names(gene_lists) <- set_names
  }

  p <- ggvenn(
    gene_lists,
    fill_color = unname(fill_colors),
    fill_alpha = fill_alpha,
    stroke_color = "black",
    stroke_size = stroke_size,
    set_name_color = "black",
    set_name_size = set_name_size,
    text_color = "black",
    text_size = count_size,
    show_percentage = FALSE
  ) +
    coord_equal() +
    labs(title = title) +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = title_size,
        hjust = 0.5,
        margin = margin(b = 10)
      ),
      plot.margin = margin(8, 8, 8, 8)
    )

  return(p)
}

get_downstream_TFs <- function(gene,comb_adj_mat,total_TF){
  nodes_ppi <- rownames(comb_adj_mat)
  comb_adj_lists <- build_adj_lists(comb_adj_mat)
  com_TF <- intersect(nodes_ppi, total_TF)
  d_out <- multi_source_bfs_indices(gene, comb_adj_lists$out, comb_adj_lists$idx_map,2)
  d_in <- multi_source_bfs_indices(com_TF, comb_adj_lists$'in', comb_adj_lists$idx_map,2)
  sumd <- d_out + d_in
  candidate_idx <- which(!is.infinite(sumd) & (sumd <= 2))
  candidate_nodes <- names(d_out)[candidate_idx]
  if (length(intersect(gene,com_TF))>0){
    keep_node <- sapply(1:length(candidate_nodes), function(k) {
      current_sum <- sumd[candidate_nodes[k]]
      if (current_sum != 2) {
          return(TRUE) 
      }
      in_neighs_idx <- comb_adj_lists$'in'[[candidate_nodes[k]]]
      out_neighs_idx <- comb_adj_lists$out[[candidate_nodes[k]]]
      in_neighs <- names(comb_adj_lists$idx_map)[in_neighs_idx]
      out_neighs <- names(comb_adj_lists$idx_map)[out_neighs_idx]
      mut_in_neighs <- intersect(gene,in_neighs)
      TF_out_neighs <- intersect(GRN_list$total_TF,out_neighs)
      if (length(mut_in_neighs) == 0 || length(TF_out_neighs) == 0) { #本身是TF或突变基因
        return(TRUE)
      }
      if (length(TF_out_neighs)>1 || length(mut_in_neighs)>1) {
        return(TRUE)
      } 
      if (TF_out_neighs[1]==mut_in_neighs[1]){
        return(FALSE)
      } else {
        return(TRUE)
      }
    })
    med_nodes <- candidate_nodes[keep_node]
  } else {
    med_nodes <- candidate_nodes
  }
  return(med_nodes)
}

get_subgraph <- function(sub_adj_mat,adj_mat0,dg_TF,gene){
  g <- graph_from_adjacency_matrix(sub_adj_mat, mode = "directed", weighted = TRUE, diag = FALSE)
  edge_df <- igraph::as_data_frame(g, what = "edges")
  edge_df <- edge_df %>%
    rowwise() %>%
    mutate(
      is_directed = (adj_mat0[from, to] != 0)
    ) %>%
    ungroup()
  directed_edges <- edge_df %>%
    filter(is_directed == TRUE) %>%
    select(source = from, target = to, edge_weight = weight) %>%
    mutate(edge_type = "Directed")
  undirected_edges <- edge_df %>%
    filter(is_directed == FALSE) %>%
    rowwise() %>%
    mutate(pair_id = paste(sort(c(from, to)), collapse = "_")) %>%
    group_by(pair_id) %>%
    summarise(
      source = first(from),
      target = first(to),
      edge_weight = mean(weight),
      edge_type = "Undirected",
      .groups = "drop"
    ) %>%
    select(source, target, edge_weight, edge_type)
  edge_table_clean <- bind_rows(directed_edges, undirected_edges)
  if (!dir.exists(file.path(fig_dir, "subnet"))) {
    dir.create(file.path(fig_dir, "subnet"), recursive = TRUE)
  }
  write.csv(edge_table_clean, file.path(fig_dir, "subnet", paste0(gene, "_subnetwork_edges.csv")), row.names = FALSE, quote = FALSE)
  sub_nodes <- rownames(sub_adj_mat)
  node_table <- data.frame(node_id = sub_nodes) %>%
    mutate(
      node_class = case_when(
        node_id == gene ~ paste0("Source_", gene),
        node_id %in% names(dg_TF) ~ "Highlighted_TF",
        TRUE ~ "Other_mediate_nodes"
      ),
      outdegree = ifelse(node_id %in% names(dg_TF), 
                            dg_TF[node_id], 0.5)
    )
  write.csv(node_table, file.path(fig_dir, "subnet", paste0(gene, "_subnetwork_nodes.csv")), row.names = FALSE,quote = FALSE)
  return(list(edges = edge_table_clean, nodes = node_table))
}
get_rank_based_overlap <- function(l1,l2){
  K <- min(length(l1), length(l2))
  score <- 0
  for(k in 1:K){
    overlap <- length(intersect(l1[1:k], l2[1:k]))
    score <- score + overlap / k
  }
  score <- score / K
  return(score)
}

plot_survival_curve_RBO <- function(cancer,CGC,fig_dir){
  res_file <- file.path(methods[["HyperNetWalk"]],cancer,"results.Rdata")
  load(res_file)
  ref <- get_filter_ref(objs$mut_mat,CGC,500)
  N <- ref$N_pers
  sim_mat <- matrix(0,nrow=length(objs$mut_genes_ranks),ncol=length(objs$mut_genes_ranks))
  rownames(sim_mat) <- colnames(sim_mat) <- names(objs$mut_genes_ranks)
  for (i in 1:(nrow(sim_mat)-1)){
    sample_i <- rownames(sim_mat)[i]
    genes_i <- objs$mut_genes_ranks[[sample_i]][1:min(N,length(objs$mut_genes_ranks[[sample_i]]))]
    for (j in (i+1):ncol(sim_mat)){
        sample_j <- colnames(sim_mat)[j]
        genes_j <- objs$mut_genes_ranks[[sample_j]][1:min(N,length(objs$mut_genes_ranks[[sample_j]]))]
        sim_ij <- get_rank_based_overlap(genes_i, genes_j)
        sim_mat[i,j] <- sim_ij
        sim_mat[j,i] <- sim_ij
    }
  }
  diag(sim_mat) <- 1
  dist_mat <- as.dist(1 - sim_mat)

  surv_file <- paste0("./data/rawdata/TCGA-",cancer,".survival.tsv")
  surv_raw <- read.table(surv_file,header=T)
  surv_clean <- surv_raw %>%
      mutate(
          OS.time = as.numeric(OS.time)/30.41, # 转换为月，检查你的 OS.time 是天还是月，如果已经是月就不需要除以 30.41
          OS = as.numeric(OS)
      ) %>%
      filter(OS.time > 0) %>%
      group_by(X_PATIENT) %>%
      summarise(
          final_time = max(OS.time, na.rm = TRUE),
          final_status = max(OS, na.rm = TRUE)
      ) %>%
      ungroup() %>%
      left_join(
          surv_raw %>%
          # 筛选肿瘤样本 (01~09)，你给的示例里是 01A
          filter(as.numeric(substr(sample, 14, 15)) < 10) %>% 
          select(X_PATIENT, sample_tumor = sample) %>%
          distinct(), # 去重
          by = "X_PATIENT"
      ) %>%
      select(
          sample = sample_tumor,
          OS = final_status, 
          OS.time = final_time,
          X_PATIENT
      ) %>%
      filter(!is.na(sample))
  surv_final <- surv_clean %>% 
      filter(X_PATIENT %in% colnames(sim_mat)) %>%
      select(OS,OS.time,X_PATIENT) %>%
      distinct()
  
  library(ConsensusClusterPlus)
  library(dplyr)
  library(tibble)
  if (!dir.exists(file.path(fig_dir, cancer,"RBO"))) {
    dir.create(file.path(fig_dir, cancer,"RBO"), recursive = TRUE)
  }
  cc_out <- ConsensusClusterPlus(
    d = dist_mat,
    maxK = 6,               # 尝试分成 2 到 6 类
    reps = 1000,            # 抽样 1000 次
    pItem = 0.8,            # 每次抽取 80% 的样本
    pFeature = 1,           # 抽取 100% 的特征
    title = file.path(fig_dir,cancer,"RBO"), # 输出目录
    clusterAlg = "hc",      # 层次聚类
    innerLinkage = "ward.D2", # 层次聚类的链接方法
    distance = "custom",   # 组学数据推荐 pearson
    seed = 921,
    plot = "pdf"            # 会自动生成一致性矩阵图、累计分布函数(CDF)图
  )
  plots <- list()
  for (k in 2:6){
    cluster_labels <- data.frame(
      X_PATIENT = names(cc_out[[k]]$consensusClass),
      Cluster = paste0("Subtype_", cc_out[[k]]$consensusClass)
    )
    surv_df <- surv_final %>%
      inner_join(cluster_labels, by = "X_PATIENT") %>%
      mutate(Cluster = as.factor(Cluster))
    fit <- survfit(Surv(OS.time, OS) ~ Cluster, data = surv_df)
    p <- ggsurvplot(
      fit,
      data = surv_df,
      size = 1,                 # 线条粗细
      palette = "npg",          # 使用 Nature Publishing Group 配色方案
      conf.int = TRUE,          # 显示置信区间 (可选)
      pval = TRUE,              # 显示 Log-rank test p 值
      risk.table = TRUE,        # 显示风险表 (极其重要)
      legend.title = "Subtypes",
      xlab = "Time (Months)",   # 检查你的 OS.time 是天还是月，如果是天建议除以 30.41
      title = paste("Survival Analysis for", cancer),
      ggtheme = theme_light()   # 清爽的主题
    )
    ggsave(
      filename = file.path(fig_dir, cancer,"RBO", paste0("Survival_RBO_k", k, ".png")),
      plot = p$plot,
      width = 7,
      height = 6,
      bg = "white"
    )
    plots[[paste0("K=", k)]] <- p
  }
  return(plots)
}

plot_survival_curve_TF <- function(cancer,fig_dir){
  res_file <- file.path("./results/HyperNetWalk",cancer,"results.Rdata")
  load(res_file)
  TF_mat <- objs$P_tf
  surv_file <- paste0("./data/rawdata/TCGA-",cancer,".survival.tsv")
  surv_raw <- read.table(surv_file,header=T)
  surv_clean <- surv_raw %>%
      mutate(
          OS.time = as.numeric(OS.time)/30.41, # 转换为月，检查你的 OS.time 是天还是月，如果已经是月就不需要除以 30.41
          OS = as.numeric(OS)
      ) %>%
      filter(OS.time > 0) %>%
      group_by(X_PATIENT) %>%
      summarise(
          final_time = max(OS.time, na.rm = TRUE),
          final_status = max(OS, na.rm = TRUE),
          .groups = "drop"
      ) %>%
      left_join(
          surv_raw %>%
          # 筛选肿瘤样本 (01~09)，你给的示例里是 01A
          filter(as.numeric(substr(sample, 14, 15)) < 10) %>% 
          select(X_PATIENT, sample_tumor = sample) %>%
          distinct(), # 去重
          by = "X_PATIENT"
      ) %>%
      transmute(
          sample = sample_tumor,
          OS = final_status, 
          OS.time = final_time,
          X_PATIENT = X_PATIENT
      ) %>%
      filter(!is.na(sample))
  surv_final <- surv_clean %>% 
      filter(X_PATIENT %in% colnames(TF_mat)) %>%
      select(OS,OS.time,X_PATIENT) %>%
      distinct()
  TF_mat <- TF_mat[,colSums(TF_mat)>0]
  TF_mat <- apply(TF_mat, 2, function(x) x/sum(x))
  common_samples <- intersect(surv_final$X_PATIENT, colnames(TF_mat))
  TF_mat <- TF_mat[, common_samples, drop = FALSE]
  surv_final <- surv_final %>% filter(X_PATIENT %in% common_samples)
  row_vars <- apply(TF_mat, 1, var)
  thre <- quantile(row_vars, 0.05)
  top_tfs <- names(row_vars)[row_vars > thre]
  dat <- TF_mat[top_tfs, , drop = FALSE]
  pca <- prcomp(t(dat), center = TRUE, scale. = TRUE)
  cum_var <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
  n_pcs <- which(cum_var >= 0.9)[1]
  pca_mat <- t(pca$x[, 1:n_pcs, drop = FALSE])
  if (!dir.exists(file.path(fig_dir, cancer,"TF"))) {
    dir.create(file.path(fig_dir, cancer,"TF"), recursive = TRUE)
  }
  cc_out <- ConsensusClusterPlus(
    d = pca_mat,
    maxK = 6,               # 尝试分成 2 到 6 类
    reps = 1000,            # 抽样 1000 次
    pItem = 0.8,            # 每次抽取 80% 的样本
    pFeature = 1,           # 抽取 100% 的特征
    title = file.path(fig_dir,cancer,"TF"), # 输出目录
    clusterAlg = "hc",      # 层次聚类
    distance = "pearson",   # 组学数据推荐 pearson
    seed = 921,
    plot = "pdf"            # 会自动生成一致性矩阵图、累计分布函数(CDF)图
  )
  plots <- list()
  pvals <- list()
  for (k in 2:6){
    cluster_labels <- data.frame(
      X_PATIENT = names(cc_out[[k]]$consensusClass),
      Cluster = paste0("S", cc_out[[k]]$consensusClass)
    )
    surv_df <- surv_final %>%
      inner_join(cluster_labels, by = "X_PATIENT") %>%
      mutate(
        Cluster = factor(Cluster,levels = paste0("S", sort(unique(cc_out[[k]]$consensusClass))))
      )
    fit <- survfit(Surv(OS.time, OS) ~ Cluster, data = surv_df)
    survdiff <- survdiff(Surv(OS.time, OS) ~ Cluster, data = surv_df)
    p_val <- 1 - pchisq(survdiff$chisq, length(survdiff$n) - 1)
    cluster_sizes <- table(surv_df$Cluster)
    event_counts <- tapply(surv_df$OS, surv_df$Cluster, sum)
    legend_labs <- paste0(names(cluster_sizes),
      " (n=", as.integer(cluster_sizes), ")")
    
    # axis_max <- choose_axis_max(surv_df)
    axis_mat <- max(surv_df$OS.time, na.rm = TRUE)
    axis_max <- floor(axis_mat / 50) * 50
    p_label <- paste0("Log-rank p = ", signif(p_val, 3))
    p <- ggsurvplot(
      fit, 
      data = surv_df,
      linewidth = 1.1,                 # 线条粗细
      palette = c("#FB8072", "#80B1D3", "#B3DE69", "#FDB462", "#BC80BD", "#8DD3C7")[1:length(unique(surv_df$Cluster))], 
      conf.int = FALSE,          # 显示置信区间 (可选)
      pval = TRUE,              # 显示 Log-rank test p 值
      risk.table = TRUE,        # 显示风险表 (极其重要)
      risk.table.height = 0.24, # 调整风险表高度
      risk.table.title = "Number at risk",
      risk.table.y.text = TRUE,
      risk.table.y.text.col = TRUE,
      censor.shape = "|",
      censor.size = 2.6,
      legend.title = "Cluster",
      legend.labs = legend_labs,
      xlab = "Time (Months)",   # 检查你的 OS.time 是天还是月，如果是天建议除以 30.41
      ylab = "Survival Probability",
      title = cancer,
      xlim = c(0, axis_max),
      break.time.by = 50,
      ggtheme = theme_classic(base_size = 13),
      tables.theme = theme_cleantable(base_size = 10)
    )
    p$plot <- p$plot +
      theme(
        legend.position = c(0.97, 0.97),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = scales::alpha("white",0.85), color = "grey80", linewidth = 0.3),
        legend.key = element_blank(),
        legend.title = element_text(size = 9.5),
        legend.text = element_text(size = 8.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(color = "black", size = 10)
      )
    # table中y轴text去掉括号里n=...
    p$table <- p$table +
      theme(
        legend.position = "none",
        axis.title.x = element_text(size = 11),
        axis.text.x = element_text(size = 10)
      )
    
    p_full <- arrange_ggsurvplots(
      list(p),
      print = FALSE,
      ncol = 1,
      nrow = 1,
      risk.table.height = 0.25
    )

    pdf(
      file.path(fig_dir, cancer,"TF", paste0("Survival_TF_k", k, ".pdf")),
      width = 6.4,
      height = 6.2,
      onefile = FALSE
    )
    print(p_full)
    dev.off()
    plots[[paste0("K=", k)]] <- p
    pvals[[paste0("K=", k)]] <- p_val
  }
  return(list(plots = plots, pvals = pvals))
}

choose_axis_max <- function(surv_df, cluster_col = "Cluster") {
  max_t <- max(surv_df$OS.time, na.rm = TRUE)

  candidate_times <- seq(0, ceiling(max_t / 50) * 50, by = 50)

  risk_counts <- sapply(candidate_times, function(t) {
    min(tapply(
      surv_df$OS.time >= t,
      surv_df[[cluster_col]],
      sum
    ))
  })

  valid_times <- candidate_times[risk_counts >= 5]

  if (length(valid_times) == 0) {
    return(ceiling(max_t / 50) * 50)
  }

  max(valid_times)
}

plot_visualize_subtype <- function(cancer, k, fig_dir){
  my_colors_palette <- c(
    "#FB8072", "#80B1D3", "#B3DE69", "#FDB462", 
    "#BC80BD", "#8DD3C7", "#FFFFB3", "#BEBADA", 
    "#FCCDE5", "#D9D9D9", "#CCEBC5", "#FFED6F"
  )
  res_file <- file.path("./results/HyperNetWalk",cancer,"results.Rdata")
  load(res_file)
  TF_mat <- objs$P_tf
  TF_mat <- TF_mat[,colSums(TF_mat)>0]
  TF_mat <- apply(TF_mat, 2, function(x) x/sum(x))
  row_vars <- apply(TF_mat, 1, var)
  thre <- quantile(row_vars, 0.05)
  top_tfs <- names(row_vars)[row_vars > thre]
  dat <- TF_mat[top_tfs, , drop = FALSE]
  pca <- prcomp(t(dat), center = TRUE, scale. = TRUE)
  cum_var <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
  n_pcs <- which(cum_var >= 0.9)[1]
  pca_mat <- t(pca$x[, 1:n_pcs, drop = FALSE])
  cc_out <- ConsensusClusterPlus(
    d = pca_mat,
    maxK = 6,               # 尝试分成 2 到 6 类
    reps = 1000,            # 抽样 1000 次
    pItem = 0.8,            # 每次抽取 80% 的样本
    pFeature = 1,           # 抽取 100% 的特征
    title = file.path(fig_dir,cancer,"TF"), # 输出目录
    clusterAlg = "hc",      # 层次聚类
    distance = "pearson",   # 组学数据推荐 pearson
    plot = "pdf"            # 会自动生成一致性矩阵图、累计分布函数(CDF)图
  )
  cluster_labels <- data.frame(
    X_PATIENT = names(cc_out[[k]]$consensusClass),
    Cluster = paste0("S", cc_out[[k]]$consensusClass)
  )
  # 用UMAP可视化,并且按照Cluster上色
  umap_config <- umap::umap.defaults
  umap_config$n_neighbors <- 30
  umap_config$min_dist <- 0.3
  umap_fit <- umap::umap(t(pca_mat), config = umap_config)
  umap_df <- data.frame(
    X_PATIENT = colnames(pca_mat),
    UMAP1 = umap_fit$layout[, 1],
    UMAP2 = umap_fit$layout[, 2]
  ) %>%
    left_join(cluster_labels, by = "X_PATIENT")
  library(TCGAbiolinks)
  subtypes_all <- PanCancerAtlas_subtypes()
  subtypes_c <- subtypes_all %>%
    filter(cancer.type == cancer) %>%
    mutate(X_PATIENT = substr(pan.samplesID, 1, 12))
  if (length(intersect(umap_df$X_PATIENT, subtypes_c$X_PATIENT)) > 0) {
    subtype_c_clean <- subtypes_c %>%
      filter(!is.na(Subtype_mRNA) & Subtype_mRNA != "NA") %>%
      select(X_PATIENT, Subtype_mRNA) %>%
      distinct(X_PATIENT, .keep_all = TRUE)
    
    umap_df_sub <- umap_df %>%
      left_join(subtype_c_clean, by = "X_PATIENT") %>%
      filter(!is.na(Subtype_mRNA))
    p <- ggplot(umap_df_sub, aes(x = UMAP1, y = UMAP2)) +
      geom_point(aes(fill = Cluster), shape = 21, color = "grey95", size = 2.5, alpha = 0.85, stroke = 0.15) +
      scale_fill_manual(values = my_colors_palette) +
      labs(title = paste0("TF-Layer-Derived Subtypes (K=", k, ")"), x = "UMAP 1", y = "UMAP 2", fill = "Subtype") +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey60", linewidth = 0.5),
        axis.text = element_text(color = "black"),
        legend.position = "right",
        legend.key = element_blank(),
        legend.text = element_text(size = 9),
        aspect.ratio = 1
      ) +
      guides(fill = guide_legend(override.aes = list(size = 4, alpha = 1)))
    p_sub <- ggplot(umap_df_sub, aes(x = UMAP1, y = UMAP2)) +
      geom_point(aes(fill = Subtype_mRNA), shape = 21, color = "grey95", size = 2.5, alpha = 0.85, stroke = 0.15) +
      scale_fill_manual(values = rev(my_colors_palette)) +
      labs(
        title = "TCGA-Defined Subtypes",
        x = "UMAP 1",
        y = "UMAP 2",
        fill = "Subtype"
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey60", linewidth = 0.5),
        axis.text = element_text(color = "black"),
        legend.position = "right",
        legend.key = element_blank(),
        legend.text = element_text(size = 9),
        aspect.ratio = 1
      ) +
      guides(fill = guide_legend(override.aes = list(size = 4, alpha = 1)))
    # 将两个图并排显示,用cancer作为图片大标题
    library(patchwork)
    p_combined <- p + p_sub + plot_annotation(tag_levels = 'A', title = cancer) &
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        plot.tag = element_text(face = "bold", size = 12)
      )
    ggsave(
      filename = file.path(fig_dir, paste0(cancer, "_Subtype_Visualization_K", k, ".pdf")),
      plot = p_combined,
      width = 12.5,
      height = 5.5,
      bg = "white"
    )
    return(p_combined)
  } else {
    p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
      geom_point(aes(fill = Cluster), shape = 21, color = "grey95", size = 2.5, alpha = 0.85, stroke = 0.15) +
      scale_fill_manual(values = my_colors_palette) +
      labs(title = paste0(cancer, " HyperNetWalk Subtypes (K=", k, ")"), x = "UMAP 1", y = "UMAP 2", fill = "Subtype") +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey60", linewidth = 0.5),
        axis.text = element_text(color = "black"),
        legend.position = "right",
        legend.key = element_blank(),
        legend.text = element_text(size = 9),
        aspect.ratio = 1
      ) +
      guides(fill = guide_legend(override.aes = list(size = 4, alpha = 1)))
    ggsave(
      filename = file.path(fig_dir, paste0(cancer, "_Subtype_Visualization_K", k, ".pdf")),
      plot = p,
      width = 6.4,
      height = 6.2,
      bg = "white"
    )
    return(p)
  }
}
