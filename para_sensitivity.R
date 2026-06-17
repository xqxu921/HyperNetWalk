# =========================
# Parameter sensitivity analysis
# =========================

source("model.R")
source("plot_formal.R")

library(dplyr)
library(ggplot2)
library(patchwork)

# -------------------------
# Cancer types
# -------------------------
cancers <- c(
    "BRCA", "COAD", "HNSC", "KIRC", "KIRP", "LIHC",
    "LUAD", "LUSC", "PRAD", "STAD", "THCA", "UCEC"
)

# -------------------------
# Parameter settings
# Default values are added later when collecting results.
# -------------------------
max_degs_list <- c(100, 250, 750, 1000)

# alpha is constrained in (0, 1), so we avoid degenerate values 0 and 1.
alpha_list <- c(0.2, 0.8)

# sigma values are absolute values, following the default scale used in the method.
sigma_list <- c(0.05, 0.3, 0.5, 1, 2)

# r and gamma include small, moderate, and large restart probabilities.
# 0.15 is the PageRank-style default and is added later.
r_list <- c(0.05, 0.3, 0.5, 0.8)
gamma_list <- c(0.05, 0.3, 0.5, 0.8)

# -------------------------
# Run HyperNetWalk under alternative parameter settings
# -------------------------
for (c in seq_along(cancers)) {
    cancer_type <- cancers[c]

    mut_data_file <- file.path("./data/processed_data", cancer_type, "mut_data.tsv")
    exp_data_file <- file.path("./data/processed_data", cancer_type, "exp_tpm_data.tsv")

    message("Running sensitivity analysis for ", cancer_type)

    for (max_degs in max_degs_list) {
        output_dir <- file.path("./results/max_degs", as.character(max_degs))
        tryCatch({
            run_hypernetwalk(
                cancer_type,
                mut_data_file,
                exp_data_file,
                output_dir,
                max_degs = max_degs
            )
        }, error = function(e) {
            message("Error in ", cancer_type, ", max_degs = ", max_degs, ": ", e$message)
        })
    }

    for (alpha in alpha_list) {
        output_dir <- file.path("./results/alpha", as.character(alpha))
        tryCatch({
            run_hypernetwalk(
                cancer_type,
                mut_data_file,
                exp_data_file,
                output_dir,
                alpha = alpha
            )
        }, error = function(e) {
            message("Error in ", cancer_type, ", alpha = ", alpha, ": ", e$message)
        })
    }

    for (sigma in sigma_list) {
        output_dir <- file.path("./results/sigma", as.character(sigma))
        tryCatch({
            run_hypernetwalk(
                cancer_type,
                mut_data_file,
                exp_data_file,
                output_dir,
                sigma = sigma
            )
        }, error = function(e) {
            message("Error in ", cancer_type, ", sigma = ", sigma, ": ", e$message)
        })
    }

    for (r in r_list) {
        output_dir <- file.path("./results/r", as.character(r))
        tryCatch({
            run_hypernetwalk(
                cancer_type,
                mut_data_file,
                exp_data_file,
                output_dir,
                r = r
            )
        }, error = function(e) {
            message("Error in ", cancer_type, ", r = ", r, ": ", e$message)
        })
    }

    for (gamma in gamma_list) {
        output_dir <- file.path("./results/gamma", as.character(gamma))
        tryCatch({
            run_hypernetwalk(
                cancer_type,
                mut_data_file,
                exp_data_file,
                output_dir,
                gamma = gamma
            )
        }, error = function(e) {
            message("Error in ", cancer_type, ", gamma = ", gamma, ": ", e$message)
        })
    }
}

# =========================
# Collect and plot sensitivity results
# =========================

fig_dir <- "./figs/para_sensitivity/"
if (!dir.exists(fig_dir)) {
    dir.create(fig_dir, recursive = TRUE)
}

CGC <- read.delim(
    "./reference_dg/CGC_Tier1.tsv",
    header = TRUE,
    as.is = TRUE
)[, 1]

# -------------------------
# Default parameters used in the main analysis
# -------------------------
default_params <- list(
    max_degs = 500,
    alpha = 0.5,
    sigma = 0.1,
    r = 0.15,
    gamma = 0.15
)

# Add default values to tested parameter values.
param_values <- list(
    max_degs = sort(unique(c(max_degs_list, default_params$max_degs))),
    alpha    = sort(unique(c(alpha_list, default_params$alpha))),
    sigma    = sort(unique(c(sigma_list, default_params$sigma))),
    r        = sort(unique(c(r_list, default_params$r))),
    gamma    = sort(unique(c(gamma_list, default_params$gamma)))
)

# Display names for parameters in figures.
param_display <- c(
    max_degs = "num_degs",
    alpha = expression(alpha),
    sigma = expression(sigma),
    r = "r",
    gamma = expression(gamma)
)

# -------------------------
# Reference cache
# -------------------------
ref_cache <- list()

for (cancer in cancers) {
    mut_data <- get_mut_data(
        file.path("./data/processed_data", cancer, "mut_data.tsv"),
        1
    )

    ref_cache[[cancer]] <- get_filter_ref(
        mut_data,
        CGC,
        N_coh = 500
    )
}

# -------------------------
# Read personalized Precision@6
# -------------------------
read_precision_pers <- function(res_file, ref_obj, k = 6) {
    if (!file.exists(res_file)) {
        return(NA_real_)
    }

    res <- tryCatch(
        read.table(res_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE),
        error = function(e) NULL
    )

    if (is.null(res) || nrow(res) == 0 || ref_obj$N_pers < k) {
        return(NA_real_)
    }

    prediction <- list()

    for (j in seq_len(nrow(res))) {
        sample_j <- substr(res[j, 1], 1, 12)
        genes_ranking_j <- unlist(strsplit(res[j, 2], ","))
        prediction[[sample_j]] <- genes_ranking_j
    }

    PRF <- compute_PRF(
        prediction,
        ref_obj$selected_samples,
        ref_obj$reference_genes_pers,
        ref_obj$mut_genes_pers,
        ref_obj$N_pers
    )

    as.numeric(PRF["Precision", k])
}

# -------------------------
# Read cohort-level Precision@100
# -------------------------
read_precision_coh <- function(res_file, ref_obj, k = 100) {
    if (!file.exists(res_file)) {
        return(NA_real_)
    }

    res <- tryCatch(
        read.table(res_file, row.names = 1),
        error = function(e) NULL
    )

    if (is.null(res) || nrow(res) == 0 || ref_obj$N_coh < k) {
        return(NA_real_)
    }

    mut_genes_sup <- setdiff(ref_obj$mut_genes_coh, rownames(res))

    if (length(mut_genes_sup) > 0) {
        padding <- matrix(0, nrow = length(mut_genes_sup), ncol = ncol(res))
        rownames(padding) <- mut_genes_sup
        colnames(padding) <- colnames(res)
        res <- rbind(res, padding)
    }

    PRF <- compute_PRF(
        res,
        ref_obj$selected_samples,
        ref_obj$reference_genes_coh,
        ref_obj$mut_genes_coh,
        ref_obj$N_coh
    )

    as.numeric(PRF["Precision", k])
}

# -------------------------
# Formatting helper for x-axis labels
# -------------------------
format_param_value <- function(x) {
    if (all(abs(x - round(x)) < 1e-10)) {
        return(as.character(round(x)))
    }

    out <- format(x, trim = TRUE, scientific = FALSE)
    out <- sub("0+$", "", out)
    out <- sub("\\.$", "", out)
    out
}

# -------------------------
# Collect all sensitivity results
# -------------------------
collect_sensitivity_data <- function() {
    all_rows <- list()
    row_id <- 1

    for (param_name in names(param_values)) {
        default_value <- default_params[[param_name]]

        for (value in param_values[[param_name]]) {
            is_default <- isTRUE(all.equal(as.numeric(value), as.numeric(default_value)))

            result_root <- if (is_default) {
                "./results/HyperNetWalk"
            } else {
                file.path("./results", param_name, as.character(value))
            }

            for (cancer in cancers) {
                ref_obj <- ref_cache[[cancer]]

                pers_file <- file.path(result_root, cancer, "genes_ranking.txt")
                coh_file  <- file.path(result_root, cancer, "genes_ranking_cohort.txt")

                p6 <- read_precision_pers(pers_file, ref_obj, k = 6)
                p100 <- read_precision_coh(coh_file, ref_obj, k = 100)

                all_rows[[row_id]] <- data.frame(
                    parameter = param_name,
                    value_num = as.numeric(value),
                    value_label = format_param_value(as.numeric(value)),
                    is_default = is_default,
                    cancer_type = cancer,
                    metric = "Precision@6",
                    score = p6,
                    stringsAsFactors = FALSE
                )
                row_id <- row_id + 1

                all_rows[[row_id]] <- data.frame(
                    parameter = param_name,
                    value_num = as.numeric(value),
                    value_label = format_param_value(as.numeric(value)),
                    is_default = is_default,
                    cancer_type = cancer,
                    metric = "Precision@100",
                    score = p100,
                    stringsAsFactors = FALSE
                )
                row_id <- row_id + 1
            }
        }
    }

    bind_rows(all_rows)
}

sensitivity_df <- collect_sensitivity_data()

write.csv(
    sensitivity_df,
    file.path(fig_dir, "sensitivity_metrics_core.csv"),
    row.names = FALSE
)

# -------------------------
# Summarize median and IQR for reporting
# -------------------------
sensitivity_summary <- sensitivity_df %>%
    filter(!is.na(score)) %>%
    group_by(parameter, value_num, value_label, is_default, metric) %>%
    summarize(
        median = median(score, na.rm = TRUE),
        q1 = quantile(score, 0.25, na.rm = TRUE),
        q3 = quantile(score, 0.75, na.rm = TRUE),
        mean = mean(score, na.rm = TRUE),
        sd = sd(score, na.rm = TRUE),
        n_cancers = n(),
        .groups = "drop"
    )

write.csv(
    sensitivity_summary,
    file.path(fig_dir, "sensitivity_metrics_summary.csv"),
    row.names = FALSE
)

# -------------------------
# Compute shared y-axis limits for each metric
# -------------------------
metric_limits <- sensitivity_df %>%
    filter(!is.na(score)) %>%
    group_by(metric) %>%
    summarize(
        ymin = max(0, min(score, na.rm = TRUE) - 0.04),
        ymax = min(1, max(score, na.rm = TRUE) + 0.04),
        .groups = "drop"
    )

get_metric_limits <- function(metric_name) {
    row <- metric_limits %>% filter(metric == metric_name)

    if (nrow(row) == 0) {
        return(NULL)
    }

    c(row$ymin[1], row$ymax[1])
}

# -------------------------
# Main plotting function
# x-axis is treated as ordered discrete settings, so values are equally spaced.
# -------------------------
make_sensitivity_plot <- function(df_all, param_name, metric_name, default_value) {
    df_plot <- df_all %>%
        filter(
            .data$parameter == param_name,
            .data$metric == metric_name,
            !is.na(.data$score)
        )

    if (nrow(df_plot) == 0) {
        return(
            ggplot() +
                theme_void() +
                labs(title = paste0(param_name, " (", metric_name, "): No data"))
        )
    }

    value_table <- df_plot %>%
        distinct(value_num, value_label) %>%
        arrange(value_num)

    value_levels <- value_table$value_label

    df_plot <- df_plot %>%
        mutate(
            value_label = factor(.data$value_label, levels = value_levels)
        )

    summary_df <- df_plot %>%
        group_by(.data$value_num, .data$value_label) %>%
        summarize(
            med = median(.data$score, na.rm = TRUE),
            q1 = quantile(.data$score, 0.25, na.rm = TRUE),
            q3 = quantile(.data$score, 0.75, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        arrange(.data$value_num) %>%
        mutate(
            value_label = factor(.data$value_label, levels = value_levels)
        )

    default_label <- value_table$value_label[
        which.min(abs(value_table$value_num - as.numeric(default_value)))
    ]

    default_x <- match(default_label, value_levels)

    y_limits <- get_metric_limits(metric_name)

    p <- ggplot(
        df_plot,
        aes(
            x = .data$value_label,
            y = .data$score,
            group = .data$cancer_type
        )
    ) +
        geom_line(
            color = "#BDBDBD",
            linewidth = 0.45,
            alpha = 0.85
        ) +
        geom_errorbar(
            data = summary_df,
            aes(
                x = .data$value_label,
                ymin = .data$q1,
                ymax = .data$q3
            ),
            inherit.aes = FALSE,
            width = 0.18,
            linewidth = 0.75,
            color = "black"
        ) +
        geom_line(
            data = summary_df,
            aes(
                x = .data$value_label,
                y = .data$med,
                group = 1
            ),
            inherit.aes = FALSE,
            color = "#D94841",
            linewidth = 1.15
        ) +
        geom_point(
            data = summary_df,
            aes(
                x = .data$value_label,
                y = .data$med
            ),
            inherit.aes = FALSE,
            shape = 21,
            fill = "#D94841",
            color = "#D94841",
            size = 2.1,
            stroke = 0.25
        ) +
        geom_vline(
            xintercept = default_x,
            linetype = "dashed",
            linewidth = 0.7,
            color = "black"
        ) +
        theme_classic(base_size = 13) +
        theme(
            axis.text.x = element_text(color = "black", size = 10),
            axis.text.y = element_text(color = "black", size = 10),
            axis.title = element_text(face = "bold", size = 11),
            axis.line = element_line(linewidth = 0.7, color = "black"),
            axis.ticks = element_line(linewidth = 0.7, color = "black"),
            plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
            plot.tag = element_text(face = "bold", size = 14),
            panel.grid.major.y = element_line(color = "grey90", linewidth = 0.35),
            panel.grid.minor = element_blank()
        ) +
        labs(
            x = param_display[[param_name]],
            y = metric_name,
            title = param_name
        )

    if (!is.null(y_limits)) {
        p <- p + coord_cartesian(ylim = y_limits)
    }

    p
}

# -------------------------
# Generate 5 x 2 panel figure
# -------------------------
param_order <- c("max_degs", "alpha", "sigma", "r", "gamma")
metric_order <- c("Precision@6", "Precision@100")

plot_list <- list()
idx <- 1

for (param_name in param_order) {
    for (metric_name in metric_order) {
        plot_list[[idx]] <- make_sensitivity_plot(
            sensitivity_df,
            param_name,
            metric_name,
            default_params[[param_name]]
        )
        idx <- idx + 1
    }
}

p_main <- wrap_plots(
    plotlist = plot_list,
    ncol = 2,
    byrow = TRUE
) +
    plot_annotation(
        title = "Parameter sensitivity analysis",
        tag_levels = "A",
        theme = theme(
            plot.title = element_text(
                face = "bold",
                size = 16,
                hjust = 0.5
            )
        )
    )

# -------------------------
# Save figure
# -------------------------
ggsave(
    filename = file.path(fig_dir, "main_sensitivity_5x2_equal_spacing.pdf"),
    plot = p_main,
    width = 12,
    height = 16,
    device = "pdf",
    useDingbats = FALSE
)
