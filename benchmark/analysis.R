#!/usr/bin/Rscript

library(tidyverse)
library(RColorBrewer)
library(gridExtra)

data <- read.csv("benchmark_summary_final.csv")
dir.create("plots")
plot_dir <- "plots"

data <- data %>%
  mutate(
    data_type = str_extract(run, "^counts_[^_]+"),
    data_type = str_replace(data_type, "counts_", ""),
    mutation_level = str_extract(run, "_(\\d+)_k", 1),
    mutation_level = factor(as.numeric(mutation_level), levels = c(0, 1, 5, 10)),
    kmer = str_extract(run, "k(\\d+)", 1),
    kmer = factor(as.numeric(kmer)),
    gmer = str_extract(run, "g(\\d+)$", 1),
    gmer = factor(as.numeric(gmer))
  ) %>%
  filter(kmer %in% c(19, 23, 27, 31), gmer %in% c(5, 8, 11, 14)) %>%
  droplevels()

data <- data %>%
  mutate(
    assignment_rate_tx = assigned_tx / total_reads,
    assignment_rate_gene = assigned_gene / total_reads,
    unassigned_rate = (ambig_reads + short_reads) / total_reads
  )

# 1) Metric Heatmaps at 1% Mutation Rate (F1, Pearson) - ggplot version with separate color scales
plot_paper_metric_heatmaps <- function(data, data_type_filter, metric_prefix, metric_title, metric_label) {
  filtered <- data %>%
    filter(data_type == data_type_filter, mutation_level == 1)
  
  tx_col <- paste0(metric_prefix, "_tx")
  gene_col <- paste0(metric_prefix, "_gene")
  
  # Transcript heatmap
  heatmap_tx <- filtered %>%
    select(kmer, gmer, all_of(tx_col)) %>%
    rename(value = all_of(tx_col)) %>%
    mutate(
      text_color = ifelse(value > median(value), "white", "black")
    )
  
  p_tx <- ggplot(heatmap_tx, aes(x = gmer, y = kmer, fill = value)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.3f", value), color = text_color), size = 5, fontface = "bold") +
    scale_color_identity() +
    guides(color = "none") +
    scale_fill_gradientn(
      colors = colorRampPalette(c("#A8D8EA", "#9ECAE1", "#2E5090"))(10),
      breaks = scales::pretty_breaks(n = 5)
    ) +
    scale_x_discrete(name = "g-mer (minimizer length)") +
    scale_y_discrete(name = "k-mer size", limits = rev(levels(heatmap_tx$kmer))) +
    theme_bw(base_size = 11) +
    labs(
      title = paste(metric_title, " - Transcript"),
      fill = metric_label
    ) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11, face = "bold")
    )
  
  # Gene heatmap
  heatmap_gene <- filtered %>%
    select(kmer, gmer, all_of(gene_col)) %>%
    rename(value = all_of(gene_col)) %>%
    mutate(
      text_color = ifelse(value > median(value), "white", "black")
    )
  
  p_gene <- ggplot(heatmap_gene, aes(x = gmer, y = kmer, fill = value)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.3f", value), color = text_color), size = 5, fontface = "bold") +
    scale_color_identity() +
    guides(color = "none") +
    scale_fill_gradientn(
      colors = colorRampPalette(c("#A8D8EA", "#9ECAE1", "#2E5090"))(10),
      breaks = scales::pretty_breaks(n = 5)
    ) +
    scale_x_discrete(name = "g-mer (minimizer length)") +
    scale_y_discrete(name = "k-mer size", limits = rev(levels(heatmap_gene$kmer))) +
    theme_bw(base_size = 11) +
    labs(
      title = paste(metric_title, " - Gene"),
      fill = metric_label
    ) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11, face = "bold")
    )
  
  # Combine plots with separate legends
  combined <- gridExtra::grid.arrange(p_tx, p_gene, nrow = 1)
  return(combined)
}

# 1b) Differential Transcript Recall Heatmap
plot_diff_tx_recall_heatmap <- function(data) {
  diff_tx_num <- data$diff_tx_num[1]  # depends on ground truth so always the same
  filtered <- data %>%
    filter(data_type == "all", mutation_level == 1) %>%
    select(kmer, gmer, diff_tx_recall) %>%
    rename(value = diff_tx_recall) %>%
    mutate(
      text_color = ifelse(value > median(value, na.rm = TRUE), "white", "black")
    )
  
  p <- ggplot(filtered, aes(x = gmer, y = kmer, fill = value)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.3f", value), color = text_color), size = 5, fontface = "bold") +
    scale_color_identity() +
    guides(color = "none") +
    scale_fill_gradientn(
      colors = colorRampPalette(c("#A8D8EA", "#9ECAE1", "#2E5090"))(10),
      breaks = scales::pretty_breaks(n = 5),
      na.value = "lightgrey"
    ) +
    scale_x_discrete(name = "g-mer (minimizer length)") +
    scale_y_discrete(name = "k-mer size", limits = rev(levels(filtered$kmer))) +
    theme_bw(base_size = 11) +
    labs(
      title = paste0("Differential Transcript Recall (n = ", diff_tx_num, ")"),
      fill = "Recall"
    ) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11, face = "bold")
    )
  
  return(p)
}


# 2) Phase metrics: Heuristic vs Scoring reads across mutation rates (for "all" biotype)
plot_phase_metrics <- function(data) {
  filtered <- data %>%
    filter(data_type == "all") %>%
    na.omit()
  
  plot_data_tx <- filtered %>%
    group_by(mutation_level, kmer, gmer) %>%
    summarise(
      heuristic = mean(reads_heuristic_tx, na.rm = TRUE),
      scoring = mean(reads_scoring_tx, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = c("heuristic", "scoring"),
      names_to = "phase",
      values_to = "reads"
    ) %>%
    mutate(level_type = "Transcript")
  
  plot_data_gene <- filtered %>%
    group_by(mutation_level, kmer, gmer) %>%
    summarise(
      heuristic = mean(reads_heuristic_gene, na.rm = TRUE),
      scoring = mean(reads_scoring_gene, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = c("heuristic", "scoring"),
      names_to = "phase",
      values_to = "reads"
    ) %>%
    mutate(level_type = "Gene")
  
  plot_data <- bind_rows(plot_data_tx, plot_data_gene) %>%
    mutate(
      phase = factor(phase, levels = c("heuristic", "scoring")),
      level_type = factor(level_type, levels = c("Transcript", "Gene"))
    )
  
  p <- ggplot(plot_data, aes(x = mutation_level, y = reads, fill = level_type, alpha = phase)) +
    geom_col(position = "dodge", width = 0.7) +
    facet_grid(kmer ~ gmer, labeller = labeller(
      kmer = function(x) paste("k =", x),
      gmer = function(x) paste("g =", x)
    )) +
    scale_fill_manual(
      values = c("Transcript" = "#2166AC", "Gene" = "#1B9E77"),
      name = "Level"
    ) +
    scale_alpha_manual(
      values = c("heuristic" = 0.6, "scoring" = 1.0),
      labels = c("heuristic" = "Phase 1 (Heuristic)", "scoring" = "Phase 2 (Scoring)"),
      name = "Phase"
    ) +
    theme_bw(base_size = 10) +
    labs(
      title = "Reads by Detection Phase - all",
      x = "Mutation Rate (%)",
      y = "Number of Reads"
    ) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank(),
      strip.background = element_blank()
    )
  
  return(p)
}

# 4) Read distribution: Assigned, Ambiguous, Short across mutation rates
plot_read_distribution <- function(data, data_type_filter) {
  filtered <- data %>%
    filter(data_type == data_type_filter) %>%
    na.omit()
  
  plot_data <- filtered %>%
    group_by(mutation_level, kmer, gmer) %>%
    summarise(
      assigned_tx = mean(assigned_tx, na.rm = TRUE),
      assigned_gene = mean(assigned_gene, na.rm = TRUE),
      ambig_reads = mean(ambig_reads, na.rm = TRUE),
      short_reads = mean(short_reads, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = c("assigned_tx", "assigned_gene", "ambig_reads", "short_reads"),
      names_to = "category",
      values_to = "reads"
    ) %>%
    mutate(
      category = factor(category, levels = c("assigned_tx", "assigned_gene", "ambig_reads", "short_reads"))
    )
  
  p <- ggplot(plot_data, aes(x = mutation_level, y = reads, fill = category)) +
    geom_col(position = "stack", width = 0.7) +
    facet_grid(kmer ~ gmer, labeller = labeller(
      kmer = function(x) paste("k =", x),
      gmer = function(x) paste("g =", x)
    )) +
    scale_fill_manual(
      values = c(
        "assigned_tx" = "#2166AC",
        "assigned_gene" = "#1B9E77",
        "ambig_reads" = "#D8B365",
        "short_reads" = "#F4A460"
      ),
      labels = c(
        "assigned_tx" = "Assigned (Tx)",
        "assigned_gene" = "Assigned (Gene)",
        "ambig_reads" = "Ambiguous",
        "short_reads" = "Short"
      )
    ) +
    theme_bw(base_size = 10) +
    labs(
      title = paste("Read Distribution across Mutation Rates -", data_type_filter),
      x = "Mutation Rate (%)",
      y = "Number of Reads",
      fill = "Category"
    ) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank(),
      strip.background = element_blank()
    )
  
  return(p)
}

# 5) Focused Read Distribution (k=23, g=11)
plot_read_distribution_focused <- function(data, data_type_filter) {
  filtered <- data %>%
    filter(data_type == data_type_filter, kmer == 23, gmer == 11) %>%
    na.omit()
  
  plot_data <- filtered %>%
    select(mutation_level, assigned_tx, assigned_gene, ambig_reads, short_reads) %>%
    pivot_longer(
      cols = c("assigned_tx", "assigned_gene", "ambig_reads", "short_reads"),
      names_to = "category",
      values_to = "reads"
    ) %>%
    mutate(
      category = factor(category, levels = c("assigned_tx", "assigned_gene", "ambig_reads", "short_reads"))
    )
  
  p <- ggplot(plot_data, aes(x = mutation_level, y = reads, fill = category)) +
    geom_col(position = "stack", width = 0.6) +
    scale_fill_manual(
      values = c(
        "assigned_tx" = "#2166AC",
        "assigned_gene" = "#1B9E77",
        "ambig_reads" = "#D8B365",
        "short_reads" = "#F4A460"
      ),
      labels = c(
        "assigned_tx" = "Assigned (Tx)",
        "assigned_gene" = "Assigned (Gene)",
        "ambig_reads" = "Ambiguous",
        "short_reads" = "Short"
      )
    ) +
    theme_bw(base_size = 11) +
    labs(
      title = "(a) Read Distribution (k=23, g=11)",
      x = "Mutation Rate (%)",
      y = "Number of Reads",
      fill = "Category"
    ) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      axis.text.x = element_text(size = 13),
      legend.position = "bottom"
    )
  
  return(p)
}

plot_phase_metrics_focused <- function(data) {
  filtered <- data %>%
    filter(data_type == "all", kmer == 23, gmer == 11) %>%
    na.omit()
  
  plot_data_tx <- filtered %>%
    select(mutation_level, reads_heuristic_tx, reads_scoring_tx) %>%
    pivot_longer(
      cols = c("reads_heuristic_tx", "reads_scoring_tx"),
      names_to = "phase",
      values_to = "reads"
    ) %>%
    mutate(
      phase = ifelse(phase == "reads_heuristic_tx", "Phase 1", "Phase 2"),
      level_type = "Transcript"
    )
  
  plot_data_gene <- filtered %>%
    select(mutation_level, reads_heuristic_gene, reads_scoring_gene) %>%
    pivot_longer(
      cols = c("reads_heuristic_gene", "reads_scoring_gene"),
      names_to = "phase",
      values_to = "reads"
    ) %>%
    mutate(
      phase = ifelse(phase == "reads_heuristic_gene", "Phase 1", "Phase 2"),
      level_type = "Gene"
    )
  
  plot_data <- bind_rows(plot_data_tx, plot_data_gene) %>%
    mutate(
      phase = factor(phase, levels = c("Phase 1", "Phase 2")),
      level_type = factor(level_type, levels = c("Transcript", "Gene"))
    )
  
  p <- ggplot(plot_data, aes(x = mutation_level, y = reads, fill = level_type, alpha = phase)) +
    geom_col(position = "dodge", width = 0.7) +
    scale_fill_manual(
      values = c("Transcript" = "#2166AC", "Gene" = "#1B9E77"),
      name = "Level"
    ) +
    scale_alpha_manual(
      values = c("Phase 1" = 0.6, "Phase 2" = 1.0),
      name = "Phase"
    ) +
    theme_bw(base_size = 11) +
    labs(
      title = "(b) Phase Metrics (k=23, g=11)",
      x = "Mutation Rate (%)",
      y = "Number of Reads"
    ) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      axis.text.x = element_text(size = 13),
      legend.position = "bottom"
    )
  
  return(p)
}

# 8) Combined focused plots (read distribution and phase metrics stacked vertically)
plot_combined_focused <- function(data) {
  p_dist <- plot_read_distribution_focused(data, "all") +
    theme(
      legend.box = "horizontal",
      legend.box.margin = margin(1, 1, 1, 1),
      legend.margin = margin(0, 0, 0, 0),
      legend.spacing = unit(0.05, "cm"),
      legend.spacing.x = unit(0.1, "cm"),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.key.size = unit(0.4, "cm")
    )
  
  p_phase <- plot_phase_metrics_focused(data) +
    theme(
      legend.box = "horizontal",
      legend.box.margin = margin(1, 1, 1, 1),
      legend.margin = margin(0, 0, 0, 0),
      legend.spacing = unit(0.05, "cm"),
      legend.spacing.x = unit(0.1, "cm"),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.key.size = unit(0.4, "cm")
    )
  
  combined <- gridExtra::grid.arrange(p_dist, p_phase, nrow = 2, heights = c(1, 1.1), padding = unit(0.3, "cm"))
  
  return(combined)
}

# 7) Read Distribution at 1% Mutation (k-g combinations on x-axis)
plot_read_distribution_mutation1 <- function(data) {
  filtered <- data %>%
    filter(data_type == "all", mutation_level == 1) %>%
    na.omit()
  
  plot_data <- filtered %>%
    group_by(kmer, gmer) %>%
    summarise(
      assigned_tx = mean(assigned_tx, na.rm = TRUE),
      assigned_gene = mean(assigned_gene, na.rm = TRUE),
      ambig_reads = mean(ambig_reads, na.rm = TRUE),
      short_reads = mean(short_reads, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      kg_combination = paste0("k=", kmer, ", g=", gmer)
    ) %>%
    mutate(
      k_val = as.numeric(as.character(kmer)),
      g_val = as.numeric(as.character(gmer))
    ) %>%
    arrange(k_val, g_val) %>%
    mutate(
      kg_combination = factor(kg_combination, levels = unique(kg_combination))
    ) %>%
    select(-k_val, -g_val) %>%
    pivot_longer(
      cols = c("assigned_tx", "assigned_gene", "ambig_reads", "short_reads"),
      names_to = "category",
      values_to = "reads"
    ) %>%
    mutate(
      category = factor(category, levels = c("assigned_tx", "assigned_gene", "ambig_reads", "short_reads"))
    )
  
  p <- ggplot(plot_data, aes(x = kg_combination, y = reads, fill = category)) +
    geom_col(position = "stack", width = 0.7) +
    scale_fill_manual(
      values = c(
        "assigned_tx" = "#2166AC",
        "assigned_gene" = "#1B9E77",
        "ambig_reads" = "#D8B365",
        "short_reads" = "#F4A460"
      ),
      labels = c(
        "assigned_tx" = "Assigned (Tx)",
        "assigned_gene" = "Assigned (Gene)",
        "ambig_reads" = "Ambiguous",
        "short_reads" = "Short"
      )
    ) +
    theme_bw(base_size = 10) +
    labs(
      title = "Read Distribution at 1% Mutation Rate",
      x = NULL,
      y = "Number of Reads",
      fill = "Category"
    ) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

# 7b) Phase Metrics at 1% Mutation (k-g combinations on x-axis)
plot_phase_metrics_mutation1 <- function(data) {
  filtered <- data %>%
    filter(data_type == "all", mutation_level == 1) %>%
    na.omit()
  
  plot_data_tx <- filtered %>%
    group_by(kmer, gmer) %>%
    summarise(
      heuristic = mean(reads_heuristic_tx, na.rm = TRUE),
      scoring = mean(reads_scoring_tx, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      kg_combination = paste0("k=", kmer, ", g=", gmer),
      level_type = "Transcript"
    ) %>%
    pivot_longer(
      cols = c("heuristic", "scoring"),
      names_to = "phase",
      values_to = "reads"
    ) %>%
    mutate(
      phase = ifelse(phase == "heuristic", "Phase 1", "Phase 2")
    )
  
  plot_data_gene <- filtered %>%
    group_by(kmer, gmer) %>%
    summarise(
      heuristic = mean(reads_heuristic_gene, na.rm = TRUE),
      scoring = mean(reads_scoring_gene, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      kg_combination = paste0("k=", kmer, ", g=", gmer),
      level_type = "Gene"
    ) %>%
    pivot_longer(
      cols = c("heuristic", "scoring"),
      names_to = "phase",
      values_to = "reads"
    ) %>%
    mutate(
      phase = ifelse(phase == "heuristic", "Phase 1", "Phase 2")
    )
  
  plot_data <- bind_rows(plot_data_tx, plot_data_gene) %>%
    mutate(
      k_val = as.numeric(str_extract(kg_combination, "k=(\\d+)", 1)),
      g_val = as.numeric(str_extract(kg_combination, "g=(\\d+)", 1))
    ) %>%
    arrange(k_val, g_val) %>%
    mutate(
      kg_combination = factor(kg_combination, levels = unique(kg_combination)),
      phase = factor(phase, levels = c("Phase 1", "Phase 2")),
      level_type = factor(level_type, levels = c("Transcript", "Gene"))
    ) %>%
    select(-k_val, -g_val)
  
  p <- ggplot(plot_data, aes(x = kg_combination, y = reads, fill = level_type, alpha = phase)) +
    geom_col(position = "dodge", width = 0.7) +
    scale_fill_manual(
      values = c("Transcript" = "#2166AC", "Gene" = "#1B9E77"),
      name = "Level"
    ) +
    scale_alpha_manual(
      values = c("Phase 1" = 0.6, "Phase 2" = 1.0),
      name = "Phase"
    ) +
    theme_bw(base_size = 10) +
    labs(
      title = "Phase Metrics at 1% Mutation Rate",
      x = NULL,
      y = "Number of Reads"
    ) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

generate_all_plots <- function(data) {
  data_types <- unique(data$data_type)
  mutation_levels <- sort(unique(as.numeric(as.character(data$mutation_level))))
  kmers <- sort(unique(data$kmer))
  gmers <- sort(unique(data$gmer))
  
  # Create output directory for plots
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir)
  }
  
  for (dtype in data_types) {
    p <- plot_paper_metric_heatmaps(data, dtype, "pearson", "Pearson Correlation", "Pearson r")
    ggsave(sprintf("%s/pearson_heatmaps_%s_mutation1.png", plot_dir, dtype), p, width = 14, height = 6, dpi = 100)
    
    p <- plot_paper_metric_heatmaps(data, dtype, "f1", "F1 Score", "F1 Score")
    ggsave(sprintf("%s/f1_heatmaps_%s_mutation1.png", plot_dir, dtype), p, width = 14, height = 6, dpi = 100)
    
    p <- plot_read_distribution(data, dtype)
    ggsave(sprintf("%s/read_distribution_%s.png", plot_dir, dtype), p, width = 14, height = 8, dpi = 100)
    
    if (dtype == "all") {
      p <- plot_phase_metrics(data)
      ggsave(sprintf("%s/phase_metrics_%s.png", plot_dir, dtype), p, width = 14, height = 8, dpi = 100)
      
      p <- plot_diff_tx_recall_heatmap(data)
      ggsave(sprintf("%s/heatmap_diff_tx_recall.png", plot_dir), p, width = 8, height = 6, dpi = 100)
      
      p_combined <- plot_combined_focused(data)
      ggsave(sprintf("%s/combi_distro_phases_%s.png", plot_dir, dtype), p_combined, width = 5.5, height = 10, dpi = 100)
      
      p <- plot_read_distribution_mutation1(data)
      ggsave(sprintf("%s/read_distribution_mutation1_%s.png", plot_dir, dtype), p, width = 9, height = 6, dpi = 100)
      
      p <- plot_phase_metrics_mutation1(data)
      ggsave(sprintf("%s/phase_metrics_mutation1_%s.png", plot_dir, dtype), p, width = 9, height = 6, dpi = 100)
    }
  }
  
  cat("\nAll plots generated successfully!\n")
}

generate_all_plots(data)

# Print summary statistics
cat("\n=== SUMMARY STATISTICS ===\n")
summary_stats <- data %>%
  group_by(data_type, mutation_level) %>%
  summarise(
    mean_f1_tx = mean(f1_tx, na.rm = TRUE),
    max_f1_tx = max(f1_tx, na.rm = TRUE),
    mean_f1_gene = mean(f1_gene, na.rm = TRUE),
    max_f1_gene = max(f1_gene, na.rm = TRUE),
    mean_diff_tx_recall = mean(diff_tx_recall, na.rm = TRUE),
    max_diff_tx_recall = max(diff_tx_recall, na.rm = TRUE),
    .groups = "drop"
  )
print(summary_stats)

# Find best configurations
cat("\n=== BEST CONFIGURATIONS BY DATA TYPE AND MUTATION LEVEL ===\n")
options(tibble.width = Inf, tibble.print_max = 50)
best_configs <- data %>%
  group_by(data_type, mutation_level) %>%
  slice_max(f1_tx, n = 1) %>%
  select(kmer, gmer, data_type, mutation_level, f1_tx, f1_gene, precision_tx, recall_tx, diff_tx_recall) %>%
  arrange(data_type, mutation_level, -f1_tx)
print(best_configs)
