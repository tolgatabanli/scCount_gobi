#!/usr/bin/Rscript
library(tidyverse)
library(cowplot)
library(igraph)

star_genes <- read_tsv('/mnt/raidbio2/extstud/praktikum/genprakt-ws25/gruppe_a/simulated_reads/counts_all_1/star_gene_counts.txt', skip=1, col_names = T) %>%
    select(c(1,7)) %>%
    rename(gene_id = 1, star_count = 2)
star_transcripts <- read_tsv('/mnt/raidbio2/extstud/praktikum/genprakt-ws25/gruppe_a/simulated_reads/counts_all_1/star_transcript_counts.txt', skip=1, col_names = T) %>%
    select(c(1,7)) %>%
    rename(transcript_id = 1, star_count = 2)

miniqut3_results <- read_tsv('out/counts_all_1_k23_g8/counts.tsv', col_names = c('id', 'miniqut3_count'))

# Load paralog information
paralogs_raw <- read_tsv('/mnt/raidbio2/extstud/praktikum/genprakt-ws25/gruppe_a/paralogs.tsv', 
                         col_names = c("gene", "paralogs"), show_col_types = FALSE)

# Create gene to group mapping using igraph
edges <- paralogs_raw %>%
    separate_rows(paralogs, sep = ",") %>%
    filter(paralogs != "")

g <- graph_from_data_frame(edges, directed = FALSE)
clusters <- components(g)

gene_to_group <- data.frame(
    gene = names(clusters$membership),
    group_id = paste0("grp_", clusters$membership)
) %>% as_tibble()

# compare gene-level
gene_data <- star_genes %>%
    inner_join(miniqut3_results, by = c('gene_id' = 'id')) %>%
    mutate(across(ends_with('count'), ~log2(. + 1)))

# INTERESTING:
# gene_data %>% filter(gene_id %in% c("ENSSSCG00000056624", "ENSSSCG00000041275", "ENSSSCG00000050853", "ENSSSCG00000050067", "ENSSSCG00000018834", "ENSSSCG00000058119", "ENSSSCG00000019719"))

gene_cor <- cor.test(gene_data$star_count, gene_data$miniqut3_count, method = 'pearson')
gene_r <- round(gene_cor$estimate, 3)
n_genes <- nrow(gene_data)
star_total_counts <- sum(star_genes$star_count)
miniqut3_total_counts <- sum(miniqut3_results$miniqut3_count)
gene_y_range <- max(gene_data$miniqut3_count) - min(gene_data$miniqut3_count)
gene_x_range <- max(gene_data$star_count) - min(gene_data$star_count)

# compare transcript-level
tx_data <- star_transcripts %>%
    inner_join(miniqut3_results, by = c('transcript_id' = 'id')) %>%
    mutate(across(ends_with('count'), ~log2(. + 1)))

tx_cor <- cor.test(tx_data$star_count, tx_data$miniqut3_count, method = 'pearson')
tx_r <- round(tx_cor$estimate, 3)
n_tx <- nrow(tx_data)
star_total_counts_tx <- sum(star_transcripts$star_count)
miniqut3_total_counts_tx <- sum(miniqut3_results$miniqut3_count)
tx_y_range <- max(tx_data$miniqut3_count) - min(tx_data$miniqut3_count)
tx_x_range <- max(tx_data$star_count) - min(tx_data$star_count)

# Set common scale limits for both plots
all_data <- bind_rows(gene_data, tx_data)
y_min <- min(all_data$miniqut3_count)
y_max <- max(all_data$miniqut3_count)
x_min <- min(all_data$star_count)
x_max <- max(all_data$star_count)

comparison_gene <- gene_data %>%
    ggplot(aes(x = star_count, y = miniqut3_count)) +
    geom_point(size = 0.5, alpha = 0.6, color = '#1f77b4') +
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = '#d62728', linewidth = 0.5) +
    scale_x_continuous(limits = c(x_min, x_max), 
                       breaks = scales::pretty_breaks(n = 10),
                       labels = scales::math_format(2^.x)) +
    scale_y_continuous(limits = c(y_min, y_max),
                       breaks = scales::pretty_breaks(n = 5),
                       labels = scales::math_format(2^.x)) +
    geom_segment(data = data.frame(x = x_min + (x_max - x_min) * 0.87, xend = x_min + (x_max - x_min) * 0.93,
                                    y = y_min + (y_max - y_min) * 0.15, yend = y_min + (y_max - y_min) * 0.15),
                 aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE,
                 linetype = 'dashed', color = '#d62728', linewidth = 0.75) +
    geom_label(data = data.frame(x = x_min + (x_max - x_min) * 0.95, y = y_min + (y_max - y_min) * 0.05,
                                  label = paste0('r=', gene_r, '\nSTAR: ', format(star_total_counts, big.mark=','), '\nminiQuT3: ', format(miniqut3_total_counts, big.mark=','))),
               aes(x = x, y = y, label = label), inherit.aes = FALSE,
               size = 3.5, hjust = 1, vjust = 0, fill = 'white', alpha = 0.8, color = 'black', label.padding = unit(0.3, 'lines')) +
    labs(title = 'Gene-Level', x = 'STAR 70% Overlap count', y = 'MiniQuT3 count') +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'))

comparison_tx <- tx_data %>%
    ggplot(aes(x = star_count, y = miniqut3_count)) +
    geom_point(size = 0.5, alpha = 0.6, color = '#1f77b4') +
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = '#d62728', linewidth = 0.5) +
    scale_x_continuous(limits = c(x_min, x_max),
                       breaks = scales::pretty_breaks(n = 10),
                       labels = scales::math_format(2^.x)) +
    scale_y_continuous(limits = c(y_min, y_max),
                       breaks = scales::pretty_breaks(n = 5),
                       labels = scales::math_format(2^.x)) +
    geom_segment(data = data.frame(x = x_min + (x_max - x_min) * 0.87, xend = x_min + (x_max - x_min) * 0.93,
                                    y = y_min + (y_max - y_min) * 0.15, yend = y_min + (y_max - y_min) * 0.15),
                 aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE,
                 linetype = 'dashed', color = '#d62728', linewidth = 0.75) +
    geom_label(data = data.frame(x = x_min + (x_max - x_min) * 0.95, y = y_min + (y_max - y_min) * 0.05,
                                  label = paste0('r=', tx_r, '\nSTAR: ', format(star_total_counts_tx, big.mark=','), '\nminiQuT3: ', format(miniqut3_total_counts_tx, big.mark=','))),
               aes(x = x, y = y, label = label), inherit.aes = FALSE,
               size = 3.5, hjust = 1, vjust = 0, fill = 'white', alpha = 0.8, color = 'black', label.padding = unit(0.3, 'lines')) +
    labs(title = 'Transcript-Level', x = 'STAR 70% Overlap count', y = 'MiniQuT3 count') +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'))

combined_plot <- plot_grid(comparison_gene, comparison_tx, ncol = 2, align = 'h', axis = 'tb')

ggsave('star_comparison.png', combined_plot, width = 12, height = 5.5, dpi = 300)

# ============================================================================
# Grouped by paralogs
# ============================================================================

# Group genes by paralogs
gene_data_grouped <- star_genes %>%
    left_join(gene_to_group, by = c('gene_id' = 'gene')) %>%
    mutate(group_id = ifelse(is.na(group_id), gene_id, group_id)) %>%
    group_by(group_id) %>%
    summarise(star_count = sum(star_count), .groups = 'drop') %>%
    rename(id = group_id) %>%
    inner_join(miniqut3_results, by = 'id') %>%
    mutate(across(ends_with('count'), ~log2(. + 1)))

# Group transcripts by paralogs (use gene mapping)
tx_data_grouped <- star_transcripts %>%
    left_join(gene_to_group, by = c('transcript_id' = 'gene')) %>%
    mutate(group_id = ifelse(is.na(group_id), transcript_id, group_id)) %>%
    group_by(group_id) %>%
    summarise(star_count = sum(star_count), .groups = 'drop') %>%
    rename(id = group_id) %>%
    inner_join(miniqut3_results, by = 'id') %>%
    mutate(across(ends_with('count'), ~log2(. + 1)))

# Calculate stats for grouped data
gene_grouped_cor <- cor.test(gene_data_grouped$star_count, gene_data_grouped$miniqut3_count, method = 'pearson')
gene_grouped_r <- round(gene_grouped_cor$estimate, 3)
n_groups_gene <- nrow(gene_data_grouped)
star_total_counts_grouped <- sum(star_genes$star_count)
miniqut3_total_counts_grouped <- sum(miniqut3_results$miniqut3_count)
gene_grouped_y_range <- max(gene_data_grouped$miniqut3_count) - min(gene_data_grouped$miniqut3_count)
gene_grouped_x_range <- max(gene_data_grouped$star_count) - min(gene_data_grouped$star_count)

tx_grouped_cor <- cor.test(tx_data_grouped$star_count, tx_data_grouped$miniqut3_count, method = 'pearson')
tx_grouped_r <- round(tx_grouped_cor$estimate, 3)
n_groups_tx <- nrow(tx_data_grouped)
star_total_counts_grouped_tx <- sum(star_transcripts$star_count)
miniqut3_total_counts_grouped_tx <- sum(miniqut3_results$miniqut3_count)
tx_grouped_y_range <- max(tx_data_grouped$miniqut3_count) - min(tx_data_grouped$miniqut3_count)
tx_grouped_x_range <- max(tx_data_grouped$star_count) - min(tx_data_grouped$star_count)

# Set common scale limits for grouped plots
all_data_grouped <- bind_rows(gene_data_grouped, tx_data_grouped)
y_min_grouped <- min(all_data_grouped$miniqut3_count)
y_max_grouped <- max(all_data_grouped$miniqut3_count)
x_min_grouped <- min(all_data_grouped$star_count)
x_max_grouped <- max(all_data_grouped$star_count)

comparison_gene_grouped <- gene_data_grouped %>%
    ggplot(aes(x = star_count, y = miniqut3_count)) +
    geom_point(size = 0.5, alpha = 0.6, color = '#1f77b4') +
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = '#d62728', linewidth = 0.5) +
    scale_x_continuous(limits = c(x_min_grouped, x_max_grouped),
                       breaks = scales::pretty_breaks(n = 10),
                       labels = scales::math_format(2^.x)) +
    scale_y_continuous(limits = c(y_min_grouped, y_max_grouped),
                       breaks = scales::pretty_breaks(n = 5),
                       labels = scales::math_format(2^.x)) +
    geom_segment(data = data.frame(x = x_min_grouped + (x_max_grouped - x_min_grouped) * 0.87, xend = x_min_grouped + (x_max_grouped - x_min_grouped) * 0.93,
                                    y = y_min_grouped + (y_max_grouped - y_min_grouped) * 0.15, yend = y_min_grouped + (y_max_grouped - y_min_grouped) * 0.15),
                 aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE,
                 linetype = 'dashed', color = '#d62728', linewidth = 0.75) +
    geom_label(data = data.frame(x = x_min_grouped + (x_max_grouped - x_min_grouped) * 0.95, y = y_min_grouped + (y_max_grouped - y_min_grouped) * 0.05,
                                  label = paste0('r=', gene_grouped_r, '\nSTAR: ', format(star_total_counts_grouped, big.mark=','), '\nminiQuT3: ', format(miniqut3_total_counts_grouped, big.mark=','))),
               aes(x = x, y = y, label = label), inherit.aes = FALSE,
               size = 3.5, hjust = 1, vjust = 0, fill = 'white', alpha = 0.8, color = 'black', label.padding = unit(0.3, 'lines')) +
    labs(title = 'Gene-Level (Grouped by Paralogs)', x = 'STAR 70% Overlap count', y = 'MiniQuT3 count') +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'))

comparison_tx_grouped <- tx_data_grouped %>%
    ggplot(aes(x = star_count, y = miniqut3_count)) +
    geom_point(size = 0.5, alpha = 0.6, color = '#1f77b4') +
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = '#d62728', linewidth = 0.5) +
    scale_x_continuous(limits = c(x_min_grouped, x_max_grouped),
                       breaks = scales::pretty_breaks(n = 10),
                       labels = scales::math_format(2^.x)) +
    scale_y_continuous(limits = c(y_min_grouped, y_max_grouped),
                       breaks = scales::pretty_breaks(n = 5),
                       labels = scales::math_format(2^.x)) +
    geom_segment(data = data.frame(x = x_min_grouped + (x_max_grouped - x_min_grouped) * 0.87, xend = x_min_grouped + (x_max_grouped - x_min_grouped) * 0.93,
                                    y = y_min_grouped + (y_max_grouped - y_min_grouped) * 0.15, yend = y_min_grouped + (y_max_grouped - y_min_grouped) * 0.15),
                 aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE,
                 linetype = 'dashed', color = '#d62728', linewidth = 0.75) +
    geom_label(data = data.frame(x = x_min_grouped + (x_max_grouped - x_min_grouped) * 0.95, y = y_min_grouped + (y_max_grouped - y_min_grouped) * 0.05,
                                  label = paste0('r=', tx_grouped_r, '\nSTAR: ', format(star_total_counts_grouped_tx, big.mark=','), '\nminiQuT3: ', format(miniqut3_total_counts_grouped_tx, big.mark=','))),
               aes(x = x, y = y, label = label), inherit.aes = FALSE,
               size = 3.5, hjust = 1, vjust = 0, fill = 'white', alpha = 0.8, color = 'black', label.padding = unit(0.3, 'lines')) +
    labs(title = 'Transcript-Level (Grouped by Paralogs)', x = 'STAR 70% Overlap count', y = 'MiniQuT3 count') +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'))

combined_plot_grouped <- plot_grid(comparison_gene_grouped, comparison_tx_grouped, ncol = 2, align = 'h', axis = 'tb')

ggsave('star_comparison_grouped_paralogs.png', combined_plot_grouped, width = 12, height = 5.5, dpi = 300)

# ============================================================================
# STAR vs Ground Truth
# ============================================================================

# Load ground truth data
ground_truth <- read_tsv('../simulation/counts_all.tsv', col_names = T) %>%
    rename(true_count = count)

# Gene-level comparison: STAR vs Ground Truth
ground_truth_gene <- ground_truth %>%
    group_by(gene) %>%
    summarise(true_count = sum(true_count), .groups = 'drop')

gene_truth_data <- star_genes %>%
    inner_join(ground_truth_gene, by = c('gene_id' = 'gene')) %>%
    mutate(across(ends_with('count'), ~log2(. + 1)))

gene_truth_cor <- cor.test(gene_truth_data$star_count, gene_truth_data$true_count, method = 'pearson')
gene_truth_r <- round(gene_truth_cor$estimate, 3)
star_total_truth <- sum(ground_truth$true_count)

# Transcript-level comparison: STAR vs Ground Truth
tx_truth_data <- star_transcripts %>%
    inner_join(ground_truth, by = c('transcript_id' = 'transcript')) %>%
    mutate(across(ends_with('count'), ~log2(. + 1)))

tx_truth_cor <- cor.test(tx_truth_data$star_count, tx_truth_data$true_count, method = 'pearson')
tx_truth_r <- round(tx_truth_cor$estimate, 3)

# Set common scale limits for STAR vs Truth plots
all_data_truth <- bind_rows(gene_truth_data, tx_truth_data)
y_min_truth <- min(all_data_truth$true_count)
y_max_truth <- max(all_data_truth$true_count)
x_min_truth <- min(all_data_truth$star_count)
x_max_truth <- max(all_data_truth$star_count)

comparison_gene_truth <- gene_truth_data %>%
    ggplot(aes(x = star_count, y = true_count)) +
    geom_point(size = 0.5, alpha = 0.6, color = '#1f77b4') +
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = '#d62728', linewidth = 0.5) +
    scale_x_continuous(limits = c(x_min_truth, x_max_truth), 
                       breaks = scales::pretty_breaks(n = 10),
                       labels = scales::math_format(2^.x)) +
    scale_y_continuous(limits = c(y_min_truth, y_max_truth),
                       breaks = scales::pretty_breaks(n = 5),
                       labels = scales::math_format(2^.x)) +
    geom_segment(data = data.frame(x = x_min_truth + (x_max_truth - x_min_truth) * 0.87, xend = x_min_truth + (x_max_truth - x_min_truth) * 0.93,
                                    y = y_min_truth + (y_max_truth - y_min_truth) * 0.15, yend = y_min_truth + (y_max_truth - y_min_truth) * 0.15),
                 aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE,
                 linetype = 'dashed', color = '#d62728', linewidth = 0.75) +
    geom_label(data = data.frame(x = x_min_truth + (x_max_truth - x_min_truth) * 0.95, y = y_min_truth + (y_max_truth - y_min_truth) * 0.05,
                                  label = paste0('r=', gene_truth_r, '\nSTAR: ', format(sum(star_genes$star_count), big.mark=','), '\nGround Truth: ', format(star_total_truth, big.mark=','))),
               aes(x = x, y = y, label = label), inherit.aes = FALSE,
               size = 3.5, hjust = 1, vjust = 0, fill = 'white', alpha = 0.8, color = 'black', label.padding = unit(0.3, 'lines')) +
    labs(title = 'Gene-Level', x = 'STAR 70% Overlap count', y = 'Ground Truth count') +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'))

comparison_tx_truth <- tx_truth_data %>%
    ggplot(aes(x = star_count, y = true_count)) +
    geom_point(size = 0.5, alpha = 0.6, color = '#1f77b4') +
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = '#d62728', linewidth = 0.5) +
    scale_x_continuous(limits = c(x_min_truth, x_max_truth),
                       breaks = scales::pretty_breaks(n = 10),
                       labels = scales::math_format(2^.x)) +
    scale_y_continuous(limits = c(y_min_truth, y_max_truth),
                       breaks = scales::pretty_breaks(n = 5),
                       labels = scales::math_format(2^.x)) +
    geom_segment(data = data.frame(x = x_min_truth + (x_max_truth - x_min_truth) * 0.87, xend = x_min_truth + (x_max_truth - x_min_truth) * 0.93,
                                    y = y_min_truth + (y_max_truth - y_min_truth) * 0.15, yend = y_min_truth + (y_max_truth - y_min_truth) * 0.15),
                 aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE,
                 linetype = 'dashed', color = '#d62728', linewidth = 0.75) +
    geom_label(data = data.frame(x = x_min_truth + (x_max_truth - x_min_truth) * 0.95, y = y_min_truth + (y_max_truth - y_min_truth) * 0.05,
                                  label = paste0('r=', tx_truth_r, '\nSTAR: ', format(sum(star_transcripts$star_count), big.mark=','), '\nGround Truth: ', format(star_total_truth, big.mark=','))),
               aes(x = x, y = y, label = label), inherit.aes = FALSE,
               size = 3.5, hjust = 1, vjust = 0, fill = 'white', alpha = 0.8, color = 'black', label.padding = unit(0.3, 'lines')) +
    labs(title = 'Transcript-Level', x = 'STAR 70% Overlap count', y = 'Ground Truth count') +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'))

combined_plot_truth <- plot_grid(comparison_gene_truth, comparison_tx_truth, ncol = 2, align = 'h', axis = 'tb')

ggsave('star_vs_groundtruth.png', combined_plot_truth, width = 12, height = 5.5, dpi = 300)

## Pig Data
miniqut3_pig <- read_tsv('/home/t/tabanli/Desktop/scCount/out/pig_counts/counts.tsv', col_names = c('id', 'miniqut3_count'))
star_pig_genes <- read_tsv('/mnt/raidbio2/extstud/praktikum/genprakt-ws25/gruppe_a/data/feature_Counts/star_gene_70p.txt', skip=1, col_names = T) %>%
    select(c(1,7)) %>%
    rename(gene_id = 1, star_count = 2)
star_pig_transcripts <- read_tsv('/mnt/raidbio2/extstud/praktikum/genprakt-ws25/gruppe_a/data/feature_Counts/star_70p.txt', skip=1, col_names = T) %>%
    select(c(1,7)) %>%
    rename(transcript_id = 1, star_count = 2)

# Compare pig gene-level
pig_gene_data <- star_pig_genes %>%
    inner_join(miniqut3_pig, by = c('gene_id' = 'id')) %>%
    mutate(across(ends_with('count'), ~log2(. + 1)))

pig_gene_cor <- cor.test(pig_gene_data$star_count, pig_gene_data$miniqut3_count, method = 'pearson')
pig_gene_r <- round(pig_gene_cor$estimate, 3)
pig_star_total_counts <- sum(star_pig_genes$star_count)
pig_miniqut3_total_counts <- sum(miniqut3_pig$miniqut3_count)

# Compare pig transcript-level
pig_tx_data <- star_pig_transcripts %>%
    inner_join(miniqut3_pig, by = c('transcript_id' = 'id')) %>%
    mutate(across(ends_with('count'), ~log2(. + 1)))

pig_tx_cor <- cor.test(pig_tx_data$star_count, pig_tx_data$miniqut3_count, method = 'pearson')
pig_tx_r <- round(pig_tx_cor$estimate, 3)
pig_star_total_counts_tx <- sum(star_pig_transcripts$star_count)
pig_miniqut3_total_counts_tx <- sum(miniqut3_pig$miniqut3_count)

# Set common scale limits for pig plots
pig_all_data <- bind_rows(pig_gene_data, pig_tx_data)
pig_y_min <- min(pig_all_data$miniqut3_count)
pig_y_max <- max(pig_all_data$miniqut3_count)
pig_x_min <- min(pig_all_data$star_count)
pig_x_max <- max(pig_all_data$star_count)

comparison_pig_gene <- pig_gene_data %>%
    ggplot(aes(x = star_count, y = miniqut3_count)) +
    geom_point(size = 0.5, alpha = 0.6, color = '#1f77b4') +
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = '#d62728', linewidth = 0.5) +
    scale_x_continuous(limits = c(pig_x_min, pig_x_max), 
                       breaks = scales::pretty_breaks(n = 10),
                       labels = scales::math_format(2^.x)) +
    scale_y_continuous(limits = c(pig_y_min, pig_y_max),
                       breaks = scales::pretty_breaks(n = 5),
                       labels = scales::math_format(2^.x)) +
    geom_segment(data = data.frame(x = pig_x_min + (pig_x_max - pig_x_min) * 0.87, xend = pig_x_min + (pig_x_max - pig_x_min) * 0.93,
                                    y = pig_y_min + (pig_y_max - pig_y_min) * 0.15, yend = pig_y_min + (pig_y_max - pig_y_min) * 0.15),
                 aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE,
                 linetype = 'dashed', color = '#d62728', linewidth = 0.75) +
    geom_label(data = data.frame(x = pig_x_min + (pig_x_max - pig_x_min) * 0.95, y = pig_y_min + (pig_y_max - pig_y_min) * 0.05,
                                  label = paste0('r=', pig_gene_r, '\nSTAR: ', format(pig_star_total_counts, big.mark=','), '\nminiQuT3: ', format(pig_miniqut3_total_counts, big.mark=','))),
               aes(x = x, y = y, label = label), inherit.aes = FALSE,
               size = 3.5, hjust = 1, vjust = 0, fill = 'white', alpha = 0.8, color = 'black', label.padding = unit(0.3, 'lines')) +
    labs(title = 'Gene-Level (Pig)', x = 'STAR 70% Overlap count', y = 'MiniQuT3 count') +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'))

comparison_pig_tx <- pig_tx_data %>%
    ggplot(aes(x = star_count, y = miniqut3_count)) +
    geom_point(size = 0.5, alpha = 0.6, color = '#1f77b4') +
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = '#d62728', linewidth = 0.5) +
    scale_x_continuous(limits = c(pig_x_min, pig_x_max),
                       breaks = scales::pretty_breaks(n = 10),
                       labels = scales::math_format(2^.x)) +
    scale_y_continuous(limits = c(pig_y_min, pig_y_max),
                       breaks = scales::pretty_breaks(n = 5),
                       labels = scales::math_format(2^.x)) +
    geom_segment(data = data.frame(x = pig_x_min + (pig_x_max - pig_x_min) * 0.87, xend = pig_x_min + (pig_x_max - pig_x_min) * 0.93,
                                    y = pig_y_min + (pig_y_max - pig_y_min) * 0.15, yend = pig_y_min + (pig_y_max - pig_y_min) * 0.15),
                 aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE,
                 linetype = 'dashed', color = '#d62728', linewidth = 0.75) +
    geom_label(data = data.frame(x = pig_x_min + (pig_x_max - pig_x_min) * 0.95, y = pig_y_min + (pig_y_max - pig_y_min) * 0.05,
                                  label = paste0('r=', pig_tx_r, '\nSTAR: ', format(pig_star_total_counts_tx, big.mark=','), '\nminiQuT3: ', format(pig_miniqut3_total_counts_tx, big.mark=','))),
               aes(x = x, y = y, label = label), inherit.aes = FALSE,
               size = 3.5, hjust = 1, vjust = 0, fill = 'white', alpha = 0.8, color = 'black', label.padding = unit(0.3, 'lines')) +
    labs(title = 'Transcript-Level (Pig)', x = 'STAR 70% Overlap count', y = 'MiniQuT3 count') +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'))

combined_plot_pig <- plot_grid(comparison_pig_gene, comparison_pig_tx, ncol = 2, align = 'h', axis = 'tb')

ggsave('pig_star_comparison.png', combined_plot_pig, width = 12, height = 5.5, dpi = 300)

