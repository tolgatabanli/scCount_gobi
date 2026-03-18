#!/usr/bin/Rscript
library(tidyverse)
library(igraph)

## ignore
genes_star <- read_tsv('?', skip=1) %>%
    select("Geneid", "Aligned.sortedByCoord.out.bam") %>%
    rename(gene = Geneid, count = Aligned.sortedByCoord.out.bam)

p <- genes_star %>% inner_join(obs_genes, by = c("gene" = "id"), suffix = c("_star", "_obs")) %>%
    ggplot(aes(x = count_star, y = count_obs)) +
    geom_point()
ggsave(filename = "simul_noMut_gamma_withStar.png", p)

# star against truth
p <- genes_star %>% inner_join(genes_truth, by = c("gene" = "gene"), suffix = c("_star", "_truth")) %>%
    ggplot(aes(x = count_star, y = count_truth)) +
    geom_point()
ggsave(filename = "simul_noMut_gamma_star_truth.png", p)
## ignore end


paralogs_raw <- read_tsv('/mnt/raidbio2/extstud/praktikum/genprakt-ws25/gruppe_a/paralogs.tsv', col_names = c("gene", "paralogs"), show_col_types = FALSE)

# Create an edge list for the network
edges <- paralogs_raw %>%
    separate_rows(paralogs, sep = ",") %>%
    filter(paralogs != "") # Clean up empty strings

# Use igraph to find connected components (clusters of paralogs)
g <- graph_from_data_frame(edges, directed = FALSE)
clusters <- components(g)

gene_to_group <- data.frame(
    gene = names(clusters$membership),
    group_id = paste0("grp_", clusters$membership)
) %>% as_tibble()


create_scatter <- function(df_truth, df_obs, label, is_grouped) {
  merged <- df_truth %>% 
    inner_join(df_obs, by = "id", suffix = c("_truth", "_obs"))
  
  correlation <- cor(merged$count_truth, merged$count_obs, method = "pearson")
  
  # Calculate position for annotation (Top Left)
  # We take the minimum X and maximum Y to keep it in the corner
  x_pos <- min(merged$count_truth, na.rm = TRUE)
  y_pos <- max(merged$count_obs, na.rm = TRUE)
  
  ggplot(merged, aes(x = count_truth, y = count_obs)) +
    geom_point(alpha = 0.4, size = 1) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    # Use actual data points for positioning
    annotate("text", x = x_pos, y = y_pos, 
             label = paste0("R = ", round(correlation, 4)), 
             hjust = 0, vjust = 1, fontface = "bold", size = 5) +
    labs(title = paste(label, "Comparison"),
         subtitle = ifelse(is_grouped, "Grouped by Paralogs", "Individual Genes"),
         x = "Truth Counts",
         y = "Observed Counts") +
    scale_x_log10() +
    scale_y_log10() +
    theme_minimal()
}

#' Compare Truth vs Observed Counts with Paralog Grouping
#' 
#' @param truth_path Path to the truth .tsv file
#' @param obs_path Path to the observed .tsv file (id, count)
#' @param gene_to_group Data frame mapping genes to paralog groups (columns: gene, group_id)
#' @param output_prefix String for output naming
plot_counts_comparison <- function(truth_path, obs_path, gene_to_group = NULL, output_prefix = "comparison") {
  
  # 1. Load Data
  truth <- read_tsv(truth_path, show_col_types = FALSE)
  obs <- read_tsv(obs_path, col_names = c("id", "count"), show_col_types = FALSE)
  
  # 3. Process Truth Data (Gene Level)
  genes_truth <- truth %>% 
    group_by(gene) %>% 
    summarise(count = sum(count), .groups = 'drop')
  
  # 4. Apply Paralog Grouping to Genes
  if (!is.null(gene_to_group)) {
    # Left join to keep all genes; genes without paralogs use their own ID as group_id
    genes_truth <- genes_truth %>%
      left_join(gene_to_group, by = "gene") %>%
      mutate(group_id = ifelse(is.na(group_id), gene, group_id)) %>%
      group_by(group_id) %>%
      summarise(count = sum(count), .groups = 'drop') %>%
      rename(id = group_id)
    
    obs_grouped <- obs %>%
      left_join(gene_to_group, by = c("id" = "gene")) %>%
      mutate(group_id = ifelse(is.na(group_id), id, group_id)) %>%
      group_by(group_id) %>%
      summarise(count = sum(count), .groups = 'drop') %>%
      rename(id = group_id)
  } else {
    # Standard behavior if no paralogs file
    genes_truth <- genes_truth %>% rename(id = gene)
    obs_grouped <- obs
  }  
  

  p_gene <- create_scatter(genes_truth, obs_grouped, "Gene/Group", !is.null(gene_to_group))
  
  ggsave(filename = paste0(output_prefix, "_scatter.png"), plot = p_gene, width = 8, height = 7)
  
  message("Plot saved: ", output_prefix, "_scatter.png")
}


plot_counts_comparison('/mnt/cip/home/t/tabanli/Desktop/scCount/simulation/simcounts_1_5000_protcoding.tsv',
            '/mnt/cip/home/t/tabanli/Desktop/scCount/out/count_simul_noMut_gamma_protcoding/counts.tsv',
            gene_to_group,
            "simul_noMut_gamma_protcoding")
plot_counts_comparison('/mnt/cip/home/t/tabanli/Desktop/scCount/simulation/simcounts_1_5000.tsv',
            '/mnt/cip/home/t/tabanli/Desktop/scCount/out/count_simul_noMut_gamma_all/counts.tsv',
            gene_to_group,
            "simul_noMut_gamma_all")
plot_counts_comparison('/mnt/cip/home/t/tabanli/Desktop/scCount/simulation/simcounts_1_5000_pseudogene.tsv',
            '/mnt/cip/home/t/tabanli/Desktop/scCount/out/count_simul_noMut_gamma_pseudogene/counts.tsv',
            gene_to_group,
            "simul_noMut_gamma_pseudogene")
plot_counts_comparison('/mnt/cip/home/t/tabanli/Desktop/scCount/simulation/simcounts_1_5000.tsv',
            '/mnt/cip/home/t/tabanli/Desktop/scCount/out/count_simul_1mut_gamma_all/counts.tsv',
            gene_to_group,
            "simul_1mut_gamma_all")

# triplets
plot_counts_comparison('/mnt/cip/home/t/tabanli/Desktop/scCount/simulation/simcounts_1_5000.tsv',
            '/home/t/tabanli/Desktop/scCount/out/count_triplet_simcounts_1_5000/counts.tsv',
            NULL,
            "simul_noMut_gamma_triplet")
plot_counts_comparison('/mnt/cip/home/t/tabanli/Desktop/scCount/simulation/simcounts_1_5000.tsv',
            '/home/t/tabanli/Desktop/scCount/out/count_triplet_simcounts_1_5000_1/counts.tsv',
            NULL,
            "simul_noMut_gamma_triplet_1mut")
plot_counts_comparison('/mnt/cip/home/t/tabanli/Desktop/scCount/simulation/simcounts_1_5000.tsv',
            '/home/t/tabanli/Desktop/scCount/out/count_triplet_simcounts_1_5000/counts.tsv',
            gene_to_group,
            "simul_noMut_gamma_triplet_grouped")

plot_counts_comparison('/mnt/cip/home/t/tabanli/Desktop/scCount/simulation/counts_all.tsv',
            '/mnt/cip/home/t/tabanli/Desktop/scCount/benchmark/out/counts_all_0_k15_g5/counts.tsv',
            NULL,
            "tmp0")
