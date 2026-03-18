#!/usr/bin/Rscript
library(tidyverse)
library(data.table)
library(stringr)
library(parallel)

# --- Configuration ---
gtf_path <- "/mnt/raidbio2/extstud/praktikum/genprakt-ws25/gruppe_a/data/Sus_scrofa.Sscrofa11.1.gtf"
sim_dir  <- "../simulation"
out_base <- "out"
results_dir <- "analysis_results"
k_values <- c(23)
g_values <- c(11)
biotype_values <- c("all")

dir.create(results_dir)

# 1) GTF: transcript -> gene mapping with biotype
gtf_raw <- read_lines(gtf_path)
# Extract lines containing transcript_id and gene_id
gtf_lines <- gtf_raw[str_detect(gtf_raw, 'transcript_id "') & str_detect(gtf_raw, 'gene_id "')]
tx_to_gene <- data.frame(
  transcript_id = str_extract(gtf_lines, 'transcript_id "([^"]+)"', group = 1),
  gene_id = str_extract(gtf_lines, 'gene_id "([^"]+)"', group = 1),
  biotype = str_extract(gtf_lines, 'gene_biotype "([^"]+)"', group = 1)
) %>% distinct()

# 2. Identify all benchmark runs
all_run_dirs <- list.dirs(out_base, full.names = TRUE, recursive = FALSE)

# Filter by biotype, k and g values if specified
run_dirs <- all_run_dirs
if (!is.null(biotype_values) | !is.null(k_values) | !is.null(g_values)) {
  run_dirs <- sapply(all_run_dirs, function(dir) {
    dir_name <- basename(dir)
    
    # Extract biotype from directory name: counts_BIOTYPE_MUTRATE_kNUM_gNUM
    biotype_match <- str_extract(dir_name, "^counts_([^_]+)_", group = 1)
    k_match <- str_extract(dir_name, "_k(\\d+)_", group = 1)
    g_match <- str_extract(dir_name, "_g(\\d+)$", group = 1)
    
    biotype_ok <- is.null(biotype_values) || (!is.na(biotype_match) && biotype_match %in% biotype_values)
    k_ok <- is.null(k_values) || (!is.na(k_match) && as.numeric(k_match) %in% k_values)
    g_ok <- is.null(g_values) || (!is.na(g_match) && as.numeric(g_match) %in% g_values)
    
    if (biotype_ok && k_ok && g_ok) dir else NA_character_
  }, USE.NAMES = FALSE)
  run_dirs <- run_dirs[!is.na(run_dirs)]
  
  filter_msg <- c()
  if (!is.null(biotype_values)) filter_msg <- c(filter_msg, paste("biotype =", paste(biotype_values, collapse = ",")))
  if (!is.null(k_values)) filter_msg <- c(filter_msg, paste("k =", paste(k_values, collapse = ",")))
  if (!is.null(g_values)) filter_msg <- c(filter_msg, paste("g =", paste(g_values, collapse = ",")))
  
  cat("Filtered runs based on:", paste(filter_msg, collapse = ", "), "\n")
  cat("Selected", length(run_dirs), "directories for analysis\n\n")
}

# Process a single run (used for mclapply)
process_run <- function(run, sim_dir, results_dir, tx_to_gene) {
  run_name <- basename(run)
  cat("Processing:", run_name, "\n")
  
  sim_name <- str_remove(run_name, "_k[0-9]+_g[0-9]+$")
  ground_truth_file <- file.path(sim_dir, paste0(sub("_[0-9]+$", "", sim_name), ".tsv"))
  
  count_file  <- file.path(run, "counts.tsv")
  summary_file <- file.path(run, "miniqut3.summary")
  assign_file <- file.path(run, "read_assignments.tsv") # single column with gene or transcript id (whichever could be assigned) - for amgiguous, * for too short
  mapping_info <- file.path(sim_dir, sim_name, "read.mappinginfo") # truth
  
  # --- Part A: Summary Stats ---
  summary_text <- read_file(summary_file)
  stats <- list(
    run = run_name,
    total_reads = as.numeric(str_extract(summary_text, "Total reads processed: (\\d+)", group = 1)),
    assigned_tx = as.numeric(str_extract(summary_text, "Reads assigned to transcripts: (\\d+)", group = 1)),
    assigned_gene = as.numeric(str_extract(summary_text, "Reads assigned to genes: (\\d+)", group = 1)),
    ambig_reads = as.numeric(str_extract(summary_text, "Ambiguous reads: (\\d+)", group = 1)),
    short_reads = as.numeric(str_extract(summary_text, "Short reads discarded: (\\d+)", group = 1)),
    reads_heuristic_tx = as.numeric(str_extract(summary_text, "Heuristic phase - transcripts: (\\d+)", group = 1)),
    reads_heuristic_gene = as.numeric(str_extract(summary_text, "Heuristic phase - genes: (\\d+)", group = 1)),
    reads_scoring_tx = as.numeric(str_extract(summary_text, "Scoring phase - transcripts: (\\d+)", group = 1)),
    reads_scoring_gene = as.numeric(str_extract(summary_text, "Scoring phase - genes: (\\d+)", group = 1))
  )

  # --- Part B: Correlation and MSE ---
  truth <- fread(ground_truth_file, col.names = c("gene", "transcript", "true_count"))
  tool_counts <- fread(count_file, col.names = c("id", "est_count"))
  
  # -- Tx-level --
  comparison_tx <- tool_counts %>%
    inner_join(truth, by = c("id" = "transcript")) %>%
    mutate(true_count = replace_na(true_count, 0)) %>%
    mutate(across(ends_with('count'), ~log2(. + 1)))
  
  # Pearson correlation on raw counts
  pearson_tx <- cor(comparison_tx$est_count, comparison_tx$true_count, method = "pearson")
  
  # -- Gene-level (use gene IDs directly from count file, already aggregated) --
  tool_counts_gene <- tool_counts %>%
    filter(id %in% tx_to_gene$gene_id) %>%
    rename(gene_id = id)
  
  truth_gene <- truth %>%
    group_by(gene) %>%
    summarise(true_count = sum(true_count), .groups = 'drop')
  
  comparison_gene <- tool_counts_gene %>%
    left_join(truth_gene, by = c("gene_id" = "gene")) %>%
    mutate(true_count = replace_na(true_count, 0)) %>%
    mutate(across(ends_with('count'), ~log2(. + 1)))
  
  # Pearson correlation on raw counts
  pearson_gene <- cor(comparison_gene$est_count, comparison_gene$true_count, method = "pearson")
  
  # Placeholder storage for metrics to be calculated from per-read data (below)
  tx_precision <- NA
  tx_recall <- NA
  tx_f1 <- NA
  gene_precision <- NA
  gene_recall <- NA
  gene_f1 <- NA
  diff_tx_recall <- NA
  diff_tx_num <- NA

  # --- Part C: Per-Read Assignment Accuracy ---
  cat("  Calculating per-read accuracy...\n")
  assigns <- fread(assign_file, header = FALSE, col.names = "assigned_id")
  mapping <- fread(mapping_info, select = 3:4) # 3rd col Gene, 4th col Tx
  colnames(mapping) <- c("true_gene", "true_tx")
  
  # Combine assignment data with ground truth
  acc_data <- data.frame(
    assigned = assigns$assigned_id,
    true_gene = mapping$true_gene,
    true_tx = mapping$true_tx
  )
  
  # --- Transcript-level metrics ---
  # TP: assigned matches true transcript
  # FP: assigned doesn't match true transcript AND is not unassigned
  # FN: unassigned reads ('-' or '*')
  
  acc_data$is_unassigned_tx <- acc_data$assigned %in% c('-', '*')
  acc_data$is_correct_tx <- acc_data$assigned == acc_data$true_tx & !acc_data$is_unassigned_tx
  acc_data$is_wrong_tx <- acc_data$assigned != acc_data$true_tx & !acc_data$is_unassigned_tx
  
  tx_tp <- sum(acc_data$is_correct_tx)
  tx_fp <- sum(acc_data$is_wrong_tx)
  tx_fn <- sum(acc_data$is_wrong_tx) + sum(acc_data$is_unassigned_tx)
  
  tx_precision <- ifelse((tx_tp + tx_fp) > 0, tx_tp / (tx_tp + tx_fp), 0)
  tx_recall <- ifelse((tx_tp + tx_fn) > 0, tx_tp / (tx_tp + tx_fn), 0)
  tx_f1 <- ifelse((tx_precision + tx_recall) > 0, 2 * (tx_precision * tx_recall) / (tx_precision + tx_recall), 0)
  
  # --- Gene-level metrics ---
  # Map true transcripts to genes
  acc_data <- acc_data %>%
    left_join(tx_to_gene, by = c("true_tx" = "transcript_id")) %>%
    rename(true_gene_mapped = gene_id) %>%
    mutate(true_gene_mapped = coalesce(true_gene_mapped, true_gene))
  
  # Map assigned IDs to genes
  acc_data <- acc_data %>%
    mutate(
      assigned_gene = case_when(
        assigned %in% c('-', '*') ~ assigned,
        assigned %in% tx_to_gene$transcript_id ~ tx_to_gene$gene_id[match(assigned, tx_to_gene$transcript_id)],
        TRUE ~ assigned  # Already a gene ID
      )
    )
  
  acc_data$is_unassigned_gene <- acc_data$assigned_gene %in% c('-', '*')
  acc_data$is_correct_gene <- acc_data$assigned_gene == acc_data$true_gene_mapped & !acc_data$is_unassigned_gene
  acc_data$is_wrong_gene <- acc_data$assigned_gene != acc_data$true_gene_mapped & !acc_data$is_unassigned_gene
  
  gene_tp <- sum(acc_data$is_correct_gene)
  gene_fp <- sum(acc_data$is_wrong_gene)
  gene_fn <- sum(acc_data$is_wrong_gene) + sum(acc_data$is_unassigned_gene)
  
  gene_precision <- ifelse((gene_tp + gene_fp) > 0, gene_tp / (gene_tp + gene_fp), 0)
  gene_recall <- ifelse((gene_tp + gene_fn) > 0, gene_tp / (gene_tp + gene_fn), 0)
  gene_f1 <- ifelse((gene_precision + gene_recall) > 0, 2 * (gene_precision * gene_recall) / (gene_precision + gene_recall), 0)
  
  # --- Differential Transcript Recall (for isoform genes) ---
  # Identify genes with 2+ isoforms in both GTF and ground truth
  isoforms_in_gtf <- tx_to_gene %>%
    group_by(gene_id) %>%
    summarise(n_isoforms = n_distinct(transcript_id), .groups = 'drop') %>%
    filter(n_isoforms >= 2)
  
  isoforms_in_truth <- truth %>%
    group_by(gene) %>%
    summarise(n_isoforms = n_distinct(transcript), .groups = 'drop') %>%
    filter(n_isoforms >= 2)
  
  isoform_genes <- intersect(isoforms_in_gtf$gene_id, isoforms_in_truth$gene)
  diff_tx_num <- length(isoform_genes)
  if (diff_tx_num > 0) {
    # Filter acc_data to isoform genes
    isoform_data <- acc_data %>%
      filter(true_gene_mapped %in% isoform_genes)
    
    # Calculate per-gene transcript recall for isoform genes
    per_gene_recalls <- isoform_data %>%
      group_by(true_gene_mapped) %>%
      summarise(
        correct = sum(is_correct_tx),
        total = n(),
        .groups = 'drop'
      ) %>%
      mutate(gene_tx_recall = correct / total) %>%
      pull(gene_tx_recall)
    
    # Average across isoform genes
    diff_tx_recall <- mean(per_gene_recalls, na.rm = TRUE)
  } else {
    diff_tx_recall <- NA
  }
  
  # Simple plot for gene-level comparison with F1-score
  p <- ggplot(comparison_gene, aes(x = true_count, y = est_count)) +
    geom_point(alpha = 0.6, size = 0.6, color = "#2076b4") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10), labels = scales::math_format(2^.x), expand = c(0.05, 0)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5), labels = scales::math_format(2^.x), expand = c(0.05, 0)) +
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = '#d62728', linewidth = 0.5) +
    labs(title = paste(run_name, "(Gene-level)"), subtitle = paste("Pearson:", round(pearson_gene, 4), "Precision:", round(gene_precision, 4), "Recall:", round(gene_recall, 4), "F1:", round(gene_f1, 4)), x = "True Counts", y = "MiniQuT3 Counts") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'))
  
  ggsave(file.path(results_dir, paste0(run_name, "_corr_gene.png")), p)
  
  # Plot for transcript-level comparison with F1-score
  p_tx <- ggplot(comparison_tx, aes(x = true_count, y = est_count)) +
    geom_point(alpha = 0.6, size = 0.6, color = "#2076b4") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10), labels = scales::math_format(2^.x), expand = c(0.05, 0)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5), labels = scales::math_format(2^.x), expand = c(0.05, 0)) +
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = '#d62728', linewidth = 0.5) +
    labs(title = paste(run_name, "(Transcript-level)"), subtitle = paste("Pearson:", round(pearson_tx, 4), "Precision:", round(tx_precision, 4), "Recall:", round(tx_recall, 4), "F1:", round(tx_f1, 4)), x = "True Counts", y = "MiniQuT3 Counts") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'))
  
  ggsave(file.path(results_dir, paste0(run_name, "_corr_tx.png")), p_tx)

  # Store results - transcript-level metrics
  stats$pearson_tx <- pearson_tx
  stats$precision_tx <- tx_precision
  stats$recall_tx <- tx_recall
  stats$f1_tx <- tx_f1
  
  # Store results - gene-level metrics
  stats$pearson_gene <- pearson_gene
  stats$precision_gene <- gene_precision
  stats$recall_gene <- gene_recall
  stats$f1_gene <- gene_f1
  
  # Store results - differential transcript recall
  stats$diff_tx_recall <- diff_tx_recall
  stats$diff_tx_num <- diff_tx_num
  
  list(run_name = run_name, stats = as.data.frame(stats))
}

# Apply process_run in parallel using mclapply
results <- mclapply(run_dirs, process_run, sim_dir = sim_dir, results_dir = results_dir, 
                    tx_to_gene = tx_to_gene, mc.cores = detectCores())

# Combine results
all_metrics <- setNames(lapply(results, function(x) x$stats), 
                        lapply(results, function(x) x$run_name))

# Combine all results and save
final_table <- bind_rows(all_metrics)
write_csv(final_table, "tmp_benchmark_summary_final_all.csv")