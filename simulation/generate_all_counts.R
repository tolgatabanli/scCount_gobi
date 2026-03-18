#!/usr/bin/env Rscript

library(rtracklayer)
library(dplyr)

# Configuration
GTF <- "/mnt/raidbio2/extstud/praktikum/genprakt-ws25/gruppe_a/data/Sus_scrofa.Sscrofa11.1.gtf"
NUM_TRANSCRIPTS <- 10000
MIN_COUNT <- 1
MAX_COUNT <- 10000
BIOTYPES <- c("all", "protein_coding", "lncRNA", "miRNA", "rRNA")

set.seed(42)

message("Reading GTF file...")
gtf_data <- import(GTF)

df_all <- as.data.frame(gtf_data) %>%
  filter(type == "transcript") %>%
  select(any_of(c("gene_id", "transcript_id", "gene_biotype", "transcript_biotype", "gene_type")))

# Determine biotype column
biotype_col <- intersect(names(df_all), c("transcript_biotype", "gene_biotype", "gene_type"))[1]
if (is.na(biotype_col)) {
  stop("Could not find a biotype column in the GTF.")
}

# Generate counts for each biotype
for (biotype in BIOTYPES) {
  
  message("Processing: ", biotype)
  
  # Filter by biotype if not "all"
  if (biotype == "all") {
    df <- df_all
  } else {
    df <- df_all %>% filter(!!sym(biotype_col) == biotype)
  }
  
  # Sample transcripts
  num_available <- nrow(df)
  num_to_sample <- min(NUM_TRANSCRIPTS, num_available)
  df <- df %>% sample_n(num_to_sample)
  
  # Generate right-skewed random counts
  n_rows <- nrow(df)
  range_val <- MAX_COUNT - MIN_COUNT
  target_scale <- range_val * 0.07
  
  raw_counts <- rgamma(n_rows, shape = 1.3, scale = target_scale)
  final_counts <- round(raw_counts + MIN_COUNT)
  final_counts <- pmin(final_counts, MAX_COUNT)
  
  # Create output dataframe
  final_df <- df %>%
    transmute(
      gene = gene_id,
      transcript = transcript_id,
      count = final_counts
    )
  
  # Write to TSV
  output_file <- paste0("counts_", biotype, ".tsv")
  write.table(final_df, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  message(paste("  -> Wrote", nrow(final_df), "transcripts to:", output_file))
}

message("All count files generated successfully!")
