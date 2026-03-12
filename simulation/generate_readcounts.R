#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(rtracklayer)
  library(dplyr)
})

# 1. Define Command Line Arguments
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL, 
              help = "Path to the input GTF file", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = "transcript_counts.tsv", 
              help = "Path to the output TSV file [default= %default]", metavar = "character"),
  make_option(c("--minCount"), type = "integer", default = 0, 
              help = "Minimum count value [default= %default]", metavar = "integer"),
  make_option(c("--maxCount"), type = "integer", default = 100000, 
              help = "Maximum count value [default= %default]", metavar = "integer"),
  make_option(c("--numTranscripts"), type = "integer", default = NULL, 
              help = "Number of transcripts to randomly sample [default= all]", metavar = "integer"),
  make_option(c("--biotype"), type = "character", default = NULL, 
              help = "Filter for biotype(s), comma-separated", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input GTF file is required.", call. = FALSE)
}

# 2. Import GTF File
message("Reading GTF file...")
gtf_data <- import(opt$input)

df <- as.data.frame(gtf_data) %>%
  filter(type == "transcript") %>%
  select(any_of(c("gene_id", "transcript_id", "gene_biotype", "transcript_biotype", "gene_type")))

# 3. Apply Biotype Filtering
if (!is.null(opt$biotype)) {
  biotypes_to_keep <- unlist(strsplit(opt$biotype, ","))
  biotype_col <- intersect(names(df), c("transcript_biotype", "gene_biotype", "gene_type"))[1]
  
  if (is.na(biotype_col)) {
    stop("Could not find a biotype column in the GTF.")
  }
  df <- df %>% filter(!!sym(biotype_col) %in% biotypes_to_keep)
}

# 4. Handle Sampling
if (!is.null(opt$numTranscripts)) {
  num_available <- nrow(df)
  num_to_sample <- min(opt$numTranscripts, num_available)
  df <- df %>% sample_n(num_to_sample)
}

# 5. Generate Skewed Random Counts
# We use a Gamma distribution (shape=2) to get that "low-but-not-zero" peak.
# The 'scale' determines where that peak sits. 
# Setting scale to 1/10th of the range puts the peak at ~10% of the max.
message("Generating right-skewed counts...")
set.seed(42)
n_rows <- nrow(df)
range_val <- opt$maxCount - opt$minCount

# rgamma(n, shape, scale)
# Mode of Gamma is (shape - 1) * scale. With shape=2, mode = scale.
# We set scale to 7% of the range to get that 500-1000 peak for a 10k range.
target_scale <- range_val * 0.07

raw_counts <- rgamma(n_rows, shape = 1.3, scale = target_scale)
final_counts <- round(raw_counts + opt$minCount)

# Ensure we don't exceed maxCount due to the long tail
final_counts <- pmin(final_counts, opt$maxCount)

final_df <- df %>%
  transmute(
    gene = gene_id,
    transcript = transcript_id,
    count = final_counts
  )

# 6. Write to TSV
write.table(final_df, file = opt$output, sep = "\t", quote = FALSE, row.names = FALSE)
message(paste("Successfully wrote", nrow(final_df), "transcripts to:", opt$output))