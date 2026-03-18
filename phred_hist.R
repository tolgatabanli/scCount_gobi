#!/usr/bin/env Rscript

# Read the k parameter from command line arguments (default to 4 if not provided)
args <- commandArgs(trailingOnly = TRUE)
k <- if (length(args) > 0) as.integer(args[1]) else 4

# Function to calculate kmer sums for a single sequence
get_kmer_sums <- function(phred_seq, k) {
  scores <- utf8ToInt(phred_seq) - 33
  n <- length(scores)
  
  # Return empty if the read is shorter than k
  if (k > n) return(numeric(0)) 
  
  cs <- cumsum(scores)
  cs[k:n] - c(0, cs[1:(n - k)])
}

# Read all data piped from standard input
# warning=FALSE prevents warnings if the pipe doesn't end with a newline
lines <- readLines(file("stdin"), warn = FALSE)

if (length(lines) < 4) {
  stop("Not enough lines provided. FASTQ requires at least 4 lines per read.")
}

# Extract the Phred quality strings (every 4th line)
idx <- seq(4, length(lines), by = 4)
qual_lines <- lines[idx]

# Apply the function to all quality lines and flatten into one numeric vector
cat(sprintf("Processing %d reads for %d-mers...\n", length(qual_lines), k), file = stderr())
all_sums <- unlist(lapply(qual_lines, get_kmer_sums, k = k))

if (length(all_sums) == 0) {
  stop("No valid k-mers found. Check your input data and k parameter.")
}

# Generate and save the plot
output_file <- paste0("phred_", k, "mer_hist.png")
png(output_file, width = 800, height = 600, res = 120)

hist(all_sums, 
     main = paste("Histogram of", k, "-mer Phred Score Sums"),
     xlab = paste("Sum of Phred Scores (k =", k, ")"),
     ylab = "Frequency",
     col = "steelblue",
     border = "white",
     breaks = 50) 

invisible(dev.off()) # Suppress the "null device" print

cat(sprintf("Success! Histogram saved to: %s\n", output_file), file = stderr())
