#!/usr/bin/Rscript

library(tidyverse)
library(ggplot2)

# Get list of .idx files from the idx/ directory
idx_dir <- "idx/"
idx_files <- list.files(idx_dir, pattern = "\\.idx$", full.names = TRUE)

# Extract file information and parse k/g parameters
idx_data <- tibble(
  filepath = idx_files,
  filename = basename(idx_files)
) %>%
  mutate(
    # Extract k and g from filename like "sccount_k15_g12.idx"
    k = as.numeric(str_extract(filename, "(?<=k)\\d+")),
    g = as.numeric(str_extract(filename, "(?<=g)\\d+"))
  ) %>%
  mutate(
    # Get file size in megabytes
    size_mb = file.size(filepath) / (1024^2)
  ) %>%
  # Convert k and g to factors with proper ordering
  mutate(
    k = factor(k, levels = sort(unique(k))),
    g = factor(g, levels = sort(unique(g))),
    text_color = if_else(
      size_mb > quantile(size_mb, 0.85, na.rm = TRUE),
      "black",
      "white"
    )
  )

# Create heatmap using ggplot2
p <- ggplot(idx_data, aes(x = g, y = k, fill = size_mb)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(size_mb, 1), color = text_color), size = 3) +
  scale_color_identity() +
  scale_fill_viridis_c(name = "Size (MB)", direction = 1) +
  labs(
    title = ".idx File Sizes Comparison",
    x = "g (gap size)",
    y = "k (kmer size)",
    subtitle = "File sizes in megabytes"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave("idx_size_comparison_heatmap.png", plot = p, width = 8, height = 6, dpi = 150)

# Print summary statistics
cat("\n=== Summary Statistics ===\n")
print(idx_data %>%
  select(filename, k, g, size_mb) %>%
  arrange(as.numeric(k), as.numeric(g)))

cat("\n")
cat("Smallest file:", min(idx_data$size_mb), "MB\n")
cat("Largest file:", max(idx_data$size_mb), "MB\n")
cat("Mean file size:", mean(idx_data$size_mb), "MB\n")
