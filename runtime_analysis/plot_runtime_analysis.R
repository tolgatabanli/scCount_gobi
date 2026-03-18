#!/usr/bin/Rscript

library(tidyverse)
library(ggplot2)

# Set output directory for plots
plot_dir <- "."

# Read data
reads_data <- read_tsv("reads_runtime_scaling.tsv")
count_data <- read_tsv("runtime_count.tsv")
index_data <- read_tsv("runtime_index.tsv")

# Plot 1: Counting runtime vs number of reads (colored by k, faceted by g)
# Convert reads to numeric (remove 'M' suffix)
reads_data_plot <- reads_data %>%
  mutate(reads_num = as.numeric(gsub("M", "", reads)) * 1e6,
         reads_label = reads) %>%
  group_by(reads_num, reads_label, k, g) %>%
  summarise(mean_time = mean(time),
            sd_time = sd(time),
            .groups = "drop") %>%
  arrange(reads_num)

# Plot with lines colored by k, faceted by g
p1 <- ggplot(reads_data_plot, aes(x = reads_num, y = mean_time, color = factor(k), group = factor(k))) +
  geom_point(size = 7) +
  geom_line(linewidth = 2.5) +
  geom_errorbar(aes(ymin = mean_time - sd_time, ymax = mean_time + sd_time),
                width = 3e5,
                linewidth = 1.2, alpha = 0.6) +
  facet_wrap(~g, labeller = labeller(g = c("5" = "g=5", "8" = "g=8", "11" = "g=11", "14" = "g=14")), scales = "free_y") +
  scale_x_continuous(labels = function(x) paste0(x/1e6, "M")) +
  scale_color_manual(values = c("#2E86AB", "#06A77D", "#A23B72", "#F18F01"),
                     name = "k-mer Size (k)") +
  labs(title = "Counting Runtime vs Number of Reads",
       x = "Number of Reads",
       y = "Wall-Clock Time (seconds)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 28, face = "bold", color = "#1a1a1a", margin = margin(b = 20)),
        axis.text = element_text(size = 26, color = "#1a1a1a"),
        axis.title = element_text(size = 30, color = "#1a1a1a", face = "bold", margin = margin(t = 10)),
        legend.title = element_text(size = 26, face = "bold", color = "#1a1a1a"),
        legend.text = element_text(size = 24, color = "#1a1a1a"),
        legend.key.size = unit(1.5, "cm"),
        strip.text = element_text(size = 32, face = "bold", color = "#1a1a1a", margin = margin(b = 10, t = 10)),
        panel.grid.major = element_line(color = "#e5e5e5", linewidth = 0.6),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 20, r = 20, b = 20, l = 20))

ggsave(file.path(plot_dir, "reads_scaling_time.png"), p1, width = 26, height = 10, dpi = 300)
print("Plot saved: reads_scaling_time.png")

combined_data <- bind_rows(
  #index_data %>% mutate(step = "Index"),
  count_data %>% mutate(step = "Count")
) %>% #mutate(step = factor(step, levels = c("Index", "Count"))) %>%
  mutate(k_g = paste0("k", k, "_g", g)) %>%
  group_by(k, g, step) %>%
  summarise(mean_time = mean(time),
            sd_time = sd(time),
            .groups = "drop")

p2 <- ggplot(combined_data, aes(x = factor(g), y = mean_time, fill = factor(k), group = factor(k))) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = mean_time - sd_time, ymax = mean_time + sd_time),
                position = position_dodge(width = 0.8),
                width = 0.2, linewidth = 0.6, alpha = 0.7) +
  #facet_wrap(~step, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("#2166ab", "#1b9e77", "#d8b465", "#d7140e"), 
                    name = "k-mer Size (k)") +
  labs(title = "Runtime Analysis: k-mer vs minimizer Length",
       x = "Minimizer Length (g)",
       y = "Wall-Clock Time (seconds)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 13, face = "bold", color = "#1a1a1a"),
        axis.text = element_text(size = 10, color = "#1a1a1a"),
        axis.title = element_text(size = 11, color = "#1a1a1a"),
        legend.title = element_text(size = 10, color = "#1a1a1a"),
        legend.text = element_text(size = 9, color = "#1a1a1a"),
        strip.text = element_text(size = 11, face = "bold", color = "#1a1a1a"),
        panel.grid.major = element_line(color = "#e5e5e5", linewidth = 0.3),
        panel.grid.minor = element_blank())

#ggsave(file.path(plot_dir, "kmer_minimizer_analysis.png"), p2, width = 12, height = 6, dpi = 300)
ggsave(file.path(plot_dir, "kmer_minimizer_analysis_count.png"), p2, width = 5.5, height = 6, dpi = 300)
