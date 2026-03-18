#!/usr/bin/Rscript

library(tidyverse)

txlengths <- read_csv("tx_tail_lengths.csv")
minimlengths <- read_tsv("transcript_minimizer_counts.tsv")

p <- txlengths %>%
    inner_join(minimlengths, by = "transcript_id") %>%
    pivot_longer(cols = ends_with("_length"), names_to = "variable", values_to = "length") %>%
    ggplot(aes(x = length, y = minimizer_count), alpha = 0.5) +
    geom_point() +
    facet_wrap(~ variable, scales = "free_x") +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = "Tail Length (bp)", y = "Minimizer Count", title = "Minimizer Count vs Tail Length") +
    theme_minimal()

ggsave("minimizer_count_vs_tail_length.png", plot = p, width = 8, height = 6)
