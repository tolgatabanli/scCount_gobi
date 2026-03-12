#!/usr/bin/Rscript
library(readr)
library(dplyr)
library(ggplot2)

truth <- read_tsv('/home/t/tabanli/Desktop/scCount/simulation/random_counts.tsv')
genes_truth <- truth %>% group_by(gene) %>% summarise(count = sum(count))
tx_truth <- truth %>% select(transcript, count)

obs_genes <- read_tsv('/mnt/cip/home/t/tabanli/Desktop/scCount/out/count_out_simul_noMut_gamma/counts.tsv', col_names = c("id", "count"))

p <- genes_truth %>% inner_join(obs_genes, by = c("gene" = "id"), suffix = c("_truth", "_obs")) %>%
    ggplot(aes(x = count_truth, y = count_obs)) +
    geom_point()
ggsave(filename = "simul_noMut_gamma_genes.png", p)

p <- tx_truth %>% inner_join(obs_genes, by = c("transcript" = "id"), suffix = c("_truth", "_obs")) %>%
    ggplot(aes(x = count_truth, y = count_obs)) +
    geom_point()
ggsave(filename = "simul_noMut_gamma_transcripts.png", p)
