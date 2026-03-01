#!/usr/bin/Rscript

library(dplyr)
library(readr)
library(ggplot2)

x <- read_tsv("read2_36_10_phred_counts.tsv", col_names=c("char", "occurrence"))
y <-  x %>% mutate(score = sapply(char, utf8ToInt) - 33)

p <- ggplot(y, aes(x = score, weight = occurrence)) + geom_histogram(binwidth = 1)
ggsave(filename="phred_distro_read2_36_10.png", p)

