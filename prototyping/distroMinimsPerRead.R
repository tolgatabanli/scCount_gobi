#!/usr/bin/Rscript

#library(readr)
#library(ggplot2)

#x <- read_tsv("minimizerCountsPerRead.txt", col_names = FALSE)
#p <- ggplot(x, aes(x = x)) + geom_histogram()

#ggsave(filename = "distroMinimsCountsPerRead.png", plot = p)

nums <- scan("minimizerCountsPerRead.txt")
png("distroMinimsCountsPerRead.png", width = 600, height = 600)
hist(nums, breaks = 30)
dev.off()
