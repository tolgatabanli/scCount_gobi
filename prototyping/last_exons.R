library(rtracklayer)
library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(GenomicFeatures)

gtf_path <- '/mnt/raidbio2/extdata/praktikum/genprakt/genprakt-ws25/Block/pig-genome/Sus_scrofa.Sscrofa11.1.115.chr.gtf.gz'
gtf <- import(gtf_path)

last_exons <- gtf %>% as.data.frame() %>%
    filter(type == "exon" | type == "three_prime_utr") %>%
    mutate(exon_number = as.integer(exon_number)) %>% # idk why casting does not work
    group_by(transcript_id, strand, type) %>%
    slice_max(exon_number) # with ties return also all 1, so there is no tie

p <- last_exons %>%
    filter(width < 70) %>%
    ggplot(aes(x = width, fill = type)) +
    geom_histogram(bins = 30)
ggsave(filename = "hist_ends.png", plot = p)

gtf %>% filter(type == "transcript") %>% 

transcripts <- makeTxDbFromGFF(gtf_path, format = 'gtf')
exons <- exonsBy(transcripts, by = "tx")
