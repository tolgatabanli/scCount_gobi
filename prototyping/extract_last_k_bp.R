library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(readr)

gtf <- import('/mnt/raidbio2/extdata/praktikum/genprakt/genprakt-ws25/Block/pig-genome/Sus_scrofa.Sscrofa11.1.115.chr.gtf.gz')
gtf <- as.data.frame(gtf)

cutoff = 100

for (cutoff in seq(100, 1000, 100)) {
    gtf %>%
        filter(type == "exon") %>%
        group_by(seqnames, transcript_id) %>%
        arrange(if_else(strand == '+', start, -start)) %>% # order exons according to their tx order
        mutate(cumulative_width = cumsum(width)) %>%
        filter(cumulative_width <= cutoff) %>%
        summarize(strand = strand,
                    start = if_else(strand == '+', min(start), min(end)),
                    end = if_else(strand == '+', max(end), max(start))) %>%
        mutate(score = ".") %>%
        select(seqnames, start, end, transcript_id, score, strand) %>%
        write_tsv(file.path("t_bed", paste0("t_coverage_last", cutoff, ".bed")), col_names = FALSE)
}



