library(rtracklayer)
library(tidyverse)

gtf <- "/mnt/raidbio2/extstud/praktikum/genprakt-ws25/gruppe_a/data/Sus_scrofa.Sscrofa11.1.gtf"
gtf <- import(gtf) %>% as.data.frame()

gtf %>% filter(type == "exon" | type == "three_prime_utr") %>%
    group_by(transcript_id) %>%
    mutate(utr_length = sum(if_else(type == "three_prime_utr", width, 0)),
            tx_length = sum(width)) %>%
    select(transcript_id, utr_length, tx_length) %>%
    mutate(tail_length = utr_length + 500) %>%
    mutate(tail_length = if_else(tail_length > tx_length, tx_length, tail_length)) %>%
    distinct() %>%
    ungroup() %>%
    write_csv("tx_tail_lengths.csv")
