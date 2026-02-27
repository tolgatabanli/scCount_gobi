cutoff <- 100

x %>%
        filter(type == "exon") %>%
        group_by(seqnames, transcript_id) %>%
        arrange(if_else(strand == '+', -start, start)) %>% # order exons according to their tx order, last exon first
        mutate(cumulative_width = cumsum(width),
                prev_cum_width = lag(cumulative_width, default = 0)) %>%
        filter(prev_cum_width < cutoff) %>%
        mutate(
                # How many bases of this exon do we actually need?
                needed = pmin(width, cutoff - prev_cum_width),
                # Adjust start or end based on strand to "trim" the 5' side of the exon
                start = case_when(
                    strand == "+" ~ end - needed + 1,
                    strand == "-" ~ start,
                    TRUE ~ start
                ),
                end = case_when(
                    strand == "+" ~ end,
                    strand == "-" ~ start + needed - 1,
                    TRUE ~ end
                )
            ) %>%
    arrange(start) %>%
    summarize(
        gene_id = first(gene_id),
        tx_signature = paste(start, end, sep = "-", collapse = ",") # represents exact genomic signature of the last X bp
    ) %>%
    group_by(tx_signature) %>%
    mutate(is_unique = n() == 1) %>%
    group_by(gene_id) %>%
    summarize(
        total_transcripts = n(),
        num_uni_identified_tx = sum(is_unique),
        uni_ratio = num_uni_identified_tx / total_transcripts
    ) %>%
    filter(total_transcripts > 1) -> z

hist(z$uni_ratio)