#!/usr/bin/Rscript
library(readr)
library(rtracklayer)
library(ggplot2)
library(parallel)
library(dplyr)

x <- as.data.frame(import("/mnt/raidbio2/extdata/praktikum/genprakt/genprakt-ws25/Block/pig-genome/Sus_scrofa.Sscrofa11.1.115.chr.gtf.gz"))

calc_uni_ratios <- function(df, cutoff) {
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
      gene_id = dplyr::first(gene_id),
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
    filter(total_transcripts > 1) %>%
    mutate(cutoff = cutoff) -> z
  return(z)
}

res <- mclapply(seq(100, 1000, 100), function(cut) calc_uni_ratios(x, cut),
                mc.cores = max(detectCores() - 1, 1))
res <- bind_rows(res)


# distribution of uniquely identifiable transcripts
p <- ggplot(res, aes(x = uni_ratio, group = cutoff, fill = cutoff)) +
  geom_histogram(bins = 10, position = "identity") +
  labs(caption = "Bins: 10, cutoffs: seq(100, 1000, 100)")
ggsave(filename = "hist_uniquely_id_tx.png", plot = p)

# bins to facets, cutoffs to x axis
res_binned <- res %>%
  mutate(ratio_bin = cut(uni_ratio,
                         breaks = seq(0, 1, length.out = 47),
                         include.lowest = TRUE,
                         labels = paste0(seq(0, 90, 2), "-", seq(10, 100, 2), "%"))) %>%
  # 2. Count how many genes fall into each bin for each cutoff
  group_by(cutoff, ratio_bin) %>%
  summarize(gene_count = n(), .groups = "drop")

# 3. Plot: Facet by the bin, show cutoff on the X-axis
p <- ggplot(res_binned, aes(x = cutoff, y = gene_count, group = ratio_bin)) +
  geom_line(color = "steelblue", size = 1) +
  geom_point() +
  facet_wrap(~ratio_bin, scales = "free_y") +
  labs(title = "Gene Count per Uni-Ratio Bin across Cutoffs",
       x = "Cutoff (bp)",
       y = "Number of Genes") +
  theme_minimal()
ggsave(filename = "facetted_hist_uniquely_id_tx.png", plot = p)