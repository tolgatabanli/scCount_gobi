library(GenomicRanges)
library(rtracklayer)

gtf <- import('data/Sus_scrofa.Sscrofa11.1.gtf')

exons <- gtf[gtf$type == 'exon']
exons_by_tx <- split(exons, exons$transcript_id)

get_last_k_spliced <- function(gr, k) {
    if (as.character(strand(gr)[1]) == '+') {
        gr <- gr[order(start(gr))]
    } else {
        gr <- gr[order(start(gr), decreasing = TRUE)]
    }

    widths <- width(gr)
    cum_widths <- cumsum(widths)
    keep_exons <- which(cum_widths > (sum(widths) - k))
    selected <- gr[keep_exons]
    excess <- sum(widths) - k 
    if (excess > 0) {
        trim_amount <- excess - sum(widths[seq_len(min(keep_exons)-1)])
        if(as.character(strand(selected[1]))== '+') {
            start(selected[1]) <- start(selected[1]) + trim_amount
        } else {
            end(selected[1]) <- end(selected[1]) - trim_amount
        }
    }
    selected
}

k <- 400

ks <- seq(100, 2000, by = 100)
last_k_list <- lapply(exons_by_tx,
                      get_last_k_spliced,
                      k = k)

last_k_grl <- GRangesList(last_k_list)
last_k_regions <- unlist(last_k_grl)

mcols(last_k_regions)$name <- rep(
    names(last_k_list),
    lengths(last_k_list)
)
export(last_k_regions,"data/last_400_bp.bed", format = "bed")

