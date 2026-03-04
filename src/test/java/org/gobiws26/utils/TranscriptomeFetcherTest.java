package org.gobiws26.utils;


import org.junit.jupiter.api.Test;


import static org.gobiws26.Test.*;
import static org.junit.jupiter.api.Assertions.*;

class TranscriptomeFetcherTest {

    @Test
    void positiveStrandIsFetchedCorrectly(){
        String s = TranscriptomeFetcher.getStringOf(
                tf.fetchTranscriptSequenceOf(transcripts.get("ENSSSCT00000081764")) // doesnt have a UTR, since pseudogene
        );
        assertEquals("CTGGAGGAGGACGCATGCAAGGGAGAAAAGAATATTTCGCTGCCTGTTTGGAAATCTGAATTGGTTATTGAATATTATGTGTCTCTGATGCACAGGTCCTCCCCCACTTATGAATGTTGA",
                s);
    }

    @Test
    void negativeStrandIsFetchedCorrectly(){
        String s = TranscriptomeFetcher.getStringOf(
                tf.fetchTranscriptSequenceOf(transcripts.get("ENSSSCT00000103324"))
        );
        assertEquals("ACGGAGACGTCCGCTCGCACGGTGGGAACCCGTGCTTCCTGCGTCTGCAGCTCAGGGTGCTGTGCGGGGGCGGGGGGTGACCCCGAGATGGCGGCGCCCTGGCACGTGATGGTGCGTCCTTGTCTGTTTTGTCTGTGTCCCTGGACCGAGTTCCAGTTCCCACTCTGCCACTTCCCGGCTGTGCCCTGGGGCAGATGT",
                s);
    }

    @Test
    void tailFetchWorks() {
        String s = TranscriptomeFetcher.getStringOf(
                tf.fetchTranscriptSequenceOf(transcripts.get("ENSSSCT00000081764"), 10)
        );
        assertEquals("TGAATGTTGA", s);
    }

    // x %>% filter(type %in% c("exon", "three_prime_utr"), gene_biotype=="protein_coding") %>% group_by(transcript_id) %>% mutate(tx_width = sum(width), num_utr = sum(type == "three_prime_utr")) %>% filter(tx_width < 500, num_utr > 0) %>% arrange(tx_width) %>% select(transcript_id, tx_width, num_utr)

    @Test
    void tailUTRFetchWorks() {
        String s = TranscriptomeFetcher.getStringOf(
                tf.fetchTranscriptSequenceOfUTRPlus(transcripts.get("ENSSSCT00000100887"), 10)
        );
        assertEquals("TATACTCTAACAAACAGTGGAAGAGACT", s);
    }

}