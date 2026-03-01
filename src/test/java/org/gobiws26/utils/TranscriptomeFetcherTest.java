package org.gobiws26.utils;


import org.junit.jupiter.api.Test;


import static org.gobiws26.utils.Test.*;
import static org.junit.jupiter.api.Assertions.*;

class TranscriptomeFetcherTest {

    @Test
    void positiveStrandIsFetchedCorrectly(){
        String s = TranscriptomeFetcher.getStringOf(
                tf.fetchTranscriptSequenceOf(transcripts.get("ENSSSCT00000081764"))
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

}