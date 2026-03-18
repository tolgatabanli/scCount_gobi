package org.gobiws26.utils;

import java.lang.Long;
import org.junit.jupiter.api.Test;
import static org.gobiws26.Test.tf;
import static org.gobiws26.Test.transcripts;

import static org.junit.jupiter.api.Assertions.*;


public class KmerIteratorTests {
    byte[] s = tf.fetchTranscriptSequenceOf(transcripts.get("ENSSSCT00000081764"));
    KmerIteratorLong kiLong = new KmerIteratorLong(s, null, 20);
    long firstKmer = kiLong.nextLong(); // "CTGGAGGAGGACGCATGCAA"


    @Test
    public void LongIteratorWorks() {
        assertEquals(Long.parseUnsignedLong("110110111100111100111100011101001011010000", 2), firstKmer);
    }

    @Test
    public void MinimizerWorks() {
        int minimizer = Minimizers.getMinimizer(firstKmer); // expected: ACGCATGC -> 0001110100101101
        assertEquals(Integer.parseInt("0001110100101101", 2), minimizer);
    }
}
