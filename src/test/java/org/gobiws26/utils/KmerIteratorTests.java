package org.gobiws26.utils;

import java.lang.Long;
import org.junit.jupiter.api.Test;
import static org.gobiws26.utils.Test.tf;
import static org.gobiws26.utils.Test.transcripts;

import static org.junit.jupiter.api.Assertions.*;


public class KmerIteratorTests {

    @Test
    public void LongIteratorWorks() {
        byte[] s = tf.fetchTranscriptSequenceOf(transcripts.get("ENSSSCT00000081764"));
        KmerIteratorLong kiLong = new KmerIteratorLong(s, 20);
        long firstKmer = kiLong.nextLong(); // "CTGGAGGAGGACGCATGCAA"
        // TODO
        assertEquals(Long.parseUnsignedLong("110110111100111100111100011101001011010000", 2), firstKmer);
    }
}
