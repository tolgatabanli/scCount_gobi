package org.gobiws26.Indexing;

import it.unimi.dsi.fastutil.shorts.ShortArrayList;
import org.junit.jupiter.api.Test;

import static org.gobiws26.Test.*;
import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;

public class IndexerTests {
    @Test
    public void basicMinimizersFoundCorrectly() {
        String txId = "ENSSSCT00000081764"; // CTGGAGGAGGACGCATGCAAGGGAGAAAAGAATATTTCGCTGCCTGTTTGGAAATCTGAATTGGTTATTGAATATTATGTGTCTCTGATGCACAGGTCCTCCCCCACTTATGAATGTTGA
        String[] expectedMinimizersString = new String[]{"ACGCATG", "AAGGGAG", "AAAAGAA", "AAAGAAT", "AATATTT", "ATATTTC", "ATTTCGC", "CCTGTTT", "AAATCTG", "AATCTGA", "AATTGGT", "AATATTA", "ATATTAT", "ATTATGT", "ATGCACA", "ACAGGTC", "ACTTATG", "AATGTTG"};

        Indexer idxer = new Indexer(transcripts, fasta);
        idxer.runIndex();
        int txInId = idxer.getTranscriptToIndex().get(txId);

        ShortArrayList txMinimizers = idxer.getTranscriptToMinimizerPath().get(txInId);
        String[] observedMinimizersString = new String[txMinimizers.size()];
        for (int i = 0; i < observedMinimizersString.length; i++) {
            observedMinimizersString[i] = shortToNucleotideString(txMinimizers.getShort(i));
        }
        assertArrayEquals(expectedMinimizersString, observedMinimizersString);
    }
}
