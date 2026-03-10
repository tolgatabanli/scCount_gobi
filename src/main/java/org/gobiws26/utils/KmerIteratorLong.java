package org.gobiws26.utils;

import it.unimi.dsi.fastutil.longs.LongIterator;
import java.util.NoSuchElementException;

public class KmerIteratorLong implements LongIterator {
    private final int K;
    private final byte[] fullSeq;
    private final byte[] phredSeq; // Can be null for transcripts
    private final long kmerSizeMask;

    // The biggest two digits serve as quality flags
    private int windowStart = 0;
    private int nCounter = 0;
    private int currentQScore = 0; // Track rolling Q-score sum
    private long currentKmer = 0;

    public static class TooShortSequenceException extends IllegalArgumentException {
        public TooShortSequenceException(String message) {
            super(message);
        }
    }

    /**
     * Constructor for sequences WITH qualities.
     * @param sequence A sequence of bytes.
     * @param phredSequence Quality scores (can be null).
     * @param k kmer size. Max allowed is 31. The highest two bits are for quality flags.
     * WARNING: If K=31, the flag hits the sign bit. ALWAYS use Long.compareUnsigned() to sort!
     */
    public KmerIteratorLong(byte[] sequence, byte[] phredSequence, int k) {
        if (k <= 0) throw new IllegalArgumentException("k must be positive");
        if (k > 31) throw new IllegalArgumentException("k cannot be larger than 31 (limit of 64-bit long with flag)");
        if (k > sequence.length) throw new TooShortSequenceException("k cannot be larger than sequence length");
        if (phredSequence != null && sequence.length != phredSequence.length) {
            throw new IllegalArgumentException("Sequence and quality arrays must be the same length");
        }

        this.K = k;
        this.fullSeq = sequence;
        this.phredSeq = phredSequence;
        this.kmerSizeMask = (k == 31) ? -1L : (1L << (k * 2 + 2)) - 1;

        // Initialize the first k-mer and rolling counters
        for (int i = 0; i < K; i++) {
            byte nuc = encode(fullSeq[i]);
            if (fullSeq[i] == 78) nCounter++; // 'N' = 78

            if (phredSeq != null) {
                currentQScore += phredSeq[i]; // Add initial Q-scores
            }

            currentKmer <<= 2;
            currentKmer |= nuc;
        }
        applyQualityFlag();
    }

    private void applyQualityFlag() {
        if (nCounter > 0) {
            currentKmer &= ~(3L << (K * 2)); // keep N flag (00)
            return;
        } else {
            currentKmer |= (3L << (K * 2));  // reset N flag (11)
        }

        currentKmer &= kmerSizeMask;

        if (phredSeq != null && currentQScore < 19 * K) {
            currentKmer &= ~(3L << (K * 2)); // clear the two bits
            currentKmer |= (1L << (K * 2));  // set 01
        }
    }

    @Override
    public boolean hasNext() {
        return windowStart + K <= fullSeq.length;
    }

    @Override
    public long nextLong() {
        if (!hasNext()) throw new NoSuchElementException();

        if (windowStart == 0) {
            windowStart++;
            return currentKmer;
        }

        // See if an N is incoming at the right edge, and add its Q-score
        byte nextNuc = fullSeq[windowStart + K - 1];
        if (nextNuc == 78) nCounter++;

        if (phredSeq != null) {
            currentQScore += phredSeq[windowStart + K - 1] - 33;
        }

        // Shift kmer, place the new nucleotide at the LSD
        currentKmer <<= 2;
        currentKmer |= encode(nextNuc);
        currentKmer &= kmerSizeMask;

        // See if an N is leaving at the left edge, and subtract its Q-score
        if (fullSeq[windowStart - 1] == 78) nCounter--;

        if (phredSeq != null) {
            currentQScore -= phredSeq[windowStart - 1] - 33;
        }

        applyQualityFlag();
        windowStart++;
        return currentKmer;
    }

    private byte encode(byte nuc) {
        byte res = nuc;
        res >>= 1;
        res &= 3;
        return res;
    }
}