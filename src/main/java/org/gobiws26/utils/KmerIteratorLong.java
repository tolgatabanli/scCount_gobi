package org.gobiws26.utils;

import java.util.Iterator;
import java.util.NoSuchElementException;

public class KmerIteratorLong implements Iterator<Long> {
    private final int K;
    private final byte[] fullSeq;
    private int windowStart = 0;

    // Holds the integer representation of the current k-mer window
    private long currentKmer = 0L;

    /**
     * @param sequence A sequence of bytes (e.g., ASCII 'A', 'C', 'G', 'T').
     * @param k kmer size. Max allowed is 32.
     */
    public KmerIteratorLong(byte[] sequence, int k) {
        if (k <= 0) throw new IllegalArgumentException("k must be positive");
        if (k > 32) throw new IllegalArgumentException("k cannot be larger than 32 (limit of 64-bit long)");
        if (k > sequence.length) throw new IllegalArgumentException("k cannot be larger than sequence length");

        this.K = k;
        this.fullSeq = sequence;

        // Initialize the first k-mer
        for (int i = 0; i < K; i++) {
            long val = encode(fullSeq[i]);
            // Shift left by 2*i because the first nucleotide goes to the right/LSD
            currentKmer |= (val << (2 * i));
        }
    }

    @Override
    public boolean hasNext() {
        return windowStart + K <= fullSeq.length;
    }

    @Override
    public Long next() {
        if (!hasNext()) throw new NoSuchElementException();

        // If this isn't the first k-mer, we slide the window to the right
        if (windowStart > 0) {
            // Drop the oldest nucleotide by shifting logical right by 2
            currentKmer >>>= 2;

            // Get the newest nucleotide coming into the window
            long newVal = encode(fullSeq[windowStart + K - 1]);

            // Place the new nucleotide at the highest valid bit position
            currentKmer |= (newVal << ((K - 1) * 2));
        }

        windowStart++;
        return currentKmer;
    }

    /**
     * Maps ASCII nucleotide bytes to their 2-bit representations.
     */
    private long encode(byte nuc) {
        switch (nuc) {
            case 'A': case 'a': return 0L; // 00
            case 'C': case 'c': return 1L; // 01
            case 'G': case 'g': return 2L; // 10
            case 'T': case 't': return 3L; // 11
            default:
                throw new IllegalArgumentException("Invalid nucleotide in sequence: " + (char) nuc);
        }
    }
}
