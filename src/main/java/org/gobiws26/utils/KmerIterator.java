package org.gobiws26.utils;

import java.util.Arrays;
import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 *
 */
public class KmerIterator {
    private final int K;
    private final byte[] fullSeq;
    private int windowStart = 0;

    // samtools return a byte array of nucleotides (chars: A:65, C:67, G:71, T:84)
    // trick is 1) one left shift 2) ANDing with 3 (11) to get
    private static long pack(byte[] bases, int start, int k) {
        long kmer = 0L; // 64 bits
        for (int i = start; i < start + k; i++) {
            kmer <<= 2;
            kmer |= ((bases[i] >> 1) & 3); // &3 for masking (11)
        }
        return kmer;
    }

    public KmerIterator(byte[] sequence, int k) {
        if (k > sequence.length) throw new IllegalArgumentException("k cannot be larger than sequence length");
        if (k <= 0) throw new IllegalArgumentException("k must be positive");

        K = k;
        fullSeq = sequence;
    }

    public boolean hasNext() {
        return windowStart + K <= fullSeq.length;
    }

    public long next() {
        if (!hasNext()) throw new NoSuchElementException();
        long kmer = pack(fullSeq, windowStart, K);
        windowStart++;
        return kmer;
    }
}
