package org.gobiws26.utils;

import java.util.NoSuchElementException;

/**
 * Takes base sequence as byte array as output from samtools and creates an iterator
 * which returns a long representing kmers of the sequence on each next() call.
 * Note: samtools return a byte array of nucleotides (chars: A:65, C:67, G:71, T:84)
 */
public class KmerIterator {
    private final int K;
    private final byte[] fullSeq;
    private int windowStart = 0;

    private long last = 0;
    private int lastHash = 0;

    private RollingHash rollingHash = new RollingHash();

    private static byte[] HASH = new byte[256];
    static {
        HASH['A'] = 0;
        HASH['C'] = 1;
        HASH['G'] = 2;
        HASH['T'] = 3;
    }

    /**
     *
     * @param bases A sequence.
     * @return a long, representing a kmer of length <= 32. First nucleotide is the biggest bit and
     * the kmer is placed in the first K least significant bits, in other words, right-aligned when one writes
     * the binary number where the smallest bit is on the right hand side.
     */
    private long convertFirst(byte[] bases) {
        long kmer = 0L; // 64 bits
        for (int i = 0; i < K; i++) {
            kmer <<= 2;
            kmer |= (HASH[bases[i]] & 3); // &3 for masking (11)
        }
        return kmer;
    }

    /**
     *
     * @param sequence A sequence.
     * @param k kmer size. Could also use Config.K. Might be removed.
     */
    public KmerIterator(byte[] sequence, int k) {
        if (k > sequence.length) throw new IllegalArgumentException("k cannot be larger than sequence length");
        if (k <= 0) throw new IllegalArgumentException("k must be positive");

        K = k;
        fullSeq = sequence;
        last = convertFirst(fullSeq);
        windowStart++;
    }

    /**
     *
     * @return if there is a kmer that can be returned.
     */
    public boolean hasNext() {
        return windowStart + K <= fullSeq.length;
    }

    /**
     *
     * @return a long representing the next kmer of the sequence.
     */
    public long next() {
        if (!hasNext()) throw new NoSuchElementException();
        long current = current << 2;
        current |= HASH[fullSeq[windowStart + K - 1]];

        if (K < 32) { // set left bits to zero
            current &= (1L << (2 * K)) - 1;
        }

        windowStart++;
        return current;
    }
}
