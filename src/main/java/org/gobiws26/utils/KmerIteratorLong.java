package org.gobiws26.utils;

import it.unimi.dsi.fastutil.longs.LongIterator;
import org.gobiws26.Config;

import java.util.NoSuchElementException;

public class KmerIteratorLong implements LongIterator {
    private final int K;
    private final byte[] fullSeq;
    private byte[] phredSeq;
    private final long kmerSizeMask;

    // The biggest two digits serve as quality flags
    private int windowStart = 0;
    private int nCounter = 0;
    private long currentKmer = 0;

    /**
     * @param sequence A sequence of bytes.
     * @param k kmer size. Max allowed is 31. The highest two bits are for quality scores (e.g. N).
     */
    public KmerIteratorLong(byte[] sequence, int k) {
        if (k <= 0) throw new IllegalArgumentException("k must be positive");
        if (k > 31) throw new IllegalArgumentException("k cannot be larger than 31 (limit of 64-bit long with flag)");
        if (k > sequence.length) throw new IllegalArgumentException("k cannot be larger than sequence length");

        this.K = k;
        this.fullSeq = sequence;

        this.kmerSizeMask = (k == 31) ? -1L : (1L << (k * 2 + 2)) - 1;

        // Initialize the first k-mer
        for (int i = 0; i < K; i++) {
            byte nuc = encode(fullSeq[i]);
            if (fullSeq[i] == 78) nCounter++; // 'N' = 78
            currentKmer <<= 2;
            currentKmer |= nuc; // slide in the new nucleotide
        }
        applyNFlag();
    }

    public KmerIteratorLong(byte[] sequence, byte[] phredSequence, int k) {
        this(sequence, k);
        this.phredSeq = phredSequence;
    }

    private void applyNFlag() {
        if (nCounter > 0) {
            currentKmer &= ~(3L << (K * 2)); // keep N flag (00)
        } else {
            currentKmer |= (3L << (K * 2));  // reset N flag (11)
        }
        currentKmer &= kmerSizeMask;
    }

    @Override
    public boolean hasNext() {
        return windowStart + K <= fullSeq.length;
    }

    /**
     * Returns the next kmer where the first nucleotide is in the biggest digit, allowing easy sorting.
     */
    @Override
    public long nextLong() {
        if (!hasNext()) throw new NoSuchElementException();

        if (windowStart == 0) {
            windowStart++;
            return currentKmer;
        }

        // See if an N is incoming at the right edge of the window
        byte nextNuc = fullSeq[windowStart + K - 1];
        if (nextNuc == 78) nCounter++;

        // Shift kmer, place the new nucleotide at the LSD
        currentKmer <<= 2;
        currentKmer |= encode(nextNuc);

        currentKmer &= kmerSizeMask;

        // See if an N is leaving at the left edge of the old window
        if (fullSeq[windowStart - 1] == 78) nCounter--;

        applyNFlag();
        windowStart++;
        return currentKmer;
    }

    private static short gmerMask = 0x3FFF; // take the rightmost 14 bits
    public static short getMinimizer(long kmer) {
        short biggest = gmerMask;
        byte qualityBits = (byte) (kmer >>> (Config.K * 2));
        for (int i = 0; i < Config.K - 7 + 1; i++) {
            short currentGmer = (short) (kmer & gmerMask);
            if (currentGmer < biggest) biggest = currentGmer;

            kmer >>>= 2;
        }
        biggest |= (short) (qualityBits << 14); // 'infect' the quality from kmer to gmer
        return biggest;
    }

    /**
     * Maps ASCII nucleotide bytes to their 2-bit representations.
     * Least modification unique encoding: A: 00, C: 01, T: 10, G: 11
     */
    private byte encode(byte nuc) {
        byte res = nuc;
        res >>= 1;
        res &= 3;
        return res;
    }
}