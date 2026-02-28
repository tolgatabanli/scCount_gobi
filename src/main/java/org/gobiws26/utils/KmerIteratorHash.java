package org.gobiws26.utils;

import java.util.NoSuchElementException;

/**
 * Takes base sequence as byte array as output from samtools and creates an iterator
 * which returns a 64-bit long representing the ntHash k-mer hash (now canonical)
 * of the sequence on each next() call.
 * TODO: do we need stranded, iow, forward-hash and reverse-hash?
 */
public class KmerIteratorHash {
    private final int K;
    private final byte[] fullSeq;
    private int windowStart = 0;

    private long currentFhVal = 0L;
    private long currentRhVal = 0L;

    // Seed tables for base nucleotides and their reverse complements
    private static final long[] HASH = new long[256];
    private static final long[] RC_HASH = new long[256];

    static {
        // Base hashes as defined in the ntHash specification
        HASH['A'] = 0x3c8bfbb395c60474L;
        HASH['C'] = 0x3193c18562a02b4cL;
        HASH['G'] = 0x20323ed082572324L;
        HASH['T'] = 0x295549f54be24456L;
        HASH['N'] = 0L;

        // Reverse complement mapping (~A = T, ~C = G, ~G = C, ~T = A)
        RC_HASH['A'] = HASH['T'];
        RC_HASH['C'] = HASH['G'];
        RC_HASH['G'] = HASH['C'];
        RC_HASH['T'] = HASH['A'];
        RC_HASH['N'] = 0L;
    }


    private static long rol(long x, int k) {
        // Java automatically applies bitwise AND 63 to the shift amount for longs
        return (x << k) | (x >>> (64 - k));
    }

    private static long ror(long x, int k) {
        return (x >>> k) | (x << (64 - k));
    }

    /**
     * @param sequence A sequence byte array.
     * @param k Window size. See Config.K
     */
    public KmerIteratorHash(byte[] sequence, int k) {
        if (k > sequence.length) throw new IllegalArgumentException("k cannot be larger than sequence length");
        if (k <= 0) throw new IllegalArgumentException("k must be positive");

        K = k;
        fullSeq = sequence;

        // Initialize the forward and reverse hashes for the very first k-mer
        for (int i = 0; i < K; i++) {
            int baseIndex = fullSeq[i] & 0xFF; // Prevent negative index from signed bytes
            currentFhVal ^= rol(HASH[baseIndex], K - 1 - i);
            currentRhVal ^= rol(RC_HASH[baseIndex], i);
        }
    }


    public boolean hasNext() {
        return windowStart + K <= fullSeq.length;
    }

    /**
     * @return a long representing the ntHash of k-mer.
     */
    public long next() {
        if (!hasNext()) throw new NoSuchElementException();

        long hash;

        if (windowStart == 0) {
            hash = Math.min(currentFhVal, currentRhVal); // first was calculated from constructor
        } else {
            int outBase = fullSeq[windowStart - 1] & 0xFF;
            int inBase = fullSeq[windowStart + K - 1] & 0xFF;

            // Forward rolling hash
            currentFhVal = rol(currentFhVal, 1) ^ rol(HASH[outBase], K) ^ HASH[inBase];

            // TODO: do we need this for example?
            currentRhVal = ror(currentRhVal, 1) ^ ror(RC_HASH[outBase], 1) ^ rol(RC_HASH[inBase], K - 1);

            // as per definition of canonical (just take min)
            hash = Math.min(currentFhVal, currentRhVal);
        }

        windowStart++;
        return hash;
    }
}