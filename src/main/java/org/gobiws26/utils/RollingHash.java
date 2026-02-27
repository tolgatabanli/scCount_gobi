package org.gobiws26.utils;

// TODO: ntHash
public class RollingHash {
    private static final long hA = 0x3c8bfbb395c60474L;
    private static final long hC = 0x3193c18562a02b4cL;
    private static final long hG = 0x20323ed082572324L;
    private static final long hT = 0x295549f54be24456L;
    private static final long hN = 0L;

    private static final long[] HASH = new long[256];
    private static final long[] HASH_REVCOMP = new long[256];

    static {
        HASH['A'] = hA;
        HASH['C'] = hC;
        HASH['G'] = hG;
        HASH['T'] = hT;

        HASH_REVCOMP['A'] = hT;
        HASH_REVCOMP['C'] = hG;
        HASH_REVCOMP['G'] = hC;
        HASH_REVCOMP['T'] = hA;
    }

    public static class HashState {
        public long fh; // forward strand
        public long rh; // reverse strand

        public HashState(long fh, long rh) {
            this.fh = fh;
            this.rh = rh;
        }

        // Returns h(s) = min(f(s), r(s)) using unsigned comparison
        public long getCanonical() {
            return Long.compareUnsigned(fh, rh) < 0 ? fh : rh; // avoid signed comparison
        }
    }

    // first kmer
    public static HashState computeFirst(byte[] bases, int start, int k) {
        long fh = 0;
        long rh = 0;
        for (int i = 0; i < k; i++) { // O(k) time
            byte b = bases[start + i];
            fh ^= Long.rotateLeft(HASH[b], k - 1 - i);
            rh ^= Long.rotateLeft(HASH_REVCOMP[b], i);
        }
        return new HashState(fh, rh);
    }


    // leftBase is leaving the window, rightBase is incoming
    public static void roll(HashState state, byte leftBase, byte rightBase, int k) {
        state.fh = Long.rotateLeft(state.fh, 1)
                ^ Long.rotateLeft(HASH[leftBase], k)
                ^ HASH[rightBase];

        state.rh = Long.rotateRight(state.rh, 1)
                ^ Long.rotateRight(HASH_REVCOMP[leftBase], 1)
                ^ Long.rotateLeft(HASH_REVCOMP[rightBase], k - 1);
    }
}
