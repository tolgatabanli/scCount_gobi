package org.gobiws26.utils;

// TODO: ntHash
public class RollingHash {
    private static final long hA = 0x3c8bfbb395c60474L;
    private static final long hC = 0x3193c18562a02b4cL;
    private static final long hG = 0x20323ed082572324L;
    private static final long hT = 0x295549f54be24456L;
    private static final long hN = 0L;

    private static final long[] HASH = new long[4];
    private static final long[] HASH_REVCOMP = new long[4];

    static {
        HASH[0] = hA;
        HASH[1] = hC;
        HASH[2] = hG;
        HASH[3] = hT;

        HASH_REVCOMP[0] = hT;
        HASH_REVCOMP[1] = hG;
        HASH_REVCOMP[2] = hC;
        HASH_REVCOMP[3] = hA;
    }

    public static class HashState {
        public long fh; // forward strand
        public long rh; // reverse strand

        public HashState(long fh, long rh) {
            this.fh = fh;
            this.rh = rh;
        }

        // Returns h(s) = min(f(s), r(s))
        public long getCanonical() {
            return Long.compareUnsigned(fh, rh) < 0 ? fh : rh; // avoid signed comparison
        }
    }

    // first kmer
    public static HashState computeFirst(long kmer, int start, int k) {
        long fh = 0;
        long rh = 0;
        for (int i = 0; i < k; i++) { // O(k) time
            byte b = kmer[start + i];
            fh ^= Long.rotateLeft(HASH[b], k - 1 - i);
            rh ^= Long.rotateLeft(HASH_REVCOMP[b], i);
        }
        return new HashState(fh, rh);
    }


    // leftBase is leaving the window, rightBase is incoming
    public static void roll(HashState state, byte leftBase, byte rightBase, int k) {
        state.fh = Long.rotateLeft(state.fh, 1)
                ^ Long.rotateLeft(HASH[leftBase & 0xFF], k)
                ^ HASH[rightBase & 0xFF];

        state.rh = Long.rotateRight(state.rh, 1)
                ^ Long.rotateRight(HASH_REVCOMP[leftBase & 0xFF], 1)
                ^ Long.rotateLeft(HASH_REVCOMP[rightBase & 0xFF], k - 1);
    }

    // INSTANCE PARTS
    private KmerIteratorByte ki;
    public RollingHash(KmerIteratorByte ki) {
        this.ki = ki;
    }

    public int getNextHash()


}
