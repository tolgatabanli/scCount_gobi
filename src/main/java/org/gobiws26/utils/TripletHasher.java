package org.gobiws26.utils;

/**
 * Hashes triplets of minimizers (ints) using xxHash32.
 * -> deterministic, high-quality, fast, non-cryptographic hash algorithm
 * -> consistent results across sessions and versions -> persistent index files.
 */
public class TripletHasher {
    
    // xxHash32 constants
    private static final int PRIME32_2 = 0x85ebca77;
    private static final int PRIME32_3 = 0xc2b2ae3d;
    private static final int PRIME32_4 = 0x27d4eb2d;
    private static final int PRIME32_5 = 0x165667b1;
    

    public static int hash(int m1, int m2, int m3) {
        int h32 = PRIME32_5;
        h32 += Integer.BYTES * 3; // 12 bytes total (3 ints × 4 bytes)
        
        // Update with each minimizer
        h32 += m1 * PRIME32_3;
        h32 = Integer.rotateLeft(h32, 17) * PRIME32_4;
        
        h32 += m2 * PRIME32_3;
        h32 = Integer.rotateLeft(h32, 17) * PRIME32_4;
        
        h32 += m3 * PRIME32_3;
        h32 = Integer.rotateLeft(h32, 17) * PRIME32_4;
        
        // Finalization
        h32 ^= h32 >>> 15;
        h32 *= PRIME32_2;
        h32 ^= h32 >>> 13;
        
        return h32;
    }
}
