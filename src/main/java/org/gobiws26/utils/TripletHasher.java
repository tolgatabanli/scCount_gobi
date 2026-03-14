package org.gobiws26.utils;

import java.util.Objects;

/**
 * Hashes triplets of minimizers (shorts) to integers for use in triplet-based de Bruijn graph.
 * Collisions are tolerable and will be handled by the graph structure.
 */
public class TripletHasher {
    
    /**
     * Hash three minimizers into a single integer.
     * Uses Objects.hash() for a simple, collision-tolerant hash.
     */
    public static int hash(short m1, short m2, short m3) {
        return Objects.hash(m1, m2, m3);
    }
}
