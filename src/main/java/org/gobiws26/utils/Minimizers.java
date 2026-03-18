package org.gobiws26.utils;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import org.gobiws26.Config;

public class Minimizers {
    public static IntArrayList of(byte[] bases, byte[] phredSeq) {
        KmerIteratorLong kiLong = new KmerIteratorLong(bases, phredSeq, Config.K);
        IntArrayList minimList = new IntArrayList();

        // For each Kmer in the read, find the minimizer
        int currentMinim = getMinimizer(kiLong.nextLong());
        int lastFoundMinim = currentMinim;

        while (kiLong.hasNext()) {
            currentMinim = getMinimizer(kiLong.nextLong());
            if (currentMinim != lastFoundMinim) {
                minimList.add(lastFoundMinim);
                lastFoundMinim = currentMinim;
            }
        }

        minimList.add(lastFoundMinim); // flush
        return minimList;
    }

    /**
     * Extracts the minimizer of length Config.minimLength from a kmer.
     * The minimizer is the lexicographically smallest window of Config.minimLength bases.
     * Each position is compared as a 2-bit encoded value (right-aligned in the kmer).
     */
    public static int getMinimizer(long kmer) {
        int minimBits = Config.minimLength * 2; // Each base is 2 bits
        long minimMask = (1L << minimBits) - 1; // Mask for minimLength-sized window
        int windows = Config.K - Config.minimLength + 1;
        if (windows < 1) windows = 1;

        long smallest = minimMask; // Maximum possible value for minimLength bits
        for (int i = 0; i < windows; i++) {
            long currentMinim = kmer & minimMask;
            if (currentMinim < smallest) {
                smallest = currentMinim;
            }
            kmer >>>= 2;
        }

        return (int) smallest;
    }
}
