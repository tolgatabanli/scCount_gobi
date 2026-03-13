package org.gobiws26.utils;

import it.unimi.dsi.fastutil.shorts.ShortArrayList;
import org.gobiws26.Config;

public class Minimizers {
    public static ShortArrayList of(byte[] bases, byte[] phredSeq) {
        KmerIteratorLong kiLong = new KmerIteratorLong(bases, phredSeq, Config.K);
        ShortArrayList minimList = new ShortArrayList();

        // For each Kmer in the read, find the minimizer
        short currentMinim = getMinimizer(kiLong.nextLong());
        short lastFoundMinim = currentMinim;

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

    private static final int gmerMask = 0xFFFF;
    public static short getMinimizer(long kmer) {
        int windows = Config.K - 8 + 1;
        if (windows < 1) windows = 1;

        int smallest = gmerMask;
        for (int i = 0; i < windows; i++) {
            int currentGmer = (int) (kmer & gmerMask);
            if (currentGmer < smallest) smallest = currentGmer;
            kmer >>>= 2;
        }

        return (short) smallest;
    }
}
