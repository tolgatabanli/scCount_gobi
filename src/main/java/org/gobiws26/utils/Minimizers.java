package org.gobiws26.utils;

import it.unimi.dsi.fastutil.shorts.ShortArrayList;
import org.gobiws26.Config;

public class Minimizers {
    public static ShortArrayList of(byte[] bases){
        KmerIteratorLong kiLong = new KmerIteratorLong(bases, Config.K);
        ShortArrayList minimSet = new ShortArrayList();

        // For each Kmer in the read, find the minimizer and store pairwise different minimizers with the worst found quality assigned from kmers
        short currentMinim = getMinimizer(kiLong.nextLong());
        byte worstQualityLastMinim = (byte) ((currentMinim >>> 14) & 0b11); // quality score ordering: 11 > 10 > 01 > 00; TODO: will decide what all may correspond to (11 is normal, no 'N' etc.)
        short lastFoundMinim = currentMinim;
        //int windowSize = Config.K - 7; // to consider "real" sequential duplicates, 7 magic number is okay since minimizers are 7mer since short.
        while(kiLong.hasNext()) {
            currentMinim = getMinimizer(kiLong.nextLong());

            byte currentQualityBits = (byte) ((currentMinim >>> 14) & 0b11); // unsigned shift!
            if (((currentMinim & 0x3FFF) != (lastFoundMinim & 0x3FFF))) { // only add when we change the minimizer sequence, so that we're sure the qualities are considered correctly // || windowSize == 0
                lastFoundMinim = (short) ((lastFoundMinim & 0x3FFF) | (worstQualityLastMinim << 14)); // 'infect' with worst kmer's quality

                minimSet.add(lastFoundMinim);

                lastFoundMinim = currentMinim;
                //windowSize = Config.K - 7;
                worstQualityLastMinim = currentQualityBits;
            } else {
                if (worstQualityLastMinim > currentQualityBits) {
                    worstQualityLastMinim = currentQualityBits;
                }
                //windowSize--;
            }
        }
        lastFoundMinim = (short) ((lastFoundMinim & 0x3FFF) | (worstQualityLastMinim << 14));

        minimSet.add(lastFoundMinim); // flush
        return minimSet;
    }

    private static final short gmerMask = 0x3FFF; // take the rightmost 14 bits
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
}
