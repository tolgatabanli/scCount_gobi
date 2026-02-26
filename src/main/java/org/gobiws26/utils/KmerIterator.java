package org.gobiws26.utils;

import java.util.Arrays;
import java.util.Iterator;
import java.util.NoSuchElementException;

public class KmerIterator implements Iterator<byte[]> {
    private final int K;
    private final byte[] fullSeq;
    private int windowStart = 0;

    public KmerIterator(byte[] sequence, int k) {
        if (k > sequence.length) throw new IllegalArgumentException("k cannot be larger than sequence length");
        if (k <= 0) throw new IllegalArgumentException("k must be positive");

        K = k;
        fullSeq = sequence;
    }

    @Override
    public boolean hasNext() {
        return windowStart + K <= fullSeq.length;
    }

    @Override
    public byte[] next() {
        if (!hasNext()) throw new NoSuchElementException();
        byte[] kmer = Arrays.copyOfRange(fullSeq, windowStart, windowStart + k);
        windowStart++;
        return kmer;
    }
}
