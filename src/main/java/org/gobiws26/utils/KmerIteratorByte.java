package org.gobiws26.utils;

import java.util.Arrays;
import java.util.NoSuchElementException;

import java.util.Arrays;
import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * FOR DEBUGGING PROBABLY
 * Takes base sequence as byte array as output from samtools and creates an iterator
 * which returns kmers of the sequence on each next() call.
 * Note: samtools return a byte array of nucleotides (chars: A:65, C:67, G:71, T:84)
 */
public class KmerIteratorByte implements Iterator<byte[]> {
    private final int K;
    private final byte[] fullSeq;
    private int windowStart = 0;

    /**
     * @param sequence A sequence.
     * @param k kmer size. Could also use Config.K. Might be removed.
     */
    public KmerIteratorByte(byte[] sequence, int k) {
        if (k <= 0) throw new IllegalArgumentException("k must be positive");
        if (k > sequence.length) throw new IllegalArgumentException("k cannot be larger than sequence length");

        this.K = k;
        this.fullSeq = sequence;
    }

    @Override
    public boolean hasNext() {
        return windowStart + K <= fullSeq.length;
    }

    @Override
    public byte[] next() {
        if (!hasNext()) throw new NoSuchElementException();
        byte[] kmer = Arrays.copyOfRange(fullSeq, windowStart, windowStart + K);
        windowStart++;
        return kmer;
    }
}
