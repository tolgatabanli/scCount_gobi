package org.gobiws26.Readers;

import htsjdk.samtools.fastq.FastqRecord;

import java.io.File;

public class FastqReader implements AutoCloseable {
    private htsjdk.samtools.fastq.FastqReader fastqReader;
    public FastqReader(File file) {
        fastqReader = new htsjdk.samtools.fastq.FastqReader(file);
    }

    private byte[] currentBases;
    private byte[] currentQs;

    public FastqRecord getNextRecord() {
        return fastqReader.next();
    }

    public boolean  hasNext() {
        return fastqReader.hasNext();
    }


    @Override
    public void close() throws Exception {
        fastqReader.close();
    }
}
