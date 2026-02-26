package org.gobiws26.Readers;

import htsjdk.samtools.fastq.FastqRecord;

import java.io.File;

public class FastqReader implements AutoCloseable {
    private htsjdk.samtools.fastq.FastqReader fastqReader;
    public FastqReader(File file) {
        fastqReader = new htsjdk.samtools.fastq.FastqReader(file);
    }

    public String getNextRead() {
        FastqRecord fr = fastqReader.next();
        return fr.getReadName() + "\n" + fr.getReadString();
    }


    @Override
    public void close() throws Exception {
        fastqReader.close();
    }
}
