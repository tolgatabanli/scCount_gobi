package org.gobiws26.utils;

import org.gobiws26.Config;
import org.gobiws26.Readers.GTFReader;
import org.gobiws26.genomicstruct.Transcript;

import java.io.IOException;
import java.util.HashMap;
import org.junit.jupiter.api.Test;
import static org.gobiws26.utils.Test.*;

public class HashingBenchmark {

    @Test
    public void KmerIterTiming() throws IOException {
        HashMap<String, Transcript> transcripts = (new GTFReader()).read(gtfFile);

        Transcript debugT = transcripts.get("ENSSSCT00000092142");

        try {
            byte[] debugSeq = fasta.getSubsequenceAt(debugT.getChr(), debugT.getStart(), debugT.getEnd()).getBases();

            // Benchmarking for KmerIteratorHash (long)
            long startHash = System.nanoTime();
            KmerIteratorHash kiHash = new KmerIteratorHash(debugSeq, Config.K);
            while (kiHash.hasNext()) {
                long kmer = kiHash.next();
                //System.out.println(kmer);
            }
            long endHash = System.nanoTime();

            // Benchmarking for KmerIteratorByte (byte[])
            long startByte = System.nanoTime();
            KmerIteratorByte kiByte = new KmerIteratorByte(debugSeq, Config.K);
            while (kiByte.hasNext()) {
                byte[] kmer = kiByte.next();
                //System.out.println(Arrays.toString(kmer));
            }
            long endByte = System.nanoTime();

            // Benchmarking for KmerIteratorLong
            long startLong = System.nanoTime();
            KmerIteratorLong kiLong = new KmerIteratorLong(debugSeq, Config.K);
            int i = 0;
            while (kiLong.hasNext()) {
                i++;
                long kmer = kiLong.next();
                //System.out.println(String.format("%40s", Long.toBinaryString(kmer)).replace(' ', '0'));
            }
            System.out.println(i);
            long endLong = System.nanoTime();

            // Results
            System.out.println("--- Benchmarking Results ---");
            System.out.println("KmerIteratorHash (long):   " + (endHash - startHash) / 1_000_000.0 + " ms");
            System.out.println("KmerIteratorByte (byte[]): " + (endByte - startByte) / 1_000_000.0 + " ms");
            System.out.println("KmerIteratorLong (long): " + (endLong - startLong) / 1_000_000.0 + " ms");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
