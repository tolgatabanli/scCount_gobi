package org.gobiws26;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import org.gobiws26.Indexing.BinaryIndexWriter;
import org.gobiws26.Indexing.Indexer;
import org.gobiws26.Readers.GTFReader;
import org.gobiws26.genomicstruct.Transcript;
import org.gobiws26.utils.TranscriptomeFetcher;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;

public class Test {
    public static File fastaFile = new File("/mnt/raidbio2/extdata/praktikum/genprakt/genprakt-ws25/Block/pig-genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz");
    public static File fastaIndex = new File("/mnt/raidbio2/extdata/praktikum/genprakt/genprakt-ws25/Block/pig-genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz.fai");
    public static File gtfFile = new File("/mnt/raidbio2/extdata/praktikum/genprakt/genprakt-ws25/Block/pig-genome/Sus_scrofa.Sscrofa11.1.115.chr.gtf.gz");
    public static File indexFile = new File("/mnt/cip/home/t/tabanli/Desktop/scCount/out/sccount.idx");
    public static HashMap<String, Transcript> transcripts;
    public static ReferenceSequenceFile fasta = ReferenceSequenceFileFactory.getReferenceSequenceFile(fastaFile);
    public static TranscriptomeFetcher tf = new TranscriptomeFetcher(fasta);

    static {
        Config.K = 20;
        try {
            transcripts = (new GTFReader()).read(gtfFile);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        Indexer indexer = new Indexer(transcripts, ReferenceSequenceFileFactory.getReferenceSequenceFile(fastaFile));
        indexer.runIndex();
        try {
            BinaryIndexWriter.write(indexFile, indexer);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }


    @BeforeAll
    public static void globalSetup() throws IOException {

    }

    @AfterAll
    static void closingRoutine() throws IOException {
        fasta.close();
    }

     // ===================== \\
    //        Utilities        \\

    public static String shortToNucleotideString(short value) {
        StringBuilder sequence = new StringBuilder();

        // The first 2 bits are quality flags (ignore!).
        for (int i = 12; i >= 0; i -= 2) {
            int bits = (value >> i) & 0b11;

            switch (bits) {
                case 0b00 -> sequence.append('A');
                case 0b01 -> sequence.append('C');
                case 0b11 -> sequence.append('G');
                case 0b10 -> sequence.append('T');
            }
        }

        return sequence.toString();
    }

    public static short nucleotideStringToShort(String sequence) {
        // init with first two bits (11000000 00000000)
        short value = (short) 0xC000;

        for (int i = 0; i < Math.min(sequence.length(), 7); i++) {
            char nucleotide = sequence.charAt(i);
            int bits = switch (nucleotide) {
                case 'A' -> 0b00;
                case 'C' -> 0b01;
                case 'T' -> 0b10;
                case 'G' -> 0b11;
                default -> throw new IllegalArgumentException("Invalid nucleotide: " + nucleotide);
            };

            int shiftAmount = 12 - (i * 2);
            value |= (bits << shiftAmount);
        }

        return value;
    }

}
