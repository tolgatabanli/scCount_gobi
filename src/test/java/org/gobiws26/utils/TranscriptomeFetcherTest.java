package org.gobiws26.utils;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import org.gobiws26.Readers.GTFReader;
import org.gobiws26.genomicstruct.Transcript;
import org.junit.jupiter.api.Test;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;

import static org.junit.jupiter.api.Assertions.*;

class TranscriptomeFetcherTest {
    static File fastaFile = new File("/mnt/raidbio2/extdata/praktikum/genprakt/genprakt-ws25/Block/pig-genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz");
    static File fastaIndex = new File("/mnt/raidbio2/extdata/praktikum/genprakt/genprakt-ws25/Block/pig-genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz.fai");
    static File gtfFile = new File("/mnt/raidbio2/extdata/praktikum/genprakt/genprakt-ws25/Block/pig-genome/Sus_scrofa.Sscrofa11.1.115.chr.gtf.gz");
    static HashMap<String, Transcript> transcripts;
    ReferenceSequenceFile fasta = ReferenceSequenceFileFactory.getReferenceSequenceFile(fastaFile);
    TranscriptomeFetcher tf = new TranscriptomeFetcher(fasta);

    static {
        try {
            transcripts = (new GTFReader()).read(gtfFile);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }


    @Test
    void positiveStrandIsFetchedCorrectly(){
        String s = TranscriptomeFetcher.getStringOf(
                tf.fetchTranscriptSequenceOf(transcripts.get("ENSSSCT00000081764"))
        );
        assertEquals("CTGGAGGAGGACGCATGCAAGGGAGAAAAGAATATTTCGCTGCCTGTTTGGAAATCTGAATTGGTTATTGAATATTATGTGTCTCTGATGCACAGGTCCTCCCCCACTTATGAATGTTGA",
                s);
    }

    void negativeStrandIsFetchedCorrectly(){
        String s = TranscriptomeFetcher.getStringOf(
                tf.fetchTranscriptSequenceOf(transcripts.get("ENSSSCT00000103324"))
        );
        assertEquals("ACGGAGACGTCCGCTCGCACGGTGGGAACCCGTGCTTCCTGCGTCTGCAGCTCAGGGTGCTGTGCGGGGGCGGGGGGTGACCCCGAGATGGCGGCGCCCTGGCACGTGATGGTGCGTCCTTGTCTGTTTTGTCTGTGTCCCTGGACCGAGTTCCAGTTCCCACTCTGCCACTTCCCGGCTGTGCCCTGGGGCAGATGT",
                s);
    }

}