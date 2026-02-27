package org.gobiws26.utils;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import org.gobiws26.Config;
import org.gobiws26.Readers.GTFReader;
import org.gobiws26.genomicstruct.Transcript;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;

public class Test {
    static File fastaFile = new File("/mnt/raidbio2/extdata/praktikum/genprakt/genprakt-ws25/Block/pig-genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz");
    static File fastaIndex = new File("/mnt/raidbio2/extdata/praktikum/genprakt/genprakt-ws25/Block/pig-genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz.fai");
    static File gtfFile = new File("/mnt/raidbio2/extdata/praktikum/genprakt/genprakt-ws25/Block/pig-genome/Sus_scrofa.Sscrofa11.1.115.chr.gtf.gz");
    static HashMap<String, Transcript> transcripts;
    static ReferenceSequenceFile fasta = ReferenceSequenceFileFactory.getReferenceSequenceFile(fastaFile);
    static TranscriptomeFetcher tf = new TranscriptomeFetcher(fasta);

    static {
        Config.K = 20;
        try {
            transcripts = (new GTFReader()).read(gtfFile);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    @BeforeAll
    static void globalSetup() {}

    @AfterAll
    static void closingRoutine() throws IOException {
        fasta.close();
    }

}
