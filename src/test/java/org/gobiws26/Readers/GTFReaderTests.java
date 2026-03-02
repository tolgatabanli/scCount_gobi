package org.gobiws26.Readers;

import org.gobiws26.genomicstruct.Exon;
import org.gobiws26.genomicstruct.Transcript;
import org.junit.jupiter.api.Test;

import static org.gobiws26.Test.transcripts;
import static org.junit.jupiter.api.Assertions.*;

public class GTFReaderTests {

    @Test
    public void GTFRecordIsCorrect() {
        String t_id = "ENSSSCT00000103363";
        Transcript t = transcripts.get(t_id);
        assertEquals("ENSSSCG00000028996", t.getGeneId());
        assertEquals("1", t.getChr());
        assertTrue(t.isNegativeStranded());
        assertTrue(t.getSortedExons().contains(new Exon(226202702, 226202806)));
    }
}
