package org.gobiws26.utils;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import org.gobiws26.genomicstruct.Exon;
import org.gobiws26.genomicstruct.Region;
import org.gobiws26.genomicstruct.Transcript;

import static htsjdk.samtools.util.SequenceUtil.reverseComplement;

public class TranscriptomeFetcher {
    private ReferenceSequenceFile fasta;

    public TranscriptomeFetcher(ReferenceSequenceFile fasta) {
        this.fasta = fasta;
    }

    public byte[] fetchTranscriptSequenceOf(Transcript t) {
        byte[] transcriptSequence = new byte[t.getTranscriptomicLength()];
        int transcriptIndex = 0;

        System.out.println(t);
        byte[] genomicSequence = fasta.getSubsequenceAt(t.getChr(), t.getStart(), t.getEnd()).getBases(); // TODO: THERE IS A BUG

        for (Exon exon : t.getSortedExons()) {
            for (int exonIndex = exon.getStart(); exonIndex <= exon.getEnd(); exonIndex++) {
                transcriptSequence[transcriptIndex++] = genomicSequence[exonIndex - t.getStart()];
            }
        }
        System.out.println(getStringOf(transcriptSequence));
        if (t.isNegativeStranded()) reverseComplement(transcriptSequence);
        return transcriptSequence;
    }

    public static String getStringOf(byte[] inBytes) {
        if (inBytes == null) return "";
        return new String(inBytes, java.nio.charset.StandardCharsets.US_ASCII);
    }
}
