package org.gobiws26.utils;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import org.gobiws26.genomicstruct.Exon;
import org.gobiws26.genomicstruct.Transcript;

import java.util.Arrays;

import static htsjdk.samtools.util.SequenceUtil.reverseComplement;

public class TranscriptomeFetcher {
    private ReferenceSequenceFile fasta;

    public TranscriptomeFetcher(ReferenceSequenceFile fasta) {
        this.fasta = fasta;
    }


    public byte[] fetchTranscriptSequenceOf(Transcript t) {
        byte[] transcriptSequence = new byte[t.getTranscriptomicLength()];
        int transcriptIndex = 0;

        byte[] genomicSequence = fasta.getSubsequenceAt(t.getChr(), t.getStart(), t.getEnd()).getBases();

        for (Exon exon : t.getSortedExons()) {
            for (int exonIndex = exon.getStart(); exonIndex <= exon.getEnd(); exonIndex++) {
                transcriptSequence[transcriptIndex++] = genomicSequence[exonIndex - t.getStart()];
            }
        }
        if (t.isNegativeStranded()) reverseComplement(transcriptSequence);
        return transcriptSequence;
    }

    public byte[] fetchTranscriptSequenceOf(Transcript t, int tailLength) {
        byte[] bases = fetchTranscriptSequenceOf(t);
        return Arrays.copyOfRange(bases, Math.max(0, bases.length - tailLength), bases.length);
    }

    /**
     *
     * @param t Transcript object, of which the sequence should be returned.
     * @param plus How longer the tail should be on top of the 3' UTR
     * @return byte array representing nucleotides.
     */
    public byte[] fetchTranscriptSequenceOfUTRPlus(Transcript t, int plus) {
        int tailLength = t.getThreePrimeUTRLength() + plus;
        System.out.println(tailLength);
        if (tailLength > t.getTranscriptomicLength()) tailLength = t.getTranscriptomicLength();
        return fetchTranscriptSequenceOf(t, tailLength);
    }

    public static String getStringOf(byte[] inBytes) {
        if (inBytes == null) return "";
        return new String(inBytes, java.nio.charset.StandardCharsets.US_ASCII);
    }
}
