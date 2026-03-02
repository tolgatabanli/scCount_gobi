package org.gobiws26.genomicstruct;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public class Transcript extends Region {
    private String chromosome;
    private String geneId;
    private boolean isNegativeStranded;
    private final List<Exon> exons = new ArrayList<>();

    private int transcriptomicLength = 0;
    private int threePrimeUTRLength = 0;

    public Transcript() {
        super();
    }

    public void setStart(int start) {
        this.start = start;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    public void setChr(String chromosome) {
        this.chromosome = chromosome;
    }

    public void setStrand(boolean isNegativeStranded) {
        this.isNegativeStranded = isNegativeStranded;
    }

    public void addExon(Exon exon) {
        this.exons.add(exon);
        transcriptomicLength += exon.length();
    }

    public void setGene(String geneId) {
        this.geneId = geneId;
    }

    public String getGeneId() {
        return geneId;
    }

    public String getChr() {
        return chromosome;
    }

    public boolean isNegativeStranded() {
        return isNegativeStranded;
    }

    public int getTranscriptomicLength() {
        return transcriptomicLength;
    }

    public List<Exon> getSortedExons() {
        exons.sort(Comparator.comparingInt(Exon::getStart));
        return Collections.unmodifiableList(exons);
    }

    public void addToUTR(int utrWidth) {
        this.threePrimeUTRLength += utrWidth;
    }

    @Override
    public String toString() {
        return "Transcript\t" + "chr:" + chromosome + "\t" + (isNegativeStranded ? "-" : "+") + "\t" + start + "-" + end;
    }
}
