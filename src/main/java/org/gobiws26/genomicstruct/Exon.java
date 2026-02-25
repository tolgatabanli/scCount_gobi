package org.gobiws26.genomicstruct;

public class Exon extends Region {

    public Exon(int start, int end) {
        super(start, end);
    }


    @Override
    public String toString() {
        return super.getStart() + "-" + super.getEnd();
    }
}
