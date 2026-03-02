package org.gobiws26.genomicstruct;

import java.util.Objects;

public class Exon extends Region {

    public Exon(int start, int end) {
        super(start, end);
    }


    @Override
    public String toString() {
        return super.getStart() + "-" + super.getEnd();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Region that = (Region) o;
        return start == that.start && end == that.end;
    }

    @Override
    public int hashCode() {
        return Objects.hash(start, end);
    }
}
