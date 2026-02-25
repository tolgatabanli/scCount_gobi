package org.gobiws26.genomicstruct;

public class Region {
    int start;
    int end;

    public Region() {
    }

    public Region(int start, int end) {
        this.start = start;
        this.end = end;
    }

    public int getStart() {
        return start;
    }
    public int getEnd() {
        return end;
    }

    public int length() {
        return end - start + 1;
    }
}
