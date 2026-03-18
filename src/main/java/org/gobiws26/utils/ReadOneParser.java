package org.gobiws26.utils;

import htsjdk.samtools.fastq.FastqRecord;

/**
 * Extracts spatial barcode (27bp) and UMI (16bp) from readOne sequences.
 * Expected format: [27bp barcode][16bp UMI][rest of sequence]
 */
public class ReadOneParser {
    public static final int BARCODE_LENGTH = 27;
    public static final int UMI_LENGTH = 16;


    /**
     * Extracts both barcode and UMI in one call for efficiency.
     * @param read the FastqRecord from readOne file
     * @return a BarcodeUMI object, or null if sequence is too short
     */
    public static BarcodeUMI extract(FastqRecord read) {
        byte[] sequence = read.getReadBases();
        if (sequence.length < BARCODE_LENGTH + UMI_LENGTH) {
            return null;
        }
        String barcode = new String(sequence, 0, BARCODE_LENGTH);
        String umi = new String(sequence, BARCODE_LENGTH, UMI_LENGTH);
        return new BarcodeUMI(barcode, umi);
    }

    /**
         * Container class for barcode and UMI pair.
         */
        public record BarcodeUMI(String barcode, String umi) {

            @Override
            public String toString() {
                return barcode + "-" + umi;
            }
        }
}
