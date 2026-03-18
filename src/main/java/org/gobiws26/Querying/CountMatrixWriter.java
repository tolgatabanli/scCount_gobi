package org.gobiws26.Querying;

import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import org.gobiws26.Indexing.IndexData;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public class CountMatrixWriter {
    private BufferedWriter bw;
    private IndexData idxData;

    Int2ObjectOpenHashMap<String> transcriptInterpreter;
    Int2ObjectOpenHashMap<String> geneInterpreter;
    Int2IntOpenHashMap geneOfTranscript;

    public CountMatrixWriter(File file, Int2ObjectOpenHashMap<String> transcriptInterpreter,
                             Int2ObjectOpenHashMap<String> geneInterpreter, Int2IntOpenHashMap geneOfTranscript) throws IOException {
        bw = new BufferedWriter(new FileWriter(file));

        this.transcriptInterpreter = transcriptInterpreter;
        this.geneInterpreter = geneInterpreter;
        this.geneOfTranscript = geneOfTranscript;
    }

    public CountMatrixWriter(File file, IndexData idxData) throws IOException {
        bw = new BufferedWriter(new FileWriter(file));
        this.idxData = idxData;
        this.transcriptInterpreter = idxData.int2TxString();
        this.geneInterpreter = idxData.int2GeneString();
        this.geneOfTranscript = idxData.txInt2GeneInt();
    }

    public void write(Int2IntOpenHashMap geneInt2Counts, Int2IntOpenHashMap txInt2Counts) throws IOException {
        // Group transcripts by gene
        Int2ObjectOpenHashMap<Int2IntOpenHashMap> geneToTxCounts = new Int2ObjectOpenHashMap<>();
        for (Int2IntMap.Entry entry : txInt2Counts.int2IntEntrySet()) {
            int txInt = entry.getIntKey();
            int count = entry.getIntValue();
            int geneInt = geneOfTranscript.get(txInt);
            geneToTxCounts.computeIfAbsent(geneInt, k -> new Int2IntOpenHashMap()).put(txInt, count);
        }

        // Iterate through genes and write gene-level counts followed by their transcripts
        for (Int2IntMap.Entry entry : geneInt2Counts.int2IntEntrySet()) {
            int geneInt = entry.getIntKey();
            int count = entry.getIntValue();
            String geneId = geneInterpreter.get(geneInt);

            if (geneId != null) {
                bw.write(geneId + "\t" + count);
                bw.newLine();

                // Write transcripts for this gene
                Int2IntOpenHashMap txCounts = geneToTxCounts.get(geneInt);
                if (txCounts != null) {
                    for (Int2IntMap.Entry txEntry : txCounts.int2IntEntrySet()) {
                        int txInt = txEntry.getIntKey();
                        int txCount = txEntry.getIntValue();
                        String txId = transcriptInterpreter.get(txInt);
                        if (txId != null) {
                            bw.write(txId + "\t" + txCount);
                            bw.newLine();
                        }
                    }
                }
            }
        }
        bw.flush();
    }

    /**
     * Write counts as a matrix with spatial barcodes as columns.
     * Rows are genes and transcripts, columns are spatial barcodes.
     * Format:
     *   gene_id\tbarcode1\tbarcode2\tbarcode3...
     */
    public void writeBarcodeMatrix(Int2ObjectOpenHashMap<Int2IntOpenHashMap> geneCountsPerBarcode,
                                   Int2ObjectOpenHashMap<Int2IntOpenHashMap> txCountsPerBarcode,
                                   BarcodeMapper barcodeMapper) throws IOException {
        // Get sorted barcode IDs
        int[] barcodeIds = barcodeMapper.getAllIds();

        // Write header: gene_id followed by all barcodes
        bw.write("gene_id");
        for (int barcodeId : barcodeIds) {
            bw.write("\t");
            bw.write(barcodeMapper.getBarcode(barcodeId));
        }
        bw.newLine();

        // Collect all gene IDs to output
        java.util.Set<Integer> geneInts = new java.util.HashSet<>();
        for (var entry : geneCountsPerBarcode.int2ObjectEntrySet()) {
            Int2IntOpenHashMap geneCounts = entry.getValue();
            for (int geneInt : geneCounts.keySet()) {
                geneInts.add(geneInt);
            }
        }

        // Write gene rows
        for (int geneInt : geneInts) {
            String geneId = geneInterpreter.get(geneInt);
            if (geneId != null) {
                bw.write(geneId);
                for (int barcodeId : barcodeIds) {
                    bw.write("\t");
                    Int2IntOpenHashMap geneCounts = geneCountsPerBarcode.get(barcodeId);
                    int count = (geneCounts != null) ? geneCounts.getOrDefault(geneInt, 0) : 0;
                    bw.write(String.valueOf(count));
                }
                bw.newLine();
            }
        }

        // Collect all transcript IDs to output
        java.util.Set<Integer> txInts = new java.util.HashSet<>();
        for (var entry : txCountsPerBarcode.int2ObjectEntrySet()) {
            Int2IntOpenHashMap txCounts = entry.getValue();
            for (int txInt : txCounts.keySet()) {
                txInts.add(txInt);
            }
        }

        // Write transcript rows
        for (int txInt : txInts) {
            String txId = transcriptInterpreter.get(txInt);
            if (txId != null) {
                bw.write(txId);
                for (int barcodeId : barcodeIds) {
                    bw.write("\t");
                    Int2IntOpenHashMap txCounts = txCountsPerBarcode.get(barcodeId);
                    int count = (txCounts != null) ? txCounts.getOrDefault(txInt, 0) : 0;
                    bw.write(String.valueOf(count));
                }
                bw.newLine();
            }
        }

        bw.flush();
    }

    /**
     * Write counts as a sparse matrix in Market Matrix format (MTX).
     * Creates three files in the output directory:
     * - matrix.mtx: Sparse coordinate matrix (industry-standard format)
     * - barcodes.tsv: Barcode identifiers (one per line)
     * - features.tsv: Gene/transcript identifiers (one per line)
     */
    public void writeBarcodeSparseMatrix(Int2ObjectOpenHashMap<Int2IntOpenHashMap> geneCountsPerBarcode,
                                         Int2ObjectOpenHashMap<Int2IntOpenHashMap> txCountsPerBarcode,
                                         BarcodeMapper barcodeMapper,
                                         File outputDir) throws IOException {
        
        long startTime = System.currentTimeMillis();
        
        // Collect all unique gene and transcript IDs (features)
        java.util.Set<Integer> featureInts = new java.util.HashSet<>();
        long totalNonZero = 0;
        
        for (var entry : geneCountsPerBarcode.int2ObjectEntrySet()) {
            Int2IntOpenHashMap geneCounts = entry.getValue();
            featureInts.addAll(geneCounts.keySet());
            totalNonZero += geneCounts.size();
        }
        
        for (var entry : txCountsPerBarcode.int2ObjectEntrySet()) {
            Int2IntOpenHashMap txCounts = entry.getValue();
            featureInts.addAll(txCounts.keySet());
            totalNonZero += txCounts.size();
        }
        
        int numBarcodes = barcodeMapper.size();
        int numFeatures = featureInts.size();
        
        System.out.printf("Writing sparse matrix: %d features × %d barcodes (%d non-zero entries)%n", 
                          numFeatures, numBarcodes, totalNonZero);
        
        // Write matrix.mtx (Market Matrix format - sparse coordinate)
        File matrixFile = new File(outputDir, "matrix.mtx");
        try (BufferedWriter mtxWriter = new BufferedWriter(new FileWriter(matrixFile))) {
            mtxWriter.write("%%MatrixMarket matrix coordinate integer general\n");
            mtxWriter.write("%\n");
            mtxWriter.write("% Sparse count matrix from scCount\n");
            mtxWriter.write("% Rows: genes/transcripts | Columns: spatial barcodes\n");
            mtxWriter.write("%\n");
            mtxWriter.write(numFeatures + " " + numBarcodes + " " + totalNonZero + "\n");
            
            // Write all non-zero entries in 1-indexed format
            // MTX format: feature_index(1-based) barcode_index(1-based) count
            int featureIdx = 1;
            java.util.Map<Integer, Integer> featureIntToIdx = new java.util.HashMap<>();
            
            for (int featureInt : featureInts) {
                featureIntToIdx.put(featureInt, featureIdx++);
            }
            
            // Write entries for genes
            for (var barcodeEntry : geneCountsPerBarcode.int2ObjectEntrySet()) {
                int barcodeInt = barcodeEntry.getIntKey();
                int barcodeIdx = barcodeInt + 1; // 1-indexed for MTX
                Int2IntOpenHashMap geneCounts = barcodeEntry.getValue();
                
                for (var geneEntry : geneCounts.int2IntEntrySet()) {
                    int geneInt = geneEntry.getIntKey();
                    int count = geneEntry.getIntValue();
                    int geneIdx = featureIntToIdx.get(geneInt);
                    mtxWriter.write(geneIdx + " " + barcodeIdx + " " + count + "\n");
                }
            }
            
            // Write entries for transcripts
            for (var barcodeEntry : txCountsPerBarcode.int2ObjectEntrySet()) {
                int barcodeInt = barcodeEntry.getIntKey();
                int barcodeIdx = barcodeInt + 1; // 1-indexed for MTX
                Int2IntOpenHashMap txCounts = barcodeEntry.getValue();
                
                for (var txEntry : txCounts.int2IntEntrySet()) {
                    int txInt = txEntry.getIntKey();
                    int count = txEntry.getIntValue();
                    int txIdx = featureIntToIdx.get(txInt);
                    mtxWriter.write(txIdx + " " + barcodeIdx + " " + count + "\n");
                }
            }
        }
        
        // Write barcodes.tsv
        File barcodesFile = new File(outputDir, "barcodes.tsv");
        try (BufferedWriter barcodesWriter = new BufferedWriter(new FileWriter(barcodesFile))) {
            int[] barcodeIds = barcodeMapper.getAllIds();
            for (int barcodeId : barcodeIds) {
                barcodesWriter.write(barcodeMapper.getBarcode(barcodeId));
                barcodesWriter.newLine();
            }
        }
        
        // Write features.tsv (gene/transcript names and type)
        File featuresFile = new File(outputDir, "features.tsv");
        try (BufferedWriter featuresWriter = new BufferedWriter(new FileWriter(featuresFile))) {
            java.util.List<Integer> sortedFeatureInts = new ArrayList<>(featureInts);
            java.util.Collections.sort(sortedFeatureInts);
            
            for (int featureInt : sortedFeatureInts) {
                String featureName = geneInterpreter.get(featureInt);
                if (featureName == null) {
                    featureName = transcriptInterpreter.get(featureInt);
                }
                
                String featureType = geneInterpreter.get(featureInt) != null ? "Gene" : "Transcript";
                
                if (featureName != null) {
                    featuresWriter.write(featureName + "\t" + featureType);
                    featuresWriter.newLine();
                }
            }
        }
        
        long elapsedTime = System.currentTimeMillis() - startTime;
        System.out.printf("Sparse matrix written to %s in %.2f seconds%n", 
                          outputDir.getAbsolutePath(), elapsedTime / 1000.0);
        System.out.println("Files created:");
        System.out.println("  - matrix.mtx (sparse count matrix)");
        System.out.println("  - barcodes.tsv (barcode identifiers)");
        System.out.println("  - features.tsv (gene/transcript identifiers)");
    }
}
