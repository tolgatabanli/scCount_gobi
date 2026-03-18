package org.gobiws26.Querying;

import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import org.gobiws26.Indexing.IndexData;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class CountMatrixWriter {
    private final BufferedWriter bw;

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
        
        // Collect unique gene and transcript IDs using fastutil IntSet for efficiency
        IntOpenHashSet geneInts = new IntOpenHashSet();
        IntOpenHashSet txInts = new IntOpenHashSet();
        long totalNonZero = 0;
        
        for (var entry : geneCountsPerBarcode.int2ObjectEntrySet()) {
            geneInts.addAll(entry.getValue().keySet());
            totalNonZero += entry.getValue().size();
        }
        
        for (var entry : txCountsPerBarcode.int2ObjectEntrySet()) {
            txInts.addAll(entry.getValue().keySet());
            totalNonZero += entry.getValue().size();
        }
        
        // Sort IDs for consistent matrix indexing
        IntArrayList sortedGeneIds = new IntArrayList(geneInts);
        IntArrayList sortedTxIds = new IntArrayList(txInts);
        sortedGeneIds.sort(null);
        sortedTxIds.sort(null);
        
        int numBarcodes = barcodeMapper.size();
        int numGenes = sortedGeneIds.size();
        int numTranscripts = sortedTxIds.size();
        int numFeatures = numGenes + numTranscripts;
        
        System.out.printf("Writing sparse matrix: %d features (%d genes + %d transcripts) × %d barcodes (%d non-zero entries)%n", 
                          numFeatures, numGenes, numTranscripts, numBarcodes, totalNonZero);
        
        // Create mappings from internal IDs to matrix indices
        Int2IntOpenHashMap geneIntToIdx = new Int2IntOpenHashMap();
        Int2IntOpenHashMap txIntToIdx = new Int2IntOpenHashMap();
        
        for (int i = 0; i < sortedGeneIds.size(); i++) {
            geneIntToIdx.put(sortedGeneIds.getInt(i), i + 1);
        }
        
        for (int i = 0; i < sortedTxIds.size(); i++) {
            txIntToIdx.put(sortedTxIds.getInt(i), numGenes + i + 1);
        }
        
        // Write matrix.mtx (Market Matrix format - sparse coordinate)
        File matrixFile = new File(outputDir, "matrix.mtx");
        try (BufferedWriter mtxWriter = new BufferedWriter(new FileWriter(matrixFile))) {
            mtxWriter.write("%%MatrixMarket matrix coordinate integer general\n");
            mtxWriter.write("%\n");
            mtxWriter.write("% Sparse count matrix from MiniQuT3\n");
            mtxWriter.write("% Rows: genes/transcripts | Columns: spatial barcodes\n");
            mtxWriter.write("%\n");
            mtxWriter.write(numFeatures + " " + numBarcodes + " " + totalNonZero + "\n");
            
            // Write entries for genes
            writeEntriesPerFeature(geneCountsPerBarcode, mtxWriter, geneIntToIdx);

            // Write entries for transcripts
            writeEntriesPerFeature(txCountsPerBarcode, mtxWriter, txIntToIdx);
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
        
        // Write features.tsv (gene/transcript names only - Ensembl IDs indicate type)
        File featuresFile = new File(outputDir, "features.tsv");
        try (BufferedWriter featuresWriter = new BufferedWriter(new FileWriter(featuresFile))) {
            // Write genes first (indices 1 to numGenes) - already sorted
            writeFeatures(sortedGeneIds, featuresWriter, geneInterpreter);

            // Write transcripts next (indices numGenes+1 onwards) - already sorted
            writeFeatures(sortedTxIds, featuresWriter, transcriptInterpreter);
        }
        
        long elapsedTime = System.currentTimeMillis() - startTime;
        System.out.printf("Sparse matrix written to %s in %.2f seconds%n", 
                          outputDir.getAbsolutePath(), elapsedTime / 1000.0);
        System.out.println("Files created:");
        System.out.println("  - matrix.mtx (sparse count matrix)");
        System.out.println("  - barcodes.tsv (barcode identifiers)");
        System.out.println("  - features.tsv (gene/transcript identifiers)");
    }

    private void writeFeatures(IntArrayList sortedFeatureIds, BufferedWriter featuresWriter, Int2ObjectOpenHashMap<String> int2FeatInterpreter) throws IOException {
        for (int i = 0; i < sortedFeatureIds.size(); i++) {
            int featureId = sortedFeatureIds.getInt(i);
            String featureName = int2FeatInterpreter.get(featureId);
            if (featureName != null) {
                featuresWriter.write(featureName);
                featuresWriter.newLine();
            }
        }
    }

    private void writeEntriesPerFeature(Int2ObjectOpenHashMap<Int2IntOpenHashMap> countsPerBarcode, BufferedWriter mtxWriter, Int2IntOpenHashMap intToIdInterpreter) throws IOException {
        for (var barcodeEntry : countsPerBarcode.int2ObjectEntrySet()) {
            int barcodeInt = barcodeEntry.getIntKey();
            int barcodeIdx = barcodeInt + 1; // 1-indexed for MTX
            Int2IntOpenHashMap geneCounts = barcodeEntry.getValue();

            for (var geneEntry : geneCounts.int2IntEntrySet()) {
                int geneInt = geneEntry.getIntKey();
                int count = geneEntry.getIntValue();
                int mappedGeneIdx = intToIdInterpreter.get(geneInt);
                mtxWriter.write(mappedGeneIdx + " " + barcodeIdx + " " + count + "\n");
            }
        }
    }
}
