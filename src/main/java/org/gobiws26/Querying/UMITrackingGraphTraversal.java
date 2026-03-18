package org.gobiws26.Querying;

import htsjdk.samtools.fastq.FastqRecord;
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;

import java.util.Set;

/**
 * Tracks counts per barcode with UMI-based deduplication using shared, thread-safe sets.
 * Each thread receives references to globally-shared deduplication sets to ensure
 * correct handling of duplicate UMIs across concurrent threads.
 */
public class UMITrackingGraphTraversal {
    private final IndexGraph g;
    private final Int2IntMap tx2geneMapping;

    // Maps: barcodeInt -> (geneInt -> count)
    private final Int2ObjectOpenHashMap<Int2IntOpenHashMap> geneCountsPerBarcode = new Int2ObjectOpenHashMap<>();
    // Maps: barcodeInt -> (txInt -> count)
    private final Int2ObjectOpenHashMap<Int2IntOpenHashMap> txCountsPerBarcode = new Int2ObjectOpenHashMap<>();

    // Shared thread-safe sets for deduplication (preventing double-counting across threads)
    // Format: "barcodeInt:geneInt:umi"
    private final Set<String> sharedSeenGeneBarcodeUMIs;
    // Format: "barcodeInt:txInt:umi"
    private final Set<String> sharedSeenTxBarcodeUMIs;

    private int totalAmbigReads = 0;
    private int totalShortReadsDiscarded = 0;
    private int totalReadsAssignedToTranscripts = 0;
    private int totalReadsAssignedToGenes = 0;

    public UMITrackingGraphTraversal(IndexGraph g, Int2IntMap tx2gene, 
                                     Set<String> sharedSeenGeneBarcodeUMIs,
                                     Set<String> sharedSeenTxBarcodeUMIs) {
        this.g = g;
        this.tx2geneMapping = tx2gene;
        this.sharedSeenGeneBarcodeUMIs = sharedSeenGeneBarcodeUMIs;
        this.sharedSeenTxBarcodeUMIs = sharedSeenTxBarcodeUMIs;
    }

    /**
     * Process a read and track its counts per barcode with UMI deduplication.
     * Uses thread-safe shared deduplication sets to prevent double-counting across threads.
     * @param read the FastqRecord to process
     * @param barcodeInt the spatial barcode as an integer identifier
     * @param umi the UMI string
     */
    public void processWithBarcode(FastqRecord read, int barcodeInt, String umi) {
        // Use a temporary traversal to get gene/transcript info
        IndexGraphTraversal temp = new IndexGraphTraversal(g, tx2geneMapping);
        temp.process(read);
        
        // Get the counts from this read
        Int2IntOpenHashMap tempGeneCounts = temp.getGeneCounts();
        Int2IntOpenHashMap tempTxCounts = temp.getTxCounts();
        
        // Transfer to barcode-aware tracking with UMI dedup using SHARED thread-safe sets
        // Set.add() returns false if already present, enabling atomic check-and-set
        for (var entry : tempGeneCounts.int2IntEntrySet()) {
            int geneInt = entry.getIntKey();
            String key = barcodeInt + ":" + geneInt + ":" + umi;
            if (sharedSeenGeneBarcodeUMIs.add(key)) {
                geneCountsPerBarcode
                    .computeIfAbsent(barcodeInt, k -> new Int2IntOpenHashMap())
                    .addTo(geneInt, 1);
            }
        }
        
        for (var entry : tempTxCounts.int2IntEntrySet()) {
            int txInt = entry.getIntKey();
            String key = barcodeInt + ":" + txInt + ":" + umi;
            if (sharedSeenTxBarcodeUMIs.add(key)) {
                txCountsPerBarcode
                    .computeIfAbsent(barcodeInt, k -> new Int2IntOpenHashMap())
                    .addTo(txInt, 1);
            }
        }
        
        totalAmbigReads += temp.getAmbigReads();
        totalShortReadsDiscarded += temp.getShortReadsDiscarded();
        totalReadsAssignedToTranscripts += temp.getReadsAssignedToTranscripts();
        totalReadsAssignedToGenes += temp.getReadsAssignedToGenes();
    }

    public Int2ObjectOpenHashMap<Int2IntOpenHashMap> getGeneCountsPerBarcode() {
        return geneCountsPerBarcode;
    }

    public Int2ObjectOpenHashMap<Int2IntOpenHashMap> getTxCountsPerBarcode() {
        return txCountsPerBarcode;
    }

    public int getTotalAmbigReads() {
        return totalAmbigReads;
    }

    public int getTotalShortReadsDiscarded() {
        return totalShortReadsDiscarded;
    }

    public int getReadsAssignedToTranscripts() {
        return totalReadsAssignedToTranscripts;
    }

    public int getReadsAssignedToGenes() {
        return totalReadsAssignedToGenes;
    }
}
