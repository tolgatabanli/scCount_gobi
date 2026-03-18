package org.gobiws26.Querying;

import htsjdk.samtools.fastq.FastqRecord;
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.booleans.BooleanArrayList;

import java.io.*;
import java.util.*;
import java.util.concurrent.*;

public class ParallelGraphQuery {
    private final IndexGraph idxGraph;
    private final Int2IntMap tx2gene;
    private final int numThreads;
    private final ExecutorService executor;

    Int2IntOpenHashMap globalGeneCounts;
    Int2IntOpenHashMap globalTxCounts;
    
    // For barcode-aware mode
    private final BarcodeMapper barcodeMapper;
    private Int2ObjectOpenHashMap<Int2IntOpenHashMap> geneCountsPerBarcode;
    private Int2ObjectOpenHashMap<Int2IntOpenHashMap> txCountsPerBarcode;
    
    private int totalAmbigReads;
    private int totalShortReadsDiscarded;
    private int totalReadsAssignedToTranscripts;
    private int totalReadsAssignedToGenes;
    
    // Separated by level (transcript vs gene) and phase
    private int numTxFoundHeuristically = 0;
    private int numGenesFoundHeuristically = 0;
    private int numTxFoundWithAlignment = 0;
    private int numGenesFoundWithAlignment = 0;
    
    // Track per-read mapping details if writeDetails is true
    // Uses int encoding to save memory: -2=short, -1=ambiguous, >=0=actual ID
    private final boolean writeDetails;
    private final IntArrayList readMappingIds;
    private final BooleanArrayList readMappingIsGene;

    public ParallelGraphQuery(IndexGraph g, Int2IntMap tx2gene, int numThreads) {
        this.idxGraph = g;
        this.tx2gene = tx2gene;
        this.numThreads = Math.min(numThreads, Runtime.getRuntime().availableProcessors());
        this.executor = Executors.newFixedThreadPool(this.numThreads);
        this.globalGeneCounts = new Int2IntOpenHashMap();
        this.globalTxCounts   = new Int2IntOpenHashMap();
        this.barcodeMapper = null;
        this.writeDetails = false;
        this.readMappingIds = null;
        this.readMappingIsGene = null;
    }

    /**
     * Constructor for barcode-aware mode with UMI deduplication.
     */
    public ParallelGraphQuery(IndexGraph g, Int2IntMap tx2gene, int numThreads, boolean withBarcodes) {
        this.idxGraph = g;
        this.tx2gene = tx2gene;
        this.numThreads = Math.min(numThreads, Runtime.getRuntime().availableProcessors());
        this.executor = Executors.newFixedThreadPool(this.numThreads);
        this.globalGeneCounts = new Int2IntOpenHashMap();
        this.globalTxCounts   = new Int2IntOpenHashMap();
        this.barcodeMapper = withBarcodes ? new BarcodeMapper() : null;
        this.writeDetails = false;
        this.readMappingIds = null;
        this.readMappingIsGene = null;
        if (withBarcodes) {
            this.geneCountsPerBarcode = new Int2ObjectOpenHashMap<>();
            this.txCountsPerBarcode = new Int2ObjectOpenHashMap<>();
        }
    }

    /**
     * Constructor with optional per-read mapping tracking.
     */
    public ParallelGraphQuery(IndexGraph g, Int2IntMap tx2gene, int numThreads, boolean withBarcodes, boolean trackDetails) {
        this.idxGraph = g;
        this.tx2gene = tx2gene;
        this.numThreads = Math.min(numThreads, Runtime.getRuntime().availableProcessors());
        this.executor = Executors.newFixedThreadPool(this.numThreads);
        this.globalGeneCounts = new Int2IntOpenHashMap();
        this.globalTxCounts   = new Int2IntOpenHashMap();
        this.barcodeMapper = withBarcodes ? new BarcodeMapper() : null;
        this.writeDetails = trackDetails;
        this.readMappingIds = trackDetails ? new IntArrayList() : null;
        this.readMappingIsGene = trackDetails ? new BooleanArrayList() : null;
        if (withBarcodes) {
            this.geneCountsPerBarcode = new Int2ObjectOpenHashMap<>();
            this.txCountsPerBarcode = new Int2ObjectOpenHashMap<>();
        }
    }

    /**
     * Process reads with barcode information and UMI deduplication.
     * @param readTwoRecords the sequencing reads
     * @param readOneBarcodes the barcode strings from readOne (must match order and count)
     * @param readOneUMIs the UMI strings from readOne (must match order and count)
     */
    public void processAllWithBarcodes(List<FastqRecord> readTwoRecords,
                                       List<String> readOneBarcodes,
                                       List<String> readOneUMIs) throws InterruptedException {
        if (barcodeMapper == null) {
            throw new IllegalStateException("Barcode-aware mode not enabled");
        }
        if (readTwoRecords.size() != readOneBarcodes.size() || readTwoRecords.size() != readOneUMIs.size()) {
            throw new IllegalArgumentException("Read counts must match: readTwo=" + readTwoRecords.size() +
                    ", barcodes=" + readOneBarcodes.size() + ", UMIs=" + readOneUMIs.size());
        }

        // Create shared thread-safe deduplication sets that ALL threads will use
        Set<String> sharedSeenGeneUMIs = Collections.newSetFromMap(new ConcurrentHashMap<>());
        Set<String> sharedSeenTxUMIs = Collections.newSetFromMap(new ConcurrentHashMap<>());

        List<Future<UMITrackingGraphTraversal>> futures = new ArrayList<>(numThreads);
        int chunkSize = (readTwoRecords.size() + numThreads - 1) / numThreads;

        for (int i = 0; i < numThreads; i++) {
            final int start = i * chunkSize;
            if (start >= readTwoRecords.size()) break;
            final int end = Math.min(start + chunkSize, readTwoRecords.size());

            futures.add(executor.submit(() -> {
                // Pass shared deduplication sets to each worker
                UMITrackingGraphTraversal worker = new UMITrackingGraphTraversal(idxGraph, tx2gene, 
                                                                                 sharedSeenGeneUMIs, 
                                                                                 sharedSeenTxUMIs);
                for (int j = start; j < end; j++) {
                    String barcode = readOneBarcodes.get(j);
                    String umi = readOneUMIs.get(j);
                    int barcodeInt = barcodeMapper.getId(barcode);
                    worker.processWithBarcode(readTwoRecords.get(j), barcodeInt, umi);
                }
                return worker;
            }));
        }

        List<UMITrackingGraphTraversal> workers = new ArrayList<>(futures.size());
        for (Future<UMITrackingGraphTraversal> f : futures) {
            try {
                workers.add(f.get());
            } catch (ExecutionException e) {
                throw new RuntimeException("Worker task failed", e.getCause());
            }
        }

        mergeBarcodeResults(workers);
        System.out.println("Total ambiguous reads: " + totalAmbigReads);
        System.out.println("Total short reads discarded (< 3 minimizers): " + totalShortReadsDiscarded);
        System.out.println("Reads assigned to transcripts: " + totalReadsAssignedToTranscripts);
        System.out.println("Reads assigned to genes: " + totalReadsAssignedToGenes);
        
        // Print memory usage for diagnostics
        Runtime runtime = Runtime.getRuntime();
        long usedMemory = (runtime.totalMemory() - runtime.freeMemory()) / (1024 * 1024);
        long totalGeneEntries = geneCountsPerBarcode.values().stream()
            .mapToLong(Int2IntOpenHashMap::size).sum();
        System.out.printf("Memory used: %d MB | Gene entries: %d | Barcodes: %d%n", 
                          usedMemory, totalGeneEntries, geneCountsPerBarcode.size());
    }

    /**
     * Process a batch of reads in parallel. Results accumulate across multiple
     * batch calls. Used for memory-efficient streaming of large read sets.
     */
    public void processBatch(List<FastqRecord> batchReads) throws InterruptedException {
        int chunkSize = (batchReads.size() + numThreads - 1) / numThreads;
        List<Future<IndexGraphTraversal>> futures = new ArrayList<>(numThreads);

        for (int i = 0; i < numThreads; i++) {
            final int start = i * chunkSize;
            if (start >= batchReads.size()) break;
            final int end = Math.min(start + chunkSize, batchReads.size());
            final List<FastqRecord> slice = batchReads.subList(start, end);

            futures.add(executor.submit(() -> {
                IndexGraphTraversal worker = new IndexGraphTraversal(idxGraph, tx2gene, writeDetails);
                for (FastqRecord fastqRecord : slice) {
                    worker.process(fastqRecord);
                }
                return worker;
            }));
        }

        List<IndexGraphTraversal> workers = new ArrayList<>(futures.size());
        for (Future<IndexGraphTraversal> f : futures) {
            try {
                workers.add(f.get());
            } catch (ExecutionException e) {
                throw new RuntimeException("Worker task failed", e.getCause());
            }
        }

        mergeResults(workers);
        
        // Collect per-read mapping results in order if tracking details
        if (writeDetails && readMappingIds != null) {
            for (IndexGraphTraversal worker : workers) {
                IntArrayList workerIds = worker.getReadMappingIds();
                BooleanArrayList workerIsGene = worker.getReadMappingIsGene();
                if (workerIds != null && workerIsGene != null) {
                    readMappingIds.addAll(workerIds);
                    readMappingIsGene.addAll(workerIsGene);
                }
            }
        }
        
        System.out.println("Total ambiguous reads: " + totalAmbigReads);
        System.out.println("Total short reads discarded (< 3 minimizers): " + totalShortReadsDiscarded);
        System.out.println("Heuristic phase - transcripts: " + numTxFoundHeuristically);
        System.out.println("Heuristic phase - genes: " + numGenesFoundHeuristically);
        System.out.println("Scoring phase - transcripts: " + numTxFoundWithAlignment);
        System.out.println("Scoring phase - genes: " + numGenesFoundWithAlignment);
    }

    public void shutdown() throws InterruptedException {
        executor.shutdown();
        executor.awaitTermination(1, TimeUnit.HOURS);
    }

    public Int2IntOpenHashMap getGlobalGeneCounts() { return globalGeneCounts; }
    public Int2IntOpenHashMap getGlobalTxCounts()   { return globalTxCounts;   }
    
    public Int2ObjectOpenHashMap<Int2IntOpenHashMap> getGeneCountsPerBarcode() { 
        return geneCountsPerBarcode; 
    }
    
    public Int2ObjectOpenHashMap<Int2IntOpenHashMap> getTxCountsPerBarcode() { 
        return txCountsPerBarcode; 
    }
    
    public BarcodeMapper getBarcodeMapper() { 
        return barcodeMapper; 
    }

    private void mergeResults(List<IndexGraphTraversal> workers) {
        for (IndexGraphTraversal worker : workers) {
            worker.getGeneCounts().int2IntEntrySet().fastForEach(e ->
                    globalGeneCounts.addTo(e.getIntKey(), e.getIntValue()));
            worker.getTxCounts().int2IntEntrySet().fastForEach(e ->
                    globalTxCounts.addTo(e.getIntKey(), e.getIntValue()));
            numTxFoundHeuristically += worker.getNumTxFoundHeuristically();
            numGenesFoundHeuristically += worker.getNumGenesFoundHeuristically();
            numTxFoundWithAlignment += worker.getNumTxFoundWithAlignment();
            numGenesFoundWithAlignment += worker.getNumGenesFoundWithAlignment();

            this.totalAmbigReads += worker.getAmbigReads();
            this.totalShortReadsDiscarded += worker.getShortReadsDiscarded();
            this.totalReadsAssignedToTranscripts += worker.getReadsAssignedToTranscripts();
            this.totalReadsAssignedToGenes += worker.getReadsAssignedToGenes();
        }
    }

    private void mergeBarcodeResults(List<UMITrackingGraphTraversal> workers) {
        for (UMITrackingGraphTraversal worker : workers) {
            // Merge per-barcode gene counts
            for (var barcodeEntry : worker.getGeneCountsPerBarcode().int2ObjectEntrySet()) {
                int barcodeInt = barcodeEntry.getIntKey();
                Int2IntOpenHashMap workerGeneCounts = barcodeEntry.getValue();
                
                Int2IntOpenHashMap barcodeGenes = geneCountsPerBarcode.computeIfAbsent(barcodeInt, 
                    k -> new Int2IntOpenHashMap());
                
                for (var geneEntry : workerGeneCounts.int2IntEntrySet()) {
                    barcodeGenes.addTo(geneEntry.getIntKey(), geneEntry.getIntValue());
                }
            }
            
            // Merge per-barcode transcript counts
            for (var barcodeEntry : worker.getTxCountsPerBarcode().int2ObjectEntrySet()) {
                int barcodeInt = barcodeEntry.getIntKey();
                Int2IntOpenHashMap workerTxCounts = barcodeEntry.getValue();
                
                Int2IntOpenHashMap barcodeTx = txCountsPerBarcode.computeIfAbsent(barcodeInt, 
                    k -> new Int2IntOpenHashMap());
                
                for (var txEntry : workerTxCounts.int2IntEntrySet()) {
                    barcodeTx.addTo(txEntry.getIntKey(), txEntry.getIntValue());
                }
            }
            
            this.totalAmbigReads += worker.getTotalAmbigReads();
            this.totalShortReadsDiscarded += worker.getTotalShortReadsDiscarded();
            this.totalReadsAssignedToTranscripts += worker.getReadsAssignedToTranscripts();
            this.totalReadsAssignedToGenes += worker.getReadsAssignedToGenes();
        }
    }


    public void writeSummary(File summaryFile) throws IOException {
        try (FileWriter writer = new FileWriter(summaryFile)) {
            writer.write("=== miniQuT3 Query Summary ===\n\n");
            writer.write("Total reads processed: " + (totalShortReadsDiscarded + totalReadsAssignedToTranscripts + totalReadsAssignedToGenes + totalAmbigReads) + "\n");
            writer.write("Reads assigned to transcripts: " + totalReadsAssignedToTranscripts + "\n");
            writer.write("Reads assigned to genes: " + totalReadsAssignedToGenes + "\n");
            writer.write("Ambiguous reads: " + totalAmbigReads + "\n");
            writer.write("Short reads discarded: " + totalShortReadsDiscarded + "\n");
            writer.write("\nHeuristic phase - transcripts: " + numTxFoundHeuristically + "\n");
            writer.write("Heuristic phase - genes: " + numGenesFoundHeuristically + "\n");
            writer.write("Scoring phase - transcripts: " + numTxFoundWithAlignment + "\n");
            writer.write("Scoring phase - genes: " + numGenesFoundWithAlignment + "\n");
        }
        System.out.println("Wrote summary to: " + summaryFile.getAbsolutePath());
    }

    /**
     * Write per-read mapping details to a file.
     * Format: one line per read, containing gene name, transcript name, "-" for ambiguous, or "*" for short reads
     */
    public void writeReadMappingDetails(File detailsFile, Int2ObjectOpenHashMap<String> int2GeneString,
                                       Int2ObjectOpenHashMap<String> int2TxString) throws IOException {
        if (readMappingIds == null || readMappingIds.isEmpty()) {
            System.out.println("No read mapping details to write (tracking was disabled or no reads processed)");
            return;
        }
        
        try (java.io.FileWriter writer = new java.io.FileWriter(detailsFile)) {
            for (int i = 0; i < readMappingIds.size(); i++) {
                int id = readMappingIds.getInt(i);
                
                if (id == -2) { // was too short
                    writer.write("*\n");
                } else if (id == -1) { // ambig
                    writer.write("-\n");
                } else { // valids
                    boolean isGene = readMappingIsGene.getBoolean(i);
                    if (isGene) {
                        // Gene ID: look up gene name
                        String geneName = int2GeneString.get(id);
                        if (geneName != null) {
                            writer.write(geneName);
                        } else {
                            writer.write(String.valueOf(id));
                        }
                    } else {
                        // Transcript ID: look up transcript name
                        String txName = int2TxString.get(id);
                        if (txName != null) {
                            writer.write(txName);
                        } else {
                            writer.write(String.valueOf(id));
                        }
                    }
                    writer.write('\n');
                }
            }
        }
        
        System.out.println("Wrote " + readMappingIds.size() + " read mapping results to " + detailsFile.getAbsolutePath());
    }
}