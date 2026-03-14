package org.gobiws26.Querying;

import htsjdk.samtools.fastq.FastqRecord;
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;

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

    public ParallelGraphQuery(IndexGraph g, Int2IntMap tx2gene, int numThreads) {
        this.idxGraph = g;
        this.tx2gene = tx2gene;
        this.numThreads = Math.min(numThreads, Runtime.getRuntime().availableProcessors());
        this.executor = Executors.newFixedThreadPool(this.numThreads);
        this.globalGeneCounts = new Int2IntOpenHashMap();
        this.globalTxCounts   = new Int2IntOpenHashMap();
        this.barcodeMapper = null;
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
        if (withBarcodes) {
            this.geneCountsPerBarcode = new Int2ObjectOpenHashMap<>();
            this.txCountsPerBarcode = new Int2ObjectOpenHashMap<>();
        }
    }

    /**
     * Process reads without barcode information (legacy mode).
     */
    public void processAll(List<FastqRecord> allReads) throws InterruptedException {
        List<Future<IndexGraphTraversal>> futures = new ArrayList<>(numThreads);
        int chunkSize = (allReads.size() + numThreads - 1) / numThreads;

        for (int i = 0; i < numThreads; i++) {
            final int start = i * chunkSize;
            if (start >= allReads.size()) break;
            final int end = Math.min(start + chunkSize, allReads.size());
            final List<FastqRecord> slice = allReads.subList(start, end);

            futures.add(executor.submit(() -> {
                IndexGraphTraversal worker = new IndexGraphTraversal(idxGraph, tx2gene);
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
        System.out.println("Total ambiguous reads: " + totalAmbigReads);
        System.out.println("Total short reads discarded (< 3 minimizers): " + totalShortReadsDiscarded);
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
        
        // Print memory usage for diagnostics
        Runtime runtime = Runtime.getRuntime();
        long usedMemory = (runtime.totalMemory() - runtime.freeMemory()) / (1024 * 1024);
        long totalGeneEntries = geneCountsPerBarcode.values().stream()
            .mapToLong(m -> m.size()).sum();
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
                IndexGraphTraversal worker = new IndexGraphTraversal(idxGraph, tx2gene);
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
        System.out.println("Total ambiguous reads: " + totalAmbigReads);
        System.out.println("Total short reads discarded (< 3 minimizers): " + totalShortReadsDiscarded);
        System.out.println("numReadsFoundHeuristically: " + numReadsFoundHeuristically);
        System.out.println("numReadsFoundWithAlignment: " + numReadsFoundWithAlignment);
        System.out.println("numReadsEmptyCandidates: " + numReadsEmptyCandidates);
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

    public int getTotalShortReadsDiscarded() {
        return totalShortReadsDiscarded;
    }

    private void mergeResults(List<IndexGraphTraversal> workers) {
        for (IndexGraphTraversal worker : workers) {
            worker.getGeneCounts().int2IntEntrySet().fastForEach(e ->
                    globalGeneCounts.addTo(e.getIntKey(), e.getIntValue()));
            worker.getTxCounts().int2IntEntrySet().fastForEach(e ->
                    globalTxCounts.addTo(e.getIntKey(), e.getIntValue()));
            numReadsFoundHeuristically += worker.getNumReadsFoundHeuristically();
            numReadsFoundWithAlignment += worker.getNumReadsFoundWithAlignment();
            numReadsEmptyCandidates += worker.getNumReadsEmptyCandidates();

            this.totalAmbigReads += worker.getAmbigReads();
            this.totalShortReadsDiscarded += worker.getShortReadsDiscarded();
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
        }
    }

    private int numReadsFoundHeuristically = 0;
    private int numReadsFoundWithAlignment = 0;
    private int numReadsEmptyCandidates = 0;
}