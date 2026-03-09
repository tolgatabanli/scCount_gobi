package org.gobiws26.Querying;

import htsjdk.samtools.fastq.FastqRecord;
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;

import java.util.List;
import java.util.ArrayList;
import java.util.concurrent.*;

public class ParallelGraphQuery {
    private final IndexGraph idxGraph;
    private final Int2IntMap tx2gene;
    private final int numThreads;
    private final ExecutorService executor;

    Int2IntOpenHashMap globalGeneCounts;
    Int2IntOpenHashMap globalTxCounts;
    private int totalAmbigReads;

    public ParallelGraphQuery(IndexGraph g, Int2IntMap tx2gene, int numThreads) {
        this.idxGraph = g;
        this.tx2gene = tx2gene;

        this.numThreads = Math.min(numThreads, Runtime.getRuntime().availableProcessors());
        this.executor = Executors.newFixedThreadPool(this.numThreads);
        this.globalGeneCounts = new Int2IntOpenHashMap();
        this.globalTxCounts   = new Int2IntOpenHashMap();
    }

    /**
     * Convenience constructor that uses all available hardware threads.
     */
    public ParallelGraphQuery(IndexGraph g, Int2IntMap tx2gene) {
        this(g, tx2gene, Runtime.getRuntime().availableProcessors());
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
                for (int j = 0; j < slice.size(); j++) {   // index loop — no Itr
                    worker.process(slice.get(j));
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

        mergeBatchResults(workers);
        System.out.println("Total ambiguous reads: " + totalAmbigReads);
    }

    /**
     * Process all reads in one call. For backward compatibility or small read sets.
     */
    public void processAll(List<FastqRecord> allReads) throws InterruptedException {
        processBatch(allReads);
        shutdown();
        System.out.println("Total ambiguous reads: " + totalAmbigReads);
    }

    /**
     * Shutdown the executor after all batches have been processed.
     */
    public void shutdown() throws InterruptedException {
        executor.shutdown();
        executor.awaitTermination(1, TimeUnit.HOURS);
    }

    public Int2IntOpenHashMap getGlobalGeneCounts() { return globalGeneCounts; }
    public Int2IntOpenHashMap getGlobalTxCounts()   { return globalTxCounts;   }

    /**
     * Merge batch results into global accumulators.
     *
     * FIX 2: Index-based loop instead of for-each — eliminates the
     * ArrayList$Itr allocation that JMC flagged as the top allocator on the
     * main thread. fastForEach() uses a primitive consumer, avoiding boxing
     * on the inner forEach as well.
     */
    private void mergeBatchResults(List<IndexGraphTraversal> workers) {
        // FIX 2: index loop — no ArrayList$Itr allocated
        for (int i = 0, n = workers.size(); i < n; i++) {
            IndexGraphTraversal worker = workers.get(i);

            // FIX 2: fastForEach uses a primitive Int2IntConsumer — no boxing
            worker.getGeneCounts().int2IntEntrySet().fastForEach(e ->
                    globalGeneCounts.addTo(e.getIntKey(), e.getIntValue()));
            worker.getTxCounts().int2IntEntrySet().fastForEach(e ->
                    globalTxCounts.addTo(e.getIntKey(), e.getIntValue()));

            this.totalAmbigReads += worker.getAmbigReads();
        }
    }

    @Deprecated
    private void mergeResults(List<IndexGraphTraversal> workers) {
        mergeBatchResults(workers);
    }
}