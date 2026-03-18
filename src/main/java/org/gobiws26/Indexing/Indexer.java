package org.gobiws26.Indexing;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import org.gobiws26.genomicstruct.Transcript;
import org.gobiws26.utils.Minimizers;
import org.gobiws26.utils.TranscriptomeFetcher;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Takes an output file and transcriptome to write an index file that maps the following information:
 *  1) gtf_tx_id -> int
 *  2) gtf_gene_id -> int
 *  3) tx -> gene
 *  4) triplet hash -> BitSet of transcript IDs containing that triplet (inverted index)
 */
public class Indexer {
    // inputs
    private final ReferenceSequenceFile refSeqFile;
    private final HashMap<String, Transcript> transcripts;

    // outputs
    private String[] txIdArray;
    private String[] geneIdArray;
    private int[] txToGeneArray;
    private Int2ObjectOpenHashMap<IntArrayList> txIdToMinimizerPath;  // txId -> minimizer sequence


    public Indexer(HashMap<String, Transcript> transcripts, ReferenceSequenceFile refSeqFile) {
        this.transcripts = transcripts;
        this.refSeqFile = refSeqFile;
    }


    public void runIndex() {
        int expectedSize = transcripts.size();

        // temporary, thread-safe structures used during parallel build
        ConcurrentHashMap<String, Integer> tmpTranscriptToIndex = new ConcurrentHashMap<>(expectedSize);
        ConcurrentHashMap<Integer, IntArrayList> tmpMinimizerPaths = new ConcurrentHashMap<>(expectedSize);
        ConcurrentHashMap<Integer, String> idToTxId = new ConcurrentHashMap<>(expectedSize);
        ConcurrentHashMap<Integer, Integer> idToGene = new ConcurrentHashMap<>(expectedSize);
        ConcurrentHashMap<String, Integer> geneToIndex = new ConcurrentHashMap<>();

        AtomicInteger internalTxId = new AtomicInteger(0);
        AtomicInteger internalGeneId = new AtomicInteger(0);

        int nThreads = Runtime.getRuntime().availableProcessors();
        try (ExecutorService pool = Executors.newFixedThreadPool(nThreads)) {

            // Submit one task per transcript
            for (Map.Entry<String, Transcript> txEntry : transcripts.entrySet()) {
                String txId = txEntry.getKey();
                Transcript tx = txEntry.getValue();

                pool.submit(new Worker(txId, tx, refSeqFile, internalTxId, geneToIndex, internalGeneId, tmpTranscriptToIndex, tmpMinimizerPaths, idToTxId, idToGene));
            }

            // wait for completion
            pool.shutdown();
            try {
                boolean ok = pool.awaitTermination(1, TimeUnit.HOURS);
                if (!ok) throw new IllegalStateException("Indexing did not complete in time");
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
                throw new RuntimeException(e);
            }
        }

        // Build final structures in deterministic id order
        int txCount = internalTxId.get();
        this.txIdArray = new String[txCount];
        this.txToGeneArray = new int[txCount];

        for (int id = 0; id < txCount; id++) {
            String tid = idToTxId.get(id);
            Integer g = idToGene.get(id);

            this.txIdArray[id] = tid;
            this.txToGeneArray[id] = (g == null) ? -1 : g;
        }

        // Store minimizer paths per transcript (compact, efficient)
        this.txIdToMinimizerPath = new Int2ObjectOpenHashMap<>(tmpMinimizerPaths);

        // Build geneId array from geneToIndex map
        this.geneIdArray = new String[geneToIndex.size()];
        geneToIndex.forEach((id, idx) -> geneIdArray[idx] = id);
    }

    public Int2ObjectOpenHashMap<IntArrayList> getTxIdToMinimizerPath() {
        if (txIdToMinimizerPath == null) throw new IllegalStateException("Indexer has not been run!");
        return txIdToMinimizerPath;
    }

    public String[] getTxIdArray() {
        return txIdArray;
    }

    public String[] getGeneIdArray () {
        return geneIdArray;
    }

    public int[] getTxToGeneArray() {
        return txToGeneArray;
    }

    // Worker performs per-transcript work: fetch sequence, compute minimizers, and record results
        private record Worker(String txId, Transcript tx, ReferenceSequenceFile refSeqFile, AtomicInteger internalTxId,
                              ConcurrentHashMap<String, Integer> geneToIndex, AtomicInteger internalGeneId,
                              ConcurrentHashMap<String, Integer> tmpTranscriptToIndex,
                              ConcurrentHashMap<Integer, IntArrayList> tmpMinimizerPaths,
                              ConcurrentHashMap<Integer, String> idToTxId,
                              ConcurrentHashMap<Integer, Integer> idToGene) implements Runnable {

        @Override
            public void run() {
                // Use a local fetcher but synchronize on the refSeqFile to be safe
                TranscriptomeFetcher tf = new TranscriptomeFetcher(refSeqFile);
                byte[] seq;
                try {
                    synchronized (refSeqFile) {
                        seq = tf.fetchTranscriptSequenceOfUTRPlus(tx, 500);
                    }
                } catch (RuntimeException e) {
                    // Fetch failure: skip this transcript
                    return;
                }

                IntArrayList minimizerPath;
                try {
                    minimizerPath = Minimizers.of(seq, null);
                } catch (IllegalArgumentException e) {
                    // invalid sequence, skip
                    return;
                }

                int id = internalTxId.getAndIncrement();

                tmpTranscriptToIndex.put(txId, id);
                idToTxId.put(id, txId);
                tmpMinimizerPaths.put(id, minimizerPath);

                String geneId = tx.getGeneId();
                int geneIdx = geneToIndex.computeIfAbsent(geneId, k -> internalGeneId.getAndIncrement());
                idToGene.put(id, geneIdx);
            }
        }
}
