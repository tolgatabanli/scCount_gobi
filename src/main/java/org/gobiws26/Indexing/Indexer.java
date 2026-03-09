package org.gobiws26.Indexing;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.shorts.*;
import org.gobiws26.genomicstruct.Transcript;
import org.gobiws26.utils.Minimizers;
import org.gobiws26.utils.TranscriptomeFetcher;
import java.util.*;

/**
 * Takes an output file and transcriptome to write an index file that maps the following information:
 *  1) gtf_tx_id -> int
 *  2) gtf_gene_id -> int
 *  3) tx -> gene
 *  4) tx (int) to an ordered list of minimizers (ShortArrayList)
 */
public class Indexer {
    // inputs
    private ReferenceSequenceFile refSeqFile;
    private HashMap<String, Transcript> transcripts;

    // outputs
    private HashMap<String, Integer> transcriptToIndex; // our internal IDs
    private String[] txIdArray;
    private String[] geneIdArray;
    private int[] txToGeneArray;
    private Int2ObjectOpenHashMap<ShortArrayList> transcriptToMinimizerPath;


    public Indexer(HashMap<String, Transcript> transcripts, ReferenceSequenceFile refSeqFile) {
        this.transcripts = transcripts;
        this.refSeqFile = refSeqFile;
    }


    public void runIndex() {
        int expectedSize = transcripts.size();

        // minimize rehashing
        this.transcriptToMinimizerPath = new Int2ObjectOpenHashMap<>(expectedSize);
        this.transcriptToIndex = new HashMap<>(expectedSize);

        HashMap<String, Integer> geneToIndex = new HashMap<>();
        List<String> validTxIds = new ArrayList<>();
        List<Integer> txToGeneList = new ArrayList<>();

        TranscriptomeFetcher tf = new TranscriptomeFetcher(refSeqFile);
        int internalId = 0;
        int geneCounter = 0;

        for (Map.Entry<String, Transcript> txEntry : transcripts.entrySet()) {
            String txId = txEntry.getKey();
            Transcript tx = txEntry.getValue();

            // 1. Validate the sequence FIRST
            byte[] seq = tf.fetchTranscriptSequenceOf(tx, 500); // TODO: Magic number
            ShortArrayList minimizerPath;
            try {
                minimizerPath = Minimizers.of(seq, null);
            } catch (IllegalArgumentException e) {
                // Log and skip entirely
                continue;
            }

            // 2. Now that we know it's valid, assign IDs
            this.transcriptToIndex.put(txId, internalId);
            this.transcriptToMinimizerPath.put(internalId, minimizerPath);
            validTxIds.add(txId);

            // 3. Handle Gene Mapping
            String geneId = tx.getGeneId();
            int finalGeneCounter = geneCounter;
            int geneIdx = geneToIndex.computeIfAbsent(geneId, k -> finalGeneCounter);
            if (geneIdx == geneCounter) geneCounter++;

            txToGeneList.add(geneIdx);

            internalId++;
        }

        this.txIdArray = validTxIds.toArray(new String[0]);
        this.txToGeneArray = txToGeneList.stream().mapToInt(i -> i).toArray();
        this.geneIdArray = new String[geneToIndex.size()];
        geneToIndex.forEach((id, idx) -> geneIdArray[idx] = id);
    }

    public HashMap<String, Integer> getTranscriptToIndex() {
        if  (transcriptToIndex == null) throw new IllegalStateException("Indexer has not been run!");
        return transcriptToIndex;
    }

    public Int2ObjectOpenHashMap<ShortArrayList> getTranscriptToMinimizerPath() {
        if (transcriptToMinimizerPath == null) throw new IllegalStateException("Indexer has not been run!");
        return transcriptToMinimizerPath;
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
}
