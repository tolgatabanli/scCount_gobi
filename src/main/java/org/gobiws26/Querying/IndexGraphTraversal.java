package org.gobiws26.Querying;

import htsjdk.samtools.fastq.FastqRecord;
import it.unimi.dsi.fastutil.ints.*;
import it.unimi.dsi.fastutil.shorts.ShortArrayList;
import org.gobiws26.utils.Minimizers;
import org.gobiws26.utils.TripletHasher;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.List;


public class IndexGraphTraversal {

    // Per-read outputs
    private final IndexGraph g;
    private final Int2IntMap tx2geneMapping;

    private final Int2IntOpenHashMap geneCounts = new Int2IntOpenHashMap();
    private final Int2IntOpenHashMap txCounts   = new Int2IntOpenHashMap();
    private int ambigReads = 0;
    private int shortReadsDiscarded = 0;

    // Getters for merging threads (unchanged public API)
    public Int2IntOpenHashMap getGeneCounts() { return geneCounts; }
    public Int2IntOpenHashMap getTxCounts()   { return txCounts;   }
    public int getAmbigReads()                { return ambigReads; }
    public int getShortReadsDiscarded()       { return shortReadsDiscarded; }

    int txCount;
    int[] scores; // should be reset
    private final BitSet candidatesBitSet = new BitSet();
    private final BitSet validatedBitSet = new BitSet();
    private final BitSet bestTxsBitSet = new BitSet();
    private final BitSet unionCandidates = new BitSet();
    
    public IndexGraphTraversal(IndexGraph g, Int2IntMap tx2gene) {
        this.g = g;
        this.tx2geneMapping = tx2gene;

        txCount = tx2geneMapping.size();
        scores   = new int[txCount];
    }

    // for benchmarking/logging
    private int numReadsFoundHeuristically = 0;

    public int getNumReadsFoundHeuristically() {
        return numReadsFoundHeuristically;
    }

    private int numReadsFoundWithAlignment = 0;

    public int getNumReadsFoundWithAlignment() {
        return numReadsFoundWithAlignment;
    }

    private int numReadsEmptyCandidates = 0;

    public int getNumReadsEmptyCandidates() {
        return numReadsEmptyCandidates;
    }


    public void process(FastqRecord read) {
        ShortArrayList minims = Minimizers.of(read.getReadBases(), read.getBaseQualities());
        
        // Discard reads with < 3 minimizers (cannot form a triplet)
        if (minims.size() < 3) {
            shortReadsDiscarded++;
            return;
        }

        // Convert minimizers to triplet nodes
        List<Integer> triplets = new ArrayList<>();
        for (int i = 0; i + 2 < minims.size(); i++) {
            short m1 = minims.getShort(i);
            short m2 = minims.getShort(i + 1);
            short m3 = minims.getShort(i + 2);
            triplets.add(TripletHasher.hash(m1, m2, m3));
        }

        // Phase 1: Heuristic - intersect transcripts from all triplet nodes
        candidatesBitSet.clear();
        boolean candidatesEmpty = false;

        for (int i = 0; i < triplets.size(); i++) {
            int tripletHash = triplets.get(i);
            IndexGraph.Node node = g.getNode(tripletHash);
            
            if (node == null) {
                candidatesBitSet.clear();
                candidatesEmpty = true;
                break;
            }

            BitSet nodeTxs = node.getTranscriptsBitSet();
            if (i == 0) {
                candidatesBitSet.or(nodeTxs);
            } else {
                candidatesBitSet.and(nodeTxs); // retainAll / intersection
            }

            if (candidatesBitSet.isEmpty()) {
                candidatesEmpty = true;
                break;
            }

            // Early stop: unique gene found
            if (getGeneIfAllIsoforms(candidatesBitSet) != -1) break;
        }

        if (!candidatesEmpty && !candidatesBitSet.isEmpty()) {
            int geneId = getGeneIfAllIsoforms(candidatesBitSet);
            if (candidatesBitSet.cardinality() == 1) {
                // Unambiguous transcript
                int txId = candidatesBitSet.nextSetBit(0);
                txCounts.addTo(txId, 1);
                geneCounts.addTo(tx2geneMapping.get(txId), 1);
                numReadsFoundHeuristically++;
                return;
            } else if (geneId != -1) {
                // Unambiguous gene
                geneCounts.addTo(geneId, 1);
                numReadsFoundHeuristically++;
                return;
            }
        } else {
            numReadsEmptyCandidates++;
        }

        // Phase 2: Union-based quasi-alignment with consecutiveness scoring
        // Compute union of all triplet nodes' transcripts
        unionCandidates.clear();
        for (int tripletHash : triplets) {
            IndexGraph.Node node = g.getNode(tripletHash);
            if (node != null) {
                unionCandidates.or(node.getTranscriptsBitSet());
            }
        }

        if (unionCandidates.isEmpty()) {
            ambigReads++;
            return;
        }

        // Score based on consecutiveness: transcript gains +1 for each consecutive triplet pair it appears in
        Arrays.fill(scores, 0);

        for (int i = 0; i + 1 < triplets.size(); i++) {
            int tripletHash = triplets.get(i);
            int nextTripletHash = triplets.get(i + 1);

            IndexGraph.Node node = g.getNode(tripletHash);
            IndexGraph.Node nextNode = g.getNode(nextTripletHash);

            if (node == null || nextNode == null) continue;

            BitSet nodeTxs = node.getTranscriptsBitSet();
            BitSet nextNodeTxs = nextNode.getTranscriptsBitSet();

            for (int txId = unionCandidates.nextSetBit(0); txId >= 0; txId = unionCandidates.nextSetBit(txId + 1)) {
                if (nodeTxs.get(txId) && nextNodeTxs.get(txId)) {
                    scores[txId]++; // Consecutive appearance: +1 point
                }
            }
        }

        // Find maximum score among union candidates
        int maxScore = Integer.MIN_VALUE;
        for (int txId = unionCandidates.nextSetBit(0); txId >= 0; txId = unionCandidates.nextSetBit(txId + 1)) {
            if (scores[txId] > maxScore) {
                maxScore = scores[txId];
            }
        }

        if (maxScore < 0) {
            ambigReads++;
            return;
        }

        // Choose transcripts with best score
        bestTxsBitSet.clear();
        for (int txId = unionCandidates.nextSetBit(0); txId >= 0; txId = unionCandidates.nextSetBit(txId + 1)) {
            if (scores[txId] == maxScore) {
                bestTxsBitSet.set(txId);
            }
        }

        int geneId = getGeneIfAllIsoforms(bestTxsBitSet);
        int bestCount = bestTxsBitSet.cardinality();
        
        if (bestCount == 1) {
            int txId = bestTxsBitSet.nextSetBit(0);
            txCounts.addTo(txId, 1);
            geneCounts.addTo(tx2geneMapping.get(txId), 1);
        } else if (geneId != -1) {
            geneCounts.addTo(geneId, 1);
        } else {
            // ambiguous across genes -> do nothing
            ambigReads++;
            return;
        }
        
        numReadsFoundWithAlignment++;
    }

    private int getGeneIfAllIsoforms(BitSet txSet) {
        if (txSet.isEmpty()) return -1;
        int geneId = -1;
        for (int tx = txSet.nextSetBit(0); tx >= 0; tx = txSet.nextSetBit(tx + 1)) {
            int g = tx2geneMapping.get(tx);
            if (geneId == -1) geneId = g;
            else if (geneId != g) return -1;
        }
        return geneId;
    }
}
