package org.gobiws26.Querying;

import htsjdk.samtools.fastq.FastqRecord;
import it.unimi.dsi.fastutil.ints.*;
import it.unimi.dsi.fastutil.booleans.BooleanArrayList;
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
    private int readsAssignedToTranscripts = 0;
    private int readsAssignedToGenes = 0;
    
    // Per-read mapping tracking: integers with companion boolean list
    // ID encoding: -1 = ambiguous, -2 = short read, >=0 = actual ID
    private final IntArrayList readMappingIds;
    private final BooleanArrayList readMappingIsGene;  // true if ID is geneId, false if txId
    private final boolean trackDetails;

    // Getters for merging threads (unchanged public API)
    public Int2IntOpenHashMap getGeneCounts() { return geneCounts; }
    public Int2IntOpenHashMap getTxCounts()   { return txCounts;   }
    public int getAmbigReads()                { return ambigReads; }
    public int getShortReadsDiscarded()       { return shortReadsDiscarded; }
    public int getReadsAssignedToTranscripts() { return readsAssignedToTranscripts; }
    public int getReadsAssignedToGenes()      { return readsAssignedToGenes; }
    public IntArrayList getReadMappingIds() { return readMappingIds; }
    public BooleanArrayList getReadMappingIsGene() { return readMappingIsGene; }

    int txCount;
    int[] scores; // should be reset
    private final BitSet candidatesBitSet = new BitSet();
    private final BitSet bestTxsBitSet = new BitSet();
    private final BitSet unionCandidates = new BitSet();
    
    public IndexGraphTraversal(IndexGraph g, Int2IntMap tx2gene) {
        this.g = g;
        this.tx2geneMapping = tx2gene;

        txCount = tx2geneMapping.size();
        scores   = new int[txCount];
        
        this.trackDetails = false;
        this.readMappingIds = null;
        this.readMappingIsGene = null;
    }
    
    /**
     * Constructor with optional per-read mapping tracking.
     */
    public IndexGraphTraversal(IndexGraph g, Int2IntMap tx2gene, boolean trackDetails) {
        this.g = g;
        this.tx2geneMapping = tx2gene;

        txCount = tx2geneMapping.size();
        scores   = new int[txCount];
        
        this.trackDetails = trackDetails;
        this.readMappingIds = trackDetails ? new IntArrayList() : null;
        this.readMappingIsGene = trackDetails ? new BooleanArrayList() : null;
    }

    // for benchmarking/logging - separated by level (transcript vs gene) and phase
    private int numTxFoundHeuristically = 0;
    private int numGenesFoundHeuristically = 0;
    private int numTxFoundWithAlignment = 0;
    private int numGenesFoundWithAlignment = 0;

    public int getNumTxFoundHeuristically() { return numTxFoundHeuristically; }
    public int getNumGenesFoundHeuristically() { return numGenesFoundHeuristically; }
    public int getNumTxFoundWithAlignment() { return numTxFoundWithAlignment; }
    public int getNumGenesFoundWithAlignment() { return numGenesFoundWithAlignment; }


    public void process(FastqRecord read) {
        IntArrayList minims = Minimizers.of(read.getReadBases(), read.getBaseQualities());
        
        // Discard reads with < 3 minimizers (cannot form a triplet)
        if (minims.size() < 3) {
            shortReadsDiscarded++;
            if (trackDetails) {
                readMappingIds.add(-2);      // -2 encodes "short read"
                readMappingIsGene.add(false); // dummy value, not used for short reads
            }
            return;
        }

        // Convert minimizers to triplet nodes
        IntArrayList triplets = new IntArrayList();
        for (int i = 0; i + 2 < minims.size(); i++) {
            int m1 = minims.getInt(i);
            int m2 = minims.getInt(i + 1);
            int m3 = minims.getInt(i + 2);
            triplets.add(TripletHasher.hash(m1, m2, m3));
        }

        // Phase 1: Heuristic - intersect transcripts from all triplet nodes
        candidatesBitSet.clear();
        boolean candidatesEmpty = false;
        int phase1GeneId = -1;

        for (int i = 0; i < triplets.size(); i++) {
            int tripletHash = triplets.getInt(i);
            IndexGraph.Node node = g.getNode(tripletHash);
            
            if (node == null) {
                candidatesBitSet.clear();
                candidatesEmpty = true;
                break;
            }

            // Convert sparse IntOpenHashSet to BitSet for bitwise operations
            IntOpenHashSet nodeTxsSet = node.getTranscripts();
            BitSet nodeTxs = new BitSet();
            for (int txId : nodeTxsSet) {
                nodeTxs.set(txId);
            }
            
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
            phase1GeneId = getGeneIfAllIsoforms(candidatesBitSet);
            if (phase1GeneId != -1) break;
        }

        if (!candidatesEmpty && !candidatesBitSet.isEmpty()) {
            if (candidatesBitSet.cardinality() == 1) {
                // Unambiguous transcript
                int txId = candidatesBitSet.nextSetBit(0);
                txCounts.addTo(txId, 1);
                geneCounts.addTo(tx2geneMapping.get(txId), 1);
                readsAssignedToTranscripts++;
                numTxFoundHeuristically++;
                if (trackDetails) {
                    readMappingIds.add(txId);
                    readMappingIsGene.add(false);
                }
                return;
            }
            // If multiple transcripts found but same gene, fall through to Phase 2
            // for consecutiveness scoring to achieve transcript-level precision
        }

        // Phase 2: Union-based quasi-alignment with consecutiveness scoring
        // Compute union of all triplet nodes' transcripts
        unionCandidates.clear();
        for (int i = 0; i < triplets.size(); i++) {
            int tripletHash = triplets.getInt(i);
            IndexGraph.Node node = g.getNode(tripletHash);
            if (node != null) {
                // Convert sparse IntOpenHashSet to BitSet for bitwise operations
                IntOpenHashSet nodeTxsSet = node.getTranscripts();
                for (int txId : nodeTxsSet) {
                    unionCandidates.set(txId);
                }
            }
        }

        if (unionCandidates.isEmpty()) {
            ambigReads++;
            if (trackDetails) {
                readMappingIds.add(-1);      // -1 encodes "ambiguous"
                readMappingIsGene.add(false); // dummy value, not used for ambiguous
            }
            return;
        }

        // Score based on consecutiveness: transcript gains +1 for each consecutive triplet pair it appears in
        Arrays.fill(scores, 0);

        for (int i = 0; i + 1 < triplets.size(); i++) {
            int tripletHash = triplets.getInt(i);
            int nextTripletHash = triplets.getInt(i + 1);

            IndexGraph.Node node = g.getNode(tripletHash);
            IndexGraph.Node nextNode = g.getNode(nextTripletHash);

            if (node == null || nextNode == null) continue;
            
            // Use sparse sets directly for membership checking
            IntOpenHashSet nodeTxsSet = node.getTranscripts();
            IntOpenHashSet nextTxsSet = nextNode.getTranscripts();

            for (int txId = unionCandidates.nextSetBit(0); txId >= 0; txId = unionCandidates.nextSetBit(txId + 1)) {
                if (nodeTxsSet.contains(txId) && nextTxsSet.contains(txId)) {
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
            if (trackDetails) {
                readMappingIds.add(-1);      // -1 encodes "ambiguous"
                readMappingIsGene.add(false); // dummy value, not used for ambiguous
            }
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
            readsAssignedToTranscripts++;
            if (phase1GeneId != -1) {
                numTxFoundHeuristically++;  // Refined from Phase 1 gene-level candidate
            } else {
                numTxFoundWithAlignment++;  // Found entirely by Phase 2
            }
            if (trackDetails) {
                readMappingIds.add(txId);
                readMappingIsGene.add(false);
            }
        } else if (geneId != -1) {
            geneCounts.addTo(geneId, 1);
            readsAssignedToGenes++;
            if (phase1GeneId != -1) {
                numGenesFoundHeuristically++;  // Confirmed by Phase 1 and Phase 2
            } else {
                numGenesFoundWithAlignment++;  // Found entirely by Phase 2
            }
            if (trackDetails) {
                readMappingIds.add(geneId);
                readMappingIsGene.add(true);
            }
        } else {
            ambigReads++;
            if (trackDetails) {
                readMappingIds.add(-1);      // -1 encodes "ambiguous"
                readMappingIsGene.add(false); // dummy value, not used for ambiguous
            }
            return;
        }
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
