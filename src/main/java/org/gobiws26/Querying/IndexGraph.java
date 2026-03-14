package org.gobiws26.Querying;

import it.unimi.dsi.fastutil.ints.*;
import it.unimi.dsi.fastutil.shorts.ShortArrayList;
import org.gobiws26.utils.TripletHasher;

import java.util.BitSet;

/**
 * Triplet-based de Bruijn graph for minimizer sequences.
 * 
 * Nodes represent triplets of consecutive minimizers [m1, m2, m3] hashed to an int.
 * Edges are implicit: triplet [a, b, c] transitions to [b, c, d] when the 4th minimizer is d.
 * Each node stores which transcripts contain that triplet sequence.
 */
public class IndexGraph {

    private int size = 0;
    private final Int2ObjectOpenHashMap<Node> tripletHash2Node = new Int2ObjectOpenHashMap<>();

    // -------------------------------------------------------------------------
    // Node: represents a unique triplet of minimizers [m1, m2, m3]
    //   Stores BitSet of transcripts containing this triplet.
    //   Edges are implicit via triplet hash computation.
    // -------------------------------------------------------------------------
    public static class Node {
        private final BitSet containingTranscripts = new BitSet();

        void recordTranscript(int txId) {
            containingTranscripts.set(txId);
        }

        public BitSet getTranscriptsBitSet() {
            return containingTranscripts;
        }
    }

    // construction
    /**
     * Add a transcript path as triplets to the graph.
     * For a minimizer list [m0, m1, m2, m3, m4, ...], creates triplets:
     * [m0, m1, m2], [m1, m2, m3], [m2, m3, m4], ...
     */
    public void addTxPath(int txId, ShortArrayList minimList) {
        // Need at least 3 minimizers to form a triplet
        if (minimList.size() < 3) {
            return;
        }

        for (int i = 0; i + 2 < minimList.size(); i++) {
            short m1 = minimList.getShort(i);
            short m2 = minimList.getShort(i + 1);
            short m3 = minimList.getShort(i + 2);

            int tripletHash = TripletHasher.hash(m1, m2, m3);
            
            Node node = tripletHash2Node.computeIfAbsent(tripletHash, k -> new Node());
            node.recordTranscript(txId);
            size++;
        }
    }

    public Node getNode(int tripletHash) {
        return tripletHash2Node.get(tripletHash);
    }

    public Int2ObjectOpenHashMap<Node> getTripletHash2Node() {
        return tripletHash2Node;
    }

    public int getSize() { 
        return size; 
    }

    public void freezeGraph() {
        // No-op in new system; kept for API compatibility
    }
}