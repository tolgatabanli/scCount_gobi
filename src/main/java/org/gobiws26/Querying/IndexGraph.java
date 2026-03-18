package org.gobiws26.Querying;

import it.unimi.dsi.fastutil.ints.*;
import org.gobiws26.utils.TripletHasher;
import java.util.BitSet;

/**
 * Triplet-based de Bruijn graph for minimizer sequences.
 * Nodes represent triplets of consecutive minimizers [m1, m2, m3] hashed to an int.
 * Each node stores which transcripts contain that triplet sequence.
 */
public class IndexGraph {

    private int size = 0;
    private final Int2ObjectOpenHashMap<Node> tripletHash2Node = new Int2ObjectOpenHashMap<>();

    /**
     * Node: represents a unique triplet of minimizers [m1, m2, m3]
     * Stores transcript IDs in sparse IntOpenHashSet (memory efficient, only stores actual IDs).
     */
    public static class Node {
        // Use IntOpenHashSet (sparse) instead of BitSet (dense) to avoid memory exhaustion
        // with thousands of transcripts
        private final IntOpenHashSet containingTranscripts = new IntOpenHashSet();

        void recordTranscript(int txId) {
            containingTranscripts.add(txId);
        }

        /**
         * Returns the set of transcript IDs for this node.
         * Caller is responsible for converting to BitSet if bitwise operations needed.
         */
        public IntOpenHashSet getTranscripts() {
            return containingTranscripts;
        }
    }

    /**
     * Add a transcript path as triplets to the graph.
     * For a minimizer list [m0, m1, m2, m3, m4, ...], creates triplets:
     * [m0, m1, m2], [m1, m2, m3], [m2, m3, m4], ...
     */
    public void addTxPath(int txId, IntArrayList minimList) {
        // Need at least 3 minimizers to form a triplet
        if (minimList.size() < 3) {
            return;
        }

        for (int i = 0; i + 2 < minimList.size(); i++) {
            int m1 = minimList.getInt(i);
            int m2 = minimList.getInt(i + 1);
            int m3 = minimList.getInt(i + 2);

            int tripletHash = TripletHasher.hash(m1, m2, m3);
            
            Node node = tripletHash2Node.computeIfAbsent(tripletHash, k -> new Node());
            node.recordTranscript(txId);
            size++;
        }
    }

    public Node getNode(int tripletHash) {
        return tripletHash2Node.get(tripletHash);
    }

    public int getSize() { 
        return size; 
    }
}