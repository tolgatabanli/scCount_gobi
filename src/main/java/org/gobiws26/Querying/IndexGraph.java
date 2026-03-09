package org.gobiws26.Querying;

import it.unimi.dsi.fastutil.ints.*;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.fastutil.shorts.*;

import java.util.*;

/**
 * Graph where:
 *   - Nodes are minimizers (de-duplicated across all transcripts)
 *   - Edges are directed (prevMinim -> nextMinim) and carry the set of
 *     (txId, positionInTx) pairs that use this transition.
 * The position disambiguates which occurrence
 *   of the transition is being referred to, enabling exact path reconstruction.
 */
public class IndexGraph {

    private int size = 0;
    private final Short2ObjectOpenHashMap<Node> minim2Node = new Short2ObjectOpenHashMap<>();
    private Int2IntMap tx2minimCount;


    // Edge: directed transition between two minimizers
    //   Carries every (txId, position) pair that traverses this edge.
    //   position = index in the transcript of the SOURCE minimizer of this step.
    public static class Edge {
        private final Node target;

        // Packed array: upper 32 bits = txId, lower 32 bits = position
        private final LongArrayList occurrences = new LongArrayList();
        private final BitSet txIdIndex = new BitSet();

        private boolean isFrozen = false;

        Edge(Node target) {
            this.target = target;
        }

        void addOccurrence(int txId, int positionOfSource) {
            // txId is in first 32 bits
            long packed = ((long) txId << 32) | (positionOfSource & 0xFFFFFFFFL);

            occurrences.add(packed);
            txIdIndex.set(txId);
        }

        public BitSet getTxIdIndex() { return txIdIndex; }


        /**
         * Called after construction, before query phase.
         * Sorts the occurrences for binary search.
         */
        public void freeze() {
            if (!isFrozen) {
                // Because txId is in the upper 32 bits, sorting this array automatically
                // sorts primarily by txId, and secondarily by position.
                Arrays.sort(occurrences.elements(), 0, occurrences.size());

                // Optional: Free up unused backing array memory now that size is final
                occurrences.trim();
                isFrozen = true;
            }
        }

        /**
         * Binary search for a (txId, position) pair.
         */
        public boolean hasExactOccurrence(int targetTx, int targetPos) {
            if (!isFrozen) {
                throw new IllegalStateException("Graph must be frozen before searching!");
            }

            // Re-pack the search target into a long
            long searchKey = ((long) targetTx << 32) | (targetPos & 0xFFFFFFFFL);

            // Native binary search directly on the primitive backing array
            int index = Arrays.binarySearch(occurrences.elements(), 0, occurrences.size(), searchKey);

            return index >= 0; // If index is >= 0, the exact (txId, pos) exists
        }

        /**
         * For Step 2 of Read Mapping: Gets all starting positions for a specific transcript.
         * Extracts the lower 32 bits (positions) for all matching txIds.
         */
        public IntArrayList getPositionsForTx(int txId) {
            if (!isFrozen) {
                throw new IllegalStateException("Graph must be frozen before querying!");
            }

            IntArrayList positions = new IntArrayList();

            // Create a search key using the target txId and position 0
            long searchKey = ((long) txId << 32);

            // Binary search to find the start of the txId block
            int index = Arrays.binarySearch(occurrences.elements(), 0, occurrences.size(), searchKey);

            // If exact match for pos 0 isn't found, Arrays.binarySearch returns (-(insertion point) - 1).
            if (index < 0) {
                index = -index - 1;
            }

            // Iterate forward to extract all positions belonging to this txId
            while (index < occurrences.size()) {
                long packed = occurrences.getLong(index);

                // Shift right by 32 to unpack the current txId
                int currentTx = (int) (packed >>> 32);

                // If we've moved past our target transcript block, stop searching
                if (currentTx != txId) {
                    break;
                }

                // Cast to int to extract the lower 32 bits (the position)
                positions.add((int) packed);
                index++;
            }

            return positions;
        }
    }

    // -------------------------------------------------------------------------
    // Node: one unique minimizer value
    //   Stores transcripts (for fast node-level filtering).
    //   Outgoing edges are keyed by (nextMinim, txId) to support the case where
    //   the same next-minimizer is reached via different transcripts.
    //   But two occurrences in the SAME transcript are differentiated by position.
    // -------------------------------------------------------------------------
    public static class Node {
        private final IntOpenHashSet txIds = new IntOpenHashSet();

        // Edge has (txId, position)
        private final Short2ObjectOpenHashMap<Edge> outEdges = new Short2ObjectOpenHashMap<>();

        void recordTranscript(int txId) {
            txIds.add(txId);
        }

        /**
         * Add a directed edge to nextNode (keyed by nextMinim) and record the
         * occurrence (txId, positionOfThisNode) on that edge.
         */
        void addOutEdge(short nextMinim, Node nextNode, int txId, int positionOfSource) {
            Edge edge = outEdges.computeIfAbsent(nextMinim, k -> new Edge(nextNode));
            edge.addOccurrence(txId, positionOfSource);
        }

        public IntSet getTranscripts() {
            return IntSets.unmodifiable(txIds);
        }

        public Short2ObjectOpenHashMap<Edge> getOutEdges() {
            return outEdges;
        }

        public Edge getEdgeTo(short nextMinim) {
            return outEdges.get(nextMinim);
        }
    }

    // construction
    public void addTxPath(int txId, ShortArrayList minimList) {
        Node lastNode = null;

        for (int pos = 0; pos < minimList.size(); pos++) {
            short minim = minimList.getShort(pos);

            Node currentNode = minim2Node.computeIfAbsent(minim, k -> {
                size++;
                return new Node();
            });
            currentNode.recordTranscript(txId);

            if (lastNode != null) {
                // pos-1 is the position of the source node of this edge
                lastNode.addOutEdge(minim, currentNode, txId, pos - 1);
            }

            lastNode = currentNode;
        }
    }

    public Node getNode(short minimizer) {
        return minim2Node.get(minimizer);
    }

    public int getSize() { return size; }

    public void setTx2minimCount(Int2IntMap tx2minimCount) {
        this.tx2minimCount = tx2minimCount;
    }

    public int getMinimCountOfTranscript(int txId) {
        return tx2minimCount.get(txId);
    }

    public void freezeGraph() {
        for (Node node : minim2Node.values()) {
            for (Edge edge : node.getOutEdges().values()) {
                edge.freeze();
            }
        }
    }
}


// Read Mapping:
//
// Given a read encoded in minimizers `[m0, m1, m2, ...]`, do:
//  1. Look up Node(m0)
//  2. Get edge Node(m0) -> Node(m1), check its occurrences for candidate (txId, pos)
//  3. For each candidate, verify the next edge's occurrences contain (txId, pos+1)
//  4. Continue chaining — the transcript(s) where ALL match consecutively are hits