package org.gobiws26.Querying;

import it.unimi.dsi.fastutil.ints.*;
import it.unimi.dsi.fastutil.shorts.*;

import java.util.*;

/**
 * Graph:
 *   - Nodes -> minimizers (de-duplicated across all transcripts)
 *   - Edges: directed (prevMinim -> nextMinim) and carry the set of
 *     (txId, positionInTx) pairs that use this transition.
 * The position disambiguates which occurrence of the transition is being referred to.
 */
public class IndexGraph {

    private int size = 0;
    private final Short2ObjectOpenHashMap<Node> minim2Node = new Short2ObjectOpenHashMap<>();

    // Edge: directed transition between two minimizers
    //   Carries every (txId, position) pair that traverses this edge.
    public static class Edge {
        private final Node target;
        private final BitSet txIdIndex = new BitSet();

        private final Int2ObjectOpenHashMap<ShortArrayList> positionsByTx = new Int2ObjectOpenHashMap<>();

        Edge(Node target) { this.target = target; }

        void addOccurrence(int txId, int pos) {
            positionsByTx.computeIfAbsent(txId, k -> new ShortArrayList()).add((short) pos);
            txIdIndex.set(txId);
        }

        public void freeze() {
            positionsByTx.values().forEach(ShortArrayList::trim);
        }

        public short[] getPositionsForTx(int txId) {
            ShortArrayList list = positionsByTx.get(txId);
            return list == null ? new short[0] : list.elements();
        }

        public BitSet getTxIdIndex() { return txIdIndex; }
    }

    // -------------------------------------------------------------------------
    // Node: one unique minimizer value
    //   Stores transcripts (for fast node-level filtering).
    //   Outgoing edges are keyed by (nextMinim, txId)
    // -------------------------------------------------------------------------
    public static class Node {
        private final BitSet containingTranscripts = new BitSet();

        // Edge has (txId, position)
        private final Short2ObjectOpenHashMap<Edge> outEdges = new Short2ObjectOpenHashMap<>();

        void recordTranscript(int txId) {
            containingTranscripts.set(txId);
        }


        void addOutEdge(short nextMinim, Node nextNode, int txId, int positionOfSource) {
            Edge edge = outEdges.computeIfAbsent(nextMinim, k -> new Edge(nextNode));
            edge.addOccurrence(txId, positionOfSource);
        }

        public BitSet getTranscriptsBitSet() {
            return containingTranscripts;
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
                // pos-1 is the position of the source node of the edge
                lastNode.addOutEdge(minim, currentNode, txId, pos - 1);
            }

            lastNode = currentNode;
        }
    }

    public Node getNode(short minimizer) {
        return minim2Node.get(minimizer);
    }

    public int getSize() { return size; }

    public void freezeGraph() {
        for (Node node : minim2Node.values()) {
            for (Edge edge : node.outEdges.values()) {
                edge.freeze();
            }
        }
    }
}