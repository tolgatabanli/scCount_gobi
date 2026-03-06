package org.gobiws26.Querying;

import it.unimi.dsi.fastutil.ints.*;
import it.unimi.dsi.fastutil.shorts.*;
import java.util.Collection;

// Nodes are minimizers, keep transcripts that contain them.
// Edges are directed, pointing to the next minimizer that come after in a transcript, and are labeled with next minimizer.
public class IndexGraph {
    private int size = 0;
    private final Short2ObjectOpenHashMap<Node> minim2Node = new Short2ObjectOpenHashMap<>();

    public IndexGraph() {}

    public static class Node {
        private final IntOpenHashSet txIds = new IntOpenHashSet();
        private final Short2ObjectOpenHashMap<Node> nextNodes = new Short2ObjectOpenHashMap<>();

        private boolean addTranscript(int txId) {
            return txIds.add(txId);
        }

        private void putNextNode(short nextMinim, Node nextNode) {
            nextNodes.put(nextMinim, nextNode);
        }

        public IntSet getTranscripts() {
            return IntSets.unmodifiable(txIds);
        }

        public Node getNextNeighbor(short minim) {
            return nextNodes.get(minim);
        }

        public Collection<Node> getAllNextNeighbors() {
            return nextNodes.values();
        }
    }

    public Node getNode(short minimizer) {
        return minim2Node.get(minimizer);
    }

    public void addTxPath(int txId, ShortArrayList minimList) {
        Node lastNode = null;
        for (short minim : minimList) {
            Node currentNode = minim2Node.computeIfAbsent(minim, k -> {
                size++;
                return new Node();
            });

            currentNode.addTranscript(txId);
            if (lastNode != null) {
                lastNode.putNextNode(minim, currentNode);
            }
            lastNode = currentNode;
        }
    }

    public int getSize() { return size; }
}
