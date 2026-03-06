package org.gobiws26.Querying;

import it.unimi.dsi.fastutil.ints.*;
import it.unimi.dsi.fastutil.shorts.*;
import org.gobiws26.utils.Minimizers;

public class IndexGraphTraversal {
    private final IndexGraph g;

    private Int2IntMap int2TxMapping;
    private Int2IntMap int2GeneMapping;
    private Int2IntMap tx2geneMapping;

    private final Int2IntOpenHashMap geneCounts = new Int2IntOpenHashMap();
    private final Int2IntOpenHashMap txCounts = new Int2IntOpenHashMap();
    private int ambigReads = 0;

    // getters for merging threads
    public Int2IntOpenHashMap getGeneCounts() { return geneCounts; }
    public int getAmbigReads() { return ambigReads; }

    public IndexGraphTraversal(IndexGraph g, Int2IntMap tx2gene) {
        this.g = g;
        this.tx2geneMapping = tx2gene;
    }

    // reusables, thread-safety only if IndexGraphTraversal is worker-specific
    private final IntOpenHashSet candidateTranscripts = new IntOpenHashSet();
    private final IntOpenHashSet rescueBuffer = new IntOpenHashSet();

    public void process(byte[] read) {
        ShortArrayList minims = Minimizers.of(read);
        candidateTranscripts.clear();

        boolean rescueMode = false;
        IndexGraph.Node lastNode = null;

        for (int i = 0; i < minims.size(); i++) {
            short minim = minims.get(i);

            // try rescuing (skip low quality etc)
            if (rescueMode && getQualityFlag(minim) != 3) {
                if (lastNode != null) {
                    rescueBuffer.clear();
                    for (IndexGraph.Node n : lastNode.getAllNextNeighbors()) {
                        rescueBuffer.addAll(n.getTranscripts());
                    }
                    candidateTranscripts.retainAll(rescueBuffer);
                }
                continue;
            }

            IndexGraph.Node currentNode;
            if (lastNode != null) {
                currentNode = lastNode.getNextNeighbor(minim);
                if (currentNode != null) {
                    candidateTranscripts.retainAll(currentNode.getTranscripts());
                }
            } else {
                currentNode = g.getNode(minim);
                if (currentNode != null) {
                    candidateTranscripts.addAll(currentNode.getTranscripts());
                }
            }

            // if no evidence found, reset and enter rescueMode
            if (candidateTranscripts.isEmpty() || currentNode == null) {
                if (!rescueMode) {
                    i = -1; // Restart loop once
                    rescueMode = true;
                    lastNode = null;
                    candidateTranscripts.clear();
                    continue;
                } else {
                    break; // rescue failed, give up!
                }
            }

            lastNode = currentNode;
        }

        int geneId = getGeneIfAllIsoforms(candidateTranscripts);
        if (geneId > -1) {
            geneCounts.addTo(geneId, 1);
        } else {
            ambigReads++;
        }
    }

    private byte getQualityFlag(short minim) {
        return (byte) ((minim >>> 14) & 3);
    }

    private int getGeneIfAllIsoforms(IntOpenHashSet txSet) {
        if (txSet.isEmpty()) return -1;

        int geneId = -1;
        for (int tx : txSet) {
            int currentGene = tx2geneMapping.get(tx);
            if (geneId == -1) {
                geneId = currentGene;
            } else if (geneId != currentGene) {
                return -1;
            }
        }
        return geneId;
    }
}
