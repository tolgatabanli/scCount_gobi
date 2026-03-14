package org.gobiws26.Querying;

import htsjdk.samtools.fastq.FastqRecord;
import it.unimi.dsi.fastutil.ints.*;
import it.unimi.dsi.fastutil.shorts.ShortArrayList;
import org.gobiws26.utils.Minimizers;

import java.util.Arrays;
import java.util.BitSet;


public class IndexGraphTraversal {

    // Per-read outputs
    private final IndexGraph g;
    private final Int2IntMap tx2geneMapping;

    private final Int2IntOpenHashMap geneCounts = new Int2IntOpenHashMap();
    private final Int2IntOpenHashMap txCounts   = new Int2IntOpenHashMap();
    private int ambigReads = 0;

    // Getters for merging threads (unchanged public API)
    public Int2IntOpenHashMap getGeneCounts() { return geneCounts; }
    public Int2IntOpenHashMap getTxCounts()   { return txCounts;   }
    public int getAmbigReads()                { return ambigReads; }

    int txCount;
    int[] scores; // should be reset
    int[] expected; // should also be reset
    private final BitSet candidatesBitSet = new BitSet();
    private final BitSet validatedBitSet = new BitSet();
    private final BitSet bestTxsBitSet = new BitSet();
    private final float[] normScores;
    public IndexGraphTraversal(IndexGraph g, Int2IntMap tx2gene) {
        this.g = g;
        this.tx2geneMapping = tx2gene;

        txCount = tx2geneMapping.size();
        scores   = new int[txCount];
        expected = new int[txCount];
        normScores = new float[txCount];
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
        if (minims.isEmpty()) return;

        // 1) Heuristic/ideal case: intersect transcripts of all minimizers -> get transcripts -> see if single transcript of gene
        candidatesBitSet.clear();
        int exploredCount = 0;
        boolean candidatesEmpty = false;

        for (int i = 0; i < minims.size(); i++) {
            IndexGraph.Node node = g.getNode(minims.getShort(i));
            if (node == null) { candidatesBitSet.clear(); candidatesEmpty = true; break; }

            BitSet nodeTxs = node.getTranscriptsBitSet();
            if (exploredCount == 0) {
                candidatesBitSet.or(nodeTxs);
            } else {
                candidatesBitSet.and(nodeTxs); // retainAll / intersection
            }
            exploredCount = i + 1;

            if (candidatesBitSet.isEmpty()) { candidatesEmpty = true; break; }

            // Early stop: unique gene found
            if (getGeneIfAllIsoforms(candidatesBitSet) != -1) break;
        }

        if (!candidatesEmpty && !candidatesBitSet.isEmpty()) {
            // Verify order-preserving path for explored minimizers on each candidate
            verifyOrderedPath(minims, exploredCount, candidatesBitSet, validatedBitSet);

            if (!validatedBitSet.isEmpty()) {
                int geneId = getGeneIfAllIsoforms(validatedBitSet);
                if (validatedBitSet.cardinality() == 1) {
                    // Unambiguous transcript
                    int txId = validatedBitSet.nextSetBit(0);
                    txCounts.merge(txId, 1, Integer::sum);
                    geneCounts.merge(tx2geneMapping.get(txId), 1, Integer::sum);
                    numReadsFoundHeuristically++;
                    return;
                } else if (geneId != -1) {
                    // Unambiguous gene
                    geneCounts.merge(geneId, 1, Integer::sum);
                    numReadsFoundHeuristically++;
                    return;
                }
            }
        } else numReadsEmptyCandidates++;

        // 2) quasi-alignment
        Arrays.fill(scores, 0);
        Arrays.fill(expected, -1);

        int firstScore = 0, secondScore = 0;
        int halfMinims = (minims.size() - 1) / 2;

        for (int i = 0; i + 1 < minims.size(); i++) {
            short curMinim  = minims.getShort(i);
            short nextMinim = minims.getShort(i + 1);

            IndexGraph.Node curNode = g.getNode(curMinim);
            if (curNode == null) continue;

            IndexGraph.Edge edge = curNode.getEdgeTo(nextMinim);
            if (edge == null) continue;

            for (int txId = edge.getTxIdIndex().nextSetBit(0); txId >= 0; txId = edge.getTxIdIndex().nextSetBit(txId + 1)) {
                IntArrayList positions = edge.getPositionsForTx(txId);
                for (int p = 0; p < positions.size(); p++) {
                    int pos = positions.getInt(p);
                    if (expected[txId] == -1 || expected[txId] == pos) {
                        int s = ++scores[txId];
                        if (s > firstScore)       { secondScore = firstScore; firstScore = s; }
                        else if (s > secondScore) { secondScore = s; }
                        expected[txId] = pos + 1;
                        break;
                    }
                }
            }

            // Early termination: dominant winner after at least half the minimizers
            if (i >= halfMinims && firstScore >= secondScore * 2) break;
        }

        // Normalize scores by transcript minimizer count, keep as float
        int maxRaw = 0;
        for (int s : scores) if (s > maxRaw) maxRaw = s;
        if (maxRaw == 0) { ambigReads++; return; }

        Arrays.fill(normScores, 0f);
        float maxNorm = 0f;
        for (int txId = 0; txId < txCount; txId++) {
            if (scores[txId] == 0) continue;
            int txLen = g.getMinimCountOfTranscript(txId);
            float norm = txLen > 0 ? (float) scores[txId] / txLen : 0f;
            normScores[txId] = norm;
            if (norm > maxNorm) maxNorm = norm;
        }

        // Pick best-scoring transcript(s) within a small tolerance for float comparison
        bestTxsBitSet.clear();
        float threshold = maxNorm * 0.999f;
        for (int txId = 0; txId < txCount; txId++) {
            if (normScores[txId] >= threshold) bestTxsBitSet.set(txId);
        }

        int geneId = getGeneIfAllIsoforms(bestTxsBitSet);
        int bestCount = bestTxsBitSet.cardinality();
        if (bestCount == 1) {
            int txId = bestTxsBitSet.nextSetBit(0);
            txCounts.addTo(txId, 1);
            geneCounts.addTo(geneId, 1);
        } else if (geneId != -1) {
            geneCounts.addTo(geneId, 1);
        } else {
            ambigReads++; return;
        }
        numReadsFoundWithAlignment++;
    }


    private void verifyOrderedPath(ShortArrayList minims, int exploredCount, BitSet candidates, BitSet outValid) {
        outValid.clear();
        outValid.or(candidates);

        for (int i = 0; i + 1 < exploredCount && !outValid.isEmpty(); i++) {
            IndexGraph.Node curNode = g.getNode(minims.getShort(i));
            if (curNode == null) { outValid.clear(); return; }

            IndexGraph.Edge edge = curNode.getEdgeTo(minims.getShort(i + 1));
            if (edge == null) { outValid.clear(); return; }

            for (int txId = outValid.nextSetBit(0); txId >= 0; txId = outValid.nextSetBit(txId + 1)) {
                IntArrayList positions = edge.getPositionsForTx(txId);
                // Position at step i must be exactly i (0-indexed consecutive path)
                if (!positions.contains(i)) outValid.clear(txId);
            }
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