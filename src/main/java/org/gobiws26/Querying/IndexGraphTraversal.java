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

    public IndexGraphTraversal(IndexGraph g, Int2IntMap tx2gene) {
        this.g = g;
        this.tx2geneMapping = tx2gene;
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

    public void process(FastqRecord read) {
        ShortArrayList minims = Minimizers.of(read.getReadBases(), read.getBaseQualities());
        if (minims.isEmpty()) return;

        // 1) Heuristic/ideal case: intersect transcripts of all minimizers -> get transcripts -> see if single transcript of gene
        BitSet candidates = null;
        int exploredCount = 0;

        for (int i = 0; i < minims.size(); i++) {
            IndexGraph.Node node = g.getNode(minims.getShort(i));
            if (node == null) { candidates = new BitSet(); break; }

            BitSet nodeTxs = node.getTranscriptsBitSet();
            if (candidates == null) {
                candidates = (BitSet) nodeTxs.clone();
            } else {
                candidates.and(nodeTxs); // retainAll / intersection
            }
            exploredCount = i + 1;

            if (candidates.isEmpty()) break;

            // Early stop: unique gene found
            if (getGeneIfAllIsoforms(candidates) != -1) break;
        }

        if (!candidates.isEmpty()) {
            // Verify order-preserving path for explored minimizers on each candidate
            BitSet validated = verifyOrderedPath(minims, exploredCount, candidates);

            if (!validated.isEmpty()) {
                int geneId = getGeneIfAllIsoforms(validated);
                if (validated.size() == 1) {
                    // Unambiguous transcript
                    int txId = validated.nextSetBit(0);
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
        }

        // 2) quasi-alignment (optimized): reuse candidate transcripts if available,
        //    traverse the graph once and score only transcripts that remain candidates.
        BitSet validated = null;
        // If we produced a `candidates` BitSet earlier and it had members, verifyOrderedPath
        // may already have produced a validated set in the heuristic block. Re-run verify
        // only if candidates exists and is non-empty.
        if (candidates != null && !candidates.isEmpty()) {
            validated = verifyOrderedPath(minims, exploredCount, candidates);
            if (validated != null && !validated.isEmpty()) {
                int geneId = getGeneIfAllIsoforms(validated);
                if (validated.size() == 1) {
                    int txId = validated.nextSetBit(0);
                    txCounts.merge(txId, 1, Integer::sum);
                    geneCounts.merge(tx2geneMapping.get(txId), 1, Integer::sum);
                    numReadsFoundHeuristically++;
                    return;
                } else if (geneId != -1) {
                    geneCounts.merge(geneId, 1, Integer::sum);
                    numReadsFoundHeuristically++;
                    return;
                }
            }
        }

        int txCount = tx2geneMapping.size();
        int[] scores   = new int[txCount];
        int[] expected = new int[txCount];
        Arrays.fill(expected, -1);

        int firstScore = 0, secondScore = 0;
        int halfMinims = (minims.size() - 1) / 2;

        // Use a BitSet to restrict scoring to candidates when available. If `validated` is null
        // or empty, we score all transcripts.
        BitSet txsToScore = (validated != null && !validated.isEmpty()) ? (BitSet) validated.clone() : null;

        for (int i = 0; i + 1 < minims.size(); i++) {
            short curMinim  = minims.getShort(i);
            short nextMinim = minims.getShort(i + 1);

            IndexGraph.Node curNode = g.getNode(curMinim);
            if (curNode == null) continue;

            IndexGraph.Edge edge = curNode.getEdgeTo(nextMinim);
            if (edge == null) continue;

            // iterate over transcripts that have this edge; skip those not in txsToScore when present
            for (int txId = edge.getTxIdIndex().nextSetBit(0); txId >= 0; txId = edge.getTxIdIndex().nextSetBit(txId + 1)) {
                if (txsToScore != null && !txsToScore.get(txId)) continue;
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

        // Find max raw and quit fast if nothing matched
        int maxRaw = 0;
        for (int s : scores) if (s > maxRaw) maxRaw = s;
        if (maxRaw == 0) { ambigReads++; return; }

        // Normalize scores by transcript minimizer count, keep as float
        float[] normScores = new float[txCount];
        float maxNorm = 0f;
        for (int txId = 0; txId < txCount; txId++) {
            if (scores[txId] == 0) continue;
            if (txsToScore != null && !txsToScore.get(txId)) continue;
            int txLen = g.getMinimCountOfTranscript(txId);
            float norm = txLen > 0 ? (float) scores[txId] / txLen : 0f;
            normScores[txId] = norm;
            if (norm > maxNorm) maxNorm = norm;
        }

        if (maxNorm == 0f) { ambigReads++; return; }

        // Pick best-scoring transcript(s) within a small tolerance for float comparison
        BitSet bestTxs = new BitSet();
        float threshold = maxNorm * 0.999f;
        for (int txId = 0; txId < txCount; txId++) {
            if (normScores[txId] >= threshold) bestTxs.set(txId);
        }

        int geneId = getGeneIfAllIsoforms(bestTxs);
        if (bestTxs.cardinality() == 1) {
            int txId = bestTxs.nextSetBit(0);
            txCounts.merge(txId, 1, Integer::sum);
            geneCounts.merge(tx2geneMapping.get(txId), 1, Integer::sum);
        } else if (geneId != -1) {
            geneCounts.merge(geneId, 1, Integer::sum);
        } else {
            ambigReads++; return;
        }
        numReadsFoundWithAlignment++;
    }


    private BitSet verifyOrderedPath(ShortArrayList minims, int exploredCount, BitSet candidates) {
        BitSet valid = (BitSet) candidates.clone();

        for (int i = 0; i + 1 < exploredCount && !valid.isEmpty(); i++) {
            IndexGraph.Node curNode = g.getNode(minims.getShort(i));
            if (curNode == null) return new BitSet();

            IndexGraph.Edge edge = curNode.getEdgeTo(minims.getShort(i + 1));
            if (edge == null) return new BitSet();

            for (int txId = valid.nextSetBit(0); txId >= 0; txId = valid.nextSetBit(txId + 1)) {
                IntArrayList positions = edge.getPositionsForTx(txId);
                // Position at step i must be exactly i (0-indexed consecutive path)
                if (!positions.contains(i)) valid.clear(txId);
            }
        }
        return valid;
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

    private int getGeneIfAllIsoforms(IntOpenHashSet txSet) {
        if (txSet.isEmpty()) return -1;
        int geneId = -1;
        for (int tx : txSet) {
            int g = tx2geneMapping.get(tx);
            if (geneId == -1) geneId = g;
            else if (geneId != g) return -1;
        }
        return geneId;
    }

    private short[] maskQualities(ShortArrayList minims) {
        short[] res = new short[minims.size()];
        for (int i = 0; i < res.length; i++)
            res[i] = (short) (minims.getShort(i) | 0xC000);
        return res;
    }

    private byte getQualityFlag(short minim) {
        return (byte) ((minim >>> 14) & 3);
    }
}