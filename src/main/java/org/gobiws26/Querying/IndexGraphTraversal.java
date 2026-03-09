package org.gobiws26.Querying;

import htsjdk.samtools.fastq.FastqRecord;
import it.unimi.dsi.fastutil.ints.*;
import it.unimi.dsi.fastutil.shorts.ShortArrayList;
import org.gobiws26.utils.Minimizers;

import java.util.Arrays;


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
    private int readsFoundHeuristically = 0;
    private int readsFoundWithAlignment = 0;
    public void process(FastqRecord read) {
        ShortArrayList minims = Minimizers.of(read.getReadBases(), read.getBaseQualities());
        if (minims.isEmpty()) return;

        // 1) Heuristic: intersect transcripts of all minimizers -> get transcripts -> see if single transcript of gene
        IntOpenHashSet candidates = null;
        int exploredCount = 0;

        for (int i = 0; i < minims.size(); i++) {
            IndexGraph.Node node = g.getNode(minims.getShort(i));
            if (node == null) { candidates = new IntOpenHashSet(); break; }

            IntSet nodeTxs = node.getTranscripts();
            if (candidates == null) {
                candidates = new IntOpenHashSet(nodeTxs);
            } else {
                candidates.retainAll(nodeTxs);
            }
            exploredCount = i + 1;

            if (candidates.isEmpty()) break;

            // Early stop: unique gene found
            if (getGeneIfAllIsoforms(candidates) != -1) break;
        }

        if (!candidates.isEmpty()) {
            // Verify order-preserving path for explored minimizers on each candidate
            IntOpenHashSet validated = verifyOrderedPath(minims, exploredCount, candidates);

            if (!validated.isEmpty()) {
                int geneId = getGeneIfAllIsoforms(validated);
                if (validated.size() == 1) {
                    // Unambiguous transcript
                    int txId = validated.iterator().nextInt();
                    txCounts.merge(txId, 1, Integer::sum);
                    geneCounts.merge(tx2geneMapping.get(txId), 1, Integer::sum);
                    readsFoundHeuristically++;
                    return;
                } else if (geneId != -1) {
                    // Unambiguous gene
                    geneCounts.merge(geneId, 1, Integer::sum);
                    readsFoundHeuristically++;
                    return;
                }
            }
        }
        readsFoundWithAlignment++;

        // 2) quasi-alignment
        int txCount = tx2geneMapping.size();
        int[] scores   = new int[txCount];
        int[] expected = new int[txCount];
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

        float[] normScores = new float[txCount];
        float maxNorm = 0f;
        for (int txId = 0; txId < txCount; txId++) {
            if (scores[txId] == 0) continue;
            int txLen = g.getMinimCountOfTranscript(txId);
            float norm = txLen > 0 ? (float) scores[txId] / txLen : 0f;
            normScores[txId] = norm;
            if (norm > maxNorm) maxNorm = norm;
        }

        // Pick best-scoring transcript(s) within a small tolerance for float comparison
        IntOpenHashSet bestTxs = new IntOpenHashSet();
        float threshold = maxNorm * 0.999f;
        for (int txId = 0; txId < txCount; txId++) {
            if (normScores[txId] >= threshold) bestTxs.add(txId);
        }

        int geneId = getGeneIfAllIsoforms(bestTxs);
        if (bestTxs.size() == 1) {
            int txId = bestTxs.iterator().nextInt();
            txCounts.merge(txId, 1, Integer::sum);
            geneCounts.merge(tx2geneMapping.get(txId), 1, Integer::sum);
        } else if (geneId != -1) {
            geneCounts.merge(geneId, 1, Integer::sum);
        } else {
            ambigReads++;
        }
    }


    private IntOpenHashSet verifyOrderedPath(ShortArrayList minims, int exploredCount, IntOpenHashSet candidates) {
        IntOpenHashSet valid = new IntOpenHashSet(candidates);

        for (int i = 0; i + 1 < exploredCount && !valid.isEmpty(); i++) {
            IndexGraph.Node curNode = g.getNode(minims.getShort(i));
            if (curNode == null) return new IntOpenHashSet();

            IndexGraph.Edge edge = curNode.getEdgeTo(minims.getShort(i + 1));
            if (edge == null) return new IntOpenHashSet();

            IntIterator it = valid.iterator();
            while (it.hasNext()) {
                int txId = it.nextInt();
                IntArrayList positions = edge.getPositionsForTx(txId);
                // Position at step i must be exactly i (0-indexed consecutive path)
                if (!positions.contains(i)) it.remove();
            }
        }
        return valid;
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