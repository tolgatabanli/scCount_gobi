package org.gobiws26.Querying;

import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import org.gobiws26.Indexing.IndexData;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class CountMatrixWriter {
    private BufferedWriter bw;
    private IndexData idxData;

    public CountMatrixWriter(File file, IndexData idxData) throws IOException {
        bw = new BufferedWriter(new FileWriter(file));
        this.idxData = idxData;
    }

    public void write(Int2IntOpenHashMap geneInt2Counts, Int2IntOpenHashMap txInt2Counts) throws IOException {
        Int2ObjectOpenHashMap<String> transcriptInterpreter = idxData.int2TxString;
        Int2ObjectOpenHashMap<String> geneInterpreter = idxData.int2GeneString;
        Int2IntOpenHashMap geneOfTranscript = idxData.txInt2GeneInt;

        // Group transcripts by gene
        Int2ObjectOpenHashMap<Int2IntOpenHashMap> geneToTxCounts = new Int2ObjectOpenHashMap<>();
        for (Int2IntMap.Entry entry : txInt2Counts.int2IntEntrySet()) {
            int txInt = entry.getIntKey();
            int count = entry.getIntValue();
            int geneInt = geneOfTranscript.get(txInt);
            geneToTxCounts.computeIfAbsent(geneInt, k -> new Int2IntOpenHashMap()).put(txInt, count);
        }

        // Iterate through genes and write gene-level counts followed by their transcripts
        for (Int2IntMap.Entry entry : geneInt2Counts.int2IntEntrySet()) {
            int geneInt = entry.getIntKey();
            int count = entry.getIntValue();
            String geneId = geneInterpreter.get(geneInt);

            if (geneId != null) {
                bw.write(geneId + "\t" + count);
                bw.newLine();

                // Write transcripts for this gene
                Int2IntOpenHashMap txCounts = geneToTxCounts.get(geneInt);
                if (txCounts != null) {
                    for (Int2IntMap.Entry txEntry : txCounts.int2IntEntrySet()) {
                        int txInt = txEntry.getIntKey();
                        int txCount = txEntry.getIntValue();
                        String txId = transcriptInterpreter.get(txInt);
                        if (txId != null) {
                            bw.write(txId + "\t" + txCount);
                            bw.newLine();
                        }
                    }
                }
            }
        }
        bw.flush();
    }
}
