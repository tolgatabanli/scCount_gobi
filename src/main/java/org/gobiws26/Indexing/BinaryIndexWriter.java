package org.gobiws26.Indexing;

import org.gobiws26.Config;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.shorts.ShortArrayList;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.*;

import static org.gobiws26.Main.PROGRAM_IDENTIFIER;

public class BinaryIndexWriter {
    public static final int INDEX_VERSION = 1;

    public static void write(File filePath, Indexer indexer) throws IOException {

        Int2ObjectOpenHashMap<ShortArrayList> paths = indexer.getTranscriptToMinimizerPath();
        String[] geneIdArray = indexer.getGeneIdArray();
        String[] txIdArray = indexer.getTxIdArray();
        int[] txToGeneArray = indexer.getTxToGeneArray();

        try (DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(filePath)))) {
            // Header
            dos.writeInt(PROGRAM_IDENTIFIER);
            dos.writeInt(INDEX_VERSION); // Version
            dos.writeInt(Config.K); // Kmer length

            // Write genes: [gene_count + geneId]
            dos.writeInt(geneIdArray.length);
            for (String gene : geneIdArray) {
                writeString(dos, gene);
            }

            // Write transcripts: [tx_count + txId + geneId]
            dos.writeInt(txIdArray.length);
            for (int i = 0; i < txIdArray.length; i++) {
                writeString(dos, txIdArray[i]);
                dos.writeInt(txToGeneArray[i]);
            }

            // Write ordered minimizers pro Tx: [tx_count + [ path_size + minimizerShort ] for each tx ]
            dos.writeInt(txIdArray.length); // Should match tx count
            for (int txInternalId = 0; txInternalId < txIdArray.length; txInternalId++) {
                ShortArrayList path = paths.get(txInternalId);

                if (path == null || path.isEmpty()) { // should hopefully not occur
                    dos.writeInt(0);
                } else {
                    dos.writeInt(path.size());
                    for (int p = 0; p < path.size(); p++) {
                        dos.writeShort(path.getShort(p));
                    }
                }
            }
        }
    }

    // writes str_length + str
    private static void writeString(DataOutputStream dos, String str) throws IOException {
        byte[] bytes = str.getBytes(StandardCharsets.UTF_8);
        dos.writeShort(bytes.length);
        dos.write(bytes);
    }
}
