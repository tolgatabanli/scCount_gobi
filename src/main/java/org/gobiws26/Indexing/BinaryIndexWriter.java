package org.gobiws26.Indexing;

import org.gobiws26.Config;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.*;

import static org.gobiws26.Main.PROGRAM_IDENTIFIER;

public class BinaryIndexWriter {
    public static final int INDEX_VERSION = 4;


    public static void write(File filePath, Indexer indexer) throws IOException {

        Int2ObjectOpenHashMap<IntArrayList> txIdToMinimizerPath = indexer.getTxIdToMinimizerPath();
        String[] geneIdArray = indexer.getGeneIdArray();
        String[] txIdArray = indexer.getTxIdArray();
        int[] txToGeneArray = indexer.getTxToGeneArray();

        try (DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(filePath)))) {
            // Header
            dos.writeInt(PROGRAM_IDENTIFIER);
            dos.writeInt(INDEX_VERSION); // Version
            dos.writeInt(Config.K); // Kmer length
            dos.writeInt(Config.minimLength); // Minimizer length

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

            // Write minimizer paths: [tx_count + for each tx: minimizer_count + minimizers]
            dos.writeInt(txIdArray.length);
            for (int txId = 0; txId < txIdArray.length; txId++) {
                IntArrayList minimizerPath = txIdToMinimizerPath.get(txId);
                
                if (minimizerPath == null || minimizerPath.isEmpty()) {
                    dos.writeInt(0);
                } else {
                    dos.writeInt(minimizerPath.size());
                    for (int i = 0; i < minimizerPath.size(); i++) {
                        dos.writeInt(minimizerPath.getInt(i));
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
