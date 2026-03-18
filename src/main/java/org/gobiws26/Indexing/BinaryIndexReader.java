package org.gobiws26.Indexing;

import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import org.gobiws26.Config;
import org.gobiws26.Querying.IndexGraph;

import java.io.*;
import java.nio.charset.StandardCharsets;

import static org.gobiws26.Indexing.BinaryIndexWriter.INDEX_VERSION;
import static org.gobiws26.Main.PROGRAM_IDENTIFIER;

public class BinaryIndexReader {

    public static IndexData read(File indexFile) throws IOException {
        try (DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(indexFile)))) {

            // Verify Header
            int magic = dis.readInt();
            if (magic != PROGRAM_IDENTIFIER) throw new IOException("Invalid file format. Not written by the same program.");
            int version = dis.readInt();
            if (version != INDEX_VERSION) throw new IOException("Invalid file format. Not written by the same version: " + version + " vs " + INDEX_VERSION);
            Config.K = dis.readInt();
            Config.minimLength = dis.readInt();

            // 1) read genes
            int numGenes = dis.readInt();
            Int2ObjectOpenHashMap<String> int2GeneString = new Int2ObjectOpenHashMap<>(numGenes);
            for (int i = 0; i < numGenes; i++) {
                int2GeneString.put(i, readString(dis));
            }

            // 2) Read Tx & tx -> gene map
            int numTx = dis.readInt();
            Int2ObjectOpenHashMap<String> int2TxString = new Int2ObjectOpenHashMap<>(numTx);
            Int2IntOpenHashMap txInt2GeneInt = new Int2IntOpenHashMap(numTx);

            for (int i = 0; i < numTx; i++) {
                int2TxString.put(i, readString(dis));
                txInt2GeneInt.put(i, dis.readInt());
            }

            // 3) Read minimizer paths and build triplet graph
            int numPaths = dis.readInt();
            if (numPaths != numTx) {
                throw new IOException("Mismatch between transcript count and path count.");
            }

            IndexGraph graph = new IndexGraph();
            for (int txId = 0; txId < numPaths; txId++) {
                int pathLength = dis.readInt();
                if (pathLength > 0) {
                    IntArrayList minimizers = new IntArrayList(pathLength);
                    for (int i = 0; i < pathLength; i++) {
                        minimizers.add(dis.readInt());
                    }
                    // Build triplet index from minimizers
                    graph.addTxPath(txId, minimizers);
                }
            }

            return new IndexData(graph, int2TxString, int2GeneString, txInt2GeneInt);
        }
    }

    // reads [field_length + field_content]
    // returns field_content
    private static String readString(DataInputStream dis) throws IOException {
        int length = dis.readShort();
        byte[] bytes = new byte[length];
        dis.readFully(bytes);
        return new String(bytes, StandardCharsets.UTF_8);
    }
}
