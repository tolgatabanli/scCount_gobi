package org.gobiws26.Indexing;

import org.junit.jupiter.api.Test;
import static org.gobiws26.Test.*;

import org.gobiws26.Config;
import org.gobiws26.Main;
import org.gobiws26.genomicstruct.Transcript;

import static org.junit.jupiter.api.Assertions.*;

import java.io.*;
import java.nio.charset.StandardCharsets;

public class IndexFileTest {
    // public static HashMap<String, Transcript> transcripts; // imported from org.gobiws26.Test

    @Test
    public void fileStructureIsValid() throws IOException {
        try (DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(indexFile)))) {
            // 1. Check Header
            assertEquals(Main.PROGRAM_IDENTIFIER, dis.readInt(), "Program identifier mismatch");
            assertEquals(BinaryIndexWriter.INDEX_VERSION, dis.readInt(), "Index version mismatch");
            int k = dis.readInt();
            assertEquals(Config.K, k, "Kmer length mismatch");

            // 2. Check Genes
            int numGenes = dis.readInt();
            assertTrue(numGenes >= 0, "Gene count cannot be negative");
            for (int i = 0; i < numGenes; i++) {
                assertNotNull(readString(dis), "Gene ID should not be null");
            }

            // 3. Check Transcripts & tx2gene mapping
            int numTx = dis.readInt();
            assertTrue(numTx >= 0, "Transcript count cannot be negative");
            for (int i = 0; i < numTx; i++) {
                assertNotNull(readString(dis), "Transcript ID should not be null");
                int geneIdx = dis.readInt();
                assertTrue(geneIdx >= 0 && geneIdx < numGenes, "Gene index out of bounds");
            }

            // 4. Check Ordered Minimizers
            int pathTxCount = dis.readInt();
            assertEquals(numTx, pathTxCount, "Transcript count for paths does not match main transcript count");

            for (int i = 0; i < pathTxCount; i++) {
                int pathSize = dis.readInt();
                assertTrue(pathSize >= 0, "Path size cannot be negative");
                for (int p = 0; p < pathSize; p++) {
                    dis.readShort(); // Verify we can read the correct number of shorts without EOFException
                }
            }

            // Ensure we reached the exact end of the file
            assertEquals(0, dis.available(), "File has trailing unexpected bytes");
        }
    }

    @Test
    public void transcriptToGeneMappingIsCorrect() throws IOException {
        String targetTx = "ENSSSCT00000081764";
        String foundGene = null;

        try (DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(indexFile)))) {
            // Skip Header
            dis.readInt(); dis.readInt(); dis.readInt();

            // Read Genes into an array for index lookup
            int numGenes = dis.readInt();
            String[] genes = new String[numGenes];
            for (int i = 0; i < numGenes; i++) {
                genes[i] = readString(dis);
            }

            // Scan Transcripts
            int numTx = dis.readInt();
            for (int i = 0; i < numTx; i++) {
                String txId = readString(dis);
                int geneIdx = dis.readInt();
                if (txId.equals(targetTx)) {
                    foundGene = genes[geneIdx];
                    break;
                }
            }
        }

        assertNotNull(foundGene, "Transcript " + targetTx + " was not found in the index file");

        // Check against the GTF ground truth imported from your test setup
        // Note: adjust 'getGeneId()' to whatever method your Transcript class uses
        Transcript expectedTranscript = transcripts.get(targetTx);
        assertNotNull(expectedTranscript, "Target transcript missing from GTF transcripts map");
        assertEquals(expectedTranscript.getGeneId(), foundGene, "Transcript to Gene mapping does not match GTF");
    }

    @Test
    public void minimizerPathIsCorrect() throws IOException {
        String targetTx = "ENSSSCT00000081764";
        String[] expectedMinimizersString = new String[]{"ACGCATG", "AAGGGAG", "AAAAGAA", "AAAGAAT", "AATATTT", "ATATTTC", "ATTTCGC", "CCTGTTT", "AAATCTG", "AATCTGA", "AATTGGT", "AATATTA", "ATATTAT", "ATTATGT", "ATGCACA", "ACAGGTC", "ACTTATG", "AATGTTG"};

        String[] foundPathStrings = null;

        try (DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(indexFile)))) {
            // Skip Header
            dis.readInt(); dis.readInt(); dis.readInt();

            // Skip Genes
            int numGenes = dis.readInt();
            for (int i = 0; i < numGenes; i++) {
                readString(dis);
            }

            // Scan Transcripts to find the internal ID (index) for targetTx
            int numTx = dis.readInt();
            int targetInternalId = -1;
            for (int i = 0; i < numTx; i++) {
                String txId = readString(dis);
                dis.readInt(); // Skip gene mapping
                if (txId.equals(targetTx)) {
                    targetInternalId = i;
                }
            }

            assertTrue(targetInternalId >= 0, "Target transcript not found in index");

            // Scan Paths
            int pathCount = dis.readInt();
            for (int i = 0; i < pathCount; i++) {
                int pathSize = dis.readInt();
                if (i == targetInternalId) {
                    foundPathStrings = new String[pathSize];
                    for (int p = 0; p < pathSize; p++) {
                        foundPathStrings[p] = shortToNucleotideString(dis.readShort());
                    }
                    break; // Found our transcript's path, exit early
                } else {
                    // Skip over the paths we don't care about
                    for (int p = 0; p < pathSize; p++) {
                        dis.readShort();
                    }
                }
            }
        }

        assertNotNull(foundPathStrings, "Path for target transcript was not found");
        assertArrayEquals(expectedMinimizersString, foundPathStrings, "Minimizer path does not match expected sequences");
    }


    // Helpers

    /**
     * Reads a string encoded identically to BinaryIndexWriter.writeString()
     */
    private String readString(DataInputStream dis) throws IOException {
        short length = dis.readShort();
        byte[] bytes = new byte[length];
        dis.readFully(bytes);
        return new String(bytes, StandardCharsets.UTF_8);
    }
}
