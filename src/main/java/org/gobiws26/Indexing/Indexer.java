package org.gobiws26.Indexing;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntArrays;
import it.unimi.dsi.fastutil.shorts.*;
import org.gobiws26.genomicstruct.Transcript;
import org.gobiws26.utils.Minimizers;
import org.gobiws26.utils.TranscriptomeFetcher;

import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

/**
 * Takes an output file and transcriptome to write an index file that maps the following information:
 *  1) transcript_id's, their gene_id's and a unique ID (int)
 *  2) map of minimizers to transcript (the above unique ID)
 */
public class Indexer {
    // inputs
    private ReferenceSequenceFile refSeqFile;
    private HashMap<String, Transcript> transcripts;

    // outputs
    private HashMap<String, Integer> transcriptToIndex = new HashMap<>();

    // Map minimizer (short) -> transcripts (int array)
    // TODO: use synchronized version or the concurrent wrapper for multi-thread
    private ConcurrentHashMap<Short, IntArrayList> minimizer2Transcripts;


    public Indexer(HashMap<String, Transcript> transcripts, ArrayList<String> transcriptIDArray, ReferenceSequenceFile refSeqFile) {
        this.transcriptIDArray = transcriptIDArray;
        this.transcripts = transcripts;
        this.refSeqFile = refSeqFile;

    }

    // TODO: REIHENFOLGE IST AUCH WICHTIG VON DEN MINIMIZERN!!
    public void runIndex() {
        minimizer2Transcripts = new ConcurrentHashMap<>();

        int txCounter = 0;
        TranscriptomeFetcher tf = new TranscriptomeFetcher(refSeqFile);
        for (Map.Entry<String, Transcript> txEntry : transcripts.entrySet()) {
            Transcript tx = txEntry.getValue();
            String txId = txEntry.getKey();
            transcriptToIndex.put(txId, txCounter);

            byte[] seq = tf.fetchTranscriptSequenceOf(tx);
            ShortArrayList minimSet = Minimizers.of(seq);
            for (short minimizer : minimSet) {
                IntArrayList txList = minimizer2Transcripts.computeIfAbsent(minimizer, k -> new IntArrayList());
                txList.add(txCounter);
            }

            txCounter++;
        }
    }
    // TODO: add thread
    // TODO: can give tail length as second arg

    public void writeToFile(BufferedWriter bw) {
        if (transcriptToIndex == null) throw new IllegalStateException("Indexer has not been run!");

    }
}
