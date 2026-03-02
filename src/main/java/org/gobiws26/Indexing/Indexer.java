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
 *  1) lines of transcript_id's and their gene_id's as they appeared in the GTF, where line number corresponds to their 'index'
 *  2) map of minimizers to transcript index (the above line index)
 */
public class Indexer {
    private ReferenceSequenceFile refSeqFile;
    private HashMap<String, Transcript> transcripts;
    private ArrayList<String> transcriptIDArray;

    // Map minimizer (short) -> transcripts (int array)
    //private Short2IntArrayMap minimizer2Transcripts; // TODO: use synchronized version or the concurrent wrapper for multi-thread
    private ConcurrentHashMap<Short, IntArrayList> minimizer2Transcripts;


    public Indexer(HashMap<String, Transcript> transcripts, ArrayList<String> transcriptIDArray, ReferenceSequenceFile refSeqFile) {
        this.transcriptIDArray = transcriptIDArray;
        this.transcripts = transcripts;
        this.refSeqFile = refSeqFile;

    }

    public void runIndex() {
        TranscriptomeFetcher tf = new TranscriptomeFetcher(refSeqFile);
        for (Map.Entry<String, Transcript> txEntry : transcripts.entrySet()) {
            Transcript tx = txEntry.getValue();
            String txId = txEntry.getKey();

            byte[] seq = tf.fetchTranscriptSequenceOf(tx);
            ShortSet minimSet = Minimizers.of(seq);
            for (short minimizer : minimSet) {
                IntArrayList txList = minimizer2Transcripts.computeIfAbsent(minimizer, k -> new IntArrayList());
                txList.add(transcriptIDArray.indexOf(txId));
            }
        }
    }
    // TODO: add thread
    // TODO: can give tail length as second arg

    public void writeToFile(BufferedWriter bw) {

    }
}
