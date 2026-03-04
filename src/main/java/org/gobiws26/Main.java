package org.gobiws26;


import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;


import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.fastq.FastqReader;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.shorts.ShortArrayList;
import org.gobiws26.Readers.GTFReader;
import org.gobiws26.genomicstruct.Transcript;
import org.gobiws26.utils.Minimizers;


public class Main {
    public static File fastaRef = null;
    public static File fastaRefIdx = null;
    public static File gtfFile = null;
    public static File readTwoFile = null;

    public static void main(String[] args) throws IOException {
        //if (args[0].equals("index")) argParserIndex(args);
        //else if (args[0].equals("count")) argParserCount(args);
        //else {
        // System.err.println("Could not identify command: " + args[0])
        // printHelp()
        //}
        argParser(args);

        // GTF gives a number to transcripts (the order in the below list)
        ArrayList<String> int2Transcript = new ArrayList<>();
        HashMap<String, Transcript> transcripts = (new GTFReader(int2Transcript)).read(gtfFile);

        //Indexer indexer = new Indexer(int2Transcript, transcripts);

        Transcript debugT = transcripts.get("ENSSSCT00000092142");

        try (ReferenceSequenceFile fasta = ReferenceSequenceFileFactory.getReferenceSequenceFile(fastaRef)) {
            byte[] debugSeq = fasta.getSubsequenceAt(debugT.getChr(), debugT.getStart(), debugT.getEnd()).getBases();
        }

        long startTime = System.nanoTime();
        int readCounter = 0;
        int readCounterM = 0;
        IntArrayList minimizerCountsPerRead = new IntArrayList();
        try (FastqReader readTwoReader = new FastqReader(readTwoFile)) {

            // Iterate through the whole read2 file
            while (readTwoReader.hasNext()) {
                FastqRecord theRead = readTwoReader.next();
                ShortArrayList minimSet = Minimizers.of(theRead.getReadBases()); // TODO: qualities
                readCounter++;

                if (readCounter % 1_000_000 == 0) {
                    readCounterM++;
                    System.out.println("Reads processed: " + readCounterM + "M");
                }

                minimizerCountsPerRead.add(minimSet.size());
            }

        } catch (Exception e) {
            throw new RuntimeException("Readers.FastqReader error:\n" + e);
        }
        long endTime = System.nanoTime();
        System.out.println("Time: " + (endTime - startTime) / 1_000_000);

        // for distribution of minimizer number (different qualities are taken as different)
        BufferedWriter bw = new BufferedWriter(new FileWriter("/mnt/cip/home/t/tabanli/Desktop/scCount/prototyping/minimizerCountsPerRead.txt"));
        for (int i : minimizerCountsPerRead) {
            bw.write(String.valueOf(i));
            bw.newLine();
        }
        bw.flush();
        bw.close();

        // TODO create index file:
        //  1. Write transcript id and gene id in the order imposed by int2Transcript// maybe find better way to store this order
        //  2. Map minimizers to transcripts
    }



    private static void argParser(String[] args) {
        if (args.length == 0 || args[0].equals("-h") || args[0].equals("--help")) {
            printHelp();
            System.exit(0);
        }
        for (int i = 0; i < args.length; i++) {
            switch (args[i]) {
                case "-f":
                case "--fasta":
                    if (i + 1 < args.length) {
                        fastaRef = new File(args[++i]);
                    } else {
                        System.err.println("Error: [-f | --fasta] Please specify a reference FASTA!");
                        System.exit(1);
                    }
                    break;

                // FAI file | optional since samtools searches it in the same dir as fa by default
                case "-fidx":
                case "--fastaIdx":
                    if (i + 1 < args.length) {
                        fastaRefIdx = new File(args[++i]);
                    }
                    break;

                case "-g":
                case "-gtf":
                    if (i + 1 < args.length) {
                        gtfFile = new File(args[++i]);
                    } else {
                        System.err.println("Error: [-g | --gtf] Please specify a GTF!");
                        System.exit(1);
                    }
                    break;

                case "-r2":
                    if (i + 1 < args.length) {
                        readTwoFile = new File(args[++i]);
                    } else {
                        System.err.println("Error: [-r2] Please specify a read file in FASTQ format!");
                        System.exit(1);
                    }
                    break;

                case "-k":
                    if (i + 1 < args.length) {
                        Config.K = Integer.parseInt(args[++i]);
                    } else {
                        System.err.println("Error: [-k] Please specify a K!");
                        System.exit(1);
                    }
                    break;

                default:
                    System.err.println("Unknown argument: " + args[i]);
                    System.err.println("Use -h or --help to see the parameter list.");
                    System.exit(1);
            }
        }

        // TODO: check required args
        if (fastaRef == null) {
            System.err.println("Error: Please specify a fasta file!");
        }

    }

    public static void printHelp() {
        System.out.println("Usage: java -jar goenrich.jar [options]");
        System.out.println("Input files:");
        System.out.printf("  %-25s %s%n", "-f <file>", "Path to the reference FASTA file.");
        System.out.printf("  %-25s %s%n", "-fidx <file>", "Path to the reference FASTA index file.");

        System.out.println();

        System.out.println("Output:");
        System.out.println();
    }
}