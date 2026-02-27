package org.gobiws26;


import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import org.gobiws26.Readers.FastqReader;
import org.gobiws26.Readers.GTFReader;
import org.gobiws26.genomicstruct.Transcript;
import org.gobiws26.utils.KmerIteratorByte;
import org.gobiws26.utils.KmerIteratorHash;
import org.gobiws26.utils.KmerIteratorLong;
import org.gobiws26.utils.TranscriptomeFetcher;

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
        HashMap<String, Transcript> transcripts = (new GTFReader()).read(gtfFile);

        Transcript debugT = transcripts.get("ENSSSCT00000092142");

        try (ReferenceSequenceFile fasta = ReferenceSequenceFileFactory.getReferenceSequenceFile(fastaRef)) {
            byte[] debugSeq = fasta.getSubsequenceAt(debugT.getChr(), debugT.getStart(), debugT.getEnd()).getBases();


        }

        try (FastqReader readTwoReader = new FastqReader(readTwoFile)) {
            //System.out.println(readTwoReader.getNextRead());

        } catch (Exception e) {
            throw new RuntimeException("Readers.FastqReader error:\n" + e);
        }

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