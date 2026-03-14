package org.gobiws26;


import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;


import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.fastq.FastqReader;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import org.gobiws26.Indexing.*;
import org.gobiws26.Querying.CountMatrixWriter;
import org.gobiws26.Querying.ParallelChainQuery;
import org.gobiws26.Querying.ParallelGraphQuery;
import org.gobiws26.Readers.GTFReader;
import org.gobiws26.genomicstruct.Transcript;
import org.gobiws26.utils.ReadOneParser;


public class Main {
    public static File fastaRef = null;
    public static File gtfFile = null;

    public static File readOneFile = null;
    public static File readTwoFile = null;

    private static File indexFile = null;
    private static File outputDir = null;

    private static int threads = 1;
    private static int batch_size = 1_000_000;

    public static final int PROGRAM_IDENTIFIER = 0x474F4249; // "GOBI"

    public static void main(String[] args) throws IOException, InterruptedException {
        long startTime = System.nanoTime();
        if (args[0].equals("index")) {
            argParserIndex(args);

            HashMap<String, Transcript> transcripts = (new GTFReader()).read(gtfFile);

            Indexer indexer = new Indexer(transcripts, ReferenceSequenceFileFactory.getReferenceSequenceFile(fastaRef));
            indexer.runIndex();
            BinaryIndexWriter.write(indexFile, indexer);
        }
        else if (args[0].equals("count")) {
            argParserCount(args);

            if (!outputDir.exists() && !outputDir.mkdirs()) {
                throw new IOException("Unable to create the output directory: " + outputDir);
            }
            File countMatrixFile = new File(outputDir, "counts.tsv");

            IndexData idxData = BinaryIndexReader.read(indexFile);
            //IndexChainData idxData = BinaryIndexChainReader.read(indexFile);

            ParallelGraphQuery pgq = new ParallelGraphQuery(idxData.graph, idxData.txInt2GeneInt, threads, 
                                                             readOneFile != null); // Enable barcode mode if r1 provided
            //ParallelChainQuery pgq = new ParallelChainQuery(idxData, threads);

            // Process reads in batches to maintain constant memory usage
            final int BATCH_SIZE = batch_size;
            List<FastqRecord> batch = new ArrayList<>(BATCH_SIZE);
            List<String> batchBarcodes = readOneFile != null ? new ArrayList<>(BATCH_SIZE) : null;
            List<String> batchUMIs = readOneFile != null ? new ArrayList<>(BATCH_SIZE) : null;

            int reportingStep = BATCH_SIZE;
            System.out.println("Starting with threads: " + threads);
            System.out.println("Batch size and reporting step is every " + reportingStep + " reads.");
            System.out.println("Barcode-aware mode: " + (readOneFile != null ? "enabled" : "disabled"));

            try (FastqReader readTwoReader = new FastqReader(readTwoFile);
                 FastqReader readOneReader = readOneFile != null ? new FastqReader(readOneFile) : null) {
                
                long totalProcessed = 0;  // Tracks ACTUALLY PROCESSED and COUNTED reads
                long lastReportedMilestone = 0;

                while (readTwoReader.hasNext()) {
                    //if (totalProcessed >= 100) break;
                    FastqRecord theRead = readTwoReader.next();
                    batch.add(theRead);

                    // Extract barcode and UMI if in barcode-aware mode
                    if (readOneFile != null && readOneReader != null) {
                        if (readOneReader.hasNext()) {
                            FastqRecord readOneRecord = readOneReader.next();
                            ReadOneParser.BarcodeUMI extracted = ReadOneParser.extract(readOneRecord);
                            if (extracted != null) {
                                batchBarcodes.add(extracted.barcode);
                                batchUMIs.add(extracted.umi);
                            } else {
                                System.err.println("Warning: Could not extract barcode/UMI from read: " + readOneRecord.getReadName());
                                batchBarcodes.add("UNKNOWN");
                                batchUMIs.add("UNKNOWN");
                            }
                        } else {
                            throw new IOException("ReadOne file has fewer reads than ReadTwo file");
                        }
                    }

                    // Process batch when it reaches target size
                    if (batch.size() >= BATCH_SIZE) {
                        if (readOneFile != null) {
                            pgq.processAllWithBarcodes(batch, batchBarcodes, batchUMIs);
                            batchBarcodes.clear();
                            batchUMIs.clear();
                        } else {
                            pgq.processBatch(batch);
                        }
                        totalProcessed += batch.size();
                        batch.clear();

                        // Report progress only after actual processing completes
                        while (totalProcessed >= lastReportedMilestone + reportingStep) {
                            lastReportedMilestone += reportingStep;
                            System.out.println("Reads mapped and counted: " + (lastReportedMilestone / reportingStep) + " batch");
                        }
                    }
                }

                // Process remaining reads
                if (!batch.isEmpty()) {
                    if (readOneFile != null) {
                        pgq.processAllWithBarcodes(batch, batchBarcodes, batchUMIs);
                    } else {
                        pgq.processBatch(batch);
                    }
                    totalProcessed += batch.size();
                }

                // Final progress report
                System.out.println("All " + totalProcessed + " reads mapped and counted");

                // Shutdown executor after all batches
                pgq.shutdown();
            } catch (Exception e) {
                throw new RuntimeException(e);
            }

            // Write output
            if (readOneFile != null) {
                // Barcode-aware mode: output as sparse matrix
                CountMatrixWriter cmw = new CountMatrixWriter(countMatrixFile, idxData.getInt2TxString(), 
                                                              idxData.getInt2GeneString(), idxData.getTxInt2GeneInt());
                cmw.writeBarcodeSparseMatrix(pgq.getGeneCountsPerBarcode(), pgq.getTxCountsPerBarcode(), 
                                           pgq.getBarcodeMapper(), outputDir);
                System.out.println("Wrote sparse matrix (Market Matrix format) to: " + outputDir);
            } else {
                // Legacy mode: output counts
                Int2IntOpenHashMap geneInt2Counts = pgq.getGlobalGeneCounts();
                Int2IntOpenHashMap txInt2Counts = pgq.getGlobalTxCounts();
                CountMatrixWriter cmw = new CountMatrixWriter(countMatrixFile, idxData.getInt2TxString(), 
                                                              idxData.getInt2GeneString(), idxData.getTxInt2GeneInt());
                cmw.write(geneInt2Counts, txInt2Counts);
                System.out.println("Wrote counts to: " + countMatrixFile);
            }
        }
        else {
            System.err.println("Could not identify command: " + args[0]);
            printHelp();
            System.exit(1);
        }

        long endTime = System.nanoTime();
        System.out.printf("Runtime: %.3f secs%n", (double) (endTime - startTime) / 1_000_000_000);
    }



    private static void argParserIndex(String[] args) {
        if (args.length == 1 || args[1].equals("-h") || args[1].equals("--help")) {
            printHelp();
            System.exit(0);
        }
        for (int i = 1; i < args.length; i++) {
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

                case "-g":
                case "-gtf":
                    if (i + 1 < args.length) {
                        gtfFile = new File(args[++i]);
                    } else {
                        System.err.println("Error: [-g | --gtf] Please specify a GTF!");
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

                case "-idx":
                    if (i + 1 < args.length) {
                        indexFile = new File(args[++i]);
                    } else {
                        System.err.println("Error: [-idx] Please specify an index file to be written!");
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

    public static void argParserCount(String[] args) {
        if (args.length == 1 || args[1].equals("-h") || args[1].equals("--help")) {
            printHelp();
            System.exit(0);
        }
        for (int i = 1; i < args.length; i++) {
            switch (args[i]) {
                case "-r1":
                    if (i + 1 < args.length) {
                        readOneFile = new File(args[++i]);
                    } else {
                        System.err.println("Error: [-r1] Please specify a read file in FASTQ format!");
                        System.exit(1);
                    }
                    break;
                case "-r2":
                    if (i + 1 < args.length) {
                        readTwoFile = new File(args[++i]);
                    } else {
                        //System.err.println("Error: [-r2] Please specify a read file in FASTQ format!");
                        //System.exit(1);
                    }
                    break;
                case "-o": // directory
                    if (i + 1 < args.length) {
                        outputDir = new File(args[++i]);
                    } else {
                        System.err.println("Error: [-o] Please specify an output directory!");
                    }
                    break;
                case "-idx":
                    if (i + 1 < args.length) {
                        indexFile = new File(args[++i]);
                    } else {
                        System.err.println("Error: [-idx] Please specify an index file!");
                    }
                    break;
                case "-threads":
                    if (i + 1 < args.length) {
                        threads = Integer.parseInt(args[++i]);
                    } else {
                        System.err.println("Error: [-threads] Please specify a number of threads!");
                    }
                    break;
                case "-batchSize":
                    if (i + 1 < args.length) {
                        batch_size = Integer.parseInt(args[++i]);
                    } else {
                        System.err.println("Error: [-batchSize] Please specify a batch size!");
                    }
                    break;

                default:
                    System.err.println("Unknown argument: " + args[i]);
                    System.err.println("Use -h or --help to see the parameter list.");
                    System.exit(1);
            }
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