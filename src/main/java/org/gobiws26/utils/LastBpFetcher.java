package org.gobiws26.utils;

import org.gobiws26.Readers.GTFReader;
import org.gobiws26.genomicstruct.Exon;
import org.gobiws26.genomicstruct.Transcript;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map;

public class LastBpFetcher {
    public static void main(String[] args) {
        final Path gtfPath = args.length > 0
                ? Path.of(args[0])
                : Path.of("data", "10x.genes.gtf");
        final Path outBedPath = Path.of("data", "human_bedfiles");

        try {
            Map<String, Transcript> transcripts = new GTFReader().read(gtfPath.toFile());
            for(int k = 100; k <= 3000; k += 150) {
                Path outBedPathForK = Path.of(outBedPath.toString(), "last_" + k + "_bp.bed");
                new LastBpFetcher().writeLastKBpsAsBed(transcripts, k, outBedPathForK);
                System.out.println("Wrote BED: " + outBedPathForK.toAbsolutePath());
            }
        } catch (IOException e) {
            throw new RuntimeException("Failed to create BED from GTF: " + gtfPath.toAbsolutePath(), e);
        }
    }

    public record BedEntry(String chrom, int chromStart, int chromEnd, String name, int score, char strand) {
        public String toBedLine() {
            return chrom + "\t" + chromStart + "\t" + chromEnd + "\t" + name + "\t" + score + "\t" + strand;
        }
    }

    public List<Exon> getLastKBpsOfTranscript(Transcript transcript, int k) {
        if (k <= 0) return List.of();

        List<Exon> exonsInTxOrder = new ArrayList<>(transcript.getSortedExons());
        if (transcript.isNegativeStranded()) {
            exonsInTxOrder = exonsInTxOrder.reversed();
        }

        int totalWidth = exonsInTxOrder.stream().mapToInt(Exon::length).sum();

        int excess = totalWidth - k;
        int cumulative = 0;
        int firstKeepIndex = -1;
        int widthBeforeFirstKept = 0;

        for (int i = 0; i < exonsInTxOrder.size(); i++) {
            cumulative += exonsInTxOrder.get(i).length();
            if (cumulative > excess) {
                firstKeepIndex = i;
                widthBeforeFirstKept = cumulative - exonsInTxOrder.get(i).length();
                break;
            }
        }

        if (firstKeepIndex < 0) return List.of();

        List<Exon> selected = new ArrayList<>();
        for (int i = firstKeepIndex; i < exonsInTxOrder.size(); i++) {
            Exon exon = exonsInTxOrder.get(i);
            selected.add(new Exon(exon.getStart(), exon.getEnd()));
        }

        int trimAmount = excess - widthBeforeFirstKept;
        if (trimAmount > 0) {
            Exon boundaryExon = selected.getFirst();
            if (transcript.isNegativeStranded()) {
                selected.set(0, new Exon(boundaryExon.getStart(), boundaryExon.getEnd() - trimAmount));
            } else {
                selected.set(0, new Exon(boundaryExon.getStart() + trimAmount, boundaryExon.getEnd()));
            }
        }

        return selected;
    }

    public List<BedEntry> toBedEntries(Transcript transcript, String transcriptId, int k) {
        List<Exon> lastK = new ArrayList<>(getLastKBpsOfTranscript(transcript, k));
        lastK.sort(Comparator.comparingInt(Exon::getStart));

        char strand = transcript.isNegativeStranded() ? '-' : '+';
        List<BedEntry> bedEntries = new ArrayList<>(lastK.size());
        for (Exon exon : lastK) {
            // BED is 0-based start, half-open end.
            bedEntries.add(new BedEntry(
                    transcript.getChr(),
                    exon.getStart() - 1,
                    exon.getEnd(),
                    transcriptId,
                    0,
                    strand
            ));
        }
        return bedEntries;
    }

    public void writeLastKBpsAsBed(Map<String, Transcript> transcripts, int k, Path outputBedFile) throws IOException {
        try (BufferedWriter writer = Files.newBufferedWriter(outputBedFile)) {
            transcripts.entrySet().stream()
                    .sorted(Map.Entry.comparingByKey())
                    .forEach(entry -> {
                        List<BedEntry> bedEntries = toBedEntries(entry.getValue(), entry.getKey(), k);
                        for (BedEntry bedEntry : bedEntries) {
                            try {
                                writer.write(bedEntry.toBedLine());
                                writer.newLine();
                            } catch (IOException e) {
                                throw new RuntimeException(e);
                            }
                        }
                    });
        } catch (RuntimeException e) {
            if (e.getCause() instanceof IOException ioException) {
                throw ioException;
            }
            throw e;
        }
    }

}
