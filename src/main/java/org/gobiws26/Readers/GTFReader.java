package org.gobiws26.Readers;

import org.gobiws26.genomicstruct.Exon;
import org.gobiws26.genomicstruct.Transcript;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

public class GTFReader {
    static final Pattern patternTranscriptID = Pattern.compile("transcript_id\\s+\"([^\"]+)\"");
    static final Pattern patternGeneID = Pattern.compile("gene_id\\s+\"([^\"]+)\"");
    private final ArrayList<String> int2Transcript;

    public GTFReader() {
        this.int2Transcript = new ArrayList<>(); // might add getter if created this way
    }

    public GTFReader(ArrayList<String> int2Transcript) {
        this.int2Transcript = int2Transcript;
    }

    public HashMap<String, Transcript> read(File gtfFile) throws IOException {
        HashMap<String, Transcript> transcripts = new HashMap<>();
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(openGtfStream(gtfFile)))) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("#")) continue;

                String[] fields = line.split("\t");
                String chr = fields[0];
                String feature = fields[2];
                int start = Integer.parseInt(fields[3]);
                int end = Integer.parseInt(fields[4]);
                boolean isNegative = fields[6].equals("-"); // defaults to positive if dot

                // TODO decide features better?
                if (!feature.equals("transcript") && !feature.equals("exon") && !feature.equals("three_prime_utr")) continue; // TODO: integrate 3' UTR

                Matcher transcriptMatcher = patternTranscriptID.matcher(fields[8]);
                Matcher geneMatcher = patternGeneID.matcher(fields[8]);
                String transcriptId;
                String geneId;

                Transcript t;
                if (transcriptMatcher.find()) {
                    transcriptId = transcriptMatcher.group(1);
                    t = transcripts.computeIfAbsent(transcriptId, k -> {
                            Transcript k = new Transcript();
                            int2Transcript.add(transcriptId);
                            return k;
                    }); // do not set info if exon comes before
                } else continue;
                if (geneMatcher.find()) {
                    geneId = geneMatcher.group(1);
                    t.setGene(geneId);
                }

                // set coordinates
                if (feature.equals("transcript")) {
                    t.setStart(start);
                    t.setEnd(end);
                    t.setChr(chr);
                    t.setStrand(isNegative);
                } else if (feature.equals("exon")) {
                    t.addExon(new Exon(start, end));
                } else { // UTR
                    t.addToUTR(end - start + 1); // closed interval
                }
            }
        }
        return transcripts;
    }

    private InputStream openGtfStream(File gtfFile) throws IOException {
        FileInputStream fis = new FileInputStream(gtfFile);
        if (gtfFile.getName().endsWith(".gz")) {
            return new GZIPInputStream(fis);
        }
        return fis;
    }
}
