package org.gobiws26.Readers;

import org.gobiws26.genomicstruct.Exon;
import org.gobiws26.genomicstruct.Transcript;

import java.io.*;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

public class GTFReader {
    static final Pattern patternTranscriptID = Pattern.compile("transcript_id\\s+\"([^\"]+)\"");

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
                if (!feature.equals("transcript") && !feature.equals("exon")) continue;

                Matcher matcher = patternTranscriptID.matcher(fields[8]);
                String transcriptId;



                Transcript t;
                if (matcher.find()) {
                    transcriptId = matcher.group(1);
                    t = transcripts.computeIfAbsent(transcriptId, k -> new Transcript()); // do not set info if exon comes before
                } else continue;

                // set coordinates
                if (feature.equals("transcript")) {
                    t.setStart(start);
                    t.setEnd(end);
                    t.setChr(chr);
                    t.setStrand(isNegative);
                } else {
                    t.addExon(new Exon(start, end));
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
