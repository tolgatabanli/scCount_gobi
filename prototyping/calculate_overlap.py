#!/usr/bin/env python3
import pysam 
import argparse
from collections import defaultdict

def parse_gtf(gtf_file):

    transcripts = defaultdict(list)

    with open(gtf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")

            if fields[2] != "exon":
                continue

            chrom = fields[0]
            start = int(fields[3]) - 1
            end = int(fields[4])

            attrs = fields[8]

            tid = None
            for attr in attrs.split(";"):
                attr = attr.strip()
                if attr.startswith("transcript_id"):
                    tid = attr.split('"')[1]

            if tid:
                transcripts[tid].append((chrom, start, end))

    return transcripts
#!/usr/bin/env python3

import pysam
import argparse
from collections import defaultdict


# ------------------------------
# Read exon coordinates per transcript
# ------------------------------
def parse_gtf(gtf_file):

    transcripts = defaultdict(list)

    with open(gtf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")

            if fields[2] != "exon":
                continue

            chrom = fields[0]
            start = int(fields[3]) - 1
            end = int(fields[4])

            attrs = fields[8]

            tid = None
            for attr in attrs.split(";"):
                attr = attr.strip()
                if attr.startswith("transcript_id"):
                    tid = attr.split('"')[1]

            if tid:
                transcripts[tid].append((chrom, start, end))

    return transcripts


# ------------------------------
# overlap helper
# ------------------------------
def read_overlap_fraction(read, exons):

    overlap = 0
    read_len = read.reference_length

    for bstart, bend in read.get_blocks():
        for chrom, estart, eend in exons:
            if chrom != read.reference_name:
                continue

            overlap += max(0, min(bend, eend) - max(bstart, estart))

    if read_len == 0:
        return 0

    return overlap / read_len


# ------------------------------
# read encoder
# ------------------------------
class ReadEncoder:

    def __init__(self):
        self.index = {}
        self.next_id = 0

    def encode(self, name):

        if name not in self.index:
            self.index[name] = self.next_id
            self.next_id += 1

        return self.index[name]


# ------------------------------
# extract reads per transcript
# ------------------------------
def reads_per_transcript(bam_path, transcripts, encoder):

    bam = pysam.AlignmentFile(bam_path)

    tx_reads = defaultdict(set)

    for tid, exons in transcripts.items():

        chrom = exons[0][0]

        start = min(e[1] for e in exons)
        end = max(e[2] for e in exons)

        for read in bam.fetch(chrom, start, end):

            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            
            rid = encoder.encode(read.query_name)
            tx_reads[tid].add(rid)

    bam.close()

    return tx_reads


# ------------------------------
# compare transcripts
# ------------------------------
def compare(b1, b2):

    results = []

    all_tx = set(b1) | set(b2)

    for tx in all_tx:

        r1 = b1.get(tx, set())
        r2 = b2.get(tx, set())

        shared = len(r1 & r2)
        unique1 = len(r1 - r2)
        unique2 = len(r2 - r1)

        results.append((tx, shared, unique1, unique2))

    return results


# ------------------------------
# main
# ------------------------------
def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("bam1")
    parser.add_argument("bam2")
    parser.add_argument("gtf")
    parser.add_argument("output")

    args = parser.parse_args()

    print("Loading transcripts from GTF...")
    transcripts = parse_gtf(args.gtf)

    encoder = ReadEncoder()

    print("Scanning BAM1...")
    bam1_reads = reads_per_transcript(
        args.bam1,
        transcripts,
        encoder
    )

    print("Scanning BAM2...")
    bam2_reads = reads_per_transcript(
        args.bam2,
        transcripts,
        encoder
    )

    print("Comparing transcripts...")
    results = compare(bam1_reads, bam2_reads)

    print("Writing output...")

    with open(args.output, "w") as out:

        out.write(
            "transcript_id\tshared_reads\tunique_bam1\tunique_bam2\n"
        )

        for tx, s, u1, u2 in results:
            out.write(f"{tx}\t{s}\t{u1}\t{u2}\n")


if __name__ == "__main__":
    main()