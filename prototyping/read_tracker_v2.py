

import sys
import pysam
import argparse
import heapq
import itertools
import pyranges as pr
from bx.intervals.intersection import IntervalTree, Interval
import pandas as pd


interval_trees = {}
padding = int
mapq = int
tracked_reads = set()
seen_in_a = set()
seen_in_b = set()
not_present_in_a = set()
not_present_in_b = set()
paralogous = set()
repeat = set()
filtered_by_quality_in_a = set()
filtered_by_quality_in_b = set()
read_status = {}
gtf_annotation = pr.PyRanges()
repeat_annotation = pr.PyRanges()
priority_counter = itertools.count()
skip_intergenic = False
paralog_table: dict[str, set[str]] = {}

class Locus:
    chromosome: str
    start: int
    end: int
    interval: Interval
    reads_in_a: set[str]
    reads_in_b: set[str]
    associated_genes: set[str]
    associated_transcripts: set[str]
    read_count: int

    def __init__(self, contig, start, end):
        self.chromosome = contig
        self.start = start
        self.end = end
        self.interval = Interval(chrom=contig, start=start, end=end)
        global gtf_annotation
        genes = query_annotation(gtf_annotation, contig, start, end, "gene_id")
        transcripts = query_annotation(gtf_annotation, contig, start, end, "transcript_id")
        self.associated_genes = genes
        self.associated_transcripts = transcripts
        self.reads_in_a = set()
        self.reads_in_b = set()

    def add_read(self, read: pysam.AlignedSegment, source_side):
        if source_side == "a":
            reads = self.reads_in_a
        else:
            reads = self.reads_in_b
        if read.query_name in reads:
            return False
        reads.add(read.query_name)
        self.start = min(self.start, read.reference_start)
        self.end = max(self.end, read.reference_end)
        return True

    #TODO: make find/has_association method


def parse_args():
    parser = argparse.ArgumentParser(description="Read tracker")
    parser.add_argument("-a", "--bamA", help="BAM file A", required=True)
    parser.add_argument("-b", "--bamB", help="BAM file B", required=True)
    parser.add_argument("-s", "--start_range", help="Comma-separated list of gene IDs (e.g. ENSG00000001) or transcript IDs (e.g. ENST00000001) to use as traversal seed loci", required=True)
    parser.add_argument("-g", "--gtf", help="Path to GTF annotation file", required=True)
    parser.add_argument("-r", "--repeats", help="Path to RepeatMasker .out file", required=True)
    parser.add_argument("--padding", help="Padding around loci for considering reads as mapped to the locus", type=int, default=500)
    parser.add_argument("--quality_filter", help="Minimum mapping quality for considering a read as mapped (default: 255, i.e. only uniquely mapped reads)", type=int, default=255)
    parser.add_argument("--skip_inter_genic", help="Whether to skip loci that do not overlap any gene annotation", action="store_true")
    parser.add_argument("--paralog_table", help="Path to paralog lookup table TSV (gene_id<TAB>paralog1,paralog2,...) produced by extract_paralogs.sh", required=True)
    return parser.parse_args()


def load_repeats(repeats_file):
    global repeat_annotation
    df = pd.read_csv(
        repeats_file,
        sep=r"\s+",
        skiprows=3,
        header=None,
        engine="python"
    ).rename(columns={
        4: "Chromosome",
        5: "Start",
        6: "End",
        9: "repeat_name",
        10: "repeat_class"
    })[["Chromosome", "Start", "End", "repeat_name", "repeat_class"]]
    df["Start"] -= 1
    df["Feature"] = "repeat"
    df["Strand"] = "."
    repeat_annotation = pr.PyRanges(df)


def load_paralog_table(filepath):
    global paralog_table
    paralog_table = {}
    with open(filepath) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            paralog_table[parts[0]] = set(parts[1].split(","))
    print(f"Loaded paralog table: {len(paralog_table)} genes")


def is_paralogs(source_locus, target_locus):
    for gene in source_locus.associated_genes:
        if paralog_table.get(gene, set()) & target_locus.associated_genes:
            return True
    return False


def is_repeat(source_locus, target_locus):
    global repeat_annotation
    source_hits = query_annotation(
        repeat_annotation,
        source_locus.chromosome,
        source_locus.start,
        source_locus.end,
        "repeat_name")
    target_hits = query_annotation(
        repeat_annotation,
        target_locus.chromosome,
        target_locus.start,
        target_locus.end,
        "repeat_name")
    return bool(source_hits & target_hits)

#return only locis where reads are mapped in on side but not the other and not explained


def mark_mapped_elsewhere(read_id, source_side):
    read_status[read_id] = "mapped_elsewhere"
    filtered_by_quality_in_a.discard(read_id)
    filtered_by_quality_in_b.discard(read_id)
    if source_side == "a":
        not_present_in_a.discard(read_id)
        seen_in_a.add(read_id)
    else:
        not_present_in_b.discard(read_id)
        seen_in_b.add(read_id)


def mark_filtered_in_other(read_id, source_side):
    if read_status.get(read_id) == "mapped_elsewhere":
        return
    read_status[read_id] = "filtered_by_quality"
    if source_side == "a":
        filtered_by_quality_in_a.add(read_id)
    else:
        filtered_by_quality_in_b.add(read_id)


def mark_not_present_in_other(read_id, source_side):
    if read_status.get(read_id) == "mapped_elsewhere":
        return
    read_status[read_id] = "not_present_in_other_bam"
    filtered_by_quality.discard(read_id)
    if source_side == "a":
        not_present_in_a.add(read_id)
    else:
        not_present_in_b.add(read_id)


def get_mapping_diff(locus: Locus, bam_a, bam_b, indices):
    # maybe add quality filter ?
    # consider strand specificity ?
    bam_a_reads = bam_a.fetch(locus.chromosome, locus.start, locus.end)
    bam_b_reads = bam_b.fetch(locus.chromosome, locus.start, locus.end)
    bam_a_ids = set()
    bam_b_ids = set()
    for read_a in  bam_a_reads:
        if quality_filter(read_a):
            locus.add_read(read_a, "a")
            bam_a_ids.add(read_a.query_name)
    for read_b in bam_b_reads:
        if quality_filter(read_b):
            locus.add_read(read_b, "b")
            bam_b_ids.add(read_b.query_name)
    #reads only in a or b that maps to the current locus look in the other for them
    only_in_a = bam_a_ids - bam_b_ids
    only_in_b = bam_b_ids - bam_a_ids
    tracked_reads.update(only_in_a)
    tracked_reads.update(only_in_b)
    #locus where reads are mapped to in a that appeared 'here' b and vice versa
    target_loci_a = get_loci(only_in_b, "a", indices, locus)
    target_loci_b = get_loci(only_in_a, "b", indices, locus)
    return target_loci_a, target_loci_b

# given a set of read ids, get the loci they map to in the other bam file and return those loci


def get_loci(read_ids, source_side, indices, source_locus):
    loci = set()
    seen_in_other = seen_in_b if source_side == "a" else seen_in_a
    for read_id in read_ids:
        if read_id in seen_in_other:
            continue
        index = indices[source_side]
        try:
            read_iter = index.find(read_id)
        except KeyError:
            mark_not_present_in_other(read_id, source_side)
            continue
        found_passing_alignment = False
        for read in read_iter:
            if not quality_filter(read):
                continue
            found_passing_alignment = True
            mark_mapped_elsewhere(read_id, source_side)
            locus = get_locus(read, source_side, indices)
            loci.add(locus)
        if not found_passing_alignment:
            mark_filtered_in_other(read_id, source_side)
    to_be_removed = set()
    for locus in loci:
        locus_reads = locus.reads_in_a if source_side == "a" else locus.reads_in_b
        tracked_locus_reads = locus_reads & read_ids
        if is_paralogs(source_locus, locus):
            paralogous.update(tracked_locus_reads)
            to_be_removed.add(locus)
        if is_repeat(source_locus, locus):
            repeat.update(tracked_locus_reads)
            to_be_removed.add(locus)
    loci = loci - to_be_removed

    return loci

# given read, query the interval trees to find the locus it maps to, if any. If it does not map to any existing
# locus, create a new one and add it to the tree


def get_locus(read, source_side, indices):
    contig = read.reference_name
    tree = interval_trees.get(contig, None)
    if tree is None:
        tree = IntervalTree()
        interval_trees[contig] = tree
    candidates_interval = tree.find(read.reference_start - padding, read.reference_end + padding)
    if candidates_interval:
        #don't know if thats necessary
        best_interval = max(
            candidates_interval,
            key=lambda interval: min(interval.end, read.reference_end) - max(interval.start, read.reference_start)
        )
        locus = best_interval
    else:
        locus = Locus(contig, read.reference_start, read.reference_end)
    locus.add_read(read, source_side)
    add_locus_to_tree(contig, locus)
    return locus


def quality_filter(read):
    #make this as argument to be specified later
    if read.mapping_quality < mapq:
        return False
    if read.is_secondary or read.is_supplementary:
        return False
    if read.is_unmapped:
        return False
    if read.cigartuples:
        ##check if read contains a junction that is to big
        for op, length in read.cigartuples:
            if op == 3 and length > 1.5e4:
                return False
    return True


def query_annotation(annotation, contig, start, end, feature):
    if contig not in annotation.Chromosome.unique():

        return set()
    query = pr.from_dict({
        "Chromosome": [contig],
        "Start": [start],
        "End": [end]
    })
    hits = annotation.overlap(query).df
    if feature not in hits.columns:
        return set()
    return set(hits[feature].dropna())

def add_locus_to_tree(contig, locus:Locus):
    if contig not in interval_trees:
        interval_trees[contig] = IntervalTree()
    interval_trees[contig].insert(locus.start, locus.end, locus)


def load_parameters(args):
    """Load BAM files, annotations, and build indices"""
    bam_a = pysam.AlignmentFile(args.bamA, "rb", require_index=True)
    bam_b = pysam.AlignmentFile(args.bamB, "rb", require_index=True)

    print("loading gtf annotations...")
    global gtf_annotation
    gtf = pr.read_gtf(args.gtf)
    gtf_annotation = gtf[gtf.Feature.isin(["gene", "transcript", "exon"])]
    print("loading repeats....")
    load_repeats(args.repeats)

    print("building indices...")
    name_index_a = pysam.IndexedReads(bam_a, multiple_iterators=False)
    name_index_a.build()
    print("index for bamA built")
    name_index_b = pysam.IndexedReads(bam_b, multiple_iterators=False)
    name_index_b.build()
    print("index for bamB built")

    indices = {
        "a": name_index_a,
        "b": name_index_b
    }
    global padding
    padding = args.padding
    global mapq
    mapq = args.quality_filter
    if args.skip_inter_genic:
        global skip_intergenic
        skip_intergenic = True
    print("loading paralog table...")
    load_paralog_table(args.paralog_table)
    return bam_a, bam_b, indices

def make_priority_entry(locus, source_side):
    read_count = len(locus.reads_in_a) if source_side == "a" else len(locus.reads_in_b)
    return (-read_count, next(priority_counter), locus)

def push_to_queue(loci_in_a, loci_in_b, priority_queue):
    queue_a = []
    #pushes loci in two heaps based on number of reads in the respective side
    #and merges into one heap later
    for locus in loci_in_a:
        heapq.heappush(queue_a, make_priority_entry(locus, "a"))
    queue_b = []
    for locus in loci_in_b:
        heapq.heappush(queue_b, make_priority_entry(locus, "b"))
    while queue_a:
        heapq.heappush(priority_queue, heapq.heappop(queue_a))
    while queue_b:
        heapq.heappush(priority_queue, heapq.heappop(queue_b))


def resolve_ids_to_loci(id_list):
    global gtf_annotation
    intervals = []
    missing_ids = set()
    df = gtf_annotation.df
    for id in id_list:
        row = None
        if "transcript_id" in df.columns and "Feature" in df.columns:
            match = df[(df["transcript_id"] == id) & (df["Feature"] == "transcript")]
            if not match.empty:
                row = match.iloc[0]
        if row is None and "gene_id" in df.columns and "Feature" in df.columns:
            match = df[(df["gene_id"] == id) & (df["Feature"] == "gene")]
            if not match.empty:
                row = match.iloc[0]
        if row is None:
            missing_ids.add(id)
            continue
        intervals.append((row["Chromosome"], int(row["Start"]), int(row["End"])))
    return intervals, missing_ids


def main():
    args = parse_args()
    bam_a, bam_b, indices= load_parameters(args)
    global tracked_reads, seen_in_a, seen_in_b, not_present_in_a, not_present_in_b, paralogous, repeat, filtered_by_quality_in_a, filtered_by_quality_in_b, read_status
    tracked_reads = set()
    seen_in_a = set()
    seen_in_b = set()
    not_present_in_a = set()
    not_present_in_b = set()
    paralogous = set()
    repeat = set()
    filtered_by_quality_in_a = set()
    filtered_by_quality_in_b = set()
    read_status = {}
    ids = [s.strip() for s in args.start_range.split(",")]
    intervals, missing = resolve_ids_to_loci(ids)
    if len(missing) == len(ids):
        print(f"Error: none of the {len(ids)} provided IDs were found in the annotation: {missing}", file=sys.stderr)
        sys.exit(1)
    if missing:
        print(f"Warning: {len(missing)}/{len(ids)} IDs not found in annotation: {missing}")
    seeds = [Locus(chrom, start, end) for chrom, start, end in intervals]
    priority_queue = []
    for seed in seeds:
        add_locus_to_tree(seed.chromosome, seed)
        heapq.heappush(priority_queue, make_priority_entry(seed, "a"))
    visited = set()
    while priority_queue:
        neg_count, _, locus = heapq.heappop(priority_queue)
        if locus in visited:
            continue
        visited.add(locus)
        print('visiting locus:', locus.chromosome, locus.start, locus.end, f'{locus.associated_genes}', f'{locus.associated_transcripts}')
        loci_in_a, loci_in_b = get_mapping_diff(locus, bam_a, bam_b, indices)
        push_to_queue(loci_in_a, loci_in_b, priority_queue)

    total_tracked_reads = len(tracked_reads)
    mapped_elsewhere = {read_id for read_id, status in read_status.items() if status == "mapped_elsewhere"}
    missing_in_other = {read_id for read_id, status in read_status.items() if status == "not_present_in_other_bam"}
    filtered_reads = filtered_by_quality_in_a | filtered_by_quality_in_b
    classified_reads = mapped_elsewhere | missing_in_other | filtered_reads
    unclassified_reads = tracked_reads - classified_reads

    print("total visited loci:", len(visited))
    print("tracked reads:", total_tracked_reads)
    print("mapped else where in a:", len(seen_in_a))
    print("mapped else where in b:", len(seen_in_b))
    if total_tracked_reads:
        if unclassified_reads:
            print(f"warning: {len(unclassified_reads)} tracked reads were not assigned a terminal status")
        print(f"mapped else where {len(mapped_elsewhere)} ({len(mapped_elsewhere)/total_tracked_reads*100:.2f}%)")
        print(f"not present in other bam {len(missing_in_other)} ({len(missing_in_other)/total_tracked_reads*100:.2f}%)")
        print(f"filtered by quality (total) {len(filtered_reads)} ({len(filtered_reads)/total_tracked_reads*100:.2f}%)")
        print(f"  filtered in bam A {len(filtered_by_quality_in_a)} ({len(filtered_by_quality_in_a)/total_tracked_reads*100:.2f}%)")
        print(f"  filtered in bam B {len(filtered_by_quality_in_b)} ({len(filtered_by_quality_in_b)/total_tracked_reads*100:.2f}%)")
        print(f"paralogous {len(paralogous)} ({len(paralogous)/total_tracked_reads*100:.2f}%)")
        print(f"repeat {len(repeat)} ({len(repeat)/total_tracked_reads*100:.2f}%)")
    else:
        print("mapped else where 0 (0.00%)")
        print("not present in other bam 0 (0.00%)")
        print("filtered by quality (total) 0 (0.00%)")
        print("  filtered in bam A 0 (0.00%)")
        print("  filtered in bam B 0 (0.00%)")
        print("paralogous 0 (0.00%)")
        print("repeat 0 (0.00%)")
if __name__ == "__main__":
    main()
