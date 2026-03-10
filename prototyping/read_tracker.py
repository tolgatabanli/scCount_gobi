#./usr/bin/env python
import pysam
import argparse
import heapq
from dataclasses import dataclass, field
from typing import Dict, List, Set, Tuple

import pandas as pd
import pyranges as pr
from bx.intervals.intersection import IntervalTree, Interval


@dataclass
class Locus:
    chromosome: str
    start: int
    end: int
    reads: List[pysam.AlignedSegment] = field(default_factory=list)
    associated_genes: Set[str] = field(default_factory=set)
    associated_transcripts: Set[str] = field(default_factory=set)
    data: Dict[str, object] = field(
        default_factory=lambda: {"read_count": 0, "is_repeat": False, "is_paralogue": False}
    )

    def add_read(self, read: pysam.AlignedSegment):
        if read not in self.reads:
            self.reads.append(read)
        self.start = min(self.start, read.reference_start)
        self.end = max(self.end, read.reference_end)
        self.data["read_count"] = len(self.reads)

    def set_associations(self, genes: Set[str], transcripts: Set[str]):
        self.associated_genes = set(genes)
        self.associated_transcripts = set(transcripts)

    def has_association_overlap(self, feature_start: int, feature_end: int) -> bool:
        # Internal coordinates are 0-based, end-exclusive.
        return not (feature_end <= self.start or feature_start >= self.end)


def load_gtf_annotations(gtf_path: str):
    gtf_df = pr.read_gtf(gtf_path).df.copy()
    gtf_df["Chromosome"] = gtf_df["Chromosome"].astype(str).str.replace("^chr", "", regex=True)

    gene_features = gtf_df[gtf_df["gene_id"].notna()][["Chromosome", "Start", "End", "gene_id"]]
    transcript_features = gtf_df[gtf_df["transcript_id"].notna()][
        ["Chromosome", "Start", "End", "transcript_id"]
    ]

    genes_gr = pr.PyRanges(gene_features) if not gene_features.empty else None
    transcripts_gr = pr.PyRanges(transcript_features) if not transcript_features.empty else None
    return genes_gr, transcripts_gr


def build_locus_ranges(loci_by_chr: Dict[str, Dict[Tuple[int, int], Locus]]) -> pr.PyRanges:
    rows = []
    for chromosome, chr_payload in loci_by_chr.items():
        for locus in chr_payload.get("loci", {}).values():
            rows.append(
                {
                    "Chromosome": chromosome,
                    "Start": locus.start,
                    "End": locus.end,
                    "locus_key": f"{chromosome}:{locus.start}-{locus.end}",
                }
            )

    if not rows:
        return pr.from_dict({"Chromosome": [], "Start": [], "End": [], "locus_key": []})

    return pr.from_dict(
        {
            "Chromosome": [row["Chromosome"] for row in rows],
            "Start": [row["Start"] for row in rows],
            "End": [row["End"] for row in rows],
            "locus_key": [row["locus_key"] for row in rows],
        }
    )


def associate_loci_with_annotations(
    loci_by_chr: Dict[str, Dict[Tuple[int, int], Locus]],
    genes_gr,
    transcripts_gr,
):
    loci_gr = build_locus_ranges(loci_by_chr)
    if loci_gr.df.empty:
        return

    locus_by_key = {}
    for chromosome, chr_payload in loci_by_chr.items():
        for locus in chr_payload.get("loci", {}).values():
            locus_by_key[f"{chromosome}:{locus.start}-{locus.end}"] = locus

    gene_map: Dict[str, Set[str]] = {}
    transcript_map: Dict[str, Set[str]] = {}

    if genes_gr is not None:
        gene_overlap_df = loci_gr.join(genes_gr).df
        if not gene_overlap_df.empty:
            for row in gene_overlap_df.itertuples(index=False):
                locus_key = row.locus_key
                locus = locus_by_key[locus_key]
                if locus.has_association_overlap(row.Start_b, row.End_b):
                    gene_map.setdefault(locus_key, set()).add(row.gene_id)

    if transcripts_gr is not None:
        transcript_overlap_df = loci_gr.join(transcripts_gr).df
        if not transcript_overlap_df.empty:
            for row in transcript_overlap_df.itertuples(index=False):
                locus_key = row.locus_key
                locus = locus_by_key[locus_key]
                if locus.has_association_overlap(row.Start_b, row.End_b):
                    transcript_map.setdefault(locus_key, set()).add(row.transcript_id)

    for chromosome, chr_payload in loci_by_chr.items():
        updated_loci = {}
        for locus in chr_payload.get("loci", {}).values():
            key_str = f"{chromosome}:{locus.start}-{locus.end}"
            locus.set_associations(gene_map.get(key_str, set()), transcript_map.get(key_str, set()))
            updated_loci[(locus.start, locus.end)] = locus
        chr_payload["loci"] = updated_loci


def parse_args():
    parser = argparse.ArgumentParser(description="Read tracker")
    parser.add_argument("-a", "--bamA", help = "BAM file A", required = True)
    parser.add_argument("-b", "--bamB", help = "BAM file B", required = True)
    parser.add_argument("-s", "--start_range", help = "Start range for tracking", required = True)
    parser.add_argument("-g", "--gtf", help="Path to GTF annotation file", required=True)
    parser.add_argument("-r", "--repeatmasker", help="Path to RepeatMasker .out file", required=False)
    return parser.parse_args()


def load_repeatmasker_annotations(repeatmasker_path: str):
    repeats_df = pd.read_csv(
        repeatmasker_path,
        sep=r"\s+",
        header=None,
        skiprows=3,
        usecols=[4, 5, 6, 9, 10],
        names=["chromosome", "start", "end", "repeat_name", "repeat_class_family"],
        engine="python",
        on_bad_lines="skip",
    )

    if repeats_df.empty:
        return pd.DataFrame(
            columns=["chromosome", "start", "end", "repeat_name", "repeat_class_family"]
        )

    repeats_df["start"] = pd.to_numeric(repeats_df["start"], errors="coerce")
    repeats_df["end"] = pd.to_numeric(repeats_df["end"], errors="coerce")
    repeats_df = repeats_df.dropna(subset=["start", "end"])

    repeats_df["chromosome"] = repeats_df["chromosome"].astype(str).str.replace("^chr", "", regex=True)
    repeats_df["start"] = repeats_df["start"].astype(int) - 1
    repeats_df["end"] = repeats_df["end"].astype(int)

    repeats_index = {}
    for row in repeats_df.itertuples(index=False):
        chromosome = str(row.chromosome)
        if chromosome not in repeats_index:
            repeats_index[chromosome] = IntervalTree()
        repeats_index[chromosome].add_interval(
            Interval(row.start, row.end, (row.repeat_name, row.repeat_class_family))
        )

    return repeats_index

def get_mapping_diff(chr, start, end, bamA, bamB):
    #maybe add quality filter ? 
    #consider strand specificity ?
    bamA_reads = bamA.fetch(chr, start, end)
    bamB_reads = bamB.fetch(chr, start, end)

    bamA_ids = set(read.query_name for read in bamA_reads if quality_filter(read))
    bamB_ids = set(read.query_name for read in bamB_reads if quality_filter(read))
    only_in_A = bamA_ids - bamB_ids
    only_in_B = bamB_ids - bamA_ids
    return only_in_A, only_in_B

def quality_filter(read):
    #make this as argument to be specified later 
    if read.mapping_quality < 255:
        return False
    if read.is_secondary or read.is_supplementary:
        return False
    if read.cigartuples:
        ##check if read contains a junction that is to big 
        for op, length in read.cigartuples:
            if op == 3 and length > 1.5e4:
                return False
    return True

def add_read_to_loci(loci, chr, read_id, read, padding=500):
    if read.is_unmapped or read.reference_start is None or read.reference_end is None:
        return

    chr = read.reference_name.replace("chr", "")
    read_start = read.reference_start
    read_end = read.reference_end
    read_id = read.query_name
    if chr not in loci:
        loci[chr] = {"tree": IntervalTree(), "loci": {}}

    chromosome = loci[chr]
    tree = chromosome["tree"]
    locus_map = chromosome["loci"]

    def rebuild_tree_from_loci():
        rebuilt = IntervalTree()
        for start, end in locus_map.keys():
            rebuilt.add_interval(Interval(start, end))
        chromosome["tree"] = rebuilt
        return rebuilt

    def normalize_loci_with_padding():
        if not locus_map:
            chromosome["tree"] = IntervalTree()
            return

        sorted_loci = sorted(locus_map.values(), key=lambda locus: (locus.start, locus.end))
        merged_loci = [sorted_loci[0]]

        for locus in sorted_loci[1:]:
            current = merged_loci[-1]
            if locus.start <= current.end + padding:
                current.start = min(current.start, locus.start)
                current.end = max(current.end, locus.end)
                for entry in locus.reads:
                    current.add_read(entry)
                current.associated_genes.update(locus.associated_genes)
                current.associated_transcripts.update(locus.associated_transcripts)
            else:
                merged_loci.append(locus)

        locus_map.clear()
        for merged in merged_loci:
            merged.data["read_count"] = len(merged.reads)
            locus_map[(merged.start, merged.end)] = merged

        rebuild_tree_from_loci()

    # Query potentially matching loci in O(log n + k) instead of scanning all keys.
    candidates = tree.find(read_start - padding, read_end + padding)

    matching_interval = None
    for interval in candidates:
        # Guard against stale candidate intervals if the tree was rebuilt.
        if (interval.start, interval.end) not in locus_map:
            continue
        if read_start >= interval.start - padding and read_end <= interval.end + padding:
            matching_interval = interval
            break

    if matching_interval is None:
        new_interval = Interval(read_start, read_end)
        tree.add_interval(new_interval)
        new_locus = Locus(chromosome=chr, start=read_start, end=read_end, reads=[read])
        new_locus.data["read_count"] = 1
        locus_map[(read_start, read_end)] = new_locus
        normalize_loci_with_padding()
        return

    locus_start, locus_end = matching_interval.start, matching_interval.end
    updated_key = (min(locus_start, read_start), max(locus_end, read_end))

    existing_locus = locus_map.pop((locus_start, locus_end))
    existing_locus.add_read(read)
    existing_locus.start = updated_key[0]
    existing_locus.end = updated_key[1]

    if updated_key in locus_map:
        merged_locus = locus_map[updated_key]
        for entry in existing_locus.reads:
            merged_locus.add_read(entry)
        merged_locus.associated_genes.update(existing_locus.associated_genes)
        merged_locus.associated_transcripts.update(existing_locus.associated_transcripts)
    else:
        locus_map[updated_key] = existing_locus

    normalize_loci_with_padding()

def is_paralogous(query_locus, candidate_locus, padding=500):
    # TODO: implement sequence/annotation-based paralogue detection.
    return False

def is_repeat(query_locus, candidate_locus, repeats_index, padding=500):
    if not repeats_index:
        return False

    query_chr, query_start, query_end = query_locus
    locus_chr, locus_start, locus_end = candidate_locus

    query_tree = repeats_index.get(query_chr)
    candidate_tree = repeats_index.get(locus_chr)
    if query_tree is None or candidate_tree is None:
        return False

    query_hits = query_tree.find(query_start, query_end)
    candidate_hits = candidate_tree.find(locus_start, locus_end)
    if not query_hits or not candidate_hits:
        return False

    query_repeat_names = {interval.value[0] for interval in query_hits if interval.value}
    candidate_repeat_names = {interval.value[0] for interval in candidate_hits if interval.value}
    return bool(query_repeat_names & candidate_repeat_names)

def classify_and_queue_loci(
    locus_list,
    query_chr,
    query_start,
    query_end,
    repeats_index,
    relevant_read_ids,
    global_locus_queue,
    visited_loci,
    source_side,
    padding=500,
):
    queued = 0
    repeat_read_ids = set()
    repeat_loci = 0
    relevant_loci = 0
    query_locus = (query_chr, query_start, query_end)
    for locus_chr, chr_payload in locus_list.items():
        for locus_key, locus in chr_payload.get("loci", {}).items():
            locus_read_ids = {read.query_name for read in locus.reads}
            if not (locus_read_ids & relevant_read_ids):
                continue
            relevant_loci += 1

            locus_start, locus_end = locus_key
            candidate_locus = (locus_chr, locus_start, locus_end)
            locus.data["read_count"] = len(locus.reads)
            locus.data["is_repeat"] = is_repeat(query_locus, candidate_locus, repeats_index, padding=padding)
            if locus.data["is_repeat"]:
                repeat_loci += 1
                repeat_read_ids.update(locus_read_ids & relevant_read_ids)
            locus.data["is_paralogue"] = False
            if not locus.data["is_repeat"]:
                locus.data["is_paralogue"] = is_paralogous(query_locus, candidate_locus, padding=padding)

            if not locus.data["is_repeat"] and not locus.data["is_paralogue"]:
                visit_key = (source_side, locus_chr, locus_start, locus_end)
                if visit_key not in visited_loci:
                    # max-priority by read count via negative key.
                    heapq.heappush(
                        global_locus_queue,
                        (-locus.data["read_count"], locus_chr, locus_start, locus_end, source_side),
                    )
                    queued += 1
    return queued, repeat_read_ids, repeat_loci, relevant_loci


def annotate_locus_store(locus_list, genes_gr, transcripts_gr):
    associate_loci_with_annotations(locus_list, genes_gr, transcripts_gr)

def track_unmatched_reads(unmatched_ids, name_index, locus_list):
    mapped_elsewhere_ids = set()
    filtered_by_quality_ids = set()
    not_found_in_opposite_ids = set()

    for read_id in unmatched_ids:
        try:
            read_iter = name_index.find(read_id)
        except KeyError:
            not_found_in_opposite_ids.add(read_id)
            continue

        mapped_this_read = False
        for read in read_iter:
            if quality_filter(read):
                add_read_to_loci(locus_list, read.reference_name, read_id, read)
                mapped_elsewhere_ids.add(read_id)
                mapped_this_read = True
                break
        if not mapped_this_read:
            filtered_by_quality_ids.add(read_id)

    return mapped_elsewhere_ids, filtered_by_quality_ids, not_found_in_opposite_ids

def parse_range(range_str):
    chr_name, interval = range_str.split(":")
    chr_name = chr_name.replace("chr", "")
    start_1based, end_1based = map(int, interval.replace(",", "").split("-"))
    if end_1based < start_1based:
        raise ValueError("start_range end must be >= start")
    # Input contract is 1-based inclusive. Internal representation is 0-based end-exclusive.
    return chr_name, start_1based - 1, end_1based

def main():
    args = parse_args()
    bamA = pysam.AlignmentFile(args.bamA, "rb", require_index=True)
    bamB = pysam.AlignmentFile(args.bamB, "rb", require_index=True)
    print("loading annotations...")
    genes_gr, transcripts_gr = load_gtf_annotations(args.gtf)
    repeats_index = None
    if args.repeatmasker:
        print("loading repeat annotations")
        repeats_index = load_repeatmasker_annotations(args.repeatmasker)
    print("building indices...")
    bamA_name_index = pysam.IndexedReads(bamA)
    bamA_name_index.build()
    print("index for bamA built")
    bamB_name_index = pysam.IndexedReads(bamB)
    bamB_name_index.build()
    print("index for bamB built")
    start_chr, start, end = parse_range(args.start_range)

    global_locus_queue = [(-1, start_chr, start, end, "seed")]
    visited_loci = set()
    seen_diff_reads = {"only_in_A": set(), "only_in_B": set()}
    locus_list_by_side = {"only_in_A": {}, "only_in_B": {}}
    considered_read_ids = set()
    mapped_elsewhere_read_ids = set()
    filtered_by_quality_read_ids = set()
    not_found_in_opposite_read_ids = set()
    repeat_explained_read_ids = set()
    total_enqueued = 1
    no_diff_nodes = 0

    while global_locus_queue:
        neg_count, locus_chr, locus_start, locus_end, source_side = heapq.heappop(global_locus_queue)
        visit_key = (source_side, locus_chr, locus_start, locus_end)
        if visit_key in visited_loci:
            continue
        visited_loci.add(visit_key)
        #compute differences for the locus 
        only_in_A, only_in_B = get_mapping_diff(locus_chr, locus_start, locus_end, bamA, bamB)
        new_only_in_A = only_in_A - seen_diff_reads["only_in_A"]
        new_only_in_B = only_in_B - seen_diff_reads["only_in_B"]
        seen_diff_reads["only_in_A"].update(new_only_in_A)
        seen_diff_reads["only_in_B"].update(new_only_in_B)
        considered_read_ids.update(new_only_in_A)
        considered_read_ids.update(new_only_in_B)
 
        print(
            f"VISIT {source_side} {locus_chr}:{locus_start}-{locus_end} "
            f"priority={-neg_count} new_only_in_A={len(new_only_in_A)} new_only_in_B={len(new_only_in_B)}"
        )

        # No new differences to track: IDs were already seen elsewhere.
        if not new_only_in_A and not new_only_in_B:
            no_diff_nodes += 1
            print("No new differences found; skipping expansion.")
            continue

        prioritized_diff_sets = [("only_in_A", new_only_in_A), ("only_in_B", new_only_in_B)]
        prioritized_diff_sets.sort(key=lambda x: len(x[1]), reverse=True)

        for side_name, unmatched_ids in prioritized_diff_sets:
            if not unmatched_ids:
                continue

            if side_name == "only_in_A":
                opposite_index = bamB_name_index
                target_label = "B"
            else:
                opposite_index = bamA_name_index
                target_label = "A"

            locus_list = locus_list_by_side[side_name]
            mapped_elsewhere_ids, filtered_by_quality_ids, not_found_in_opposite_ids = track_unmatched_reads(
                unmatched_ids, opposite_index, locus_list
            )
            mapped_elsewhere_read_ids.update(mapped_elsewhere_ids)
            filtered_by_quality_read_ids.update(filtered_by_quality_ids)
            not_found_in_opposite_read_ids.update(not_found_in_opposite_ids)
            annotate_locus_store(locus_list, genes_gr, transcripts_gr)
            queued, repeat_read_ids, repeat_loci, total_loci = classify_and_queue_loci(
                locus_list,
                locus_chr,
                locus_start,
                locus_end,
                repeats_index,
                unmatched_ids,
                global_locus_queue,
                visited_loci,
                source_side=side_name,
            )
            repeat_explained_read_ids.update(repeat_read_ids)
            total_enqueued += queued

            print(
                f"  {side_name}: {len(unmatched_ids)} new reads, "
                f"mapped_elsewhere_in_{target_label}={len(mapped_elsewhere_ids)}, "
                f"repeat_loci={repeat_loci}, "
                f"filtered_by_quality={len(filtered_by_quality_ids)}, "
                f"not_found_in_{target_label}={len(not_found_in_opposite_ids)}, "
                f"unique_loci={total_loci}, queued={queued}"
            )
            '''
            for report_chr, chr_payload in locus_list.items():
                for locus_key, locus in chr_payload.get("loci", {}).items():
                    print(
                        f"  LOCUS {report_chr}:{locus_key[0]}-{locus_key[1]} "
                        f"read_count={locus.data['read_count']} "
                        f"is_repeat={locus.data['is_repeat']} "
                        f"is_paralogue={locus.data['is_paralogue']} "
                        f"genes={sorted(locus.associated_genes)} "
                        f"transcripts={sorted(locus.associated_transcripts)}"
                    )
            '''
    print(
        f"Traversal complete: visited={len(visited_loci)}, "
        f"queued_total={total_enqueued}, "
        f"no_new_difference_nodes={no_diff_nodes}"
    )

    total_considered = len(considered_read_ids)
    total_mapped_elsewhere = len(mapped_elsewhere_read_ids)
    total_filtered = len(filtered_by_quality_read_ids)
    total_not_found = len(not_found_in_opposite_read_ids)
    total_repeat_explained = len(repeat_explained_read_ids & considered_read_ids)
    explained_percent = (total_mapped_elsewhere / total_considered * 100.0) if total_considered else 0.0
    repeat_explained_percent = (total_repeat_explained / total_considered * 100.0) if total_considered else 0.0
    accounted_read_ids = (
        mapped_elsewhere_read_ids
        | filtered_by_quality_read_ids
        | not_found_in_opposite_read_ids
    )
    total_unaccounted = len(considered_read_ids - accounted_read_ids)

    print("Read statistics:")
    print(f"  different_reads_considered={total_considered}")
    print(
        f"  explained_by_mapping_elsewhere={total_mapped_elsewhere} "
        f"({explained_percent:.2f}%)"
    )
    print(
        f"  explained_through_repeats={total_repeat_explained} "
        f"({repeat_explained_percent:.2f}%)"
    )
    print(f"  filtered_out_by_quality={total_filtered}")
    print(f"  not_found_in_opposite_bam={total_not_found}")
    print(f"  unaccounted_reads={total_unaccounted}")

            
            



if __name__ == "__main__":
    main()
