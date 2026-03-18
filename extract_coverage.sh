#!/usr/bin/env bash
set -euo pipefail

BED_DIR="data/data/bedfiles_utr"
SPACERANGER_BAM_FILE="/mnt/raidbio2/extdata/praktikum/genprakt/genprakt-ws25/Block/pig-data-visium/Vis_66-10.spaceranger.bam"
HUMAN_BAM_FILE="/mnt/biocluster/praktikum/genprakt/gruppe_a/reference/5k_Human_Donor1_PBMC_3p_gem-x_GEX_fastqs/pbmc_Aligned.sortedByCoord.out.bam"
BAM_FILE="/home/t/tan/scCount_gobi/data/data/star_run/Aligned.sortedByCoord.out.bam"
OUT_DIR="data/data/coverage_spaceranger_utr"
BEDTOOLS_BIN="/mnt/raidbio2/extsoft/software/bedtools/bedtools-2.31.1/bedtools2/bin/bedtools"
MAX_JOBS=20

if ! [[ "$MAX_JOBS" =~ ^[0-9]+$ ]] || (( MAX_JOBS < 1 )); then
	echo "Error: MAX_JOBS must be a positive integer, got: $MAX_JOBS" >&2
	exit 1
fi

if [[ ! -x "$BEDTOOLS_BIN" ]]; then
	echo "Error: bedtools binary is not executable: $BEDTOOLS_BIN" >&2
	exit 1
fi

if [[ ! -f "$BAM_FILE" ]]; then
	echo "Error: BAM file not found: $BAM_FILE" >&2
	exit 1
fi

mkdir -p "$OUT_DIR"

pids=()
declare -A pid_to_out

#generating bed files 
#/usr/bin/env /usr/lib64/jvm/java-25-openjdk-25/bin/java @/tmp/cp_9bbi9y65kjtxxxmu1r9o5b7hy.argfile org.gobiws26.utils.LastBpFetcher 

for k in $(seq 100 100 2000); do
	bed_file="$BED_DIR/last_${k}_bp.bed"
	out_file="$OUT_DIR/last_${k}_bp_merged.tsv"


	if [[ ! -f "$bed_file" ]]; then
		echo "Warning: missing BED file, skipping: $bed_file" >&2
		continue
	fi

	#sorting bedfile and merging overlapping intervals frist 

	(
		"$BEDTOOLS_BIN" sort -i "$bed_file" | "$BEDTOOLS_BIN" merge -i - > "${bed_file}.merged" 2>&1 && \
		"$BEDTOOLS_BIN" multicov -bed "${bed_file}.merged" -bams "$SPACERANGER_BAM_FILE" -f 0.70 -q 30 > "$out_file" 2>&1
	) &
	pid=$!
	pids+=("$pid")
	pid_to_out["$pid"]="$out_file"

	echo "Started k=$k -> $out_file"

	while [[ $(jobs -rp | wc -l) -ge $MAX_JOBS ]]; do
		sleep 0.2
	done
done

if [[ ${#pids[@]} -eq 0 ]]; then
	echo "Error: no coverage files were generated." >&2
	exit 1
fi

failed=0
for pid in "${pids[@]}"; do
	if ! wait "$pid"; then
		failed=1
		failed_out="${pid_to_out[$pid]}"
		echo "Debug: bedtools job failed for output: $failed_out" >&2
		if [[ -f "$failed_out" ]]; then
			echo "Debug: last 20 lines of $failed_out" >&2
			tail -n 20 "$failed_out" >&2
		fi
	fi
done

if [[ $failed -ne 0 ]]; then
	echo "Error: one or more bedtools jobs failed." >&2
	exit 1
fi

echo "Wrote per-k TSV files to $OUT_DIR"
