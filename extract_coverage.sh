#!/usr/bin/env bash
set -euo pipefail

BED_DIR="data/bedfiles"
BAM_FILE="data/Aligned.sortedByCoord.out.bam"
OUT_DIR="data/coverage_tsv"
BEDTOOLS_BIN="/mnt/raidbio2/extsoft/software/bedtools/bedtools-2.31.1/bedtools2/bin/bedtools"

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

for k in $(seq 100 100 2000); do
	bed_file="$BED_DIR/last_${k}_bp.bed"
	out_file="$OUT_DIR/last_${k}_bp.tsv"

	if [[ ! -f "$bed_file" ]]; then
		echo "Warning: missing BED file, skipping: $bed_file" >&2
		continue
	fi

	(
		{
			printf 'transcript_id\tk_count\n'
			"$BEDTOOLS_BIN" multicov -bed "$bed_file" -bams "$BAM_FILE" 2>&1 | awk '{print $4 "\t" $NF}'
		} > "$out_file" 2>&1
	) &
	pids+=("$!")

	echo "Started k=$k -> $out_file"
done

if [[ ${#pids[@]} -eq 0 ]]; then
	echo "Error: no coverage files were generated." >&2
	exit 1
fi

failed=0
for pid in "${pids[@]}"; do
	if ! wait "$pid"; then
		failed=1
	fi
done

if [[ $failed -ne 0 ]]; then
	echo "Error: one or more bedtools jobs failed." >&2
	exit 1
fi

echo "Wrote per-k TSV files to $OUT_DIR"
