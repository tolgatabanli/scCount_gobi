#!/usr/bin/env bash
# extract_paralogs.sh — Build a gene-to-paralogs lookup table from ISAR model homology files.
#
# Iterates linearly through INDEX_FILE, derives the homology file path for each
# ISAR model directly from the IsarId, reads it if present (skips otherwise),
# and extracts all same-species (taxid 9823) paralog relationships.
# Both directions of each pair are emitted so either gene can be the lookup key.
#
# Usage:
#   extract_paralogs.sh INDEX_FILE BASE_DIR OUTPUT_FILE
#
#   INDEX_FILE  – isar-model-index.tsv
#   BASE_DIR    – root of the ISAR model tree
#                 (e.g. /mnt/raidbio2/extdata/projekte/isar/isar-models/sus_scrofa)
#   OUTPUT_FILE – path to write the resulting TSV
#                 (e.g. data/paralogs.tsv)
#
# Output format (tab-separated, one line per gene):
#   gene_id<TAB>paralog1,paralog2,...

set -euo pipefail

if [ "$#" -lt 3 ]; then
    echo "Usage: $0 INDEX_FILE BASE_DIR OUTPUT_FILE" >&2
    exit 1
fi

INDEX_FILE="$1"
BASE_DIR="$2"
OUTPUT_FILE="$3"

if [ ! -f "$INDEX_FILE" ]; then
    echo "Error: index file not found: $INDEX_FILE" >&2
    exit 1
fi

echo "Scanning ISAR models for paralog relationships..."

# Single awk pass over the index file.
# For each model, the homology file path is derived directly from the IsarId:
#   ISAR0000008482x9823 → BASE/0000/0084/ISAR0000008482x9823/...tsv
# getline opens each file inline — no subprocesses spawned.
# Files that don't exist produce no getline lines and are silently skipped.
awk -v base="$BASE_DIR" -F'\t' '
/^#/ { next }
{
    isar_id = $5
    gene_id = $8

    num = isar_id
    sub(/^ISAR/, "", num)
    sub(/x[0-9]+$/, "", num)
    layer1 = substr(num, 1, 4)
    layer2 = substr(num, 5, 4)

    hom_file = base "/" layer1 "/" layer2 "/" isar_id \
               "/" isar_id ".homology.ensembl-compara_115.tsv"

    first = 1
    while ((getline line < hom_file) > 0) {
        if (first) { first = 0; continue }   # skip header row
        n = split(line, f, "\t")
        type  = f[2]
        stab  = f[4]
        taxid = f[5]
        if (type ~ /paralog/ && taxid == "9823" && stab != gene_id) {
            print gene_id "\t" stab
            print stab    "\t" gene_id
        }
    }
    close(hom_file)
}
' "$INDEX_FILE" | sort -u |

# Aggregate per-gene: consecutive same-key lines (input is sorted) are folded
# into a comma-separated paralog list.
awk -F'\t' '
{
    if ($1 != prev) {
        if (prev != "") print prev "\t" paralogs
        prev     = $1
        paralogs = $2
    } else {
        paralogs = paralogs "," $2
    }
}
END {
    if (prev != "") print prev "\t" paralogs
}
' > "$OUTPUT_FILE"

n_genes=$(wc -l < "$OUTPUT_FILE")
echo "Done. Paralog table written to $OUTPUT_FILE ($n_genes genes with at least one paralog)."

