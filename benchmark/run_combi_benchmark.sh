#!/bin/bash

set -e  # Exit on any error

# --- Argument Parsing ---
MODE="${1:-both}"  # Default to 'both' if no argument provided

# Validate MODE
if [[ ! "$MODE" =~ ^(index|count|both)$ ]]; then
    echo "ERROR: Invalid mode '$MODE'"
    echo "Usage: $0 [index|count|both]"
    echo "  index - only run indexing step"
    echo "  count - only run counting step"
    echo "  both  - run both steps (default)"
    exit 1
fi

# --- Configuration ---
PROG="../out/artifacts/miniqut3_jar/miniqut3.jar"
FASTA="/mnt/raidbio2/extdata/praktikum/genprakt/genprakt-ws25/Block/pig-genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz"
GTF="/mnt/raidbio2/extdata/praktikum/genprakt/genprakt-ws25/Block/pig-genome/Sus_scrofa.Sscrofa11.1.115.chr.gtf.gz"

# Parameter Arrays
#KS=(19 23 27 31)
KS=(15 19 27)
MINIMS=(11 12 14)

# Directory Setup
SIM_BASE="../simulation"
INDEX_BASE="idx"
OUT_BASE="out"

# --- Validation ---
echo "Validating script environment..."

# Check if we're in the right directory
if [[ ! -d "$SIM_BASE" ]]; then
    echo "ERROR: Simulation directory not found at $SIM_BASE"
    echo "This script must be run from the benchmark/ directory"
    exit 1
fi

# Check if program JAR exists
if [[ ! -f "$PROG" ]]; then
    echo "ERROR: Program JAR not found at $PROG"
    exit 1
fi

# Check if required data files exist
if [[ ! -f "$FASTA" ]]; then
    echo "ERROR: FASTA file not found at $FASTA"
    exit 1
fi

if [[ ! -f "$GTF" ]]; then
    echo "ERROR: GTF file not found at $GTF"
    exit 1
fi

# Detect available resources
NUM_THREADS=$(nproc)
THREADS_PER_JOB=$((NUM_THREADS > 30 ? 30 : NUM_THREADS - 2))
echo "Detected $NUM_THREADS CPU cores; using $THREADS_PER_JOB threads per job"

mkdir -p "$INDEX_BASE"
mkdir -p "$OUT_BASE"

# --- Main Execution Loops ---

for K in "${KS[@]}"; do
    for MINIM_LENGTH in "${MINIMS[@]}"; do
        
        INDEX="$INDEX_BASE/sccount_k${K}_g${MINIM_LENGTH}.idx"

        # 1. Indexing Step
        # We only run this if the index file doesn't exist to save time
        if [[ "$MODE" =~ ^(index|both)$ ]]; then
            if [[ ! -f "$INDEX" ]]; then
                echo ">>> Building Index: K=$K, Minim=$MINIM_LENGTH"
                java -Xmx50g -jar "$PROG" index \
                    -f "$FASTA" \
                    -g "$GTF" \
                    -k "$K" -minimLength "$MINIM_LENGTH" \
                    -idx "$INDEX"
                
                if [[ ! -f "$INDEX" ]]; then
                    echo "ERROR: Index creation failed for K=$K, Minim=$MINIM_LENGTH"
                    exit 1
                fi
            else
                echo ">>> Index $INDEX already exists, skipping to counting."
            fi
        fi

        # 2. Counting Step
        # Loop through every directory in simulation/
        if [[ "$MODE" =~ ^(count|both)$ ]]; then
            for SIM_PATH in "$SIM_BASE"/counts_*; do
                if [[ -d "$SIM_PATH" ]]; then
                    INPUT_DIR_NAME=$(basename "$SIM_PATH")
                    READ2="$SIM_PATH/read2.fastq"
                    
                    # Validate that read file exists
                    if [[ ! -f "$READ2" ]]; then
                        echo "    -> WARNING: Skipping $INPUT_DIR_NAME (read2.fastq not found)"
                        continue
                    fi
                    
                    # Construct output directory name
                    # e.g., benchmark/out/counts_all_0_k15_g5
                    OUTPUT_DIR="$OUT_BASE/${INPUT_DIR_NAME}_k${K}_g${MINIM_LENGTH}"
                    
                    # Skip if output directory exists and is non-empty
                    if [[ -d "$OUTPUT_DIR" ]] && [[ -n "$(ls -A "$OUTPUT_DIR")" ]]; then
                        echo "    -> Skipping: $OUTPUT_DIR (already has results)"
                        continue
                    fi
                    
                    echo "    -> Counting: $INPUT_DIR_NAME (K=$K, G=$MINIM_LENGTH)"
                    
                    java -Xmx90g -jar "$PROG" count \
                        -o "$OUTPUT_DIR" \
                        -idx "$INDEX" \
                        -r2 "$READ2" \
                        -details \
                        -batchSize 500000 \
                        -threads "$THREADS_PER_JOB"
                    
                    if [[ ! -d "$OUTPUT_DIR" ]] || [[ -z "$(ls -A "$OUTPUT_DIR")" ]]; then
                        echo "    -> ERROR: Counting failed for $INPUT_DIR_NAME"
                        exit 1
                    fi
                fi
            done
        fi
    done
done

echo "Benchmark complete."
