#!/bin/bash

# Runtime analysis for counting step with varying read counts
# Creates subset files from read2_Vis with 0.5M, 1M, 2M, 4M, 8M reads
# Then runs counting analysis with all k and minimizer length combinations

set -e

PROG="../out/artifacts/miniqut3_jar/miniqut3.jar"
READ_FILE="/mnt/raidbio2/extdata/omics/SeqReads/pig_aorta/visium/read2_Vis_36-15.fastq.gz"

OUT_DIR="out_time_subset"
SUBSET_DIR="read_subsets"

# Define parameter ranges
KS=(19 23 27 31)
MINIM_LENGTHS=(5 8 11 14)

# Define read counts and corresponding line counts (4 lines per read in fastq)
READ_COUNTS=(500000 1000000 2000000 4000000 8000000)
LINE_COUNTS=(2000000 4000000 8000000 16000000 32000000)
READ_LABELS=("0.5M" "1M" "2M" "4M" "8M")

# Create output directories
mkdir -p "$OUT_DIR" "$SUBSET_DIR"

# Check if results already exist
RESULTS_FILE="reads_runtime_scaling.tsv"

if [ -f "$RESULTS_FILE" ]; then
    echo "Results file $RESULTS_FILE already exists. Skipping analysis..."
    echo ""
    echo "=== Existing Results ==="
    cat "$RESULTS_FILE"
    exit 0
fi

echo "=== Creating read subsets from $READ_FILE ==="

# Create subset files
for i in "${!READ_COUNTS[@]}"; do
    READ_COUNT=${READ_COUNTS[$i]}
    LINE_COUNT=${LINE_COUNTS[$i]}
    READ_LABEL=${READ_LABELS[$i]}
    SUBSET_FILE="$SUBSET_DIR/read2_Vis_${READ_LABEL}.fastq.gz"
    
    if [ -f "$SUBSET_FILE" ]; then
        echo "Subset file $SUBSET_FILE already exists, skipping creation..."
    else
        echo "Creating subset with $READ_LABEL reads ($LINE_COUNT lines)..."
        zcat "$READ_FILE" | head -n $LINE_COUNT | gzip > "$SUBSET_FILE"
        echo "Created: $SUBSET_FILE"
    fi
done

echo ""
echo "=== Running counting runtime analysis ==="

# Initialize output file with header
RESULTS_FILE="reads_runtime_scaling.tsv"
echo -e "reads\trep\tk\tg\ttime\tcpu_time" > "$RESULTS_FILE"

# Run counting benchmarks for each k, g combination with each read file and 3 reps
for K in "${KS[@]}"; do
    for MINIM_LENGTH in "${MINIM_LENGTHS[@]}"; do
        INDEX_FILE="k${K}_g${MINIM_LENGTH}.idx"
        
        # Check if index file exists
        if [ ! -f "$INDEX_FILE" ]; then
            echo "Error: Index file $INDEX_FILE not found!"
            exit 1
        fi
        
        for i in "${!READ_COUNTS[@]}"; do
            READ_LABEL=${READ_LABELS[$i]}
            SUBSET_FILE="$SUBSET_DIR/read2_Vis_${READ_LABEL}.fastq.gz"
            
            for rep in 1 2 3; do
                echo "Running counting for $READ_LABEL reads (K=$K, g=$MINIM_LENGTH, rep $rep)..."
                
                COUNT_XMX=40
                max_attempts=5
                attempt=1
                
                while [ $attempt -le $max_attempts ]; do
                    start=$(date +%s%N)
                    /usr/bin/time -f "%U %S" java -Xmx${COUNT_XMX}G -jar $PROG count \
                                -o $OUT_DIR \
                                -idx $INDEX_FILE \
                                -r2 "$SUBSET_FILE" \
                                -batchSize 1000000 \
                                -threads 24 > /dev/null 2> temp_time.txt
                    exit_code=$?
                    end=$(date +%s%N)
                    
                    if [ $exit_code -eq 0 ]; then
                        # Success
                        elapsed=$(echo "scale=3; ($end - $start) / 1000000000" | bc)
                        read user sys < temp_time.txt
                        cpu_time=$(echo "scale=3; $user + $sys" | bc)
                        echo -e "$READ_LABEL\t$rep\t$K\t$MINIM_LENGTH\t$elapsed\t$cpu_time" >> "$RESULTS_FILE"
                        echo "  Success: time=${elapsed}s, cpu_time=${cpu_time}s"
                        break
                    else
                        # Check if it's likely a heap space error
                        if grep -qi "OutOfMemoryError\|heap" temp_time.txt 2>/dev/null; then
                            COUNT_XMX=$((COUNT_XMX + 8))
                            echo "  Memory error. Retrying with Xmx=${COUNT_XMX}G (attempt $attempt/$max_attempts)..."
                            attempt=$((attempt + 1))
                        else
                            echo "  Counting failed with non-memory error."
                            exit 1
                        fi
                    fi
                done
                
                if [ $attempt -gt $max_attempts ]; then
                    echo "Error: Max retry attempts exceeded for $READ_LABEL (K=$K, g=$MINIM_LENGTH, rep $rep)"
                    exit 1
                fi
            done
        done
    done
done

# Clean up temporary files
rm -f temp_time.txt

echo ""
echo "=== Results ==="
cat reads_runtime_scaling.tsv
