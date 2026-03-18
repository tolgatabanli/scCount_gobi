# --- Configuration ---

PROG="../out/artifacts/miniqut3_jar/miniqut3.jar"
FASTA="/mnt/raidbio2/extdata/praktikum/genprakt/genprakt-ws25/Block/pig-genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz"
GTF="/mnt/raidbio2/extdata/praktikum/genprakt/genprakt-ws25/Block/pig-genome/Sus_scrofa.Sscrofa11.1.115.chr.gtf.gz"

K=23
MINIM_LENGTH=8

INDEX_FILE="k23_g8.idx"
OUT_DIR="out"

# Binary search parameters (in MB)
RESOLUTION=256
INDEX_START=$((10 * 1024))    # 10GB in MB
COUNT_START=$((32 * 1024))    # 32GB in MB

# Convert MB to Java heap format (G or M)
format_heap() {
    local mb=$1
    if (( mb >= 1024 )); then
        echo "$((mb / 1024))G"
    else
        echo "${mb}M"
    fi
}

# Run indexing with given heap size, return 0 on success, 1 on failure
run_index() {
    local heap_mb=$1
    local heap_str=$(format_heap $heap_mb)
    
    rm -f "$INDEX_FILE"
    if java -Xmx${heap_str} -jar $PROG index \
            -f $FASTA \
            -g $GTF \
            -k 23 -minimLength 8 \
            -idx $INDEX_FILE >> index_output.log 2>&1; then
        return 0
    else
        return 1
    fi
}

# Run counting with given heap size, return 0 on success, 1 on failure
run_count() {
    local heap_mb=$1
    local heap_str=$(format_heap $heap_mb)
    
    if java -Xmx${heap_str} -jar $PROG count \
            -o $OUT_DIR \
            -idx $INDEX_FILE \
            -r2 /mnt/raidbio2/extdata/omics/SeqReads/pig_aorta/visium/read2_Vis_36-15.fastq.gz \
            -batchSize 1000000 \
            -threads 24 >> count_output.log 2>&1; then
        return 0
    else
        return 1
    fi
}

# Binary search for minimum heap size
binary_search() {
    local operation=$1
    local start=$2
    local min_working=0
    
    if [ "$operation" = "index" ]; then
        run_func="run_index"
    else
        run_func="run_count"
    fi
    
    # Start binary search from the given start value
    local low=$((start - RESOLUTION))
    local high=$start
    
    if [ $low -lt $RESOLUTION ]; then
        low=$RESOLUTION
    fi
    
    while [ $((high - low)) -gt $RESOLUTION ]; do
        local mid=$(( (low + high) / 2 ))
        # Round to nearest multiple of RESOLUTION
        mid=$(( (mid / RESOLUTION) * RESOLUTION ))
        
        if $run_func $mid; then
            high=$mid
        else
            low=$((mid + RESOLUTION))
        fi
        local heap_str=$(format_heap $high)
        echo "[${operation^^}] Range: $low - $high MB (-Xmx${heap_str})"
    done
    
    # Verify the minimum
    if $run_func $high; then
        local heap_str=$(format_heap $high)
        echo "[${operation^^}] ========================================="
        echo "[${operation^^}] MINIMUM HEAP SIZE: -Xmx${heap_str} ($high MB)"
        echo "[${operation^^}] ========================================="
        return 0
    else
        echo "[${operation^^}] ERROR: Could not verify minimum heap size"
        return 1
    fi
}

# Clear previous logs
> index_output.log
> count_output.log

# Run binary search for indexing
echo "=========================================="
echo "SEARCHING FOR MINIMAL INDEX HEAP SIZE"
echo "=========================================="
binary_search "index" $INDEX_START

# Run binary search for counting
echo ""
echo "=========================================="
echo "SEARCHING FOR MINIMAL COUNT HEAP SIZE"
echo "=========================================="
binary_search "count" $COUNT_START

