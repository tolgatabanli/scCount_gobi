#!/bin/bash

PROG="../out/artifacts/miniqut3_jar/miniqut3.jar"
FASTA="/mnt/raidbio2/extdata/praktikum/genprakt/genprakt-ws25/Block/pig-genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz"
GTF="/mnt/raidbio2/extdata/praktikum/genprakt/genprakt-ws25/Block/pig-genome/Sus_scrofa.Sscrofa11.1.115.chr.gtf.gz"

STEP="${1:-both}"

if [[ "$STEP" != "index" && "$STEP" != "count" && "$STEP" != "both" ]]; then
    echo "Usage: $0 [index|count|both]"
    echo "  index - run indexing step only"
    echo "  count - run counting step only"
    echo "  both  - run both steps (default)"
    exit 1
fi

# Define parameter ranges
KS=(19 23 27 31)
MINIM_LENGTHS=(5 8 11 14)

OUT_DIR="out_time"

# Initialize output files with headers
if [[ "$STEP" == "index" || "$STEP" == "both" ]]; then
    echo -e "k\tg\trep\ttime\tcpu_time\txmx" > runtime_index.tsv
fi
if [[ "$STEP" == "count" || "$STEP" == "both" ]]; then
    echo -e "k\tg\trep\ttime\tcpu_time\txmx" > runtime_count.tsv
fi

# Run all combinations of K and MINIM_LENGTH with 3 repetitions
for K in "${KS[@]}"; do
    for MINIM_LENGTH in "${MINIM_LENGTHS[@]}"; do
        INDEX_FILE="k${K}_g${MINIM_LENGTH}.idx"
        
        for rep in 1 2 3; do
            echo "=== K=$K, MINIM_LENGTH=$MINIM_LENGTH, Repetition $rep ==="
            
            # Remove previous index file
            rm -f "$INDEX_FILE"
            
            # Indexing with timing and retry logic
            if [[ "$STEP" == "index" || "$STEP" == "both" ]]; then
                echo "Running indexing (K=$K, g=$MINIM_LENGTH, rep $rep)..."
                INDEX_XMX=12
                while true; do
                    start=$(date +%s%N)
                    /usr/bin/time -f "%U %S" java -Xmx${INDEX_XMX}G -jar $PROG index \
                                -f $FASTA \
                                -g $GTF \
                                -k $K -minimLength $MINIM_LENGTH \
                                -idx $INDEX_FILE > /dev/null 2> temp_time.txt
                    exit_code=$?
                    end=$(date +%s%N)
                    
                    if [ $exit_code -eq 0 ]; then
                        # Success
                        elapsed=$(echo "scale=3; ($end - $start) / 1000000000" | bc)
                        read user sys < temp_time.txt
                        cpu_time=$(echo "scale=3; $user + $sys" | bc)
                        echo -e "$K\t$MINIM_LENGTH\t$rep\t$elapsed\t$cpu_time\t${INDEX_XMX}G" >> runtime_index.tsv
                        break
                    else
                        # Check if it's likely a heap space error
                        if grep -qi "OutOfMemoryError\|heap" temp_time.txt 2>/dev/null || [ $exit_code -ne 0 ]; then
                            INDEX_XMX=$((INDEX_XMX + 4))
                            echo "Indexing failed with Xmx=${INDEX_XMX}G. Retrying with Xmx=${INDEX_XMX}G..."
                            rm -f "$INDEX_FILE"
                        else
                            echo "Indexing failed with non-memory error. Giving up."
                            exit 1
                        fi
                    fi
                done
            fi
            
            # Counting with timing and retry logic
            if [[ "$STEP" == "count" || "$STEP" == "both" ]]; then
                echo "Running counting (K=$K, g=$MINIM_LENGTH, rep $rep)..."
                COUNT_XMX=40
                while true; do
                    start=$(date +%s%N)
                    /usr/bin/time -f "%U %S" java -Xmx${COUNT_XMX}G -jar $PROG count \
                                -o $OUT_DIR \
                                -idx $INDEX_FILE \
                                -r2 /mnt/raidbio2/extdata/omics/SeqReads/pig_aorta/visium/read2_Vis_36-15.fastq.gz \
                                -batchSize 1000000 \
                                -threads 24 > /dev/null 2> temp_time.txt
                    exit_code=$?
                    end=$(date +%s%N)
                    
                    if [ $exit_code -eq 0 ]; then
                        # Success
                        elapsed=$(echo "scale=3; ($end - $start) / 1000000000" | bc)
                        read user sys < temp_time.txt
                        cpu_time=$(echo "scale=3; $user + $sys" | bc)
                        echo -e "$K\t$MINIM_LENGTH\t$rep\t$elapsed\t$cpu_time\t${COUNT_XMX}G" >> runtime_count.tsv
                        break
                    else
                        # Check if it's likely a heap space error
                        if grep -qi "OutOfMemoryError\|heap" temp_time.txt 2>/dev/null || [ $exit_code -ne 0 ]; then
                            COUNT_XMX=$((COUNT_XMX + 8))
                            echo "Counting failed with Xmx=${COUNT_XMX}G. Retrying with Xmx=${COUNT_XMX}G..."
                        else
                            echo "Counting failed with non-memory error. Giving up."
                            exit 1
                        fi
                    fi
                done
            fi
        done
    done
done

# Clean up temporary files
rm -f temp_time.txt

echo ""
echo "=== Results ==="
if [[ "$STEP" == "index" || "$STEP" == "both" ]]; then
    echo "Indexing runtime:"
    cat runtime_index.tsv
fi
if [[ "$STEP" == "count" || "$STEP" == "both" ]]; then
    echo ""
    echo "Counting runtime:"
    cat runtime_count.tsv
fi