#!/bin/bash

# Configuration
NUM_TRANSCRIPTS=10000

# Arrays for the combinations
BIOTYPES=("all" "protein_coding" "lncRNA" "miRNA" "rRNA")
MUTATION_RATES=(0 1 5 10)

echo "Starting simulation grid..."

# 1. Generate ALL count files in a single R session
echo "Generating all count files using the R script..."
./generate_all_counts.R

# 2. Run simulations for each biotype and mutation rate (dynamic parallelization)
echo "Starting simulations..."
echo ""

# Set maximum number of parallel jobs (use number of CPU cores)
MAX_JOBS=$(nproc)
echo "Running up to $MAX_JOBS jobs in parallel..."
echo ""

for BIOTYPE in "${BIOTYPES[@]}"; do
    COUNT_FILE="counts_${BIOTYPE}.tsv"
    
    for MUTRATE in "${MUTATION_RATES[@]}"; do
        # Wait if we already have MAX_JOBS running
        while [ $(jobs -r | wc -l) -ge $MAX_JOBS ]; do
            sleep 0.5
        done
        
        echo "Launching: Biotype=$BIOTYPE | Mutation Rate=$MUTRATE"
        
        # Run simulation in background
        # This will create a directory named counts_${BIOTYPE}_${MUTRATE}
        ./run_simul.sh "$COUNT_FILE" "$MUTRATE" &
    done
done

# Wait for all remaining jobs to complete
echo ""
echo "All jobs submitted. Waiting for completion..."
wait

echo "All simulations finished."
