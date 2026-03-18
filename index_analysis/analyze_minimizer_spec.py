#!/usr/bin/python3

import sys
import struct
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np

def main(filepath):
    """
    Analyzes index density by reading the minimizer->transcripts mapping
    and calculating how often each minimizer appears across transcripts.
    """
    ENDIAN = ">"  # Big-Endian (same as Java)
    
    with open(filepath, "rb") as f:
        def read_int():
            b = f.read(4)
            return struct.unpack(ENDIAN + "i", b)[0] if b else None
        
        def read_short():
            b = f.read(2)
            return struct.unpack(ENDIAN + "h", b)[0] if b else None
        
        def read_string(length):
            b = f.read(length)
            return b.decode("utf-8") if b else ""
        
        try:
            # --- [ SKIP HEADER ] ---
            prog_id = read_int()
            index_version = read_int()
            k = read_int()
            minim_length = read_int()
            
            # --- [ SKIP GENE TABLE ] ---
            gene_count = read_int()
            for i in range(gene_count):
                name_len = read_short()
                gene_name = read_string(name_len)
            
            # --- [ SKIP TRANSCRIPT TABLE ] ---
            tx_count = read_int()
            for i in range(tx_count):
                name_len = read_short()
                tx_name = read_string(name_len)
                gene_idx = read_int()
            
            # --- [ READ MINIMIZER PATHS & INVERT MAPPING ] ---
            path_tx_count = read_int()
            minimizer_to_transcripts = defaultdict(list)
            
            for tx_id in range(path_tx_count):
                path_size = read_int()
                for _ in range(path_size):
                    minimizer = read_int()
                    minimizer_to_transcripts[minimizer].append(tx_id)
            
            # --- [ CALCULATE STATISTICS ] ---
            occurrences = [len(tx_list) for tx_list in minimizer_to_transcripts.values()]
            
            if not occurrences:
                print("No minimizers found in index.")
                return
            
            avg_occurrence = np.mean(occurrences)
            median_occurrence = np.median(occurrences)
            std_occurrence = np.std(occurrences)
            
            print("=" * 70)
            print("[ MINIMIZER SPECIFICITY ANALYSIS ]")
            print("=" * 70)
            print(f"Total unique minimizers:     {len(minimizer_to_transcripts):,}")
            print(f"Total transcripts:           {path_tx_count:,}")
            print(f"Total minimizer occurrences: {sum(occurrences):,}")
            print()
            print("Minimizer Occurrence Statistics:")
            print(f"  Average (mean):  {avg_occurrence:.3f}")
            print(f"  Median:          {median_occurrence:.1f}")
            print(f"  Std Dev:         {std_occurrence:.3f}")
            print(f"  Min:             {min(occurrences)}")
            print(f"  Max:             {max(occurrences)}")
            print("=" * 70)
            print()
            
            # --- [ PLOT HISTOGRAM ] ---
            plt.figure(figsize=(8, 7))
            counts, bins, patches = plt.hist(occurrences, bins=50, edgecolor='black', alpha=0.7, color='#008740')
            plt.xlabel('Number of Transcripts', fontsize=12)
            plt.ylabel('Frequency (# of minimizers)', fontsize=12)
            plt.title('Minimizer Specificity', fontsize=14, fontweight='bold')
            plt.yscale('log')  # Log-transform y-axis for better visibility of distribution
            plt.legend([f'Mean: {avg_occurrence:.2f}', f'Median: {median_occurrence:.1f}'], fontsize=11)
            plt.grid(True, alpha=0.3, linestyle='--', which='both')
            plt.tight_layout()
            plt.savefig('minimizer_density_histogram.png', dpi=150, bbox_inches='tight')
            print("Histogram saved to: minimizer_density_histogram.png (log scale)")
            
        except struct.error as e:
            print(f"\n[!] Error reading binary structure: {e}")
        except Exception as e:
            print(f"\n[!] An unexpected error occurred: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: python3 {sys.argv[0]} <binary_index_file>")
        print(f"Example: python3 {sys.argv[0]} runtime_analysis/k19_g5.idx")
        sys.exit(1)
    
    main(sys.argv[1])
