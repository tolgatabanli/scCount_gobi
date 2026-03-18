#!/usr/bin/python3

import sys
import struct

def main(filepath, output_filepath=None):
    """
    Extracts minimizer counts for each transcript from a binary index file.
    Writes a TSV with: transcript_id, minimizer_count
    """
    if output_filepath is None:
        output_filepath = "transcript_minimizer_counts.tsv"
    
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
            
            # --- [ READ TRANSCRIPT TABLE ] ---
            tx_count = read_int()
            transcript_names = []
            for i in range(tx_count):
                name_len = read_short()
                tx_name = read_string(name_len)
                gene_idx = read_int()
                transcript_names.append(tx_name)
            
            # --- [ READ MINIMIZER PATHS & COUNT ] ---
            path_tx_count = read_int()
            transcript_counts = []
            
            for tx_id in range(path_tx_count):
                path_size = read_int()
                tx_name = transcript_names[tx_id] if tx_id < len(transcript_names) else f"Unknown_Tx_{tx_id}"
                transcript_counts.append((tx_name, path_size))
                
                # Skip over the minimizers
                for _ in range(path_size):
                    read_int()  # Read minimizer value
            
            # --- [ WRITE TO TSV ] ---
            with open(output_filepath, 'w') as out:
                out.write("transcript_id\tminimizer_count\n")
                for tx_name, count in transcript_counts:
                    out.write(f"{tx_name}\t{count}\n")
            
            print(f"Successfully wrote {len(transcript_counts)} transcripts to: {output_filepath}")
            
        except struct.error as e:
            print(f"\n[!] Error reading binary structure: {e}")
        except Exception as e:
            print(f"\n[!] An unexpected error occurred: {e}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"Usage: python3 {sys.argv[0]} <binary_index_file> [output_file.tsv]")
        print(f"Example: python3 {sys.argv[0]} ../runtime_analysis/k19_g5.idx")
        print(f"         python3 {sys.argv[0]} ../runtime_analysis/k19_g5.idx custom_output.tsv")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else "transcript_minimizer_counts.tsv"
    
    main(input_file, output_file)
