#!/usr/bin/python3

import sys
import struct

def short_to_nucleotide(value):
    """Converts a 16-bit short into a 7-character nucleotide sequence."""
    sequence = []
    # Loop from 12 down to 0, stepping by -2 (ignores the first 2 bits)
    for i in range(14, -1, -2):
        bits = (value >> i) & 0b11

        if bits == 0b00:
            sequence.append('A')
        elif bits == 0b01:
            sequence.append('C')
        elif bits == 0b11:
            sequence.append('G')
        elif bits == 0b10:
            sequence.append('T')

    return "".join(sequence)

def main(filepath):
    # ">" specifies Big-Endian (standard for Java binary output).
    # If the file was generated in C++ on an x86 machine, change this to "<" (Little-Endian).
    ENDIAN = ">"

    with open(filepath, "rb") as f:
        # Helper functions to read specific byte sizes
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
            # --- [ HEADER ] ---
            prog_id = read_int()
            if prog_id is None:
                print("Error: Empty file.")
                return

            index_version = read_int()
            k = read_int()

            print("[ HEADER ]")
            print(prog_id)
            print(index_version)
            print(k)
            print()

            # --- [ GENE TABLE ] ---
            gene_count = read_int()
            print("[ GENE TABLE ]")
            for i in range(gene_count):
                name_len = read_short()
                gene_name = read_string(name_len)
                print(f"{i}\t{gene_name}")
            print()

            # --- [ TRANSCRIPT TABLE ] ---
            tx_count = read_int()
            transcripts = [] # Storing names to use later in PATH DATA
            print("[ TRANSCRIPT TABLE ]")
            for i in range(tx_count):
                name_len = read_short()
                tx_name = read_string(name_len)
                gene_idx = read_int()

                transcripts.append(tx_name)
                print(f"{i}\t{tx_name}\t{gene_idx}")
            print()

            # --- [ PATH DATA (MINIMIZERS) ] ---
            path_tx_count = read_int()
            print("[ PATH DATA (MINIMIZERS) ]")
            for i in range(path_tx_count):
                path_size = read_int()
                minimizers = []

                for _ in range(path_size):
                    m_val = read_short()
                    minimizers.append(short_to_nucleotide(m_val))

                # Retrieve the transcript name mapped to this index
                tx_name = transcripts[i] if i < len(transcripts) else f"Unknown_Tx_{i}"

                # Formatting the list of minimizers as a comma-separated string
                min_str = ", ".join(minimizers)
                print(f"{tx_name}\t{i}\t{min_str}")

        except struct.error as e:
            print(f"\n[!] Error reading binary structure: {e}")
            print("[!] Make sure the ENDIAN setting matches how the file was encoded.")
        except Exception as e:
            print(f"\n[!] An unexpected error occurred: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: python3 {sys.argv[0]} <binary_file>")
        sys.exit(1)

    main(sys.argv[1])
