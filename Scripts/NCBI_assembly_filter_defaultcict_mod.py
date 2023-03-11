"""
Custom script for trimming scaffolds.fasta files at the parameters
of the minimum length at 200 bp and minimum coverage at 5

Department of Biochemistry and Molecular Biology, Dalhousie University

Based on NCBI_assembly_filter_defaultdict.by

This version is modified and implemented by Dmytro Tymoshenko (RA at mentioned department), March 10/2023,
Github <https://github.com/edgeemer>
"""

from collections import defaultdict
import sys

min_len: int = 200
min_cov: int = 5


def main():
    # Define function to parse command-line arguments and provide usage information
    def usage():
        print("Usage: python3 NCBI_assembly_filter_defaultdict.py <input.txt> <output.txt>")
        sys.exit()

    # Check that exactly two command-line arguments are provided, otherwise print usage and exit
    if len(sys.argv) != 3:
        usage()

    # Open input file and read lines into list
    with open(sys.argv[1], 'r') as handle:
        lines = handle.readlines()

    # Create defaultdict to store sequences keyed by header
    sequences = defaultdict(list)
    for line in lines:
        if line.startswith(">"):
            header = line.strip()
        else:
            sequences[header].append(line.strip())

    # Open output file and write filtered sequences
    with open(sys.argv[2], 'w') as outfile:
        for header in sequences.keys():
            if header.startswith(">"):
                # Extract length and coverage information from header
                _, _, _, length, _, cov = header.split("_")
                # Filter sequences based on length and coverage criteria
                if float(length) >= min_len and float(cov) >= min_cov:
                    outfile.write(f"{header}\n{''.join(sequences[header])}\n")


if __name__ == '__main__':
    main()
