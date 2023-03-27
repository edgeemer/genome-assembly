from collections import defaultdict
from argparse import ArgumentParser, ArgumentError

desc = "Custom script for trimming scaffolds.fasta files at the parameters of the minimum length at 200 bp and " \
       "minimum coverage at 5. In addition, summary for each assembly file is generated.\n " \
       "======\n" \
       "Usage: python3 NCBI_assembly_filter_mod.py -i /path/to/input/file/file.fasta -o /path/to/output/directory\n" \
       "======\n" \
       "Based on NCBI_assembly_filter_defaultdict.by\nDepartment of Biochemistry and Molecular Biology, Dalhousie " \
       "University\n" \
       "This version is modified and implemented by Dmytro Tymoshenko (RA at the mentioned " \
       "department), March 10/2023, Github <https://github.com/edgeemer>"

min_len: int = 200
min_cov: int = 5

parser = ArgumentParser(description=desc)
parser.add_argument('-i', '--input', type=str, help='Input file path (use * for multiple files)', required=True)
parser.add_argument('-o', '--output', type=str, help='Output directory path', required=True)


def main():
    # Get arguments and raise error in case of problems
    try:
        args = parser.parse_args()
        input_file, output_path = args.input, args.output.rstrip('/')

    except ArgumentError as e:
        print(f"Error: {e}")
        parser.print_help()
        exit(1)

    # Open input file and read lines into list
    with open(input_file, 'r') as handle:
        lines = handle.readlines()

    # Create defaultdict to store sequences keyed by header
    sequences = defaultdict(list)
    for line in lines:
        if line.startswith(">"):
            header = line.strip()
        else:
            sequences[header].append(line.strip())

    # Open output file and write filtered sequences, calculate total coverage and total length
    #                    total coverage    ((length(0) x cov(0)) + ((length(1) x cov(1)) + .... + ((length(n) x cov(n))
    # Average Coverage = --------------- = ----------------------------------------------------------------------------
    #                     total length                      length(0) + length(1) + ... + length(n)

    total_coverage, total_length = 0.0, 0.0
    with open(f"{output_path}/{'.'.join(input_file.split('.')[0:-1]).split('/')[-1]}_custom.fasta", 'w+') as outfile:
        for header in sequences.keys():
            if header.startswith(">"):
                # Extract length and coverage information from header
                _, _, _, length, _, cov = header.split("_")
                length, cov = float(length), float(cov)
                # Filter sequences based on length and coverage criteria
                if length >= min_len and cov >= min_cov:
                    outfile.write(f"{header}\n{''.join(sequences[header])}\n")
                    total_coverage = total_coverage + (length * cov)
                    total_length = total_length + length

    # Calculate average coverage into separate file
    average_coverage = (total_coverage / total_length)

    with open(f'{output_path}/average_coverage.txt', 'a+') as cov_outfile:

        # Templates for lines

        header = '#{:^98}#{:^28}#\n'.format('Assembly name', 'Average Coverage')
        empty_line_header = '#{:^98}#{:^28}#\n'.format('#' * 98, '#' * 28)
        empty_line = '|{:^98}|{:^28}|\n'.format('-' * 98, '-' * 28)
        output_line = '|{:^98}|{:^28}|\n'.format(
            f"{'.'.join(input_file.split('.')[0:-1]).split('/')[-1]}_custom.fasta", average_coverage)

        # Create header in case of new file

        if cov_outfile.tell() == 0:
            cov_outfile.write(empty_line_header)
            cov_outfile.write(header)
            cov_outfile.write(empty_line_header)

        # Add average coverage result to the relevant assembly

        cov_outfile.write(output_line)
        cov_outfile.write(empty_line)


if __name__ == '__main__':
    main()