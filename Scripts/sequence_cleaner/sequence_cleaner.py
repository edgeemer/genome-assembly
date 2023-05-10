from collections import defaultdict
from argparse import ArgumentParser, RawTextHelpFormatter, ArgumentError, ArgumentTypeError
from pathlib import Path
from re import sub
from typing import Any

desc = "\nScript is used for sequences filtration (removing symbols)" \
       "\n\n└INPUT_DIRECTORY (exists)               └OUTPUT_DIRECTORY (exists or will be created)   " \
       "\n   ├ {name_1}.fasta (exists)      =>      ├ {name_1}_clean.fasta (will be created/overwritten) " \
       "\n   ├ {name_2}.fasta (exists)      =>      ├ {name_2}_clean.fasta (will be created/overwritten) " \
       "\n   └ {name_3}.fasta (exists)              ├ {name_3}_clean.fasta (will be created/overwritten) " \
       "\n                                          └ removal_report.md    (will be created/overwritten)" \
       "\n\nANNACCYRGC????CTT => symbols: \"?NY\" => AACCRGCCTT" \
       "\n\nIn the OUTPUT DIRECTORY report file is generated." \
       "\n\nDepartment of Biochemistry and Molecular Biology, Dalhousie University\n" \
       "This version is developed and implemented by Dmytro Tymoshenko (RA at the mentioned department), " \
       "MAY 02/2023, Github <https://github.com/edgeemer>"

usage = "\npython <script_name>.py [-h] -i INPUT -o OUTPUT [-s SYMBOLS]\n" \
        "Options:\n" \
        "-i, --input             Input directory path\n" \
        "-o, --output            Output directory path\n" \
        "-s, --symbols           Symbols for removing (default: ?)\n"

parser = ArgumentParser(description=desc, formatter_class=RawTextHelpFormatter, usage=usage)
parser.add_argument('-i', '--input', type=str, help='Input directory path', required=True)
parser.add_argument('-o', '--output', type=str, help='Output directory path', required=True)
parser.add_argument('-s', '--symbols', type=str, help='Symbols for removing (default: ?)',
                    default="?")


def initialisation():
    # Get arguments and raise error in case of problems
    try:
        args = parser.parse_args()

        # Required
        input_path: Path = Path(args.input).resolve()
        output_path: Path = Path(args.output).resolve()

        # Optional
        symbols: str = args.symbols

        return input_path, output_path, symbols

    except (ArgumentError, ArgumentTypeError) as e:
        print(f"Check help! Error: {e}")
        parser.print_help()
        exit(1)


def main():
    input_path, output_path, symbols = initialisation()

    # Create {output_path} if it doesn't exist
    output_path.mkdir() if not output_path.exists() else None
    # Creating {filtration_report.md} now for possibility of adding data only from current assembly
    with open(f'{output_path}/removal_report.md', 'w') as report_outfile:
        report_outfile.write('')

    # Open input file and read lines into list
    for input_file in input_path.glob('*.fasta'):
        with open(input_file, 'r') as handle:
            lines = handle.readlines()

        # Create defaultdict to store sequences keyed by header
        sequences = defaultdict(list)
        header = ''
        for line in lines:
            header = line.strip() if line.startswith(">") else sequences[header].append(line.strip())

        # Create defaultdict to store initial sequences in str type
        str_sequences: defaultdict[Any, str] = defaultdict(str)
        # Create defaultdict to store cleaned sequences in str type
        cleaned_sequences: defaultdict[Any, str] = defaultdict(str)
        # Create defaultdict to store cleaned sequences' parameters in the nested dict
        sequence_left: defaultdict[Any, float] = defaultdict(float)

        for header, sequence in sequences.items():
            # Transform the original format of sequences (list) to the str
            str_sequences[header] = ''.join(sequences[header])
            # Remove all SYMBOLS from initial sequences and store them in a new defaultdict
            cleaned_sequences[header] = sub(f"[{symbols}]", "", str_sequences[header])
            # Calculate relevant statistics
            sequence_left[header] = round((len(cleaned_sequences[header]) / len(str_sequences[header])) * 100, 2)

        # Create/Overwrite the output file in the output directory
        with open(f"{output_path}/{input_file.stem}_cleaned.fasta", 'w') \
                as outfile:

            # Write cleaned sequences to the output file
            for header in cleaned_sequences.keys():
                if header.startswith(">"):
                    outfile.write(f"{header}\n{cleaned_sequences[header]}\n")

        # Open previously created report file
        with open(f'{output_path}/removal_report.md', 'a+') as report_outfile:

            # Templates for lines
            report_header = '|{:^98}|{:^28}|\n'.format('Sequence name', f'% after removing {symbols}')
            empty_line_header = '|:{:^96}:|:{:^26}:|\n'.format('-' * 96, '-' * 26)
            header_output_line = '|{:^98}|{:^28}|\n'.format(f"<b>{input_file.stem}_clean.fasta</b>",
                                                            f"{len(cleaned_sequences.keys())} sequences", ' ' * 26)

            # Create header in case of new file
            if report_outfile.tell() == 0:
                report_outfile.write(report_header)
                report_outfile.write(empty_line_header)

            # Add general sample name and number of left after filtration plasmids
            report_outfile.write(header_output_line)
            # Add similarity result to the relevant plasmid assembly
            for seq_name in cleaned_sequences.keys():
                report_outfile.write('|{:^98}|{:^28}|\n'.format(seq_name.lstrip(">"), sequence_left[seq_name]))

    print("Successful run!")  # For log output


if __name__ == '__main__':
    main()
