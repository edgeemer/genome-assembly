from collections import defaultdict
from argparse import ArgumentParser, RawTextHelpFormatter, ArgumentError, ArgumentTypeError
from pathlib import Path

desc = "\nScript is used for assembly filtration and quality analysis (originally for geneious consensus output)." \
       "\nIt also uses reflen parameter in the headline to calculate coverage and generates report of" \
       "the resulting filtration:\n" \
       "\n└INPUT_DIRECTORY (exists)               └OUTPUT_DIRECTORY (exists or will be created)   " \
       "\n   ├ {name_1}.fasta (exists)      =>      ├ {name_1}_filtered.fasta (will be created/overwritten) " \
       "\n   ├ {name_2}.fasta (exists)      =>      ├ {name_2}_filtered.fasta (will be created/overwritten) " \
       "\n   └ {name_3}.fasta (exists)              ├ {name_3}_filtered.fasta (will be created/overwritten) " \
       "\n                                          └ filtration_report.md    (will be created/overwritten)" \
       "\n\nIn each file, sequences are trimmed to the first non-? SYMBOL from both sides:" \
       "\n      ?????????ACCTT????GGCC?????? => ACCTT????GGCC" \
       "\nIf REFERENCE LENGTH inside HEADERS is provided ( >{your_text}_reflen_1500_{your_text} )" \
       "and PERCENTAGE is available (75 by default) sequences with lower COVERAGE will not pass to the final files." \
       " Filtration is based on COUNTED SYMBOLS (default: ?N) parsed as string. Example: -cs RYMN?" \
       "\nIn the OUTPUT DIRECTORY report file is generated." \
       "\n\nDepartment of Biochemistry and Molecular Biology, Dalhousie University\n" \
       "This version is developed and implemented by Dmytro Tymoshenko (RA at the mentioned department), " \
       "MAY 22/2023, Github <https://github.com/edgeemer>"

usage = "\npython <script_name>.py [-h] -i INPUT -o OUTPUT " \
        "[-p PERCENTAGE] [-cs COUNTED_SYMBOLS]\n" \
        "Options:\n" \
        "-i, --input             Input directory path\n" \
        "-o, --output            Output directory path\n" \
        "-p, --percentage        Percentage of known nucleotides in sequences (default: 75)\n" \
        "-cs, --counted_symbols  Symbols for percentage calculation (default: N?)\n"

parser = ArgumentParser(description=desc, formatter_class=RawTextHelpFormatter, usage=usage)
parser.add_argument('-i', '--input', type=str, help='Input directory path', required=True)
parser.add_argument('-o', '--output', type=str, help='Output directory path', required=True)
parser.add_argument('-p', '--percentage', type=int, help='Percentage of known nucleotides in sequences (default: 75)',
                    default=75)
parser.add_argument('-cs', '--counted_symbols', type=str, help='Symbols for subtracting from a sequence statistics'
                                                               ' (default: ?N)',
                    default="?N")


def initialisation():
    # Get arguments and raise error in case of problems
    try:
        args = parser.parse_args()

        # Required
        input_path: Path = Path(args.input).resolve()
        output_path: Path = Path(args.output).resolve()

        # Optional
        percentage: int = args.percentage
        counted_symbols: str = args.counted_symbols

        return input_path, output_path, percentage, ''.join(set(counted_symbols))  # Exclude repeated symbols

    except (ArgumentError, ArgumentTypeError) as e:
        print(f"Check help! Error: {e}")
        parser.print_help()
        exit(1)


def main():
    input_path, output_path, percentage, counted_symbols = initialisation()

    # Create dict with all possible variants for report statistics and coverage calculation.
    # IUPAC ambiguous nucleotide codes: R: A or G | Y: C or T | M: A or C | K: G or T | S: C or G | W: A or T |
    # H: A or C or T | B: C or G or T | V: A or C or G | D: A or G or T | N: A or C or G or T

    all_symbols = {
        "counted": counted_symbols,
        "two_nucleotides": "".join(i for i in "RYMKSW" if i not in counted_symbols),
        "three_nucleotides": "".join(i for i in "HBVD" if i not in counted_symbols),
        "N|gap": "".join(i for i in "N-" if i not in counted_symbols)
    }

    # Create {output_path} if it doesn't exist
    output_path.mkdir() if not output_path.exists() else None
    # Creating {filtration_report.md} now for possibility of adding data only from current assembly
    with open(f'{output_path}/filtration_report.md', 'w') as report_outfile:
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

        # Create defaultdict to store trimmed sequences (both sides to first non ?) in str type
        str_sequences = defaultdict(str)
        for header in sequences.keys():
            str_sequences[header] = ''.join(sequences[header]).strip("?")

        # Create/Overwrite the output file in the output directory
        with open(f"{output_path}/{input_file.stem}_filtered.fasta", 'w') \
                as outfile:

            # Create defaultdict to store similarity and nucleotide calculations for the report
            calculations = defaultdict(dict)

            """
            Extract reflen (reference length) from header (xxxxx_reflen_1500_xxxx.fasta) and calculate coverage:
                            Trimmed Sequence Length - Number of Input Symbols Inside
            Coverage = ------------------------------------------------------------  * 100%
                                             Reference Length
            if reference length is absent, reflen and coverage are equal to zero 
            """

            for header in str_sequences.keys():
                if header.startswith(">"):

                    # Calculate nucleotide types in order to categorize
                    sequence_calculations = {
                        "counted": 0,
                        "two_nucleotides": 0,
                        "three_nucleotides": 0,
                        "N|gap": 0
                    }

                    # Calculate symbol repetitions by groups
                    for symbol_list_name, symbol_list in all_symbols.items():
                        for symbol in symbol_list:
                            sequence_calculations[symbol_list_name] += str_sequences[header].count(symbol)

                    # Calculate coverage based on nucleotide types
                    header_arg_list = header.split("_")
                    # If reference length is presented
                    if "reflen" in header_arg_list:
                        reflen = int(header_arg_list[header_arg_list.index("reflen") + 1])
                        seq_length = len(str_sequences[header])
                        not_counted_nucleotides = seq_length - sequence_calculations["counted"]
                        coverage = int(round(not_counted_nucleotides / reflen * 100))
                    # # If reference length is not presented
                    else:
                        reflen, coverage, = 0, 0

                    # Calculate parameters of the sequence
                    counted = round(((sequence_calculations["counted"] / seq_length) * 100), 2)
                    two_nucleotides = round(((sequence_calculations["two_nucleotides"] / seq_length) * 100), 2)
                    three_nucleotides = round(((sequence_calculations["three_nucleotides"] / seq_length) * 100), 2)
                    N_and_gap = round(((sequence_calculations["N|gap"] / seq_length) * 100), 2)

                    # Filter sequences from str_sequences and write them to output file
                    # based on similarity criteria or reference length absence
                    if not reflen or (coverage >= percentage):
                        outfile.write(f"{header}\n{str_sequences[header]}\n")
                        calculations[header] = {"assembly_length": seq_length,
                                                "coverage": coverage,
                                                "counted": counted,
                                                "two_nucleotides": two_nucleotides,
                                                "three_nucleotides": three_nucleotides,
                                                "N|gap": N_and_gap
                                                }

        # Open previously created report file
        with open(f'{output_path}/filtration_report.md', 'a+') as report_outfile:

            # Templates for lines
            report_header = ('|{:^98}' + '|{:^28}' * 6 + '|\n').format(
                'Assembly name',
                f'Assembly length',
                f'Coverage to reference (>{percentage}%)',
                f'{all_symbols["counted"]} %',
                f'{all_symbols["two_nucleotides"]} %',
                f'{all_symbols["three_nucleotides"]} %',
                f'{all_symbols["N|gap"]} %')
            empty_line_header = ('|{:^96}' + '|{:^26}' * 6 + '|\n').format(
                ':' + '-' * 96 + ':', ':' + '-' * 26 + ':', ':' + '-' * 26 + ':', ':' + '-' * 26 + ':',
                ':' + '-' * 26 + ':', ':' + '-' * 26 + ':', ':' + '-' * 26 + ':'
            )
            plasmid_header_output_line = ('|{:^98}' + '|{:^28}' * 6 + '|\n').format(
                f"<b>{input_file.stem}_filtered.fasta</b>",
                f"<b>{len(calculations.keys())} plasmids/b>", '*' * 26, '*' * 26, '*' * 26, '*' * 26, '*' * 26, )

            # Create header in case of new file
            if report_outfile.tell() == 0:
                report_outfile.write(report_header)
                report_outfile.write(empty_line_header)

            # Add general sample name and number of left after filtration plasmids
            report_outfile.write(plasmid_header_output_line)
            # Add similarity result to the relevant plasmid assembly
            for plasmid_assembly_name in calculations.keys():
                report_outfile.write(('|{:^98}' + '|{:^28}' * 6 + '|\n').format(
                    plasmid_assembly_name.lstrip(">"),
                    calculations[plasmid_assembly_name]["assembly_length"],
                    calculations[plasmid_assembly_name]["coverage"],
                    calculations[plasmid_assembly_name]["counted"],
                    calculations[plasmid_assembly_name]["two_nucleotides"],
                    calculations[plasmid_assembly_name]["three_nucleotides"],
                    calculations[plasmid_assembly_name]["N|gap"])
                )

    print("Successful run!")  # For log


if __name__ == '__main__':
    main()
