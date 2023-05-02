from collections import defaultdict
from argparse import ArgumentParser, RawTextHelpFormatter, ArgumentError, ArgumentTypeError
from pathlib import Path

desc = "\nScript is used for assembly filtration and symbol replacing (originally geneious consensus output)." \
       "\nIt also uses reflen parameter in the headline to calculate similarity and generates report of" \
       "the resulting filtration:\n" \
       "\n└INPUT_DIRECTORY (exists)               └OUTPUT_DIRECTORY (exists or will be created)   " \
       "\n   ├ {name_1}.fasta (exists)      =>      ├ {name_1}_filtered.fasta (will be created/overwritten) " \
       "\n   ├ {name_2}.fasta (exists)      =>      ├ {name_2}_filtered.fasta (will be created/overwritten) " \
       "\n   ├ {name_3}.fasta (exists)              ├ {name_3}_filtered.fasta (will be created/overwritten) " \
       "\n\nIn each file, INPUT SYMBOL is changed to OUTPUT SYMBOL inside sequences," \
       "as well as sequences are trimmed to the first non-INPUT SYMBOL from both sides:" \
       "\n      ?????????ACCTT????GGCC??????   => INPUT SYMBOL = ? | OUTPUT SYMBOL = N => ACCTTNNNNGGCC" \
       "\nIf REFERENCE LENGTH inside HEADERS is provided ( >{your_text}_reflen_1500_{your_text} )" \
       "and SIMILARITY is availavle (75 by default) sequences with lower SIMILARITY will not pass to the final files." \
       "\nIn the OUTPUT DIRECTORY report file is generated." \
       "\n\nDepartment of Biochemistry and Molecular Biology, Dalhousie University\n" \
       "This version is developed and implemented by Dmytro Tymoshenko (RA at the mentioned department), " \
       "MAY 02/2023, Github <https://github.com/edgeemer>"

usage = "\npython <script_name>.py [-h] -i INPUT -o OUTPUT " \
        "[-p PERCENTAGE] [-is INPUT_SYMBOL] [-os OUTPUT_SYMBOL]\n" \
        "Options:\n" \
        "-i, --input            Input directory path\n" \
        "-o, --output           Output directory path\n" \
        "-p, --percentage       Percentage of known nucleotides in sequences (default: 75)\n" \
        "-is, --input_symbol    Symbol for replacement (default: ?)\n" \
        "-os, --output_symbol   Used for replacement symbol (default: N)\n"

parser = ArgumentParser(description=desc, formatter_class=RawTextHelpFormatter, usage=usage)
parser.add_argument('-i', '--input', type=str, help='Input directory path', required=True)
parser.add_argument('-o', '--output', type=str, help='Output directory path', required=True)
parser.add_argument('-p', '--percentage', type=int, help='Percentage of known nucleotides in sequences (default: 75)',
                    default=75)
parser.add_argument('-is', '--input_symbol', type=str, help='Symbol for replacement (default: ?)',
                    default="?")
parser.add_argument('-os', '--output_symbol', type=str, help='Used for replacement symbol (default: N)',
                    default="N")


def initialisation():
    # Get arguments and raise error in case of problems
    try:
        args = parser.parse_args()

        # Required
        input_path: Path = Path(args.input).resolve()
        output_path: Path = Path(args.output).resolve()

        # Optional
        percentage: int = args.percentage
        input_symbol: str = args.input_symbol
        output_symbol: str = args.output_symbol

        return input_path, output_path, percentage, input_symbol, output_symbol

    except (ArgumentError, ArgumentTypeError) as e:
        print(f"Check help! Error: {e}")
        parser.print_help()
        exit(1)


def main():
    input_path, output_path, percentage, input_symbol, output_symbol = initialisation()

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
        for line in lines:
            header = line.strip() if line.startswith(">") else sequences[header].append(line.strip())

        # Create defaultdict to store trimmed sequences (trimmed: input_symbol, both sides to first non input_symbol)
        # in str type
        str_sequences = defaultdict(str)
        for header, sequence in sequences.items():
            str_sequences[header] = ''.join(sequences[header]).strip(input_symbol)

        # Create/Overwrite the output file in the output directory
        with open(f"{output_path}/{input_file.stem}_filtered.fasta", 'w') \
                as outfile:

            # Create defaultdict to store similarity calculations for the report
            similarity_to_reference = defaultdict(int)

            """
            Extract reflen (reference length) from header (xxxxx_reflen_1500_xxxx.fasta) and calculate coverage:
                            Trimmed Sequence Length - Number of Input Symbols Inside
            Similarity = ------------------------------------------------------------  * 100%
                                             Reference Length
            if reference length is absent, reflen and similarity are equal to zero 
            """

            for header in str_sequences.keys():
                if header.startswith(">"):
                    header_arg_list = header.split("_")
                    if "reflen" in header_arg_list:
                        reflen = int(header_arg_list[header_arg_list.index("reflen") + 1])
                        non_input_symbols = len(str_sequences[header]) - str_sequences[header].count(input_symbol)
                        similarity = int(round(non_input_symbols / reflen * 100))
                    else:
                        reflen, similarity = 0, 0

                    # Filter sequences from str_sequences and write them to output file
                    # based on similarity criteria or reference length absence
                    if not reflen or (similarity >= percentage):
                        outfile.write(f"{header}\n{str_sequences[header].replace(input_symbol, output_symbol)}\n")
                        similarity_to_reference[header] = similarity

        # Open previously created report file
        with open(f'{output_path}/filtration_report.md', 'a+') as report_outfile:

            # Templates for lines
            report_header = '|{:^98}|{:^28}|\n'.format('Assembly name', f'Similarity (>{percentage}%)')
            empty_line_header = '|:{:^96}:|:{:^26}:|\n'.format('-' * 96, '-' * 26)
            plasmid_header_output_line = '|{:^98}|{:^28}|\n'.format(f"<b>{input_file.stem}_filtered.fasta</b>",
                                                                    f"{len(similarity_to_reference.keys())} plasmids")

            # Create header in case of new file
            if report_outfile.tell() == 0:
                report_outfile.write(report_header)
                report_outfile.write(empty_line_header)

            # Add general sample name and number of left after filtration plasmids
            report_outfile.write(plasmid_header_output_line)
            # Add similarity result to the relevant plasmid assembly
            for plasmid_assembly_name, similarity in similarity_to_reference.items():
                report_outfile.write(
                    '|{:^98}|{:^28}|\n'.format(plasmid_assembly_name.lstrip(">"), similarity)
                )

    print("Successful run!")  # For log output


if __name__ == '__main__':
    main()
