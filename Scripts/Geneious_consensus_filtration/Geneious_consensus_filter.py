from collections import defaultdict
from argparse import ArgumentParser, ArgumentError
import os

desc = "" \
       "======\n" \
       "Usage: python3\n" \
       "======\n" \
       "Department of Biochemistry and Molecular Biology, Dalhousie University\n" \
       "This version is modified and implemented by Dmytro Tymoshenko (RA at the mentioned " \
       "department), April 06/2023, Github <https://github.com/edgeemer>"

parser = ArgumentParser(description=desc)
parser.add_argument('-i', '--input', type=str, help='Input directory path', required=True)
parser.add_argument('-o', '--output', type=str, help='Output directory path', required=True)
parser.add_argument('-p', '--percentage', type=int, help='Percentage of known nucleotides in sequences (default: 75)',
                    default=75)
parser.add_argument('-is', '--input_symbol', type=str, help='symbol for replacement (default: ?)',
                    default="?")
parser.add_argument('-os', '--output_symbol', type=str, help='used for replacement symbol (default: N)',
                    default="N")


def initialisation():
    # Get arguments and raise error in case of problems
    try:
        args = parser.parse_args()
        input_path: str = args.input.rstrip('/')
        output_path: str = args.output.rstrip('/')
        percentage: int = args.percentage
        input_symbol: str = args.input_symbol
        output_symbol: str = args.output_symbol
        return input_path, output_path, percentage, input_symbol, output_symbol

    except ArgumentError as e:
        print(f"Error: {e}")
        parser.print_help()
        exit(1)


def main():
    input_path, output_path, percentage, input_symbol, output_symbol = initialisation()

    for input_file in os.listdir(input_path):
        if (
                os.path.isfile(os.path.join(input_path, input_file))
                and input_file.split('.')[-1] == 'fasta'
        ):
            # Open input file and read lines into list
            with open(f'{input_path}/{input_file}', 'r') as handle:
                lines = handle.readlines()
            handle.close()

            # Create defaultdict to store sequences keyed by header
            sequences = defaultdict(list)
            str_sequences = defaultdict(str)
            for line in lines:
                if line.startswith(">"):
                    header = line.strip()
                else:
                    sequences[header].append(line.strip())
                str_sequences[header] = ''.join(sequences[header]).strip(input_symbol)

            with open(f"{output_path}/{'.'.join(input_file.split('.')[0:-1]).split('/')[-1]}_filtered.fasta", 'w+') \
                    as outfile:

                for header in str_sequences.keys():

                    # Extract length and coverage information from header
                    if header.startswith(">"):
                        header_arg_list = header.split("_")
                        if "reflen" in header_arg_list:
                            reflen = int(header_arg_list.index("reflen") + 1)
                        else:
                            reflen = -1

                        # Filter sequences based on length and coverage criteria
                        if reflen == -1 or (round(len(str_sequences[header]) / reflen, 3) >= (percentage / 100)):
                            outfile.write(f"{header}\n{str_sequences[header].replace(input_symbol, output_symbol)}\n")
                outfile.close()
    print("Successful run!")


if __name__ == '__main__':
    main()
