from collections import defaultdict
from argparse import ArgumentParser, RawTextHelpFormatter, ArgumentError, ArgumentTypeError
from pathlib import Path
from copy import deepcopy
from time import sleep

desc = "Script is used for working with CONTIGS output from geneious to fill in gaps and decrease level of ambiguous" \
       "data based on the enhancer value and consensus optimization.\n" \
       "Uninterrupted gaps with length of GAP are ignored for taking from the reference.\n\n" \
       "Total = Primary + Enhancer\n" \
       "Primary = working sequence, assembly to be improved\n" \
       "Enhancer = reference sequence/reads to reference sequence\n\n" \
       "If _reflen_ is provided (XXXX_reflen_500_YYYY) it will be displayed in the output report as:\n" \
       "XXX plasmids | Reference length: 500\n\n" \
       "If there is no nucleotides with a good score, N will be placed only in case if number of N values for the" \
       "same base position >=2. If the base quality is decreased (total score), Primary value will be used"

usage = "python3 <script_name>.py [-h] -i INPUT -o OUTPUT [-g GAP]\n\n" \
        "Options:\n" \
        "-i, --input             Input directory path\n" \
        "-o, --output            Output directory path\n" \
        "-g, --gap               Gap size to be ignored (default: 150)" \

parser = ArgumentParser(description=desc, formatter_class=RawTextHelpFormatter, usage=usage)
parser.add_argument('-i', '--input', type=str, help='Input directory path', required=True)
parser.add_argument('-o', '--output', type=str, help='Output directory path', required=True)
parser.add_argument('-g', '--gap', type=int, help='Gap size in Primary sequence to be ignored by the consensus',
                    default=150)


def initialization():
    try:
        args = parser.parse_args()

        # Required
        input_path: Path = Path(args.input).resolve()
        output_path: Path = Path(args.output).resolve()

        # Optional
        gap_size: int = args.gap

        return input_path, output_path, gap_size

    except (ArgumentError, ArgumentTypeError) as e:
        print(f'Check help! Error: {e}')
        exit(1)


def output_folder_initialization(input_path, output_path):
    # Create {output_path} if it doesn't exist
    if not output_path.exists():
        output_path.mkdir()
    # Creating {{name}_filtration_report.md} now for possibility of adding data only from current assembly
    with open(f'{output_path}/{input_path.name}_consensus_generator_report.md', 'w') as report_outfile:
        report_outfile.write('')


def start_and_end_positions(sequences, primary_name):
    start, end, length = 0, 0, 0

    for header, sequence in sequences.items():
        if header == primary_name:
            length = len(sequence)
            for index, char in enumerate(sequence):
                if char != '-':
                    start = index
                    break
            for index, char in enumerate(reversed(sequence)):
                if char != '-':
                    end = index
                    break

            for index, char in enumerate(reversed(sequence)):
                if char != '-':
                    if index < end:
                        end = index
                    break

    end = length - end - 1
    return start, end


def file_to_dict(file):
    with open(file, 'r') as handle:
        lines = handle.readlines()

    # Create defaultdict to store sequences keyed by header, get header name for the reference separately
    sequences = defaultdict(list)
    primary_name, header, reference_length = 'Error', '', 'N/A'
    primary_assigned = False
    enhancer_id = 1
    for line in lines:
        if not primary_assigned:
            if line.startswith('>'):
                header_arg_list = line.strip().split("_")
                if "reflen" in header_arg_list:
                    header = '_'.join(header_arg_list[0:header_arg_list.index("reflen")])
                    reference_length = header_arg_list[header_arg_list.index("reflen") + 1]  # If ref length is present
                else:
                    header = line.strip()
                primary_name, primary_assigned = deepcopy(header), True if header else None
            else:
                sequences[header].append(line.strip())

        else:
            if line.startswith('>'):
                header = line.strip() + '_' + str(enhancer_id)
                enhancer_id += 1
            else:
                sequences[header].append(line.strip())

    # Create defaultdict to store trimmed sequences in str type

    str_sequences = defaultdict(str)
    for header in sequences.keys():
        str_sequences[header] = ''.join(sequences[header])

    return primary_name, str_sequences, reference_length


def gap_definer(sequences, primary_name, gap_size_defined, seq_start_pos, seq_end_pos):
    gap_positions = []
    for index in range(seq_start_pos, seq_end_pos + 1):
        if sequences[primary_name][index] in '-' and index not in gap_positions:
            gap_start_position, gap_last_position = index, index
            end_trigger = False
            for gap_position in range(index, seq_end_pos + 1):
                if gap_position + 1 <= seq_end_pos:
                    if sequences[primary_name][gap_position] in '-':
                        continue
                    else:
                        gap_last_position = gap_position - 1
                        break
                else:
                    end_trigger = True
                    break

            if end_trigger:
                break
            else:
                if gap_last_position - gap_start_position <= gap_size_defined:
                    for gap_index in range(gap_start_position, gap_last_position + 1):
                        gap_positions.append(gap_index)

    return gap_positions


def separate_score_calculator(base_from_sequences, primary_name):
    #              A  T  G  C
    total_score = [0, 0, 0, 0]
    primary_score = [0, 0, 0, 0]
    enhancer_score = [0, 0, 0, 0]
    primary_4nt, enhancer_4nt = 0, 0

    iupac_reverse_scores = {
        #     A  T  G  C
        'A': [1, 0, 0, 0],
        'T': [0, 1, 0, 0],
        'G': [0, 0, 1, 0],
        'C': [0, 0, 0, 1],
        'R': [1, 0, 1, 0],
        'Y': [0, 1, 0, 1],
        'S': [0, 0, 1, 1],
        'W': [1, 1, 0, 0],
        'K': [0, 1, 1, 0],
        'M': [1, 0, 0, 1],
        'B': [0, 1, 1, 1],
        'D': [1, 1, 1, 0],
        'H': [1, 1, 0, 1],
        'V': [1, 0, 1, 1],
        'N': [0, 0, 0, 0],
        '-': [0, 0, 0, 0],
        '?': [0, 0, 0, 0]
    }
    primary_value = ''
    for header in base_from_sequences.keys():
        if header == primary_name:
            primary_value = base_from_sequences[primary_name]
            if primary_value == 'N':
                primary_4nt += 1
            else:
                for i, score in enumerate(iupac_reverse_scores[primary_value]):
                    primary_score[i] += score
        else:
            for i, score in enumerate(iupac_reverse_scores[base_from_sequences[header]]):
                if base_from_sequences[header] == 'N':
                    enhancer_4nt += 1
                else:
                    enhancer_score[i] += score

        for i, score in enumerate(iupac_reverse_scores[base_from_sequences[header]]):
            total_score[i] += score

    return total_score, primary_score, enhancer_score, primary_4nt, enhancer_4nt, primary_value


def score_to_value(total_score):
    iupac_bases = {
        'A': 'A', 'T': 'T', 'C': 'C', 'G': 'G', 'AG': 'R', 'CT': 'Y', 'GC': 'S', 'AT': 'W', 'GT': 'K', 'AC': 'M',
        'CGT': 'B', 'AGT': 'D', 'ACT': 'H', 'ACG': 'V', 'ACGT': 'N'
    }

    # Base score calculation
    max_val = max(total_score)

    mask = ['A', 'T', 'G', 'C']

    max_bases = [mask[index] for index, value in enumerate(total_score) if value == max_val]
    total_val = iupac_bases.get(''.join(sorted(max_bases)), 'E')

    return total_val


def nt_score_calculator(base_from_sequences, primary_name):
    iupac = {
        'AGTC': 4, 'RYMKSW': 3, 'HBVD': 2, 'N-': 0
    }

    total_score, primary_score, enhancer_score, primary_4nt, enhancer_4nt, primary_val = \
        separate_score_calculator(base_from_sequences, primary_name)

    total_val = score_to_value(total_score)

    # Quality
    total_quality, primary_quality = 0, 0
    for bases, score in iupac.items():
        if total_val in bases:
            total_quality = score
        if primary_val in bases:
            primary_quality = score

    if total_quality == primary_quality == 0 and primary_4nt + enhancer_4nt >= 2:
        print(f"\n\n{total_val} | {primary_val}\n{total_quality} | {primary_quality}\n{primary_4nt} | {enhancer_4nt}")
        sleep(5)
        return 'N'
    elif total_quality == primary_quality == 0 and primary_4nt + enhancer_4nt < 2:
        return ''
    elif total_quality > primary_quality:
        return total_val
    elif primary_quality >= total_quality:
        return primary_val
    else:
        return 'E'


def consensus_generator(input_file, gap_size_defined):
    primary_name, sequences, reference_length = file_to_dict(input_file)

    # Get position of the first and the last non-gap nucleotide of the reference
    start_position, end_position = start_and_end_positions(sequences, primary_name)

    gap_positions = gap_definer(sequences, primary_name, gap_size_defined, start_position, end_position)

    consensus = {primary_name: []}

    for nt_index in range(start_position, end_position + 1):
        if nt_index in gap_positions:
            continue

        base_from_sequences = defaultdict(str)
        for header in sequences.keys():
            base_from_sequences[header] = sequences[header][nt_index]

        t = nt_score_calculator(base_from_sequences, primary_name)
        consensus[primary_name].append(t)

    consensus[primary_name] = ''.join(consensus[primary_name])
    consensus_statistics = {
        "Reference length": reference_length,
        "Consensus length": len(consensus[primary_name]),
        "2 nt bp": 0,
        "3 nt bp": 0,
        "4 nt bp": 0,
        "Errors": 0,
    }

    for bp in consensus[primary_name]:
        if bp in "ATGC":
            continue
        elif bp in "RYMKSW":
            consensus_statistics['2 nt bp'] += 1
        elif bp in "HBVD":
            consensus_statistics['3 nt bp'] += 1
        elif bp == "N":
            consensus_statistics['4 nt bp'] += 1
        else:
            consensus_statistics['Errors'] += 1

    return consensus, consensus_statistics, primary_name


def outfile_filler(consensus, input_file, output_path):
    with open(f"{output_path}/{input_file.stem}_filtered.fasta", 'w') \
            as outfile:
        for name, sequence in consensus.items():
            outfile.write(f'{name}\n{sequence}\n')


def statistics_filler(consensus_statistics_total, input_path, output_path):
    with open(f'{output_path}/{input_path.name}_consensus_generator_report.md', 'a+') as report_outfile:

        # Templates for lines
        report_header = ('|{:^98}' + '|{:^28}' * 6 + '|\n').format(f'{input_path.name} plasmids', 'Reference length',
                                                                   'Consensus length', '2 nt bp', '3 nt bp',
                                                                   '4 nt bp', 'Errors')

        empty_line_header = ('|{:^96}' + '|{:^26}' * 6 + '|\n').format(
            ':' + '-' * 96 + ':', ':' + '-' * 26 + ':', ':' + '-' * 26 + ':', ':' + '-' * 26 + ':',
            ':' + '-' * 26 + ':', ':' + '-' * 26 + ':', ':' + '-' * 26 + ':'
        )

        # Create header in case of new file
        if report_outfile.tell() == 0:
            report_outfile.write(report_header)
            report_outfile.write(empty_line_header)

        # Add similarity result to the relevant plasmid assembly
        for plasmid_identifier in consensus_statistics_total.keys():
            report_outfile.write(('|{:^98}' + '|{:^28}' * 6 + '|\n').format(
                plasmid_identifier.lstrip(">"),
                consensus_statistics_total[plasmid_identifier]["Reference length"],
                consensus_statistics_total[plasmid_identifier]["Consensus length"],
                consensus_statistics_total[plasmid_identifier]["2 nt bp"],
                consensus_statistics_total[plasmid_identifier]["3 nt bp"],
                consensus_statistics_total[plasmid_identifier]["4 nt bp"],
                consensus_statistics_total[plasmid_identifier]["Errors"])
            )


def main():
    input_path, output_path, gap_size_defined = initialization()

    output_folder_initialization(input_path, output_path)

    consensus_statistics_summary = defaultdict(dict)

    for input_file in input_path.glob('*.fasta'):
        consensus, consensus_statistics_temp, primary_name = consensus_generator(input_file, gap_size_defined)

        consensus_statistics_summary[primary_name] = consensus_statistics_temp

        outfile_filler(consensus, input_file, output_path)

    statistics_filler(consensus_statistics_summary, input_path, output_path)

    print("Successful run!")  # For logging purposes


if __name__ == '__main__':
    main()
