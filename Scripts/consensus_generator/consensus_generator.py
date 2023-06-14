from collections import defaultdict
from argparse import ArgumentParser, RawTextHelpFormatter, ArgumentError, ArgumentTypeError
from pathlib import Path
from copy import deepcopy
from time import sleep

desc = "Script is used for working with CONTIGS output from geneious to fill in gaps and decrease level of ambiguous" \
       "data based on the reference value and consensus optimization.\n" \
       "Uninterrupted gaps with length of GAP are ignored for taking from the reference.\n\n" \
       "Total = reference + contigs\n" \
       "Contigs = contigs\n" \
       "Reference = reference\n\n" \
       "If there is no nucleotides with a good score, N will be placed only in case if number of N values for the" \
       "same base position >=2. If the base quality is decreased (total score), contigs value will be used"

usage = "python3 <script_name>.py [-h] -i INPUT -o OUTPUT [-g GAP]\n\n" \
        "Options:\n" \
        "-i, --input             Input directory path\n" \
        "-o, --output            Output directory path\n" \
        "-g, --gap               Gap size to be ignored (default: 150)" \

parser = ArgumentParser(description=desc, formatter_class=RawTextHelpFormatter, usage=usage)
parser.add_argument('-i', '--input', type=str, help='Input directory path', required=True)
parser.add_argument('-o', '--output', type=str, help='Output directory path', required=True)
parser.add_argument('-g', '--gap', type=int, help='Gap size in contigs to be ignored by the consensus', default=150)


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


def start_and_end_positions(sequence):
    start: int = 0
    end: int = len(sequence) - 1
    for index, char in enumerate(sequence):
        if char != '-':
            start = index
            break

    for index, char in enumerate(reversed(sequence)):
        if char != '-':
            end = len(sequence) - index - 1
            break
    return start, end


def file_to_dict(file):
    with open(file, 'r') as handle:
        lines = handle.readlines()

    # Create defaultdict to store sequences keyed by header, get header name for the reference separately
    sequences = defaultdict(list)
    ref_name, header = '', ''
    ref_assigned = False
    for line in lines:
        if not ref_assigned:
            header = line.strip() if line.startswith('>') else sequences[header].append(line.strip())
            ref_name, ref_assigned = deepcopy(header), True if header else None
        if ref_assigned:
            header = line.strip() if line.startswith('>') else sequences[header].append(line.strip())

    # Create defaultdict to store trimmed sequences in str type

    str_sequences = defaultdict(str)
    for header in sequences.keys():
        str_sequences[header] = ''.join(sequences[header])

    return ref_name, str_sequences


def gap_definer(sequences, ref_name, gap_size_defined, seq_start_pos, seq_end_pos):
    gap_positions = []
    for index in range(seq_start_pos, seq_end_pos + 1):
        if sequences[ref_name][index] in '-' and index not in gap_positions:
            gap_start_position = index
            gap_end_position = index + 1
            for gap_position in range(index, len(sequences[ref_name])):
                if sequences[ref_name][gap_position] in '-':
                    continue
                else:
                    gap_end_position = gap_position
                    break

            if gap_end_position - gap_start_position <= gap_size_defined:
                for gap_index in range(gap_start_position, gap_end_position + 1):
                    gap_positions.append(gap_index)
        else:
            continue

    return gap_positions


def separate_score_calculator(base_from_sequences, ref_name):
    #              A  T  G  C
    total_score = [0, 0, 0, 0]
    ref_score = [0, 0, 0, 0]
    contig_score = [0, 0, 0, 0]
    ref_4nt, contig_4nt = 0, 0

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
    for header in base_from_sequences.keys():
        if header == ref_name:
            ref_value = base_from_sequences[ref_name]
            if ref_value == 'N':
                ref_4nt += 1
            else:
                for i, score in enumerate(iupac_reverse_scores[ref_value]):
                    ref_score[i] += score
        else:
            for i, score in enumerate(iupac_reverse_scores[base_from_sequences[header]]):
                if base_from_sequences[header] == 'N':
                    contig_4nt += 1
                else:
                    contig_score[i] += score

        for i, score in enumerate(iupac_reverse_scores[base_from_sequences[header]]):
            total_score[i] += score

    return total_score, ref_score, contig_score, ref_4nt, contig_4nt


def score_to_value(total_score, contig_score):
    iupac_bases = {
        'A': 'A', 'T': 'T', 'C': 'C', 'G':'G', 'AG': 'R', 'CT': 'Y', 'GC': 'S', 'AT': 'W', 'GT': 'K', 'AC': 'M',
        'CGT': 'B', 'AGT': 'D', 'ACT': 'H', 'ACG': 'V', 'ACGT': 'N'
    }

    # Base score calculation
    max_contig_val = max(contig_score)
    max_val = max(total_score)

    mask = ['A', 'T', 'G', 'C']
    contig_bases = [mask[index] for index, value in enumerate(contig_score) if value == max_contig_val]
    contig_val = iupac_bases.get(''.join(sorted(contig_bases)), 'E')

    max_bases = [mask[index] for index, value in enumerate(total_score) if value == max_val]
    total_val = iupac_bases.get(''.join(sorted(max_bases)), 'E')

    return total_val, contig_val


def nt_score_calculator(base_from_sequences, ref_name):
    iupac = {
        'AGTC': 4, 'RYMKSW': 3, 'HBVD': 2, 'N': 1, '-': 0
    }

    total_score, ref_score, contig_score, ref_4nt, contig_4nt = separate_score_calculator(base_from_sequences, ref_name)

    total_val, contig_val = score_to_value(total_score, contig_score)

    # Quality
    contig_quality, total_quality = 0, 0
    for bases, score in iupac.items():
        if contig_val in bases:
            contig_quality = score
        if total_val in bases:
            total_quality = score

    if (total_quality == contig_quality == 0 and
            ref_4nt + contig_4nt >= 2):
        return 'N'
    elif (total_quality == contig_quality == 0 and
          ref_4nt + contig_4nt < 2):
        return ''
    elif total_quality > contig_quality:
        return total_val
    elif contig_quality >= total_quality:
        return contig_val
    else:
        return 'E'


def consensus_generator(input_file, gap_size_defined):
    ref_name, sequences = file_to_dict(input_file)

    # Get position of the first and the last non-gap nucleotide of the reference
    start_position, end_position = start_and_end_positions(sequences[ref_name])

    gap_positions = gap_definer(sequences, ref_name, gap_size_defined, start_position, end_position)

    consensus = {ref_name: []}

    for nt_index in range(start_position, end_position):
        if nt_index in gap_positions:
            continue

        base_from_sequences = defaultdict(str)
        for header in sequences.keys():
            base_from_sequences[header] = sequences[header][nt_index]

        consensus[ref_name].append(nt_score_calculator(base_from_sequences, ref_name))

    consensus[ref_name] = ''.join(consensus[ref_name])
    consensus_statistics = {
        "Ref length (raw)": end_position - start_position,
        "Consensus length": len(consensus[ref_name]),
        "2 nt bp": 0,
        "3 nt bp": 0,
        "4 nt bp": 0,
        "Errors": 0,
    }

    for bp in consensus[ref_name]:
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

    return consensus, consensus_statistics, ref_name


def outfile_filler(consensus, input_file, output_path):
    with open(f"{output_path}/{input_file.stem}_filtered.fasta", 'w') \
            as outfile:
        for name, sequence in consensus.items():
            outfile.write(f'{name}\n{sequence}\n')


def statistics_filler(consensus_statistics_total, input_path, output_path):
    with open(f'{output_path}/{input_path.name}_consensus_generator_report.md', 'a+') as report_outfile:

        # Templates for lines
        report_header = ('|{:^98}' + '|{:^28}' * 6 + '|\n').format(f'{input_path.name} plasmid', 'Ref length (raw)',
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
                consensus_statistics_total[plasmid_identifier]["Ref length (raw)"],
                consensus_statistics_total[plasmid_identifier]["Consensus length"],
                consensus_statistics_total[plasmid_identifier]["2 nt bp"],
                consensus_statistics_total[plasmid_identifier]["3 nt bp"],
                consensus_statistics_total[plasmid_identifier]["4 nt bp"],
                consensus_statistics_total[plasmid_identifier]["Errors"])
            )


def main():
    input_path, output_path, gap_size_defined = initialization()

    output_folder_initialization(input_path, output_path)

    consensus_statistics_total = defaultdict(dict)

    for input_file in input_path.glob('*.fasta'):
        consensus, consensus_statistics_temp, ref_name = consensus_generator(input_file, gap_size_defined)

        consensus_statistics_total[ref_name] = consensus_statistics_temp

        outfile_filler(consensus, input_file, output_path)

    statistics_filler(consensus_statistics_total, input_path, output_path)

    print("Successful run!")  # For logging purposes


if __name__ == '__main__':
    main()
