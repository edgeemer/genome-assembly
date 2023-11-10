import argparse
import os
from Bio import SeqIO

# Setup argument parser
parser = argparse.ArgumentParser(description='Extract gene names from .gbk files and generate a report.')
parser.add_argument('-i', '--input', help='Input directory with .gbk files', required=True)
parser.add_argument('-o', '--output', help='Output directory for .txt files', required=False)


# Parse arguments
def initialization():
    args = parser.parse_args()
    input_path = args.input
    output_path = args.output if args.output else '.'
    return input_path, output_path


# Function to extract genes and write to a file
def extract_genes(gbk_path, txt_path):
    gene_names = []
    with open(gbk_path) as input_handle:
        for record in SeqIO.parse(input_handle, "genbank"):
            for feature in record.features:
                if feature.type in ["gene", "CDS"]:  # Check for both gene and CDS features
                    if 'gene' in feature.qualifiers:
                        gene_name = feature.qualifiers['gene'][0]
                        gene_names.append(gene_name)
                    # elif 'locus_tag' in feature.qualifiers:
                        # gene_name = feature.qualifiers['locus_tag'][0]
                        # gene_names.append(gene_name)
                        # print(f"Found locus tag: {gene_name}")
                    else:
                        print(f"Warning: No gene name or locus tag found for feature at location {feature.location}")
    with open(txt_path, 'w') as output_handle:
        for gene_name in gene_names:
            output_handle.write(gene_name + "\n")
    return len(gene_names)


def main(input_directory, output_directory):
    # Check if the output directory exists, if not, create it
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Start the report file with headers
    report_path = os.path.join(output_directory, 'report.md')

    with open(report_path, 'w') as report_file:
        report_file.write('| Input file name | Number of genes |\n')
        report_file.write('|:---------------:|:---------------:|\n')

    # Process each .gbk file in the input directory
    gbk_files = sorted([f for f in os.listdir(input_directory) if f.endswith('.gbk')])

    for file_name in gbk_files:
        gbk_file_path = os.path.join(input_directory, file_name)
        txt_file_path = os.path.join(output_directory, os.path.splitext(file_name)[0] + '.txt')
        gene_count = extract_genes(gbk_file_path, txt_file_path)
        # Append to the report file for each .gbk file
        with open(report_path, 'a') as report_file:
            report_file.write(f'| {file_name} | {gene_count} |\n')
        print(f"Processed {file_name} with {gene_count} genes.")

    print(f"Report generated at {report_path}")


if __name__ == '__main__':
    input_dir, output_dir = initialization()
    main(input_dir, output_dir)
