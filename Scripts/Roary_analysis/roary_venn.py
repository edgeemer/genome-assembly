import pandas as pd
import argparse
import venn
import matplotlib.pyplot as plt
import os


def read_samples(file_path):
    with open(file_path, 'r') as file:
        samples = [line.strip() for line in file]
    return samples


def process_gene_presence_absence(file_path, samples):
    df = pd.read_csv(file_path)
    # Extracting genes for each sample and converting them to sets
    gene_data = {sample: set(df[df[sample].notna()]['Gene']) for sample in samples}
    return gene_data


def create_venn_diagrams(gene_data, output_dir):
    if 2 <= len(gene_data) <= 6:
        venn.venn(gene_data)
        plot_file_path = os.path.join(output_dir, 'venn_diagram.png')
        plt.savefig(plot_file_path, bbox_inches='tight')
        plt.close()
    else:
        raise ValueError("Number of gene sets must be between 2 and 6.")

def main():
    parser = argparse.ArgumentParser(description="Process Roary's output.")
    parser.add_argument('-i', '--input', required=True, help='Input directory')
    parser.add_argument('-s', '--samples', required=True, help='File with sample names')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    args = parser.parse_args()

    input_dir = args.input
    samples_file = args.samples
    output_dir = args.output

    samples = read_samples(samples_file)
    gene_presence_absence_file = f"{input_dir}/gene_presence_absence.csv"

    gene_data = process_gene_presence_absence(gene_presence_absence_file, samples)
    create_venn_diagrams(gene_data, output_dir)

    df = pd.read_csv(gene_presence_absence_file)


if __name__ == "__main__":
    main()
