import os
import argparse
import venn
import matplotlib.pyplot as plt


# Function to read genes from a file
def read_genes(file_path):
    with open(file_path, 'r') as file:
        genes = file.read().splitlines()
    return set(genes)


# Function to create a Venn diagram for the given sets
def create_venn_diagram(sets, output_path):
    if 2 <= len(sets) <= 6:
        venn.venn(sets)
        plt.savefig(output_path, bbox_inches='tight')
        plt.close()
    else:
        raise ValueError("Number of gene sets must be between 2 and 6.")


# Function to process files and create Venn diagrams
def process_files(input_dir, output_dir):
    gene_sets = {}

    # Read the gene lists and store in the dictionary
    for file_name in os.listdir(input_dir):
        if file_name.endswith('.txt'):
            set_name = file_name.replace('.txt', '')
            file_path = os.path.join(input_dir, file_name)
            gene_sets[set_name] = read_genes(file_path)

    if not (2 <= len(gene_sets) <= 6):
        raise ValueError("This script supports only between 2 and 6 gene lists.")

    output_path = os.path.join(output_dir, "venn_diagram.png")
    create_venn_diagram(gene_sets, output_path)


# Main function to set up argparse and handle the script logic
def main():
    parser = argparse.ArgumentParser(description="Create Venn diagrams from gene lists.")
    parser.add_argument('-i', '--input', type=str, help="Input directory with gene lists.", required=True)
    parser.add_argument('-o', '--output', type=str, help="Output directory for Venn diagrams.", default=".")

    args = parser.parse_args()
    input_dir = args.input
    output_dir = args.output

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    process_files(input_dir, output_dir)


if __name__ == "__main__":
    main()
