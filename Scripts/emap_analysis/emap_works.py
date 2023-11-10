import os
import pandas as pd
import argparse


def parse_emapper_annotations(annotations_path):
    """
    Parses the emapper.annotations file to extract the COG category counts.
    """
    cog_categories = {}
    with open(annotations_path, 'r') as file:
        for line in file:
            if line.startswith('#') or line.strip() == '':
                continue  # Skip header and empty lines
            parts = line.strip().split('\t')
            # Extract the highest level of taxonomy (e.g., Bacteria, Spirochaetes) from the eggNOG_OGs column
            taxonomic_group = parts[4].split('@')[-1].split('|')[1] if '@' in parts[4] else None
            if taxonomic_group:
                cog_categories[taxonomic_group] = cog_categories.get(taxonomic_group, 0) + 1
    return cog_categories


def generate_output(roary_file, annotations_dir, output_file):
    # Read the Roary output file
    roary_data = pd.read_csv(roary_file)
    num_samples = len(roary_data.columns) - 14  # Adjust based on the number of metadata columns in the file

    # Get all the .annotations files from the directory
    annotations_files = [os.path.join(annotations_dir, f) for f in os.listdir(annotations_dir) if
                         f.endswith('.emapper.annotations')]

    # Initialize the report data structure
    report_data = {}

    # Parse each annotations file and accumulate the counts
    for annotations_file in annotations_files:
        sample_name = os.path.splitext(os.path.basename(annotations_file))[0]
        cog_categories = parse_emapper_annotations(annotations_file)
        report_data[sample_name] = cog_categories

    # Determine core, soft-core, shell, and cloud gene counts
    core_genes = (roary_data.iloc[:, 14:] != '').all(axis=1).sum()
    soft_core_genes = ((roary_data.iloc[:, 14:] != '').sum(axis=1) >= (0.95 * num_samples)).sum() - core_genes
    shell_genes = ((roary_data.iloc[:, 14:] != '').sum(axis=1) > 1).sum() - core_genes - soft_core_genes
    cloud_genes = ((roary_data.iloc[:, 14:] != '').sum(axis=1) == 1).sum()

    # Write the report to the output file
    with open(output_file, 'w') as f:
        for sample, categories in report_data.items():
            f.write(f"{sample}:\n")
            for category, count in categories.items():
                f.write(f"{category}: {count}\n")
            f.write(f"       Core: {core_genes}\n")
            f.write(f"       Soft-core: {soft_core_genes}\n")
            f.write(f"       Shell: {shell_genes}\n")
            f.write(f"       Cloud: {cloud_genes}\n")
            f.write("\n")


# Usage
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a report combining COG and Roary analysis.")
    parser.add_argument('-r', '--roary', required=True, help='Roary output file (gene_presence_absence.csv)')
    parser.add_argument('-i', '--input', required=True, help='Directory containing .emapper.annotations files')
    parser.add_argument('-o', '--output', required=True, help='Output file to generate the report')

    args = parser.parse_args()
    generate_output(args.roary, args.input, args.output)
