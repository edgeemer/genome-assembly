import os
import pandas as pd
import argparse


def parse_emapper_annotations(annotations_path):
    cog_categories = {}
    with open(annotations_path, 'r') as file:
        for line in file:
            if line.startswith('#') or line.strip() == '':
                continue  # Skip header and empty lines
            parts = line.strip().split('\t')
            taxonomic_group = parts[4].split('@')[-1].split('|')[1] if '@' in parts[4] else None
            if taxonomic_group:
                cog_categories[taxonomic_group] = cog_categories.get(taxonomic_group, 0) + 1
    return cog_categories


def parse_roary(roary_path):
    roary_data = pd.read_csv(roary_path)
    num_samples = len(roary_data.columns[14:])
    core_threshold = num_samples  # present in all samples
    soft_core_threshold = num_samples * 0.95  # present in 95% of samples

    core_genes = roary_data.iloc[:, 14:].dropna(thresh=num_samples, axis=0)
    soft_core_genes = roary_data.iloc[:, 14:].dropna(thresh=soft_core_threshold, axis=0).drop(core_genes.index)
    shell_genes = roary_data.iloc[:, 14:].dropna(thresh=2, axis=0).drop(core_genes.index).drop(soft_core_genes.index)
    cloud_genes = roary_data.iloc[:, 14:].dropna(thresh=1, axis=0).drop(core_genes.index).drop(
        soft_core_genes.index).drop(shell_genes.index)

    return core_genes, soft_core_genes, shell_genes, cloud_genes


def generate_output(annotations_dir, roary_file, output_file):
    core_genes, soft_core_genes, shell_genes, cloud_genes = parse_roary(roary_file)

    with open(output_file, 'w') as report:
        for filename in os.listdir(annotations_dir):
            if filename.endswith('.emapper.annotations'):
                sample_name = filename.split('.')[0]
                sample_path = os.path.join(annotations_dir, filename)
                cog_categories = parse_emapper_annotations(sample_path)

                report.write(f'{sample_name}:\n')
                for group, count in cog_categories.items():
                    report.write(f'{group}: {count}\n')

                core_count = sum(core_genes[sample_name].notnull())
                soft_core_count = sum(soft_core_genes[sample_name].notnull())
                shell_count = sum(shell_genes[sample_name].notnull())
                cloud_count = sum(cloud_genes[sample_name].notnull())

                report.write(f'       Core: {core_count}\n')
                report.write(f'       Soft-core: {soft_core_count}\n')
                report.write(f'       Shell: {shell_count}\n')
                report.write(f'       Cloud: {cloud_count}\n\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a report combining COG and Roary analysis.")
    parser.add_argument('-i', '--input', required=True, help='Directory containing .emapper.annotations files')
    parser.add_argument('-r', '--roary', required=True, help='Roary output file (gene_presence_absence.csv)')
    parser.add_argument('-o', '--output', required=True, help='Output file to generate the report')

    args = parser.parse_args()
    generate_output(args.input, args.roary, args.output)
