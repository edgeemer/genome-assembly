import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage
from pathlib import Path
import argparse

# Setup argument parser
parser = argparse.ArgumentParser(description='Parse emapper outputs and visualize COG categories')
parser.add_argument('-i', '--input', help='Input directory with .emapper.annotations files', required=True)
parser.add_argument('-o', '--output', help='Output directory for the results', default='.')
parser.add_argument('-p', '--prefix', help='Prefix for the output files', default='combined_COG')

# Parse arguments from command line
args = parser.parse_args()


# Visualization functions
def create_heatmap(df, output_dir, prefix):
    df = df[df.index.notnull()]
    plt.figure(figsize=(20, 16))
    sns.heatmap(df, annot=False, cmap='viridis', cbar_kws={'label': 'Count'})
    plt.title('Heatmap of COG Categories')
    plt.ylabel('COG Category')
    plt.xlabel('Sample')
    plt.savefig(os.path.join(output_dir, f"{prefix}_heatmap.png"), dpi=300)
    plt.close()


def create_clustered_bar_plot(df, output_dir, prefix):
    df = df[df.index.notnull()]
    df.plot(kind='bar', stacked=False, figsize=(24, 16))
    plt.title('Clustered Bar Plot of COG Categories')
    plt.ylabel('Count')
    plt.xlabel('COG Category')
    plt.savefig(os.path.join(output_dir, f"{prefix}_clustered_bar.png"), dpi=300)
    plt.close()


def create_stacked_bar_plot(df, output_dir, prefix):
    df = df[df.index.notnull()]
    df.plot(kind='bar', stacked=True, figsize=(24, 16))
    plt.title('Stacked Bar Plot of COG Categories')
    plt.ylabel('Count')
    plt.xlabel('COG Category')
    plt.savefig(os.path.join(output_dir, f"{prefix}_stacked_bar.png"), dpi=300)
    plt.close()


def create_hierarchical_clustering(df, output_dir, prefix):
    df = df[df.index.notnull()]
    linkage_matrix = linkage(df.T, 'ward')
    plt.figure(figsize=(20, 14))
    dendrogram(linkage_matrix, labels=df.columns)
    plt.title('Hierarchical Clustering Dendrogram')
    plt.xlabel('Sample')
    plt.ylabel('Distance')
    plt.savefig(os.path.join(output_dir, f"{prefix}_dendrogram.png"), dpi=300)
    plt.close()


def parse_emapper_to_cog_table(input_dir, output_dir, output_prefix):
    # Make sure the output directory exists, if not, create it
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Initialize a dictionary to collect COG categories across all samples
    cog_category_dict = {}

    # Get a list of all .emapper.annotations files in the input directory
    file_paths = list(Path(input_dir).rglob('*.emapper.annotations'))

    # Process each file to populate the dictionary
    for file_path in file_paths:
        # Extract the sample name from the file name (remove the '.emapper' and extension)
        sample_name = Path(file_path).stem.replace('.emapper', '')
        df = pd.read_csv(file_path, sep='\t', header=None, skiprows=5, usecols=[6], names=['COG'])
        df['COG'] = df['COG'].str.split('@').str[0]

        # Update the dictionary with counts for each category
        for category in df['COG']:
            if category not in cog_category_dict:
                cog_category_dict[category] = {sample_name: 1}
            else:
                cog_category_dict[category][sample_name] = cog_category_dict[category].get(sample_name, 0) + 1

    # Create a DataFrame from the dictionary
    combined_df = pd.DataFrame.from_dict(cog_category_dict, orient='index').fillna(0).astype(int)

    # Ensure the sample columns are sorted
    combined_df = combined_df.reindex(sorted(combined_df.columns), axis=1)

    # Save the combined table to a CSV file
    output_file_path = os.path.join(output_dir, output_prefix + '.csv')
    combined_df.to_csv(output_file_path, index_label='COG Category')

    print(f"Combined COG table saved to {output_file_path}")

    create_heatmap(combined_df, output_dir, output_prefix)
    create_clustered_bar_plot(combined_df, output_dir, output_prefix)
    create_stacked_bar_plot(combined_df, output_dir, output_prefix)
    create_hierarchical_clustering(combined_df, output_dir, output_prefix)


# Call the function with provided arguments
parse_emapper_to_cog_table(args.input, args.output, args.prefix)
