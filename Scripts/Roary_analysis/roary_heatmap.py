import argparse
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Set up argument parser
parser = argparse.ArgumentParser(description='Generate a heatmap for gene presence and absence.')
parser.add_argument('-i', '--input_directory', type=str, help='Input directory where gene_presence_absence.csv is located', required=True)
parser.add_argument('-o', '--output_directory', type=str, help='Output directory for the heatmap.', default='.')

# Parse arguments
args = parser.parse_args()

# Check if output directory exists, if not, create it
if not os.path.isdir(args.output_directory):
    os.makedirs(args.output_directory)

# Construct the full path to the gene_presence_absence.csv file
gene_presence_absence_file = os.path.join(args.input_directory, 'gene_presence_absence.csv')

# Load the gene_presence_absence.csv file
data = pd.read_csv(gene_presence_absence_file)

# The file contains a lot of information. For presence/absence, we focus on gene names and their presence (1) or absence (0) across genomes.
# Typically, the 'Gene' column holds gene identifiers, and the following columns hold presence (1) or absence (0) data.

# Extracting relevant columns (Adjust the range according to your specific file structure)
presence_absence = data.iloc[:, 14:].notnull().astype(int)

# Creating a heatmap
plt.figure(figsize=(20, 10))
sns.heatmap(presence_absence.transpose(), cmap="YlGnBu", yticklabels=False)  # Transpose to make the heatmap horizontal
plt.title('Gene Presence and Absence Heatmap')
plt.xlabel('Genomes')
plt.ylabel('Genes')

# Save the heatmap
heatmap_file = os.path.join(args.output_directory, 'gene_presence_absence_heatmap.png')
plt.savefig(heatmap_file)
print(f"Heatmap saved to {heatmap_file}")
plt.close()
