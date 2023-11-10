import argparse
import os
import pandas as pd
import plotly.graph_objects as go

# Set up argument parser
parser = argparse.ArgumentParser(description='Generate a pie chart for gene distribution in core, soft-core, shell, and cloud categories using Plotly.')
parser.add_argument('-i', '--input_directory', type=str, help='Input directory where Roary output files are located', required=True)
parser.add_argument('-o', '--output_directory', type=str, help='Output directory for the pie chart.', default='.')

# Parse arguments
args = parser.parse_args()

# Check if output directory exists, if not, create it
if not os.path.isdir(args.output_directory):
    os.makedirs(args.output_directory)

# Construct the full path to the gene_presence_absence.csv file
gene_presence_absence_file = os.path.join(args.input_directory, 'gene_presence_absence.csv')

# Load the gene_presence_absence.csv file
data = pd.read_csv(gene_presence_absence_file)

# Calculate the number of samples
num_samples = len(data.columns) - 14  # Adjust the index based on your data

# Calculate the distribution of genes
total_genes = len(data.index)
core_genes = data.iloc[:, 14:].count(axis=1) == num_samples
soft_core_genes = data.iloc[:, 14:].count(axis=1) >= 0.95 * num_samples
shell_genes = data.iloc[:, 14:].count(axis=1) >= 0.15 * num_samples
cloud_genes = data.iloc[:, 14:].count(axis=1) < 0.15 * num_samples

# Convert booleans to integers
core_count = core_genes.sum()
soft_core_count = soft_core_genes.sum() - core_count  # Exclude core genes
shell_count = shell_genes.sum() - soft_core_genes.sum()  # Exclude soft-core genes
cloud_count = total_genes - shell_genes.sum()  # All genes not in shell or above

# Prepare pie chart data
sizes = [core_count, soft_core_count, shell_count, cloud_count]
labels = ['Core', 'Soft-core', 'Shell', 'Cloud']

# Plot the pie chart only if there are genes in each category
if any(size > 0 for size in sizes):
    fig = go.Figure(data=[go.Pie(
        labels=labels,
        values=sizes,
        hole=.4,
        pull=[0.02, 0.02, 0.02, 0.02],  # Pulls each slice slightly outwards
    )])

    # Customize pie chart layout
    fig.update_layout(
        title_text=f'Total Genes: {total_genes}',
        annotations=[{
            'text': f'Total<br>{total_genes}',
            'x': 0.5, 'y': 0.5,
            'font_size': 20,
            'showarrow': False
        }],
        uniformtext_minsize=15,  # Increase font size
        uniformtext_mode='hide'
    )

    # Save the pie chart as an HTML file
    pie_chart_file = os.path.join(args.output_directory, 'gene_distribution_pie_chart.html')
    fig.write_html(pie_chart_file)
    print(f"Pie chart saved to {pie_chart_file}")

else:
    print("No genes found in one or more categories. Pie chart not generated.")
