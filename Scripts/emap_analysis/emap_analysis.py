import os
import sys
import pandas as pd
import plotly.express as px


def load_annotations(input_directory):
    cog_entries = []  # Use a list to collect entries
    for filename in os.listdir(input_directory):
        if filename.endswith('.annotations'):
            sample_name = filename.split('.')[0]
            file_path = os.path.join(input_directory, filename)
            with open(file_path, 'r') as file:
                for line in file:
                    if line.startswith('#'):
                        continue
                    parts = line.strip().split('\t')
                    cog_category = parts[5] if len(parts) > 5 else 'Unclassified'
                    cog_entries.append({'Sample': sample_name, 'COG_Category': cog_category})
    cog_data = pd.DataFrame(cog_entries)  # Convert the list of dictionaries to a DataFrame
    return cog_data


def plot_distributions(cog_data, output_directory):
    # Bar chart
    bar_plot_data = cog_data['COG_Category'].value_counts().reset_index()
    bar_plot_data.columns = ['COG_Category', 'Count']

    # Create a bar chart with Plotly Express
    fig = px.bar(bar_plot_data, y='Count', x='COG_Category', orientation='v',
                 title='COG Category Distribution')
    fig.update_layout(xaxis={'categoryorder': 'total descending'})
    fig.update_xaxes(title='COG Category')
    fig.update_yaxes(title='Count')

    # Save the figure to a file
    fig.write_html(os.path.join(output_directory, 'cog_bar_chart.html'))


def main(input_directory, output_directory='.'):
    cog_data = load_annotations(input_directory)
    plot_distributions(cog_data, output_directory)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python script.py <input_directory> [output_directory]")
        sys.exit(1)

    input_dir = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else '.'
    main(input_dir, output_dir)
