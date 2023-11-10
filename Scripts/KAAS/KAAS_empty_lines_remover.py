import os


# Function to remove empty lines from a file
def remove_empty_lines(filename):
    if not os.path.isfile(filename):
        print(f"{filename} does not exist.")
        return

    with open(filename, 'r') as file:
        lines = file.readlines()

    # Filter out the empty lines
    lines = [line for line in lines if line.strip()]

    with open(filename, 'w') as file:
        file.writelines(lines)


# Specify the directory you want to process
# Use '.' for the current directory
directory = '.'

# Iterate over files in the specified directory
for filename in os.listdir(directory):
    filepath = os.path.join(directory, filename)
    if os.path.isfile(filepath):
        remove_empty_lines(filepath)
        print(f"Processed {filename}")
