import sys

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from datetime import datetime

# Load CSV data
csv_file = "/home/maximilian/Downloads/OneDrive_2025-01-26/task 4/task4_quarterxg/density.csv"  # Replace with your CSV file path
output_directory = "./output/heatmap_" + datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
os.makedirs(output_directory)

if len(sys.argv) > 1:
    csv_file = sys.argv[1]
    print("using cli csv file")

data = pd.read_csv(csv_file)

# Ensure the CSV has the expected columns
if len(data.columns) < 51:
    raise ValueError("The CSV file must have at least 51 columns: 'time' and 50 data columns.")

# Extract time and data columns
time = data.iloc[:, 0]  # First column is time
data_columns = data.iloc[:, 1:51]  # Extract the next 50 columns

# Convert to NumPy array for plotting
data_array = data_columns.to_numpy()

number_of_heatmaps = data_columns.shape[0]

global_min = np.min(data_array)
global_max = np.max(data_array)

# Generate a heatmap for each line in the CSV
for i, row in enumerate(data_array):
    current_min = row.min(axis=0)
    current_avg = row.mean(axis=0)
    current_max = row.max(axis=0)

    plt.figure(figsize=(6, 6))  # Adjust figure size for a full heatmap
    ax = sns.heatmap(
        row.reshape(1, -1),  # Reshape row into a 2D array with one row
        cmap="summer",  # Color scheme
        cbar=True,  # Display color bar
        cbar_kws={'label': 'density in p/u³'},  # Add label to the color bar
        vmin = global_min,  # Set global minimum for color scale
        vmax = global_max,  # Set global maximum for color scale
        xticklabels=False,  # Disable x-axis tick labels
        yticklabels=False,  # Disable y-axis tick labels
        linewidths=0
    )

    # Add labels and titles
    ax.set_title(f"Density @ t = {time[i]:06.2f}", fontsize=16, fontweight="bold")
    ax.text(5, 1.05, f"min: {current_min:04.2f}, avg: {current_avg:04.2f}, max: {current_max:04.2f}")

    ax.set_xlabel("")
    ax.set_ylabel("")

    # Display the plot
    plt.tight_layout()
    #plt.show()
    plt.savefig(output_directory + f"/heatmap_{time[i]}.svg", format='svg')
    plt.close()

    percent = 100 * (i / float(number_of_heatmaps))
    filled_length = int(20 * i // number_of_heatmaps)
    bar = '█' * filled_length + '-' * (20 - filled_length)
    sys.stdout.write(f'\r|{bar}| {percent:.1f}% Complete')
    sys.stdout.flush()
    if i == number_of_heatmaps:
        print()  # Move to the next line when done.
