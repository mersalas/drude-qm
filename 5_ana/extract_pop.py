import re
import os
import csv

def extract_cluster_populations(log_file):
    """
    Extracts cluster populations from a GROMACS clustering log file.

    Args:
        log_file (str): Path to the GROMACS clustering log file.

    Returns:
        dict: A dictionary where keys are cluster IDs and values are population fractions.
    """
    populations = {}
    total_frames = None

    with open(log_file, 'r') as file:
        for line in file:
            # Extract total frames if not already found
            if total_frames is None:
                match_total = re.search(r"Number of structures for matrix (\d+)", line)
                if match_total:
                    total_frames = int(match_total.group(1))

            # Extract cluster information
            match_cluster = re.match(r"^\s*(\d+)\s*\|\s*(\d+)", line)
            if match_cluster:
                cluster_id = int(match_cluster.group(1))  # Cluster ID
                num_frames = int(match_cluster.group(2))  # Number of frames in this cluster
                populations[cluster_id] = num_frames

    # Normalize populations to fractions
    if total_frames is not None:
        for cluster_id in populations:
            populations[cluster_id] /= total_frames

    return populations

def process_multiple_logs(log_files, output_csv):
    """
    Processes multiple GROMACS log files and writes the cluster populations to a CSV file.

    Args:
        log_files (list): List of paths to GROMACS clustering log files.
        output_csv (str): Path to the output CSV file.
    """
    all_data = []

    for log_file in log_files:
        # Extract time interval from the filename (e.g., 0-150 from cluster1_0-150.log)
        time_interval = re.search(r"(\d+-\d+)", log_file)
        time_key = time_interval.group(1) if time_interval else "Unknown"
        populations = extract_cluster_populations(log_file)
        all_data.append((time_key, populations))

    # Write to CSV
    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Time Interval", "Cluster ID", "Population Fraction"])

        for time_key, populations in all_data:
            for cluster_id, fraction in populations.items():
                writer.writerow([time_key, cluster_id, fraction])

# Example usage
log_files = ["cluster1_0-10.log", "cluster1_10-50.log", "cluster1_50-150.log", "cluster1_150-300.log", "cluster1_300-500.log"]
output_csv = "cluster1_pop.csv"

# Process the log files and generate the CSV
process_multiple_logs(log_files, output_csv)

print(f"Cluster populations have been saved to {output_csv}")

