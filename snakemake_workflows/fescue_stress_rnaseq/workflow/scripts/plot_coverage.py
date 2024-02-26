# Logging 
sys.stderr = sys.stdout = open(snakemake.log[0], "w") # Redirect stderr to log file

# Import modules
import matplotlib.pyplot as plt
import numpy as np

# Define functions
# Parse bedgraph file
def parse_bedgraph(bedgraph_file, normalize_length=False, filter_data=False, min_len=100):
    # Initialize dictionaries to store coverage data and transcript lengths.
    coverage_data = {}
    transcript_lengths = {}

    # Open the bedgraph file for reading.
    with open(bedgraph_file, 'r') as file:
        current_transcript = None
        current_coverage = []

        # Loop through each line in the file.
        for line in file:
            # Split each line into its respective fields.
            parts = line.strip().split()
            transcript, start, end, coverage = parts[0], int(parts[1]), int(parts[2]), float(parts[3])

            # Check if we're still processing the same transcript or moving to a new one.
            if transcript != current_transcript and current_transcript is not None:
                # Update the length of the current transcript.
                transcript_lengths[current_transcript] = len(current_coverage)
                # Reset the coverage list for the new transcript.
                current_coverage = []

            # Update the coverage list with coverage values for each position.
            current_coverage.extend([coverage] * (end - start))
            # Update the current transcript name.
            current_transcript = transcript

        # Calculate the total number of transcripts processed.
        total_transcripts = len(transcript_lengths)

        # Calculate the 25th percentile length of all transcripts.
        # min_length = np.percentile(list(transcript_lengths.values()), 25)

        # Print the 25th and 75th percentile lengths
        # print(f"25th percentile length: {min_length}")
        print(f"First 10 Transcript Lengths: {list(transcript_lengths.values())[:10]}")
        
        # Reset the file pointer to the beginning to process it again.
        # This time, we filter transcripts based on length and optionally normalize their length.
    with open(bedgraph_file, 'r') as file:
        current_transcript = None
        current_coverage = []
        for line in file:
            parts = line.strip().split()
            transcript, start, end, coverage = parts[0], int(parts[1]), int(parts[2]), float(parts[3])
            if transcript != current_transcript and current_transcript is not None:
                # Check if filtering is enabled and filter based on transcript lengths
                if filter_data and len(current_coverage) < min_len:
                    current_coverage = []  # Clear the current coverage data if it's filtered out
                else:
                    if normalize_length:
                        normalized_coverage = [current_coverage[int(i * (len(current_coverage) - 1) / 99)] for i in range(100)]
                        coverage_data[current_transcript] = normalized_coverage
                    else:
                        coverage_data[current_transcript] = current_coverage
                    current_coverage = []
            current_coverage.extend([coverage] * (end - start))
            current_transcript = transcript

        if filter_data:
            # Calculate the number of transcripts that were filtered out.
            filtered_transcripts = total_transcripts - len(coverage_data)
            print(f"{filtered_transcripts} transcripts (out of {total_transcripts}) were filtered out due to being shorter than {min_len} bp length.")

    # Return the coverage data and transcript lengths.
    return coverage_data, transcript_lengths


def bedgraph_stats(coverage_data):
    for transcript, coverages in coverage_data.items():
        print(f"Transcript: {transcript}, Length: {len(coverages)}, Max coverage: {max(coverages)}, Min coverage: {min(coverages)}")

# Plot average coverage
def plot_average_coverage(coverage_data, output_file):
    all_coverages = []
    for coverages in coverage_data.values():
        while len(all_coverages) < len(coverages):
            all_coverages.append([])
        for i, cov in enumerate(coverages):
            all_coverages[i].append(cov)
    
    avg_coverages = [sum(covs)/len(covs) for covs in all_coverages]
    
    plt.plot(avg_coverages)
    plt.title("Average Coverage Profile")
    plt.xlabel("Relative Position (5' to 3')")
    plt.ylabel("Normalize Coverage")
    plt.xlim(0, 100)
    plt.savefig(output_file)
    
def separate_plots_for_length_ranges(coverage_data, transcript_lengths, output_file):
    # Get the lengths of transcripts which are more than 100 bp
    filtered_lengths = [length for transcript, length in transcript_lengths.items() if length > 100]
    
    # Define your bins based on these filtered transcript lengths
    percentiles = [25, 50, 75, 100]
    bins = [100] + [np.percentile(filtered_lengths, p) for p in percentiles]

    binned_data = {}
    for transcript, coverage in coverage_data.items():
        length = transcript_lengths[transcript]
        for i in range(len(bins) - 1):
            if bins[i] <= length < bins[i + 1]:
                if i not in binned_data:
                    binned_data[i] = []
                binned_data[i].append(coverage)

    plt.figure(figsize=(12, 8))
    for bin_idx, (bin, coverages) in enumerate(binned_data.items()):
        avg_coverage = [sum(position_coverages)/len(position_coverages) for position_coverages in zip(*coverages)]
        plt.plot(avg_coverage, label=f"{bins[bin]}-{bins[bin + 1]} bp")
        print(f"Number of transcripts in range {bins[bin]}-{bins[bin + 1]} bp: {len(coverages)}")

    plt.legend()  # Display the legend
    plt.title("Average Coverage Profile by Transcript Length")
    plt.xlabel("Relative Position (5' to 3')")
    plt.ylabel("Normalize Coverage")
    plt.xlim(0, 100)
    plt.savefig(output_file)

normalize = True if snakemake.params["normalize"] == "True" else False

bedgraph_file = snakemake.input[0]
output_file = snakemake.output[0]
coverage_data, transcript_lengths = parse_bedgraph(bedgraph_file, normalize_length=True, filter_data=True, min_len=100)

# Call bedgraph_stats to print the statistics
# bedgraph_stats(coverage_data)

plot_average_coverage(coverage_data, output_file)

# Modify output_file for separate_plots_for_length_ranges
output_file_length_ranges = output_file.replace(".png", "_length_ranges.png")

# Call separate_plots_for_length_ranges
separate_plots_for_length_ranges(coverage_data, transcript_lengths, output_file_length_ranges)