import matplotlib.pyplot as plt
import numpy as np
import subprocess
import re
import csv
import os
import seaborn as sns

# Configure settings
CONFIG = {
    "min_sigs": 2000,
    "max_sigs": 10000,
    "num_points": 9,
    "num_runs": 10,
    "output_data_file": "performance_data.csv",
    "output_graph_file": "performance_graph.png",
}

sns.set_theme(style="darkgrid")
sns.set_context("notebook", font_scale=1.1)

# Generate linearly spaced signature counts
signature_counts = np.linspace(
    CONFIG["min_sigs"], CONFIG["max_sigs"], CONFIG["num_points"], dtype=int
)

print(f"Testing signature counts: {signature_counts}")

# Initialize data storage as lists instead of dicts
all_runs_data = []

# Patterns to extract data from output
proving_pattern = re.compile(r'Dachshund Pack proving time: (\d+\.\d+)s')
verification_pattern = re.compile(r'Dachshund Pack verification time: (\d+\.\d+)s')
size_pattern = re.compile(r'Dachshund Pack size: (\d+\.\d+) KB')

def run_tests():
    print(f"Running performance tests ({CONFIG['num_runs']} runs per point)...")
    for run in range(CONFIG["num_runs"]):
        print(f"\nRun {run + 1}/{CONFIG['num_runs']}")
        for count in signature_counts:
            print(f"Testing with {count} signatures...", end=' ', flush=True)

            result = subprocess.run(['python', '../python/agg_sig.py', str(count)],
                                   capture_output=True, text=True)
            output = result.stdout

            proving_match = proving_pattern.search(output)
            verification_match = verification_pattern.search(output)
            size_match = size_pattern.search(output)

            if all([proving_match, verification_match, size_match]):
                prove_time = float(proving_match.group(1))
                verify_time = float(verification_match.group(1))
                proof_size = float(size_match.group(1))
                all_runs_data.append({
                    'signatures': count,
                    'run': run + 1,
                    'proving_time': prove_time,
                    'verification_time': verify_time,
                    'proof_size_KB': proof_size
                })
                print(f"Done (Prove: {prove_time:.2f}s, Verify: {verify_time:.2f}s, Size: {proof_size:.2f}KB)")
            else:
                print("Failed to extract data")
                all_runs_data.append({
                    'signatures': count,
                    'run': run + 1,
                    'proving_time': np.nan,
                    'verification_time': np.nan,
                    'proof_size_KB': np.nan
                })

def save_data():
    with open(CONFIG["output_data_file"], 'w', newline='') as csvfile:
        fieldnames = ['signatures', 'run', 'proving_time', 'verification_time', 'proof_size_KB']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_runs_data)
    print(f"Performance data saved to {CONFIG['output_data_file']}")

def load_data():
    if os.path.exists(CONFIG["output_data_file"]):
        loaded_data = []
        with open(CONFIG["output_data_file"], 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                loaded_data.append({
                    'signatures': int(row['signatures']),
                    'run': int(row['run']),
                    'proving_time': float(row['proving_time']),
                    'verification_time': float(row['verification_time']),
                    'proof_size_KB': float(row['proof_size_KB'])
                })

        # Verify we have the expected number of data points
        if len(loaded_data) == CONFIG["num_runs"] * CONFIG["num_points"]:
            return loaded_data
    return None

def print_markdown_table():
    # Calculate mean values for each signature count
    signatures = np.array([d['signatures'] for d in all_runs_data])
    unique_counts = np.unique(signatures)

    # Prepare table data
    table_data = []
    for count in unique_counts:
        mask = signatures == count
        prove_mean = np.nanmean([d['proving_time'] for d in all_runs_data if d['signatures'] == count])
        verify_mean = np.nanmean([d['verification_time'] for d in all_runs_data if d['signatures'] == count])
        size_mean = np.nanmean([d['proof_size_KB'] for d in all_runs_data if d['signatures'] == count])

        table_data.append({
            'Nº signatures': count,
            'Proving time': f"{prove_mean:.4f} ± {np.nanstd([d['proving_time'] for d in all_runs_data if d['signatures'] == count]):.4f}",
            'Verification time': f"{verify_mean:.4f} ± {np.nanstd([d['verification_time'] for d in all_runs_data if d['signatures'] == count]):.4f}",
            'Proof size': f"{size_mean:.2f} KB"
        })

    # Generate markdown table
    headers = table_data[0].keys()
    md_table = "| " + " | ".join(headers) + " |\n"
    md_table += "|" + "|".join(["---"] * len(headers)) + "|\n"

    for row in table_data:
        md_table += "| " + " | ".join(str(row[h]) for h in headers) + " |\n"

    print("\nMarkdown Table:\n")
    print(md_table)

    # Optional: Save to file
    with open("performance_table.md", "w") as f:
        f.write(md_table)
    print("\nTable also saved to 'performance_table.md'")

def create_graph():
    # Convert to numpy arrays for easier processing
    signatures = np.array([d['signatures'] for d in all_runs_data])
    prove_times = np.array([d['proving_time'] for d in all_runs_data])
    verify_times = np.array([d['verification_time'] for d in all_runs_data])
    proof_sizes = np.array([d['proof_size_KB'] for d in all_runs_data])

    # Calculate statistics per signature count
    unique_counts = np.unique(signatures)
    prove_means = []
    verify_means = []
    size_means = []
    prove_stds = []
    verify_stds = []

    for count in unique_counts:
        mask = signatures == count
        prove_means.append(np.nanmean(prove_times[mask]))
        verify_means.append(np.nanmean(verify_times[mask]))
        size_means.append(np.nanmean(proof_sizes[mask]))
        prove_stds.append(np.nanstd(prove_times[mask]))
        verify_stds.append(np.nanstd(verify_times[mask]))

    # Create the plot
    fig, (ax2) = plt.subplots(1, 1, figsize=(12, 12), sharex=True)

    # # Plot timing data on first subplot
    # ax1.plot(unique_counts, prove_means, 'b-o', label='Proving Time', linewidth=2, markersize=8)
    # ax1.fill_between(unique_counts,
    #                 np.array(prove_means) - np.array(prove_stds),
    #                 np.array(prove_means) + np.array(prove_stds),
    #                 color='blue', alpha=0.2)

    # ax1.plot(unique_counts, verify_means, 'r-s', label='Verification Time', linewidth=2, markersize=8)
    # ax1.fill_between(unique_counts,
    #                 np.array(verify_means) - np.array(verify_stds),
    #                 np.array(verify_means) + np.array(verify_stds),
    #                 color='red', alpha=0.2)

    # ax1.set_ylabel('Time (seconds)', fontsize=12)
    # ax1.set_xlabel('Nº Falcon signatures', fontsize=12)
    # ax1.set_title(f'Aggregate Signature Performance', fontsize=14)
    # ax1.legend(fontsize=10)
    # ax1.grid(True, linestyle='--', alpha=0.7)

    # # Plot proof size on second subplot
    ax2.plot(unique_counts, size_means, 'g-^', label='Proof Size (mean)', linewidth=2, markersize=8)
    ax2.set_xlabel('Nº Falcon signatures', fontsize=12)
    ax2.set_ylabel('Proof Size (KB)', fontsize=12)
    ax2.legend(fontsize=10)
    ax2.grid(True, linestyle='--', alpha=0.7)

    plt.tight_layout()
    # plt.savefig(CONFIG["output_graph_file"], dpi=300, bbox_inches='tight')
    print(f"Graph saved to {CONFIG['output_graph_file']}")
    plt.show()

# Main execution
if __name__ == "__main__":
    all_runs_data = load_data()

    if all_runs_data is None:
        all_runs_data = []
        run_tests()
        save_data()
    else:
        print(f"Loading existing data from {CONFIG['output_data_file']}")

    create_graph()
    print_markdown_table()
