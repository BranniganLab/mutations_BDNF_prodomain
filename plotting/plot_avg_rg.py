import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse
from pathlib import Path

"""
This script plots the average radius of gyration 
"""

def read_data(sequence):
    """
    Reads a file containing the radius of gyration data for a given sequence.
    
    Arguments
        sequence : str
            The sequence name.

    Returns
        data : pandas.DataFrame
            A DataFrame containing the data.
    """

    data = pd.read_csv(f'{input_path}/rg_{sequence}_wrapped_centered_PIF_equil_del.txt', sep=r'\s', header=None, names=['data'], engine='python')

    return data

def do_stats(data):
    """
    Compute the standard error of the mean (SEM) for the given data.
    
    Arguments
        data : pandas.DataFrame
            A DataFrame containing the data.

    Returns
        sem : float
            The standard error of the mean.
    """

    data = data['data'].to_numpy()
    std = np.std(data, ddof=1)
    sem = std / np.sqrt(n_independent_samples)

    return float(sem)

def average_data(data):
    """
    Computes the average of the given data.

    Arguments  
        data : pandas.DataFrame
            Contains the data

    Returns
        data['data'].mean() : float
            The average data.
    """
    
    return data['data'].mean()

def plot_avg_rg():
    """
    Plots the radius of gyration data as a bar chart.

    Arguments  
        None

    Returns
        None
    """

    fig, ax = plt.subplots(figsize=(8, 6))

    for sequence in sequences:
        rg_data = read_data(sequence)
        avg_rg = average_data(rg_data)

        ax.bar(sequence, avg_rg, color='gray')
        ax.errorbar(sequence, avg_rg, yerr=do_stats(rg_data), capsize=3, elinewidth=1, color='black')
        print(sequence, round(avg_rg, 1), round(do_stats(rg_data), 2))
    
    ax.set_ylabel(r'$\langle R_{g} \rangle$ (Å)', fontsize=20)

    ax.set_ylim(17.7, 19.4)
    ax.tick_params(axis='both', labelsize=18)

    plt.tight_layout()
    plt.savefig(f"{output_path}/avg_rg.pdf")
    plt.close()

def read_distance_data(input_path, distance, sequence):
    """
    Reads radius of gyration data for a sequence and distance cutoff.

    Arguments
        input_path : pathlib.Path
            Directory containing the distance-specific data directories.
        distance : int
            Minimum-distance cutoff in Angstroms.
        sequence : str
            The sequence name.

    Returns
        data : pandas.DataFrame
            A DataFrame containing the data.
    """

    data_file = input_path / f"max_dist_{distance}A"/f"rg_{sequence}_wrapped_centered_PIF_equil_del.txt"
    data = pd.read_csv(data_file, sep=r'\s', header=None, names=['data'], engine='python')

    return data

def plot_grouped_avg_rg(input_path, output_path, distances, sequences, n_independent_samples):
    """
    Plots grouped bars of average radius of gyration for each distance cutoff.

    Arguments
        input_path : pathlib.Path
            Directory containing the distance-specific data directories.
        output_path : pathlib.Path
            Directory for the output data.
        distances : list[int]
            Minimum-distance cutoffs in Angstroms.
        sequences : list[str]
            Sequence names to plot.
        n_independent_samples : int
            Number of independent samples.

    Returns
        None
    """

    figure, axis = plt.subplots(figsize=(10, 6))
    sequence_positions = np.arange(len(sequences))
    bar_width = 0.8 / len(distances)

    for distance_index, distance in enumerate(distances):
        average_rg_values = []
        sem_values = []

        for sequence in sequences:
            rg_data = read_distance_data(input_path, distance, sequence)
            average_rg = average_data(rg_data)
            sem = do_stats_for_data(rg_data, n_independent_samples)
            average_rg_values.append(average_rg)
            sem_values.append(sem)

        bar_positions = sequence_positions + (distance_index - (len(distances) - 1) / 2) * bar_width
        axis.bar(bar_positions, average_rg_values, width=bar_width, yerr=sem_values, capsize=3, label=f"{distance} Å")

    axis.set_ylabel(r'$\langle R_{g} \rangle$ (Å)', fontsize=20)
    axis.set_ylim(17.7, 19.4)
    axis.set_xticks(sequence_positions)
    axis.set_xticklabels(sequences)
    axis.tick_params(axis='both', labelsize=18)
    axis.legend(frameon=False, fontsize=15)

    output_path.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(output_path / "avg_rg_with_18_20_22A_distances.pdf")
    plt.close()

def do_stats_for_data(data, n_independent_samples):
    """
    Computes the standard error of the mean for the given data.

    Arguments
        data : pandas.DataFrame
            A DataFrame containing the data.
        n_independent_samples : int
            Number of independent samples.

    Returns
        sem : float
            The standard error of the mean.
    """

    data_values = data['data'].to_numpy()
    standard_deviation = np.std(data_values, ddof=1)
    sem = standard_deviation / np.sqrt(n_independent_samples)

    return float(sem)

DISTANCES = [18, 20, 22]

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Plot the average radius of gyration.")

    parser.add_argument("--input_path", required=True, help="Path to the input data.") 
    parser.add_argument("--output_path", required=True, help="Path to the output data.")
    parser.add_argument("--n_independent_samples", required=False, default=1088, help="Number of samples (1088).")
    parser.add_argument("--sequences", required=False, default="F66 M66 V66 L66 A66 Y66 I66", help="Sequence names.")
    parser.add_argument("--distances", nargs="+", type=int, default=DISTANCES, help="Distance cutoffs in Angstroms for the grouped plot.")

    args = parser.parse_args()

    input_path = args.input_path
    output_path = args.output_path
    n_independent_samples = int(args.n_independent_samples)
    sequences = args.sequences.split()

    # plot_avg_rg()
    plot_grouped_avg_rg(Path(input_path), Path(output_path), args.distances, sequences, n_independent_samples)

