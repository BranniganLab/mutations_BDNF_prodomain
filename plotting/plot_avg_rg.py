import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse

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

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Plot the average radius of gyration.")

    parser.add_argument("--input_path", required=True, help="Path to the input data.") 
    parser.add_argument("--output_path", required=True, help="Path to the output data.")
    parser.add_argument("--n_samples", required=False, default=1088, help="Number of samples (1088).")
    parser.add_argument("--sequences", required=False, default="F66 M66 V66 L66 A66 Y66 I66", help="Sequence names.")

    args = parser.parse_args()

    input_path = args.input_path
    output_path = args.output_path
    n_independent_samples = int(args.n_samples)
    sequences = args.sequences.split()

    plot_avg_rg()

