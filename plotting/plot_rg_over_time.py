import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse

"""
This script plots the radius of gyration over time.
"""

def read_data(sequence):
    """
    Reads a file containing the radius of gyration data for a given sequence (PIF_del_full_traj_wrapped_centered.xtc file).
    
    Arguments
        sequence : str
            The sequence name.

    Returns
        data : pandas.DataFrame
            A DataFrame containing the data.
    """

    data = pd.read_csv(f'{input_path}/rg_{sequence}_wrapped_centered_PIF_del.txt', sep=r'\s', header=None, names=['data'], engine='python')

    return data

def average_data(data):
    """
    Computes the average of the given data.

    Arguments  
        data : pandas.DataFrame
            Contains the data

    Returns
        data['data'].mean() : float
            The average of the data.
    """
    
    return data['data'].mean()

def plot_rg_line_graph():
    """
    Plots the data for a given sequence as a line graph with a rolling average applied to smooth the data.

    Arguments  
        None

    Returns
        None
    """

    fig, ax = plt.subplots(figsize=(10, 6))

    for sequence in sequences: # ['V66', 'M66', 'F66', 'A66', 'Y66', 'L66', 'I66']: # 
        rg_data = read_data(sequence)
        time_in_ps = np.arange(len(rg_data)) * ps_per_frame   # Convert frames to time in picoseconds
        time_in_ns = time_in_ps / 1000  # Convert time to nanoseconds

        smoothed_data = rg_data.rolling(window).mean()

        ax.plot(time_in_ns, smoothed_data, c=seq_colors[sequence], label=sequence, linewidth=2)

    ax.legend(markerscale=10, fontsize=15, title_fontsize=15, loc='upper right')

    ax.set_ylabel(r'$R_{g}$ (Å)', fontsize=20)
    ax.set_xlabel('Time (ns)', fontsize=20)

    ax.set_ylim(12.1, 29.7)
    ax.set_xlim(0, 2080)

    ax.tick_params(axis='both', labelsize=18)

    plt.axvspan(0, 800, color='lightgrey', alpha=0.5)  # Adjust alpha for transparency

    plt.tight_layout()
    plt.savefig(f"{output_path}/Rg_vs_time_PIF_del_full_traj_wrapped_centered.pdf")
    plt.close()

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Plot the radius of gyration over time.")

    parser.add_argument("--input_path", required=True, help="Path to the input data.") 
    parser.add_argument("--output_path", required=True, help="Path to the output data.")
    parser.add_argument("--sequences", required=False, default="F66 M66 V66 L66 A66 Y66 I66", help="Sequence names.")
    parser.add_argument("--window", required=False, default=1000, help="Size of the rolling window.")
    parser.add_argument("--ps_per_frame", required=False, default=100, help="Time in picoseconds per frame.")

    args = parser.parse_args()

    input_path = args.input_path
    output_path = args.output_path
    sequences = args.sequences.split()
    window = int(args.window)
    ps_per_frame = float(args.ps_per_frame)

    seq_colors = {
        'F66': 'gray', 
        'M66': 'deepskyblue', 
        'V66': 'tomato', 
        'L66': 'darkcyan', 
        'A66': 'orchid', 
        'Y66': 'limegreen', 
        'I66': 'darkviolet'
        }
    
    plot_rg_line_graph()

