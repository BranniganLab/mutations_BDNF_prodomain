import numpy as np
import matplotlib.pyplot as plt
import argparse

"""
This script plots the cumulative average Rg over time.
"""

def read_center_file(sequence):
    """
    Reads a file containing the radius of gyration data for a given sequence.
    
    Arguments
        None

    Returns
        data : numpy.ndarray
            An array containing the data.
    """

    data = np.loadtxt(f'{input_path}/rg_{sequence}_wrapped_centered_PIF_equil_del.txt')

    return data

def calc_cumulative_rg(sequence):
    """
    Calculates the cumulative average Rg at every frame.
    
    Arguments:
        None

    Returns:
        cumulative_contact_frequency : numpy.ndarray
            Array of cumulative contact frequencies.
    """

    rg_data = read_center_file(sequence)

    cumulative_rg = (np.cumsum(rg_data) / np.arange(1, len(rg_data) + 1))

    return cumulative_rg

def fig_gen():
    """
    Generates a scatter plot of the cumulative Rg over time.
    
    Arguments:
        None

    Returns:
        None
    """

    fig, ax = plt.subplots(figsize=(10, 6))

    for sequence in sequences:
        cumulative_rg = calc_cumulative_rg(sequence)

        time_in_ps = np.arange(len(cumulative_rg)) * ps_per_frame # Convert frames to time in picoseconds
        time_in_ns = time_in_ps / 1000  # Convert time to nanoseconds

        ax.plot(time_in_ns, cumulative_rg, c=seq_colors[sequence], label=sequence, linewidth=2)
        ax.set_ylim(14.5, 22)

        ax.legend(markerscale=10, fontsize=15, title_fontsize=15, loc='upper right', ncol=2)

        ax.set_xlabel(r"t (ns)", fontsize=20)
        ax.set_ylabel(r"$R_{g}$ ($\AA$)", fontsize=20)
        ax.tick_params(axis='both', labelsize=18)

    output_filename = (f"cumulative_rg_all_sequences.pdf")

    plt.savefig(f"{output_path}/{output_filename}",bbox_inches="tight")
    plt.close()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Plot cumulative Rg over time.")
    
    parser.add_argument("--input_path", required=True, help="Path to the input data.") 
    parser.add_argument("--output_path", required=True, help="Path to the output data.")
    # parser.add_argument("--sequence", required=False, default="M66", help="Sequence name.")
    parser.add_argument("--sequences", required=False, default="F66 M66 V66 L66 A66 Y66 I66", help="Sequence names.")
    parser.add_argument("--ps_per_frame", required=False, default=100, help="Time step in picoseconds.")
    
    args = parser.parse_args()

    input_path = args.input_path
    output_path = args.output_path
    sequences = args.sequences.split()
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

    fig_gen()


