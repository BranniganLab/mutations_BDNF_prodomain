import numpy as np
import matplotlib.pyplot as plt
import argparse

"""
This script plots the cumulative average Rg over time.
"""

def read_center_file():
    """
    Reads a file containing the radius of gyration data for a given sequence.
    
    Arguments
        None

    Returns
        data : numpy.ndarray
            An array containing the data.
    """

    data = np.loadtxt(f'{input_path}/rg_{sequence}_wrapped_centered_PIF_del.txt')

    return data

def calc_cumulative_rg():
    """
    Calculates the cumulative average Rg at every frame.
    
    Arguments:
        None

    Returns:
        cumulative_contact_frequency : numpy.ndarray
            Array of cumulative contact frequencies.
    """

    rg_data = read_center_file()

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

    cumulative_rg = calc_cumulative_rg()

    time_in_ps = np.arange(len(cumulative_rg)) * ps_per_frame # Convert frames to time in picoseconds
    time_in_ns = time_in_ps / 1000  # Convert time to nanoseconds

    ax.plot(time_in_ns, cumulative_rg, color="black", linewidth=2)
    ax.set_ylim(14.5, 22)

    ax.set_xlabel(r"t (ns)", fontsize=20)
    ax.set_ylabel(r"$R_{g}$ ($\AA$)", fontsize=20)
    ax.tick_params(axis='both', labelsize=18)

    output_filename = (f"cumulative_rg_{sequence}.pdf")

    plt.savefig(f"{output_path}/{output_filename}",bbox_inches="tight")
    plt.close()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Plot cumulative Rg over time.")
    
    parser.add_argument("--input_path", required=True, help="Path to the input data.") 
    parser.add_argument("--output_path", required=True, help="Path to the output data.")
    parser.add_argument("--sequence", required=False, default="I66", help="Sequence name.")
    parser.add_argument("--ps_per_frame", required=False, default=100, help="Time step in picoseconds.")
    
    args = parser.parse_args()

    input_path = args.input_path
    output_path = args.output_path
    sequence = args.sequence
    ps_per_frame = float(args.ps_per_frame)

    fig_gen()


