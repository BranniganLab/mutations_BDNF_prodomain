import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

"""
Calculates and plots the time autocorrelation function of block-averaged Rg.
The x-axis is the lag time between block means, expressed in microseconds.
"""

def read_rg_data(input_path, sequence):
    """
    Reads the radius-of-gyration time series for one sequence.

    Arguments
        input_path : pathlib.Path
            Directory containing the trajectory data files.
        sequence : str
            Sequence name.

    Returns
        numpy.ndarray
            Radius-of-gyration values in trajectory order.
    """

    data_file = input_path / f"rg_{sequence}_wrapped_centered_PIF_equil_del.txt"
    if not data_file.is_file():
        raise FileNotFoundError(data_file)
    
    rg_values = pd.read_csv(data_file, sep=r"\s+", header=None).iloc[:, 0].to_numpy(dtype=float)

    return rg_values


def calculate_block_means(rg_values, block_size):
    """
    Calculates mean Rg values for consecutive complete blocks.

    Arguments
        rg_values : numpy.ndarray
            Radius-of-gyration time series.
        block_size : int
            Number of frames in each block.

    Returns
        numpy.ndarray
            Mean Rg value for each complete block.
    """

    n_blocks = len(rg_values) // block_size

    return rg_values[:n_blocks * block_size].reshape(n_blocks, block_size).mean(axis=1)


def calculate_autocorrelation(block_means):
    """
    Calculates the normalized sample autocorrelation at each lag.

    Arguments
        block_means : numpy.ndarray
            Sequential mean Rg values from equal-width blocks.

    Returns
        numpy.ndarray
            Autocorrelation values from lag zero through the maximum lag.
    """

    centered = block_means - np.mean(block_means)
    correlation = np.correlate(centered, centered, mode="full")[len(centered) - 1:]
    
    return correlation / correlation[0]


def estimate_integrated_correlation_time(autocorrelation, lag_time_step):
    """
    Estimates integrated autocorrelation time by summing through the first non-positive lag.

    Arguments
        autocorrelation : numpy.ndarray
            Normalized autocorrelation values beginning at lag zero.
        lag_time_step : float
            Time between adjacent block means, in microseconds.

    Returns
        float
            Integrated autocorrelation time in microseconds.
    """

    positive_lags = autocorrelation[1:]
    first_nonpositive = np.flatnonzero(positive_lags <= 0)
    cutoff = int(first_nonpositive[0] + 1) if len(first_nonpositive) else len(autocorrelation)

    return lag_time_step * (0.5 + np.sum(autocorrelation[1:cutoff]))


def plot_autocorrelation(input_path, output_path, sequences, block_sizes, microseconds_per_frame):
    """
    Plots block-mean Rg autocorrelation against lag time.

    Arguments
        input_path : pathlib.Path
            Directory containing the trajectory data files.
        output_path : pathlib.Path
            Directory for the output PDF.
        sequences : list[str]
            Sequence names to plot.
        block_sizes : list[int]
            Block widths in frames.
        microseconds_per_frame : float
            Simulation time represented by one trajectory frame, in microseconds.

    Returns
        None
    """

    fig, axes = plt.subplots(len(block_sizes), 1, figsize=(9, 4 * len(block_sizes)), squeeze=False)
    for axis, block_size in zip(axes[:, 0], block_sizes):
        for sequence in sequences:
            rg_values = read_rg_data(input_path, sequence)
            block_means = calculate_block_means(rg_values, block_size)
            autocorrelation = calculate_autocorrelation(block_means)
            lag_time_step = block_size * microseconds_per_frame
            lag_time = np.arange(len(autocorrelation)) * lag_time_step
            correlation_time = estimate_integrated_correlation_time(autocorrelation, lag_time_step)
            axis.plot(lag_time, autocorrelation, color=SEQ_COLORS[sequence], linewidth=1.5, label=f"{sequence} (τ={correlation_time:.4g} μs)")
            axis.axvline(correlation_time, color=SEQ_COLORS[sequence], linestyle="--", linewidth=1)
        axis.axhline(0, color="black", linewidth=0.8)
        axis.set_title(f"Block width = {block_size * microseconds_per_frame:g} μs ({block_size} frames)", fontsize=13)
        axis.set_ylabel("Autocorrelation", fontsize=13)
        axis.tick_params(axis="both", labelsize=11)
        axis.grid(True, alpha=0.25)
        axis.legend(frameon=False, ncol=min(len(sequences), 4), fontsize=10)
        axis.set_xlim(0, 0.12)
    axes[-1, 0].set_xlabel("Lag time (μs)", fontsize=13)

    fig.tight_layout()
    output_path.mkdir(parents=True, exist_ok=True)
    output_file_path = output_path / OUTPUT_FILENAME
    fig.savefig(output_file_path, bbox_inches="tight")
    plt.close(fig)


def main():
    """
    Reads the arguments and generates the block-mean Rg autocorrelation plot.

    Arguments
        None

    Returns
        None
    """

    parser = argparse.ArgumentParser(description="Plot autocorrelation of block-averaged Rg for selected block widths.")
    parser.add_argument("--input_path", required=True, help="Path to the input data.")
    parser.add_argument("--output_path", required=True, help="Directory for the output PDF.")
    parser.add_argument("--sequences", default="M66", help="Space-separated sequence names.")
    parser.add_argument("--block_sizes", nargs="+", type=int, default=BLOCK_SIZES, help="Block widths in frames (default: 1 100).")
    parser.add_argument("--microseconds_per_frame", type=float, default=MICROSECONDS_PER_FRAME, help="Simulation time per trajectory frame in microseconds (0.0001 μs).")
    args = parser.parse_args()

    plot_autocorrelation(Path(args.input_path), Path(args.output_path), args.sequences.split(), args.block_sizes, args.microseconds_per_frame)


if __name__ == "__main__":

    SEQ_COLORS = {
        "F66": "gray",
        "M66": "black",
        "V66": "tomato",
        "L66": "darkcyan",
        "A66": "orchid",
        "Y66": "limegreen",
        "I66": "darkviolet",
    }

    MICROSECONDS_PER_FRAME = 1.2 / 12000
    BLOCK_SIZES = [1, 100]
    OUTPUT_FILENAME = "rg_block_autocorrelation.pdf"
    
    main()
