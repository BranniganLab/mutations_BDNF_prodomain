import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SEQ_COLORS = {"F66": "gray", "M66": "black", "V66": "tomato", "L66": "darkcyan", "A66": "orchid", "Y66": "limegreen", "I66": "darkviolet"}
MICROSECONDS_PER_FRAME = 1.2 / 12000

"""
Plots the mean Rg of each consecutive trajectory block for selected block widths.
Block widths and block locations are expressed in microseconds, using the configured frame duration.
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
    rg_values = pd.read_csv(data_file, sep=r"\s+", header=None).iloc[:, 0].to_numpy(dtype=float)
    
    return rg_values


def calculate_block_means(rg_values, block_size):
    """
    Calculates the mean Rg and midpoint frame for each complete block.

    Arguments
        rg_values : numpy.ndarray
            Radius-of-gyration time series.
        block_size : int
            Number of consecutive frames in each block.

    Returns
        tuple
            Block midpoint frames and corresponding mean Rg values.
    """

    n_blocks = len(rg_values) // block_size
    used_values = rg_values[:n_blocks * block_size].reshape(n_blocks, block_size)
    block_means = used_values.mean(axis=1)
    block_midpoints = (np.arange(n_blocks) + 0.5) * block_size
    
    return block_midpoints, block_means


def plot_block_means(input_path, output_path, sequences, block_sizes, microseconds_per_frame):
    """
    Plots block-level mean Rg values against block midpoint frame.

    Arguments
        input_path : pathlib.Path
            Directory containing the trajectory data files.
        output_path : pathlib.Path
            Directory for the output PDF.
        sequences : list[str]
            Sequence names to plot.
        block_sizes : list[int]
            Block widths in frames, with one subplot per width.
        microseconds_per_frame : float
            Simulation time represented by one trajectory frame, in microseconds.

    Returns
        None
    """

    fig, axes = plt.subplots(len(block_sizes), 1, figsize=(10, 3.5 * len(block_sizes)), squeeze=False, sharex=True)
    
    for axis, block_size in zip(axes[:, 0], block_sizes):
        
        for sequence in sequences:
            rg_values = read_rg_data(input_path, sequence)
            block_midpoints, block_means = calculate_block_means(rg_values, block_size)
            axis.plot(block_midpoints * microseconds_per_frame, block_means, marker="o", markersize=3, linewidth=1, color=SEQ_COLORS[sequence], label=sequence)
        
        axis.set_ylabel(r"Block mean $R_g$ (Å)", fontsize=13)
        axis.set_title(f"Block width = {block_size * microseconds_per_frame:g} μs ({block_size} frames)", fontsize=13)
        axis.tick_params(axis="both", labelsize=11)
        # axis.grid(True, alpha=0.25)
    
    axes[-1, 0].set_xlabel("Block midpoint (μs)", fontsize=13)
    axes[0, 0].legend(frameon=False, ncol=min(len(sequences), 4), fontsize=10)
    fig.tight_layout()
    
    output_path.mkdir(parents=True, exist_ok=True)
    
    fig.savefig(output_path / "rg_block_means_vs_block_time.pdf", bbox_inches="tight")
    
    plt.close(fig)


def main():
    """
    Reads the arguments and generates the block-mean Rg plot.

    Arguments
        None

    Returns
        None
    """

    parser = argparse.ArgumentParser(description="Plot mean Rg for consecutive trajectory blocks.")
    parser.add_argument("--input_path", required=True, help="Path to the input data.")
    parser.add_argument("--output_path", required=True, help="Directory for the output PDF.")
    parser.add_argument("--sequences", default="M66", help="Space-separated sequence names.")
    parser.add_argument("--block_sizes", nargs="+", type=int, default=[1, 10, 100, 500, 1000], help="Block widths in frames.")
    parser.add_argument("--microseconds_per_frame", type=float, default=MICROSECONDS_PER_FRAME, help="Simulation time per trajectory frame in microseconds (0.0001 μs).")
    args = parser.parse_args()

    plot_block_means(Path(args.input_path), Path(args.output_path), args.sequences.split(), args.block_sizes, args.microseconds_per_frame)


if __name__ == "__main__":
    main()
