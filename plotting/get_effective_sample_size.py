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

def generate_block_sizes(n_data, min_blocks=10, n_sizes=10):
    """
    Generates block sizes based on the trajectory length.
    
    Arguments
        n_data : int
            Number of trajectory data points.
        min_blocks : int
            Minimum number of blocks to retain for each block size.
        n_sizes : int
            Approximate number of block sizes to generate.

    Returns
        block_sizes : list
            A list of block sizes.
    """

    if n_data < min_blocks:
        raise ValueError(
            f"The trajectory has only {n_data} points; at least "
            f"{min_blocks} points are needed."
        )

    max_block_size = n_data // min_blocks

    block_sizes = np.unique(
        np.logspace(
            np.log10(1),
            np.log10(max_block_size),
            n_sizes
        ).astype(int)
    )

    return block_sizes.tolist()

def block_trajectory(data, block_size):
    """
    Blocks the trajectory data into blocks of a given size.
    
    Arguments
        data : pandas.DataFrame
            A DataFrame containing the data.
        block_size : int
            The size of the blocks.

    Returns
        blocked_data : pandas.DataFrame
            A DataFrame containing the blocked data.
    """

    n_complete = len(data) // block_size
    trimmed_data = data.iloc[:n_complete * block_size].copy()

    blocked_data = trimmed_data.groupby(np.arange(len(trimmed_data)) // block_size).mean()

    return blocked_data

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

    block_values = data["data"].to_numpy()
    n_blocks = len(block_values)

    sem_block_avg = np.std(block_values, ddof=1) / np.sqrt(n_blocks)

    return sem_block_avg

def estimate_effective_sample_size(data, block_size):
    """
    Estimate effective sample size from a selected blocking plateau.
    
    Arguments
        data : pandas.DataFrame
            A DataFrame containing the data.
        block_size : int
            The size of the blocks.

    Returns
        n_eff : float
            The estimated effective sample size.
        sem : float
            The standard error of the mean.
        n_blocks : int
            The number of blocks.
    """
    raw_values = data["data"].to_numpy()
    raw_std = np.std(raw_values, ddof=1)

    blocked_data = block_trajectory(data, block_size)
    blocked_sem = do_stats(blocked_data)

    if not np.isfinite(blocked_sem) or blocked_sem <= 0:
        return np.nan, blocked_sem, len(blocked_data)

    n_eff = (raw_std / blocked_sem) ** 2
    return n_eff, blocked_sem, len(blocked_data)


def find_plateau_block_size(
    data,
    block_sizes,
    min_blocks=20,
    window=4,
    relative_tolerance=0.05,
):
    """Find a block size in the first approximately stable SEM region.

    A plateau is defined here as ``window`` neighboring SEM values whose
    range is no more than ``relative_tolerance`` of their median SEM.
    This is a heuristic and should be checked against the plotted curve.
    """
    results = []

    for block_size in sorted(set(block_sizes)):
        n_blocks = len(data) // block_size
        if n_blocks < min_blocks:
            continue

        blocked_data = block_trajectory(data, block_size)
        sem = do_stats(blocked_data)
        if np.isfinite(sem):
            results.append((block_size, sem, n_blocks))

    if len(results) < window:
        return None, results

    for start in range(len(results) - window + 1):
        candidate = results[start:start + window]
        sems = np.array([result[1] for result in candidate])
        median_sem = np.median(sems)
        relative_range = (np.max(sems) - np.min(sems)) / median_sem

        if relative_range <= relative_tolerance:
            plateau_block_size = int(np.median(
                [result[0] for result in candidate]
            ))
            return plateau_block_size, results

    return None, results

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

def plot_rg_sem_vs_block_size(rg_data, block_sizes, reference_n_samples=1088, plateau_block_size=None,):
    """
    Plots the radius of gyration data as a bar chart.

    Arguments  
        None

    Returns
        None
    """

    fig, ax = plt.subplots(figsize=(8, 6))

    sem_values = []

    for block_size in block_sizes:
        rg_blocked_data = block_trajectory(rg_data, block_size)
        sem = do_stats(rg_blocked_data)
        sem_values.append(sem)

    ax.plot(block_sizes, sem_values, "o", markersize=8, color="black")

    raw_std = np.std(rg_data["data"].to_numpy(), ddof=1)
    reference_sem = raw_std / np.sqrt(reference_n_samples)
    ax.axhline(reference_sem, color="red", linestyle="--", label=f"Reference: n = {reference_n_samples}",)

    if plateau_block_size is not None:
        n_eff, plateau_sem, n_blocks = estimate_effective_sample_size(
            rg_data, plateau_block_size
        )
        print(
            f"Plateau block size: {plateau_block_size}; "
            f"number of blocks: {n_blocks}; "
            f"plateau SEM: {plateau_sem:.6g}; "
            f"estimated effective sample size: {n_eff:.6g}"
        )
        ax.axvline(
            plateau_block_size,
            color="blue",
            linestyle=":",
            label=f"Plateau block size = {plateau_block_size}",
        )
    
    # ax.set_ylim(0.02, 0.14)
    ax.set_ylabel(r'Standard error of the mean', fontsize=20)
    ax.set_xlabel(r'Block size', fontsize=20)
    ax.legend()
    ax.tick_params(axis='both', labelsize=18)

    plt.tight_layout()
    plt.savefig(f"{output_path}/rg_sem_vs_block_size_{sequence}.pdf")
    plt.close()


def main(args):
    """Run the blocking analysis and generate the SEM plot."""
    global input_path, output_path, sequence

    input_path = args.input_path
    output_path = args.output_path
    sequence = args.sequence
    rg_data = read_data(sequence)

    if args.block_sizes is None:
        block_sizes = generate_block_sizes(
            len(rg_data), min_blocks=20, n_sizes=30
        )
    else:
        block_sizes = [
            block_size for block_size in args.block_sizes
            if block_size >= 1 and len(rg_data) // block_size >= 20
        ]

        if not block_sizes:
            raise ValueError("No requested block size leaves at least 20 blocks.")

    print(f"Trajectory length: {len(rg_data)}")
    print(f"Block sizes: {block_sizes}")
    print(
        f"Reference sample size: {args.reference_n_samples} "
        "(its SEM is shown as a red dashed line)."
    )

    if args.plateau_block_size is not None:
        if args.plateau_block_size < 1:
            raise ValueError("--plateau_block_size must be at least 1.")
        if len(rg_data) // args.plateau_block_size < 20:
            raise ValueError(
                "The plateau block size must leave at least 20 complete blocks."
            )
        selected_plateau_block_size = args.plateau_block_size
    else:
        selected_plateau_block_size, plateau_results = find_plateau_block_size(
            rg_data,
            block_sizes,
            min_blocks=20,
            window=4,
            relative_tolerance=0.05,
        )

        print("Automatic plateau search:")
        for block_size, sem, n_blocks in plateau_results:
            print(
                f"  block size = {block_size}, "
                f"blocks = {n_blocks}, SEM = {sem:.6g}"
            )

        if selected_plateau_block_size is None:
            print(
                "No stable plateau was detected automatically. "
                "Inspect the plot and provide --plateau_block_size manually."
            )
        else:
            print(
                f"Automatically selected plateau block size: "
                f"{selected_plateau_block_size}"
            )

    plot_rg_sem_vs_block_size(
        rg_data,
        block_sizes,
        reference_n_samples=args.reference_n_samples,
        plateau_block_size=selected_plateau_block_size,
    )


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Plot the average radius of gyration.")

    parser.add_argument("--input_path", required=True, help="Path to the input data.") 
    parser.add_argument("--output_path", required=True, help="Path to the output data.")
    parser.add_argument("--sequence", required=False, default="M66", help="Sequence names.")
    parser.add_argument("--reference_n_samples", type=int, default=1088, help="Sample size used in the original SEM calculation.",)
    parser.add_argument("--plateau_block_size", type=int, default=None, help="Block size selected from the SEM plateau for estimating N_eff.",)
    parser.add_argument("--block_sizes", nargs='+', type=int, default=None, help="Optional block sizes. If omitted, sizes are generated automatically.")

    args = parser.parse_args()
    main(args)
