import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

"""
Calculates bootstrap confidence intervals for Pearson correlations between
sequence-level averages of Rg and Ree, and between Ree and N-C contact frequency.
The bootstrap resamples the sequence-level points used in the correlation.
"""


def read_average(input_path, sequence, data_type):
    """
    Reads a trajectory data file and returns its average.

    Arguments
        input_path : pathlib.Path
            Directory containing the trajectory data files.
        sequence : str
            Sequence name.
        data_type : str
            Data file prefix, such as rg, ree, or single_pair_blob_contact_2_14.

    Returns
        float
            The mean of the values in the data file.
    """

    data_file = input_path / f"{data_type}_{sequence}_wrapped_centered_PIF_equil_del.txt"
    if not data_file.is_file():
        raise FileNotFoundError(data_file)
    data = pd.read_csv(data_file, sep=r"\s+", header=None).iloc[:, 0].to_numpy(dtype=float)
    if data.size == 0 or not np.all(np.isfinite(data)):
        raise ValueError(f"Data file is empty or contains non-finite values: {data_file}")
    return float(np.mean(data))


def bootstrap_pearson_ci(x_values, y_values, n_bootstrap, seed):
    """
    Calculates a percentile bootstrap confidence interval for Pearson's r.

    Arguments
        x_values : array-like
            X coordinates of the sequence-level observations.
        y_values : array-like
            Y coordinates of the sequence-level observations.
        n_bootstrap : int
            Number of bootstrap resamples.
        seed : int
            Random seed for reproducibility.

    Returns
        tuple
            Observed Pearson r, lower confidence bound, upper confidence bound, and p-value.
    """

    x_values = np.asarray(x_values, dtype=float)
    y_values = np.asarray(y_values, dtype=float)
    if len(x_values) != len(y_values) or len(x_values) < 3:
        raise ValueError("Pearson correlation requires at least three paired sequence values.")
    if np.std(x_values) == 0 or np.std(y_values) == 0:
        raise ValueError("Pearson correlation is undefined when either variable has no variation.")

    result = stats.pearsonr(x_values, y_values)
    rng = np.random.default_rng(seed)
    sample_indices = rng.integers(0, len(x_values), size=(n_bootstrap, len(x_values)))
    bootstrap_r = []
    for indices in sample_indices:
        x_sample = x_values[indices]
        y_sample = y_values[indices]
        if np.std(x_sample) > 0 and np.std(y_sample) > 0:
            bootstrap_r.append(stats.pearsonr(x_sample, y_sample).statistic)

    if not bootstrap_r:
        raise ValueError("All bootstrap resamples had zero variation; confidence interval is undefined.")
    alpha = 0.05
    ci_lower, ci_upper = np.percentile(bootstrap_r, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return result.statistic, ci_lower, ci_upper, result.pvalue


def plot_correlations(correlations, output_path):
    """
    Plots the sequence-level observations for both correlations.

    Arguments
        correlations : list[dict]
            Correlation labels, data, axis labels, and bootstrap statistics.
        output_path : pathlib.Path
            Directory in which to save the figure.

    Returns
        None
    """

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
    for axis, correlation in zip(axes, correlations):
        axis.scatter(correlation["x_values"], correlation["y_values"], color="black", s=55)
        for sequence, x_value, y_value in zip(correlation["sequences"], correlation["x_values"], correlation["y_values"]):
            axis.annotate(sequence, (x_value, y_value), textcoords="offset points", xytext=(6, 6), fontsize=10)
        axis.set_xlabel(correlation["x_label"], fontsize=13)
        axis.set_ylabel(correlation["y_label"], fontsize=13)
        axis.tick_params(axis="both", labelsize=11)
        axis.text(0.04, 0.96, f"r = {correlation['r']:.3f}\n95% CI [{correlation['ci_lower']:.3f}, {correlation['ci_upper']:.3f}]\np = {correlation['p_value']:.3g}", transform=axis.transAxes, ha="left", va="top", fontsize=10)

    fig.tight_layout()
    fig.savefig(output_path / "rg_termini_correlation_bootstrap_ci.pdf", bbox_inches="tight")
    plt.close(fig)


def main():
    """
    Calculates and saves bootstrap confidence intervals for the two correlations.

    Arguments
        None

    Returns
        None
    """

    parser = argparse.ArgumentParser(description="Bootstrap Pearson correlation confidence intervals across sequence-level averages.")
    parser.add_argument("--input_path", required=True, help="Path to the input data.")
    parser.add_argument("--output_path", required=True, help="Directory for the correlation results.")
    parser.add_argument("--N_blob_index", required=True, type=int, help="Index of the N-terminal blob.")
    parser.add_argument("--C_blob_index", required=True, type=int, help="Index of the C-terminal blob.")
    parser.add_argument("--sequences", default="V66 M66 F66 A66 Y66 L66 I66", help="Space-separated sequence names.")
    parser.add_argument("--n_bootstrap", type=int, default=10000, help="Number of bootstrap resamples.")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility.")
    args = parser.parse_args()

    input_path = Path(args.input_path)
    output_path = Path(args.output_path)
    sequences = args.sequences.split()
    if args.n_bootstrap < 1:
        raise ValueError("n_bootstrap must be positive.")

    rg_values = [read_average(input_path, sequence, "rg") for sequence in sequences]
    ree_values = [read_average(input_path, sequence, "ree") for sequence in sequences]
    contact_type = f"single_pair_blob_contact_{args.N_blob_index}_{args.C_blob_index}"
    contact_values = [100 * read_average(input_path, sequence, contact_type) for sequence in sequences]
    correlation_inputs = [
        ("Rg_vs_Ree", ree_values, rg_values, r"Mean $R_{ee}$ (Å)", r"Mean $R_g$ (Å)"),
        ("Ree_vs_NC_contact_frequency", contact_values, ree_values, "N-C contact frequency (%)", r"Mean $R_{ee}$ (Å)"),
    ]
    correlations = []
    for index, (name, x_values, y_values, x_label, y_label) in enumerate(correlation_inputs):
        r_value, ci_lower, ci_upper, p_value = bootstrap_pearson_ci(x_values, y_values, args.n_bootstrap, args.seed + index)
        correlations.append({"name": name, "sequences": sequences, "x_values": x_values, "y_values": y_values, "x_label": x_label, "y_label": y_label, "r": r_value, "ci_lower": ci_lower, "ci_upper": ci_upper, "p_value": p_value})

    output_path.mkdir(parents=True, exist_ok=True)
    results_file = output_path / "rg_termini_correlation_bootstrap_ci.txt"
    with results_file.open("w") as results:
        results.write("Correlation\tN_sequences\tPearson_r\tCI_lower\tCI_upper\tP_value\tConfidence\tBootstrap_resamples\n")
        for correlation in correlations:
            results.write(f"{correlation['name']}\t{len(sequences)}\t{correlation['r']:.6g}\t{correlation['ci_lower']:.6g}\t{correlation['ci_upper']:.6g}\t{correlation['p_value']:.6g}\t0.95\t{args.n_bootstrap}\n")
    plot_correlations(correlations, output_path)


if __name__ == "__main__":
    main()
