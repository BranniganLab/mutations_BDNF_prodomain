import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

"""
Plots odds ratios using the original and blocking-analysis effective sample sizes.

The comparison uses the Lmin=4, H*=0.37 parameterization. Odds-ratio files
contain the contingency-table counts, which are scaled to each effective
sample size before calculating confidence intervals and Fisher-test p-values.
"""


def read_effective_sample_size(path):
    """
    Reads the integer blocking-analysis effective sample size.

    Arguments
        path : pathlib.Path
            Path to the effective sample-size result file.

    Returns
        tuple
            The blocking effective sample size and selected plateau block size.
    """

    values = {}
    with open(path) as file:
        for line in file:
            key, value = line.strip().split("=", 1)
            values[key] = value

    return int(values["n_independent_samples"]), int(values["plateau_block_size"])


def scale_contingency_table(counts, n_independent_samples):
    """
    Scales contingency-table counts to an effective sample size.

    Arguments
        counts : numpy.ndarray
            The four raw contingency-table counts.
        n_independent_samples : int
            The effective sample size.

    Returns
        numpy.ndarray
            Integer counts whose sum equals the effective sample size.
    """

    counts = np.asarray(counts, dtype=float)
    scaled = counts * (n_independent_samples / counts.sum())
    scaled_int = np.floor(scaled).astype(int)
    remainder = n_independent_samples - scaled_int.sum()
    for index in np.argsort(scaled - scaled_int)[::-1][:remainder]:
        scaled_int[index] += 1
    return scaled_int


def odds_ratio_and_ci(path, n_independent_samples):
    """
    Calculates the odds ratio, confidence interval, and p-value for one file.

    Arguments
        path : pathlib.Path
            Path to the odds-ratio data file.
        n_independent_samples : int
            Effective sample size used to scale the counts.

    Returns
        tuple
            Odds ratio, lower confidence limit, upper confidence limit, and p-value.
    """

    row = pd.read_csv(path, sep=r"\s+", header=None).iloc[0]
    counts = row.iloc[1:5].astype(float).to_numpy()
    scaled_counts = scale_contingency_table(counts, n_independent_samples)
    a2, b2, c2, d2 = scaled_counts
    odds_ratio, p_value = stats.fisher_exact([[a2, b2], [c2, d2]])
    standard_error = np.sqrt(np.sum(1 / counts))
    standard_error *= np.sqrt(counts.sum() / n_independent_samples)
    ci_lower = np.exp(np.log(odds_ratio) - 1.96 * standard_error)
    ci_upper = np.exp(np.log(odds_ratio) + 1.96 * standard_error)
    return odds_ratio, ci_lower, ci_upper, p_value


def significance_stars(p_value):
    """
    Converts a p-value to a significance-star label.

    Arguments
        p_value : float
            The p-value to format.

    Returns
        str
            A significance-star label.
    """

    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    if p_value < 0.05:
        return "*"
    return ""


def plot_comparison(input_path, output_path, original_n_samples, blocking_n_samples, plateau_block_size, sequences):
    """
    Plots odds ratios using the original and blocking effective sample sizes.

    Arguments
        input_path : pathlib.Path
            Directory containing the odds-ratio files to compare.
        output_path : pathlib.Path
            Path for the output PDF.
        original_n_samples : int
            Original effective sample size.
        blocking_n_samples : int
            Effective sample size estimated by blocking.
        plateau_block_size : int
            Block size selected by the blocking analysis.
        sequences : list[str]
            Sequence names to include.

    Returns
        None
    """

    fig, axes = plt.subplots(1, len(ASSOCIATIONS), figsize=(22, 7.2), sharey=True, gridspec_kw={"wspace": 0.06})
    sample_sizes = [original_n_samples, blocking_n_samples]
    labels = [f"n={original_n_samples}", f"n={blocking_n_samples} (block={plateau_block_size})"]
    width = 0.32
    x = np.arange(len(sequences))

    for axis, association in zip(axes, ASSOCIATIONS):
        for condition_index, n_samples in enumerate(sample_sizes):
            odds, lower, upper, p_values = [], [], [], []
            for sequence in sequences:
                filename = f"{sequence}_wrapped_centered_PIF_equil_del_odds_ratios_{FILE_MAPPINGS[association]}.txt"
                path = input_path / filename
                if not path.exists():
                    raise FileNotFoundError(path)
                value, ci_low, ci_high, p_value = odds_ratio_and_ci(path, n_samples)
                odds.append(value)
                lower.append(value - ci_low)
                upper.append(ci_high - value)
                p_values.append(p_value)

            positions = x + (condition_index - 0.5) * width
            for position, value, error, p_value, sequence in zip(positions, odds, upper, p_values, sequences):
                axis.errorbar(position, value, yerr=error, fmt="o", color=SAMPLE_SIZE_COLORS[condition_index], capsize=3, label=labels[condition_index] if sequence == sequences[0] else None)
                stars = significance_stars(p_value)
                if stars:
                    axis.text(position, value + error + 0.12, stars, ha="center", va="bottom", fontsize=8)

        axis.axhline(1, color="black", linestyle="--", linewidth=0.9)
        axis.set_title(association, fontsize=18)
        axis.set_xticks(x)
        axis.set_xticklabels(sequences, fontsize=15)
        axis.set_ylim(0, 6.8)
        axis.tick_params(axis="both", labelsize=18)

    axes[0].set_ylabel("OR", fontsize=18)
    axes[0].legend(loc="upper left", bbox_to_anchor=(0, -0.08), ncol=1, frameon=False, fontsize=12)
    
    fig.subplots_adjust(left=0.045, right=0.995, bottom=0.15, top=0.88, wspace=0.06)
    fig.savefig(output_path, bbox_inches="tight")
    
    plt.close(fig)


def main():
    """
    Reads the blocking effective sample size and generates the comparison plot.

    Arguments
        None

    Returns
        None
    """

    parser = argparse.ArgumentParser(description="Compare odds ratios using round-trip and blocking effective sample sizes.")
    
    parser.add_argument("--input_path", required=True, help="Directory containing the odds-ratio files to compare.")
    parser.add_argument("--output_path", required=True, help="Directory for the output PDF and effective sample-size file.")
    parser.add_argument("--effective_sample_size_file", required=False, default=None, help="Blocking effective sample-size result file.")
    parser.add_argument("--original_n_samples", type=int, default=1088, help="Original effective sample size.")
    parser.add_argument("--sequences", default=" ".join(SEQUENCES), help="Space-separated sequence names.")
    
    args = parser.parse_args()

    output_directory = Path(args.output_path)
    effective_sample_size_file = Path(args.effective_sample_size_file) if args.effective_sample_size_file else output_directory / "effective_sample_size_M66.txt"
    original_n_samples = args.original_n_samples
    sequences = args.sequences.split()
    blocking_n_samples, plateau_block_size = read_effective_sample_size(effective_sample_size_file)
    input_path = Path(args.input_path)

    output_filename = output_directory/"odds_ratios_effective_sample_size_comparison.pdf"

    plot_comparison(input_path, output_filename, original_n_samples, blocking_n_samples, plateau_block_size, sequences)

if __name__ == "__main__":

    ASSOCIATIONS = [
    r"$v$-$n$ $\rightarrow$ $n$-$c$",
    r"$v$-$c$ $\rightarrow$ $n$-$c$",
    r"$v$-$m$ $\rightarrow$ $m$-$n$",
    r"$m$-$n$ $\rightarrow$ $n$-$c$",
    ]
    FILE_MAPPINGS = {
        ASSOCIATIONS[0]: "8_2_to_2_16",
        ASSOCIATIONS[1]: "8_16_to_2_16",
        ASSOCIATIONS[2]: "8_12_to_12_2",
        ASSOCIATIONS[3]: "12_2_to_2_16",
    }
    SEQUENCES = ["F66", "M66", "V66", "L66", "A66", "Y66", "I66"]
    SAMPLE_SIZE_COLORS = ["black", "red"]

    main()
