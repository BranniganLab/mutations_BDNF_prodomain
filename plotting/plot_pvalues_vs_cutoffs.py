import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

"""
Plot scaled Fisher-test p-values versus the blob-contact cutoff distance.

The odds-ratio files contain the four contingency-table counts used to
calculate the odds ratio. This script rescales those counts to the requested
number of independent samples, matching plot_pvalues_vs_parameters.py, and
applies Fisher's exact test to the scaled table.
"""


def scale_contingency_table(counts, n_independent_samples):
    """
    Scales frame counts to the effective sample size and preserves the sum.

    Arguments
        counts : numpy.ndarray
            The four raw contingency-table counts.
        n_independent_samples : int
            The number of statistically independent samples.

    Returns
        scaled_int : numpy.ndarray
            Integer contingency-table counts that sum to the effective sample size.
    """

    counts = np.asarray(counts, dtype=float)
    total = counts.sum()
    scaled = counts * (n_independent_samples / total)
    scaled_int = np.floor(scaled).astype(int)
    remainder = n_independent_samples - scaled_int.sum()

    if remainder > 0:
        order = np.argsort(scaled - scaled_int)[::-1]
        for index in order[:remainder]:
            scaled_int[index] += 1

    return scaled_int


def get_fisher_pvalue(file_path, n_independent_samples):
    """
    Calculates the Fisher exact-test p-value for one odds-ratio file.

    Arguments
        file_path : pathlib.Path
            Path to the odds-ratio data file.
        n_independent_samples : int
            The number of statistically independent samples.

    Returns
        p_value : float
            The Fisher exact-test p-value calculated from scaled counts.
    """

    row = pd.read_csv(file_path, sep=r"\s+", header=None).iloc[0]
    counts = row.iloc[1:5].astype(float).to_numpy()
    a, b, c, d = scale_contingency_table(counts, n_independent_samples)
    _, p_value = stats.fisher_exact([[a, b], [c, d]])
    
    return p_value


def resolve_input_path(input_path):
    """Resolves the root containing the cutoff directories."""

    input_path = Path(input_path)
    if (input_path / "cutoff_5A").is_dir():
        return input_path.parent

    for candidate in [input_path, input_path / "replica_0"]:
        if (candidate / FIRSTDIRECTORY).is_dir():
            return candidate

    raise FileNotFoundError(f"Could not find {FIRSTDIRECTORY} below {input_path} or below {input_path / 'replica_0'}.")


def collect_pvalues(input_path, sequences, n_independent_samples):
    """
    Collects p-values for all requested sequences and cutoff values.

    Arguments
        input_path : pathlib.Path
            Root containing the Lmin4_Hstar0.37 cutoff directory.
        sequences : list[str]
            Sequence names to include.
        n_independent_samples : int
            The number of statistically independent samples.

    Returns
        dict
            P-values keyed by association, sequence, and cutoff.
    """

    values = {association: {} for association in ASSOCIATIONS}
    for cutoff, directory, mappings in DATASETS:
        for association in ASSOCIATIONS:
            for sequence in sequences:
                filename = f"{sequence}_wrapped_centered_PIF_equil_del_odds_ratios_{mappings[association]}.txt"
                path = input_path / directory / filename
                if not path.exists():
                    raise FileNotFoundError(path)

                values[association][(sequence, cutoff)] = get_fisher_pvalue(path, n_independent_samples)

    return values


def plot_pvalues_vs_cutoff(values, output_path, sequences):
    """Plots p-values versus cutoff distance for all associations."""

    cutoff_values = [dataset[0] for dataset in DATASETS]
    fig, axes = plt.subplots(1, len(ASSOCIATIONS), figsize=(13, 4.2), squeeze=False, sharey=True)

    for column, association in enumerate(ASSOCIATIONS):
        axis = axes[0, column]
        for sequence in sequences:
            p_values = [values[association][(sequence, cutoff)] for cutoff in cutoff_values]
            axis.plot(cutoff_values, p_values, marker="o", color=SEQCOLORS[sequence], linewidth=1.5, label=sequence)

        axis.axhline(0.05, color="black", linestyle="--", linewidth=3)
        axis.set_yscale("log")
        axis.set_xlabel("Cutoff (Å)", fontsize=18)
        axis.set_title(association, fontsize=12)
        axis.set_xticks(cutoff_values)
        axis.tick_params(axis="both", labelsize=12)
        axis.grid(True, axis="y", alpha=0.25)
        if column == 0:
            axis.set_ylabel("p-value", fontsize=18)
        else:
            axis.tick_params(axis="y", labelleft=False)

    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower left", ncol=4, bbox_to_anchor=(0, -0.25), frameon=False, fontsize=18)
    fig.subplots_adjust(left=0.06, right=0.99, bottom=0.18, top=0.9, wspace=0.08)
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def main():
    """Reads odds-ratio files and generates the cutoff p-value plot."""

    parser = argparse.ArgumentParser(description="Plot scaled Fisher p-values versus cutoff distance.")
    parser.add_argument("--input_path", required=True, help="Root containing Lmin4_Hstar0.37_original_params or replica_0.")
    parser.add_argument("--output_path", required=True, help="Directory for the output PDF file.")
    parser.add_argument("--n_independent_samples", type=int, default=1088, help="Effective independent sample count (1088).")
    parser.add_argument("--sequences", default=" ".join(SEQUENCES), help="Space-separated sequence names.")

    args = parser.parse_args()
    input_path = resolve_input_path(args.input_path)
    output_path = Path(args.output_path)
    output_path.mkdir(parents=True, exist_ok=True)
    sequences = args.sequences.split()

    n_independent_samples = args.n_independent_samples
    values = collect_pvalues(input_path, sequences, n_independent_samples)
    plot_pvalues_vs_cutoff(values, output_path / OUTPUTFILENAME, sequences)


if __name__ == "__main__":
    ASSOCIATIONS = [
        r"$v$-$n$ $\rightarrow$ $n$-$c$",
        r"$v$-$c$ $\rightarrow$ $n$-$c$",
        r"$v$-$m$ $\rightarrow$ $m$-$n$",
        r"$m$-$n$ $\rightarrow$ $n$-$c$",
    ]

    DATASETS = [
        (5.0, "Lmin4_Hstar0.37_original_params/cutoff_5A", {
            ASSOCIATIONS[0]: "8_2_to_2_16",
            ASSOCIATIONS[1]: "8_16_to_2_16",
            ASSOCIATIONS[2]: "8_12_to_12_2",
            ASSOCIATIONS[3]: "12_2_to_2_16",
        }),
        (5.5, "Lmin4_Hstar0.37_original_params/cutoff_5.5A_original_params", {
            ASSOCIATIONS[0]: "8_2_to_2_16",
            ASSOCIATIONS[1]: "8_16_to_2_16",
            ASSOCIATIONS[2]: "8_12_to_12_2",
            ASSOCIATIONS[3]: "12_2_to_2_16",
        }),
        (6.0, "Lmin4_Hstar0.37_original_params/cutoff_6A", {
            ASSOCIATIONS[0]: "8_2_to_2_16",
            ASSOCIATIONS[1]: "8_16_to_2_16",
            ASSOCIATIONS[2]: "8_12_to_12_2",
            ASSOCIATIONS[3]: "12_2_to_2_16",
        }),
    ]

    FIRSTDIRECTORY = DATASETS[0][1]
    SEQUENCES = ["F66", "M66", "V66", "L66", "A66", "Y66", "I66"]
    SEQCOLORS = {
        "F66": "gray",
        "M66": "deepskyblue",
        "V66": "tomato",
        "L66": "darkcyan",
        "A66": "orchid",
        "Y66": "limegreen",
        "I66": "darkviolet",
    }
    OUTPUTFILENAME = "pvalues_vs_cutoff.pdf"
    main()
