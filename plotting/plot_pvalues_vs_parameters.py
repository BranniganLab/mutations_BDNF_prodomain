import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

"""
Plot scaled Fisher-test p-values versus blobulation parameters.

The odds-ratio files contain the four contingency-table counts used to
calculate the odds ratio. This script rescales those counts to the requested
number of independent samples, matching plot_all_odds_ratios.py, and applies
Fisher's exact test to the scaled table.
"""

def scaled_contingency_table(counts, n_independent_samples):
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
    if total <= 0:
        raise ValueError("The contingency-table counts sum to zero.")

    scaled = counts * (n_independent_samples / total)
    scaled_int = np.floor(scaled).astype(int)
    remainder = n_independent_samples - scaled_int.sum()

    if remainder > 0:
        order = np.argsort(scaled - scaled_int)[::-1]
        for index in order[:remainder]:
            scaled_int[index] += 1

    return scaled_int


def fisher_pvalue(path, n_independent_samples):
    """
    Calculates the Fisher exact-test p-value for one odds-ratio file.

    Arguments
        path : pathlib.Path
            Path to the odds-ratio data file.
        n_independent_samples : int
            The number of statistically independent samples.

    Returns
        p_value : float
            The Fisher exact-test p-value calculated from scaled counts.
    """

    row = pd.read_csv(path, sep=r"\s+", header=None).iloc[0]
    counts = row.iloc[1:5].astype(float).to_numpy()
    a, b, c, d = scaled_contingency_table(counts, n_independent_samples)
    _, p_value = stats.fisher_exact([[a, b], [c, d]])
    return p_value


def resolve_input_path(input_path):
    """
    Resolves the input data directory.

    Arguments
        input_path : str or pathlib.Path
            Root directory containing the parameterization directories or a
            replica_0 subdirectory.

    Returns
        pathlib.Path
            The directory containing the parameterization directories.
    """

    input_path = Path(input_path)
    first_directory = DATASETS[0][2]
    if (input_path / first_directory).is_dir():
        return input_path

    replica_root = input_path / "replica_0"
    if (replica_root / first_directory).is_dir():
        return replica_root

    raise FileNotFoundError(f"Could not find {first_directory} below {input_path} or below {replica_root}.")


def collect_pvalues(input_path, sequences, n_independent_samples):
    """
    Collects p-values for all requested sequences and parameterizations.

    Arguments
        input_path : pathlib.Path
            Directory containing the parameterization directories.
        sequences : list[str]
            Sequence names to include.
        n_independent_samples : int
            The number of statistically independent samples.

    Returns
        dict
            P-values keyed by association, sequence, Lmin, and Hstar.
    """

    values = {association: {} for association in ASSOCIATIONS}
    for lmin, hstar, directory, mappings in DATASETS:
        for association in ASSOCIATIONS:
            for sequence in sequences:
                filename = f"{sequence}_wrapped_centered_PIF_equil_del_odds_ratios_{mappings[association]}.txt"
                path = input_path / directory / filename
                if not path.exists():
                    raise FileNotFoundError(path)

                p_value = fisher_pvalue(path, n_independent_samples)
                values[association][(sequence, lmin, hstar)] = p_value

    return values


def plot_pvalues_vs_hstar_for_lmins(values, output_path, selected_lmins, sequences):
    """
    Plots p-values versus Hstar for multiple selected Lmin values.

    Arguments
        values : dict
            P-values keyed by association, sequence, Lmin, and Hstar.
        output_path : pathlib.Path
            Path for the output PDF.
        selected_lmins : list[int]
            The Lmin values to hold constant, with one row per value.
        sequences : list[str]
            Sequence names to include in the plot.

    Returns
        None
    """

    fig, axes = plt.subplots(len(selected_lmins), 4, figsize=(13, 4.2 * len(selected_lmins)), squeeze=False, sharey=True)

    for row, selected_lmin in enumerate(selected_lmins):
        hstar_values = sorted({dataset[1] for dataset in DATASETS if dataset[0] == selected_lmin})

        for column, association in enumerate(ASSOCIATIONS):
            axis = axes[row, column]
            for sequence in sequences:
                color = SEQ_COLORS[sequence]
                p_values = []
                for hstar in hstar_values:
                    key = (sequence, selected_lmin, hstar)
                    p_values.append(values[association][key])

                axis.plot(hstar_values, p_values, marker="o", color=color, linewidth=1.5, label=sequence)

            axis.axhline(0.05, color="black", linestyle="--", linewidth=3)
            axis.set_yscale("log")
            axis.set_xlabel("H*", fontsize=18)
            axis.set_ylabel(f"p-value, Lmin = {selected_lmin}", fontsize=18)
            axis.set_title(association, fontsize=12)
            axis.set_xticks(hstar_values)
            axis.tick_params(axis="both", labelsize=12)
            axis.grid(True, axis="y", alpha=0.25)
            if column > 0:
                axis.set_ylabel("")
                axis.tick_params(axis="y", labelleft=False)

    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower left", ncol=4, bbox_to_anchor=(0, -0.05), frameon=False, fontsize=18)

    fig.subplots_adjust(left=0.06, right=0.99, bottom=0.09, top=0.96, wspace=0.08, hspace=0.3)
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def plot_pvalues_vs_lmin_for_hstars(values, output_path, selected_hstars, sequences):
    """
    Plots p-values versus Lmin for multiple selected Hstar values.

    Arguments
        values : dict
            P-values keyed by association, sequence, Lmin, and Hstar.
        output_path : pathlib.Path
            Path for the output PDF.
        selected_hstars : list[float]
            The Hstar values to hold constant, with one row per value.
        sequences : list[str]
            Sequence names to include in the plot.

    Returns
        None
    """

    fig, axes = plt.subplots(len(selected_hstars), 4, figsize=(13, 4.2 * len(selected_hstars)), squeeze=False, sharey=True)

    for row, selected_hstar in enumerate(selected_hstars):
        lmin_values = sorted({dataset[0] for dataset in DATASETS if dataset[1] == selected_hstar})
        hstar_label = f"H* = {selected_hstar}"

        for column, association in enumerate(ASSOCIATIONS):
            axis = axes[row, column]
            
            for sequence in sequences:
                color = SEQ_COLORS[sequence]
                p_values = []
                
                for lmin in lmin_values:
                    key = (sequence, lmin, selected_hstar)
                    p_values.append(values[association][key])

                axis.plot(lmin_values, p_values, marker="o", color=color, linewidth=1.5, label=sequence)

            axis.axhline(0.05, color="black", linestyle="--", linewidth=3)
            axis.set_yscale("log")
            axis.set_xlabel("Lmin", fontsize=18)
            axis.set_ylabel(f"p-value, {hstar_label}", fontsize=18)
            axis.set_title(association, fontsize=12)
            axis.set_xticks(lmin_values)
            axis.tick_params(axis="both", labelsize=12)
            axis.grid(True, axis="y", alpha=0.25)
            if column > 0:
                axis.set_ylabel("")
                axis.tick_params(axis="y", labelleft=False)

    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower left", ncol=4, bbox_to_anchor=(0, -0.05), frameon=False, fontsize=18)

    fig.subplots_adjust(left=0.06, right=0.99, bottom=0.09, top=0.96, wspace=0.08, hspace=0.3)
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def main():
    """
    Reads odds-ratio files and generates p-value comparison plots.

    Arguments
        None

    Returns
        None
    """

    parser = argparse.ArgumentParser(description="Plot scaled Fisher p-values versus Lmin and Hstar.")

    parser.add_argument("--input_path", required=True, help="Root containing replica_0/Lmin*_Hstar* directories.")
    parser.add_argument("--output_path", required=True, help="Directory for the output PDF files.")
    parser.add_argument("--n_independent_samples", type=int, default=1088, help="Effective independent sample count (1088).")
    parser.add_argument("--sequences", default=" ".join(SEQUENCES), help="Space-separated sequence names.")
    
    args = parser.parse_args()

    input_path = resolve_input_path(args.input_path)
    output_path = Path(args.output_path)
    output_path.mkdir(parents=True, exist_ok=True)
    sequences = args.sequences.split()

    values = collect_pvalues(input_path, sequences, args.n_independent_samples)
    plot_pvalues_vs_hstar_for_lmins(values, output_path / "pvalues_vs_Hstar.pdf", selected_lmins=[3, 4, 5], sequences=sequences)
    plot_pvalues_vs_lmin_for_hstars(values, output_path / "pvalues_vs_Lmin.pdf", selected_hstars=[0.33, 0.37, 0.40], sequences=sequences)


if __name__ == "__main__":
    
    ASSOCIATIONS = [
        r"$v$-$n$ $\rightarrow$ $n$-$c$",
        r"$v$-$c$ $\rightarrow$ $n$-$c$",
        r"$v$-$m$ $\rightarrow$ $m$-$n$",
        r"$m$-$n$ $\rightarrow$ $n$-$c$",
    ]

    # Each entry is (Lmin, Hstar, data directory, association-to-filename map).
    DATASETS = [
        (3, 0.33, "Lmin3_Hstar0.33", {
            ASSOCIATIONS[0]: "10_2_to_2_20",
            ASSOCIATIONS[1]: "10_20_to_2_20",
            ASSOCIATIONS[2]: "10_16_to_16_2",
            ASSOCIATIONS[3]: "16_2_to_2_20",
        }),
        (3, 0.37, "Lmin3_Hstar0.37", {
            ASSOCIATIONS[0]: "10_2_to_2_18",
            ASSOCIATIONS[1]: "10_18_to_2_18",
            ASSOCIATIONS[2]: "10_14_to_14_2",
            ASSOCIATIONS[3]: "14_2_to_2_18",
        }),
        (3, 0.40, "Lmin3_Hstar0.40", {
            ASSOCIATIONS[0]: "10_2_to_2_18",
            ASSOCIATIONS[1]: "10_18_to_2_18",
            ASSOCIATIONS[2]: "10_14_to_14_2",
            ASSOCIATIONS[3]: "14_2_to_2_18",
        }),
        (4, 0.37, "Lmin4_Hstar0.37_original_params/cutoff_5.5A_original_params", {
            ASSOCIATIONS[0]: "8_2_to_2_16",
            ASSOCIATIONS[1]: "8_16_to_2_16",
            ASSOCIATIONS[2]: "8_12_to_12_2",
            ASSOCIATIONS[3]: "12_2_to_2_16",
        }),
        (4, 0.33, "Lmin4_Hstar0.33", {
            ASSOCIATIONS[0]: "6_2_to_2_14",
            ASSOCIATIONS[1]: "6_14_to_2_14",
            ASSOCIATIONS[2]: "6_10_to_10_2",
            ASSOCIATIONS[3]: "10_2_to_2_14",
        }),
        (4, 0.40, "Lmin4_Hstar0.40", {
            ASSOCIATIONS[0]: "8_2_to_2_15",
            ASSOCIATIONS[1]: "8_15_to_2_15",
            ASSOCIATIONS[2]: "8_10_to_10_2",
            ASSOCIATIONS[3]: "10_2_to_2_15",
        }),
        (5, 0.33, "Lmin5_Hstar0.33", {
            ASSOCIATIONS[0]: "6_2_to_2_12",
            ASSOCIATIONS[1]: "6_12_to_2_12",
            ASSOCIATIONS[2]: "6_8_to_8_2",
            ASSOCIATIONS[3]: "8_2_to_2_12",
        }),
        (5, 0.37, "Lmin5_Hstar0.37", {
            ASSOCIATIONS[0]: "8_2_to_2_14",
            ASSOCIATIONS[1]: "8_14_to_2_14",
            ASSOCIATIONS[2]: "8_10_to_10_2",
            ASSOCIATIONS[3]: "10_2_to_2_14",
        }),
    ]

    SEQUENCES = ["F66", "M66", "V66", "L66", "A66", "Y66", "I66"]
    SEQ_COLORS = {
        "F66": "gray",
        "M66": "deepskyblue",
        "V66": "tomato",
        "L66": "darkcyan",
        "A66": "orchid",
        "Y66": "limegreen",
        "I66": "darkviolet",
    }
    main()
