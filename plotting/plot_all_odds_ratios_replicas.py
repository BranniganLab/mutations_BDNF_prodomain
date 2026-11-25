import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import plot_odds_ratios as odds_ratio_plot

"""
This script plots odds ratios for all sequences as a grouped bar chart for two replicas.
"""

def get_significance_stars(p_value):
    """
    Converts a p-value to a significance-star label.

    Arguments
        p_value : float
            The p-value to label.

    Returns
        stars : str
            The significance-star label.
    """

    if p_value < 0.001:
        stars = "***"
    elif p_value < 0.01:
        stars = "**"
    elif p_value < 0.05:
        stars = "*"
    else:
        stars = ""

    return stars

def plot_grouped_odds_ratios(replica_paths, sequences, output_file_path, n_independent_samples, y_maximum):
    """
    Plots grouped odds-ratio bars for all sequences and replicas.

    Arguments
        replica_paths : list of pathlib.Path
            The data directories for the replicas.
        sequences : list of str
            The sequences to plot.
        output_file_path : pathlib.Path
            The output figure path.
        n_independent_samples : int
            The number of independent samples.
        y_maximum : float
            The maximum y-axis value.

    Returns
        None
    """

    replica_labels = ["Replica 0", "Replica 9"]
    replica_colors = ["#0072B2", "#E98E17"]
    x_positions = np.arange(len(sequences))
    bar_width = 0.35
    figure, axes = plt.subplots(2, 2, figsize=(18, 12), sharey=True)

    for association_index, association_label in enumerate(ASSOCIATION_LABELS):
        axis = axes.flat[association_index]

        for replica_index, data_path in enumerate(replica_paths):
            data_association = REPLICA_ASSOCIATIONS[replica_index][association_index]
            odds_ratios = []
            ci_lowers = []
            ci_uppers = []
            p_values = []

            odds_ratio_plot.input_path = str(data_path)
            odds_ratio_plot.association = data_association
            odds_ratio_plot.n_independent_samples = n_independent_samples

            for sequence in sequences:
                data = odds_ratio_plot.read_OR_file(sequence)
                state_a = data.iloc[0, 1]
                state_b = data.iloc[0, 2]
                state_c = data.iloc[0, 3]
                state_d = data.iloc[0, 4]
                odds_ratio, p_value = odds_ratio_plot.compute_fisher_exact(state_a, state_b, state_c, state_d)
                ci_lower, ci_upper = odds_ratio_plot.calc_confidence_intervals(odds_ratio, state_a, state_b, state_c, state_d)
                odds_ratios.append(odds_ratio)
                ci_lowers.append(ci_lower)
                ci_uppers.append(ci_upper)
                p_values.append(p_value)

            odds_ratios = np.array(odds_ratios)
            ci_lowers = np.array(ci_lowers)
            ci_uppers = np.array(ci_uppers)
            positions = x_positions + (replica_index - 0.5) * bar_width
            errors = [odds_ratios - ci_lowers, ci_uppers - odds_ratios]

            axis.bar(positions, odds_ratios, width=bar_width, yerr=errors, capsize=4, color=replica_colors[replica_index], label=replica_labels[replica_index])

            for position, ci_upper, p_value in zip(positions, ci_uppers, p_values):
                stars = get_significance_stars(p_value)

                if stars:
                    axis.text(position, ci_upper, stars, ha="center", va="bottom", fontsize=15)

        axis.axhline(1.0, linestyle="--", color="black")
        axis.set_title(association_label, fontsize=18)
        axis.set_xticks(x_positions)
        axis.set_xticklabels(sequences)
        axis.set_ylim(0, y_maximum)
        axis.tick_params(axis="both", labelsize=18)

    axes[0, 0].set_ylabel("OR", fontsize=20)
    axes[1, 0].set_ylabel("OR", fontsize=20)
    
    handles, labels = axes[0, 0].get_legend_handles_labels()
    
    figure.legend(handles, labels, loc="lower left", bbox_to_anchor=(0.04, 0.015), ncol=2, fontsize=18, frameon=False)
    figure.tight_layout(rect=[0, 0.08, 1, 1])
    figure.savefig(output_file_path)
    
    plt.close(figure)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Plot grouped odds ratios for all sequences and two replicas.")
    parser.add_argument("--input_path", required=True, help="Path to the data directory containing replica_0 and replica_9.")
    parser.add_argument("--parameter_path", default="Lmin4_Hstar0.37_original_params/cutoff_5.5A_original_params", help="Shared parameter subdirectory below each replica directory.")
    parser.add_argument("--output_path", required=True, help="Directory for the output PDF.")
    parser.add_argument("--n_independent_samples", type=int, default=1088, help="Number of independent samples.")
    parser.add_argument("--sequences", default="F66 M66 V66 L66 A66 Y66 I66", help="Space-separated sequence names.")
    parser.add_argument("--y_maximum", type=float, default=6.5, help="Maximum y-axis value.")

    args = parser.parse_args()

    REPLICA_PATHS = [Path(args.input_path) / "replica_0" / args.parameter_path, Path(args.input_path) / "replica_9" / args.parameter_path]
    REPLICA_ASSOCIATIONS = [
        ["8_2_to_2_16", "8_16_to_2_16", "8_12_to_12_2", "12_2_to_2_16"],
        ["8_2_to_2_16", "8_16_to_2_16", "8_12_to_12_2", "12_2_to_2_16"],
    ]
    ASSOCIATION_LABELS = [
        r"$v$-$n$ $\rightarrow$ $n$-$c$",
        r"$v$-$c$ $\rightarrow$ $n$-$c$",
        r"$v$-$m$ $\rightarrow$ $m$-$n$",
        r"$m$-$n$ $\rightarrow$ $n$-$c$",
    ]
    SEQUENCES = args.sequences.split()
    OUTPUT_FILENAME = "odds_ratios_replicas_0_9.pdf"
    OUTPUT_PATH = Path(args.output_path) / OUTPUT_FILENAME
    N_INDEPENDENT_SAMPLES = args.n_independent_samples
    Y_MAXIMUM = args.y_maximum

    plot_grouped_odds_ratios(REPLICA_PATHS, SEQUENCES, OUTPUT_PATH, N_INDEPENDENT_SAMPLES, Y_MAXIMUM)
