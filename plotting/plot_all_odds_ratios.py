import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

"""
This script plots the odds ratios for several Lmin and Hstar parameter sets
in one figure. The blob numbers change between parameterizations, so the
associations are specified by their physical meaning rather than by one fixed
filename.
"""

def odds_ratio_and_ci(path, n_independent_samples):
    """
    Computes the odds ratio and confidence interval from an odds-ratio file.

    Arguments
        path : pathlib.Path
            Path to the odds-ratio data file.
        n_independent_samples : int
            Number of statistically independent samples used to correct the
            confidence interval.

    Returns
        odds_ratio : float
            The odds ratio calculated using the scaled contingency table.
        ci_lower : float
            The lower bound of the 95% confidence interval.
        ci_upper : float
            The upper bound of the 95% confidence interval.
        p_value : float
            The p-value from Fisher's exact test.
    """

    row = pd.read_csv(path, sep=r"\s+", header=None).iloc[0]
    a, b, c, d = map(float, row.iloc[1:5])

    counts = np.array([a, b, c, d])
    scaled = counts * (n_independent_samples / counts.sum())
    scaled_int = np.floor(scaled).astype(int)
    remainder = n_independent_samples - scaled_int.sum()
    for i in np.argsort(scaled - scaled_int)[::-1][:remainder]:
        scaled_int[i] += 1

    a2, b2, c2, d2 = scaled_int
    odds_ratio, p_value = stats.fisher_exact([[a2, b2], [c2, d2]])

    se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    se *= np.sqrt(counts.sum() / n_independent_samples)
    ci_lower = np.exp(np.log(odds_ratio) - 1.96 * se)
    ci_upper = np.exp(np.log(odds_ratio) + 1.96 * se)
    return odds_ratio, ci_lower, ci_upper, p_value


def significance_stars(p_value):
    """
    Converts a p-value to a significance-star label.

    Arguments
        p_value : float
            The p-value to format.

    Returns
        stars : str
            ``***`` for p < 0.001, ``**`` for p < 0.01, ``*`` for
            p < 0.05, or an empty string for non-significant results.
    """

    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    if p_value < 0.05:
        return "*"
    return ""


def main():
    """
    Reads the odds-ratio data and creates a combined comparison plot.

    Arguments
        None

    Returns
        None
    """

    global input_path, output_path, n_independent_samples
    global SEQUENCES, DATASETS, ASSOCIATIONS, OUTPUT_FILENAME

    # Accept either .../data or .../data/replica_0 as --input_path.
    if not (input_path / DATASETS[0][1]).is_dir():
        replica_root = input_path / "replica_0"
        if (replica_root / DATASETS[0][1]).is_dir():
            input_path = replica_root
        else:
            raise FileNotFoundError(
                f"Could not find {DATASETS[0][1]} below {input_path} "
                "or below its replica_0 subdirectory."
            )
    x = np.arange(len(SEQUENCES))
    width = 0.15
    hstar_groups = [
        (r"$\mathrm{H}^*=0.33$", [0, 4, 6]),
        (r"$\mathrm{H}^*=0.37$", [1, 3, 7]),
        (r"$\mathrm{H}^*=0.40$", [2, 5]),
    ]
    lmin_colors = {3: "#0072B2", 4: "#E98E17", 5: "#009607"}

    fig, axes = plt.subplots(len(hstar_groups), len(ASSOCIATIONS), figsize=(22, 15), sharey=True, gridspec_kw={"wspace": 0.06, "hspace": 0.24})
    
    for row_index, (hstar_label, dataset_indices) in enumerate(hstar_groups):
        
        for column_index, association in enumerate(ASSOCIATIONS):
            ax = axes[row_index, column_index]
            n_conditions = len(dataset_indices)
            
            for condition_index, dataset_index in enumerate(dataset_indices):
                label, directory, mapping = DATASETS[dataset_index]
                odds, lower, upper, p_values = [], [], [], []
                
                for sequence in SEQUENCES:
                    filename = (f"{sequence}_wrapped_centered_PIF_equil_del_"
                                f"odds_ratios_{mapping[association]}.txt")
                    path = input_path / directory / filename
                    if not path.exists():
                        raise FileNotFoundError(path)
                    value, ci_low, ci_high, p = odds_ratio_and_ci(path, n_independent_samples)
                    odds.append(value)
                    lower.append(value - ci_low)
                    upper.append(ci_high - value)
                    p_values.append(p)

                lmin = int(directory.split("Lmin")[1].split("_")[0])
                positions = x + (condition_index - (n_conditions - 1) / 2) * width
                ax.bar(positions, odds, width=width, color=lmin_colors[lmin], yerr=[lower, upper], capsize=2, label=f"$L_{{min}}={lmin}$", error_kw={"elinewidth": 0.8})

                for position, value, error, p_value in zip(positions, odds, upper, p_values):
                    stars = significance_stars(p_value)
                    if stars:
                        ax.text(position, value + error + 0.12, stars, ha="center", va="bottom", fontsize=8)

            ax.axhline(1, color="black", linestyle="--", linewidth=0.9)
            ax.set_title(association, fontsize=15)
            ax.set_xticks(x)
            ax.tick_params(axis="both", labelsize=15)
            ax.set_xticklabels(SEQUENCES, fontsize=15)
            ax.set_ylim(0, 10)
            if column_index == 0:
                ax.set_ylabel(f"OR, {hstar_label}", fontsize=18)

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower left", bbox_to_anchor=(0.06, 0.05), ncol=3, frameon=False, fontsize=18)
    fig.subplots_adjust(left=0.06, right=0.995, bottom=0.12, top=0.95, wspace=0.06, hspace=0.24)

    output_file_path = output_path / OUTPUT_FILENAME
    output_file_path.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(output_file_path, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Plot odds ratios for several Lmin and Hstar parameter sets.")
    parser.add_argument("--input_path", required=True, help="Root containing replica_0/Lmin*_Hstar* directories")
    parser.add_argument("--output_path", required=True, help="Directory in which to save the output PDF.")
    parser.add_argument("--n_independent_samples", required=False, type=int, default=1088, help="Number of independent samples.")
    parser.add_argument("--sequences", required=False, default="F66 M66 V66 L66 A66 Y66 I66", help="Space-separated sequence names.")

    args = parser.parse_args()

    input_path = Path(args.input_path)
    output_path = Path(args.output_path)
    n_independent_samples = args.n_independent_samples
    sequences = args.sequences.split()

    # Each tuple is (display label, directory name, physical association ->
    # filename association). These mappings are inferred from the blob indices
    # used to generate the files in each directory.
    datasets = [
        ("$L_{min}=3$, $H^*=0.33$", "Lmin3_Hstar0.33", {
            "$v$-$n$ $\\rightarrow$ $n$-$c$": "10_2_to_2_20",
            "$v$-$c$ $\\rightarrow$ $n$-$c$": "10_20_to_2_20",
            "$v$-$m$ $\\rightarrow$ $m$-$n$": "10_16_to_16_2",
            "$m$-$n$ $\\rightarrow$ $n$-$c$": "16_2_to_2_20",
        }),
        ("$L_{min}=3$, $H^*=0.37$", "Lmin3_Hstar0.37", {
            "$v$-$n$ $\\rightarrow$ $n$-$c$": "10_2_to_2_18",
            "$v$-$c$ $\\rightarrow$ $n$-$c$": "10_18_to_2_18",
            "$v$-$m$ $\\rightarrow$ $m$-$n$": "10_14_to_14_2",
            "$m$-$n$ $\\rightarrow$ $n$-$c$": "14_2_to_2_18",
        }),
        ("$L_{min}=3$, $H^*=0.40$", "Lmin3_Hstar0.40", {
            "$v$-$n$ $\\rightarrow$ $n$-$c$": "10_2_to_2_18",
            "$v$-$c$ $\\rightarrow$ $n$-$c$": "10_18_to_2_18",
            "$v$-$m$ $\\rightarrow$ $m$-$n$": "10_14_to_14_2",
            "$m$-$n$ $\\rightarrow$ $n$-$c$": "14_2_to_2_18",
        }),
        ("$L_{min}=4$, $H^*=0.37$",
         "Lmin4_Hstar0.37_original_params/cutoff_5.5A_original_params", {
             "$v$-$n$ $\\rightarrow$ $n$-$c$": "8_2_to_2_16",
             "$v$-$c$ $\\rightarrow$ $n$-$c$": "8_16_to_2_16",
             "$v$-$m$ $\\rightarrow$ $m$-$n$": "8_12_to_12_2",
             "$m$-$n$ $\\rightarrow$ $n$-$c$": "12_2_to_2_16",
         }),
        ("$L_{min}=4$, $H^*=0.33$", "Lmin4_Hstar0.33", {
            "$v$-$n$ $\\rightarrow$ $n$-$c$": "6_2_to_2_14",
            "$v$-$c$ $\\rightarrow$ $n$-$c$": "6_14_to_2_14",
            "$v$-$m$ $\\rightarrow$ $m$-$n$": "6_10_to_10_2",
            "$m$-$n$ $\\rightarrow$ $n$-$c$": "10_2_to_2_14",
        }),
        ("$L_{min}=4$, $H^*=0.40$", "Lmin4_Hstar0.40", {
            "$v$-$n$ $\\rightarrow$ $n$-$c$": "8_2_to_2_15",
            "$v$-$c$ $\\rightarrow$ $n$-$c$": "8_15_to_2_15",
            "$v$-$m$ $\\rightarrow$ $m$-$n$": "8_10_to_10_2",
            "$m$-$n$ $\\rightarrow$ $n$-$c$": "10_2_to_2_15",
        }),
        ("$L_{min}=5$, $H^*=0.33$", "Lmin5_Hstar0.33", {
            "$v$-$n$ $\\rightarrow$ $n$-$c$": "6_2_to_2_12",
            "$v$-$c$ $\\rightarrow$ $n$-$c$": "6_12_to_2_12",
            "$v$-$m$ $\\rightarrow$ $m$-$n$": "6_8_to_8_2",
            "$m$-$n$ $\\rightarrow$ $n$-$c$": "8_2_to_2_12",
        }),
        ("$L_{min}=5$, $H^*=0.37$", "Lmin5_Hstar0.37", {
            "$v$-$n$ $\\rightarrow$ $n$-$c$": "8_2_to_2_14",
            "$v$-$c$ $\\rightarrow$ $n$-$c$": "8_14_to_2_14",
            "$v$-$m$ $\\rightarrow$ $m$-$n$": "8_10_to_10_2",
            "$m$-$n$ $\\rightarrow$ $n$-$c$": "10_2_to_2_14",
        }),
    ]

    associations = [
        "$v$-$n$ $\\rightarrow$ $n$-$c$",
        "$v$-$c$ $\\rightarrow$ $n$-$c$",
        "$v$-$m$ $\\rightarrow$ $m$-$n$",
        "$m$-$n$ $\\rightarrow$ $n$-$c$",
    ]
    OUTPUT_FILENAME = "odds_ratios_comparison.pdf"
    DATASETS = datasets
    ASSOCIATIONS = associations
    SEQUENCES = sequences

    main()
