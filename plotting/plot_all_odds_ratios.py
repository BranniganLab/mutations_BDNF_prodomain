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

    global data_root, output_directory, n_independent_samples
    global SEQUENCES, DATASETS, ASSOCIATIONS, OUTPUT_FILENAME

    # Accept either .../data or .../data/replica_0 as --data_root.
    if not (data_root / DATASETS[0][1]).is_dir():
        replica_root = data_root / "replica_0"
        if (replica_root / DATASETS[0][1]).is_dir():
            data_root = replica_root
        else:
            raise FileNotFoundError(
                f"Could not find {DATASETS[0][1]} below {data_root} "
                "or below its replica_0 subdirectory."
            )
    n_conditions = len(DATASETS)
    x = np.arange(len(SEQUENCES))
    width = 0.15

    fig, axes = plt.subplots(1, len(ASSOCIATIONS), figsize=(22, 7.2),
                             sharey=True, gridspec_kw={"wspace": 0.06})
    for ax, association in zip(axes, ASSOCIATIONS):
        for condition_index, (label, directory, mapping) in enumerate(DATASETS):
            odds, lower, upper, p_values = [], [], [], []
            for sequence in SEQUENCES:
                filename = (f"{sequence}_wrapped_centered_PIF_equil_del_"
                            f"odds_ratios_{mapping[association]}.txt")
                path = data_root / directory / filename
                if not path.exists():
                    raise FileNotFoundError(path)
                value, ci_low, ci_high, p = odds_ratio_and_ci(
                    path, n_independent_samples)
                odds.append(value)
                lower.append(value - ci_low)
                upper.append(ci_high - value)
                p_values.append(p)

            positions = x + (condition_index - (n_conditions - 1) / 2) * width
            ax.bar(positions, odds, width=width, color=colors[condition_index],
                   yerr=[lower, upper], capsize=2, label=label,
                   error_kw={"elinewidth": 0.8})

            for position, value, error, p_value in zip(
                    positions, odds, upper, p_values):
                stars = significance_stars(p_value)
                if stars:
                    ax.text(position, value + error + 0.12, stars,
                            ha="center", va="bottom", fontsize=8)

        ax.axhline(1, color="black", linestyle="--", linewidth=0.9)
        ax.set_title(association, fontsize=11)
        ax.set_xticks(x)
        ax.set_xticklabels(SEQUENCES, fontsize=15)
        ax.set_ylim(0, 10)

    axes[0].set_ylabel("Odds ratio", fontsize=13)
    axes[0].legend(loc="upper left", bbox_to_anchor=(0, -0.16), ncol=3,frameon=False, fontsize=10)
    fig.subplots_adjust(left=0.045, right=0.995, bottom=0.15, top=0.88,
                        wspace=0.06)

    output_path = Path(output_directory) / OUTPUT_FILENAME
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight", format="pdf")
    plt.close(fig)


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Plot odds ratios for several Lmin and Hstar parameter sets.")
    parser.add_argument("--data_root", required=True, help="Root containing replica_0/Lmin*_Hstar* directories")
    parser.add_argument("--output", required=True, help="Directory in which to save the output PDF.")
    parser.add_argument("--n_independent_samples", required=False, type=int, default=1088, help="Number of independent samples.")

    args = parser.parse_args()

    data_root = Path(args.data_root)
    output_directory = Path(args.output)
    n_independent_samples = args.n_independent_samples
    sequences = ["F66", "M66", "V66", "L66", "A66", "Y66", "I66"]

    # Each tuple is (display label, directory name, physical association ->
    # filename association). These mappings are inferred from the blob indices
    # used to generate the files in each directory.
    datasets = [
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
    colors = ["#00000053", "#0072B2", "#E98E17", "#009607", "#C9468E"]
    OUTPUT_FILENAME = "odds_ratios_comparison.pdf"
    DATASETS = datasets
    ASSOCIATIONS = associations
    SEQUENCES = sequences

    main()
