import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

"""
This script calculates the statistical inefficiency and effective sample rate
from the autocorrelation of a radius-of-gyration time series.
"""


def read_data(input_path, sequence):
    """
    Reads the radius-of-gyration time series for a sequence.

    Arguments
        input_path : str
            Directory containing the input data.
        sequence : str
            Sequence name used in the input filename.

    Returns
        numpy.ndarray
            Radius-of-gyration values in trajectory order.
    """

    data_file = f"{input_path}/rg_{sequence}_wrapped_centered_PIF_equil_del.txt"
    data = pd.read_csv(data_file, sep=r"\s+", header=None).iloc[:, 0]

    return data.to_numpy(dtype=float)


def calculate_autocorrelation(data):
    """
    Calculates the normalized autocorrelation function of a time series.

    The calculation uses a Fourier transform so that the autocorrelation can
    be evaluated efficiently for the complete trajectory.

    Arguments
        data : numpy.ndarray
            Time-series values in trajectory order.

    Returns
        numpy.ndarray
            Normalized autocorrelation from lag zero onward.
    """

    centered_data = data - np.mean(data)
    data_length = len(centered_data)
    fft_length = 2 ** int(np.ceil(np.log2(2 * data_length - 1)))
    transformed_data = np.fft.rfft(centered_data, n=fft_length)
    autocorrelation = np.fft.irfft(
        transformed_data * np.conjugate(transformed_data),
        n=fft_length,
    )[:data_length]
    autocorrelation /= autocorrelation[0]

    return autocorrelation


def calculate_statistical_inefficiency(autocorrelation, maximum_lag):
    """
    Calculates the statistical inefficiency, 1 + 2 * tau_ac.

    The autocorrelation sum stops at the first non-positive lag so that the
    noisy long-lag tail is not included in the estimate.

    Arguments
        autocorrelation : numpy.ndarray
            Normalized autocorrelation from lag zero onward.
        maximum_lag : int
            Maximum lag in trajectory frames included in the sum.

    Returns
        float
            Statistical inefficiency for the specified maximum lag.
    """

    maximum_lag = min(maximum_lag, len(autocorrelation) - 1)
    positive_lags = autocorrelation[1:maximum_lag + 1]
    first_nonpositive = np.flatnonzero(positive_lags <= 0)

    if len(first_nonpositive) > 0:
        positive_lags = positive_lags[:first_nonpositive[0]]

    tau_ac = np.sum(positive_lags)
    statistical_inefficiency = 1 + 2 * tau_ac

    return float(statistical_inefficiency)


def get_block_widths(minimum_block_width_ns, maximum_block_width_ns, block_width_step_ns, nanoseconds_per_frame):
    """
    Converts block widths in nanoseconds into trajectory-frame widths.

    Arguments
        minimum_block_width_ns : float
            Smallest block width in nanoseconds.
        maximum_block_width_ns : float
            Largest block width in nanoseconds.
        block_width_step_ns : float
            Spacing between block widths in nanoseconds.
        nanoseconds_per_frame : float
            Simulation time represented by one trajectory frame.

    Returns
        tuple[numpy.ndarray, numpy.ndarray]
            Actual block widths in nanoseconds and corresponding frame counts.
    """

    block_widths_ns = np.arange(
        minimum_block_width_ns,
        maximum_block_width_ns + block_width_step_ns / 2,
        block_width_step_ns,
    )
    block_sizes = np.maximum(
        1,
        np.rint(block_widths_ns / nanoseconds_per_frame).astype(int),
    )
    actual_block_widths_ns = block_sizes * nanoseconds_per_frame

    return actual_block_widths_ns, block_sizes


def calculate_effective_samples_per_microsecond(statistical_inefficiencies, microseconds_per_frame):
    """
    Calculates effective independent samples per microsecond.

    The raw sampling rate is divided by the statistical inefficiency for each
    block width.

    Arguments
        statistical_inefficiencies : numpy.ndarray
            Values of 1 + 2 * tau_ac for each block width.
        microseconds_per_frame : float
            Simulation time represented by one trajectory frame.

    Returns
        numpy.ndarray
            Effective independent samples per microsecond.
    """

    samples_per_microsecond = 1 / microseconds_per_frame
    effective_samples_per_microsecond = samples_per_microsecond / statistical_inefficiencies

    return effective_samples_per_microsecond


def find_plateau_effective_samples(block_widths_ns, effective_samples_per_microsecond, relative_tolerance=0.1, minimum_plateau_points=3):
    """
    Automatically estimates the effective sample rate from the terminal plateau.

    The terminal quarter of the block-width range provides an initial plateau
    estimate. The plateau begins at the earliest block width for which all
    subsequent values remain within the requested relative tolerance. If no
    stable region is found, the terminal quarter is used and reported as an
    approximate plateau.

    Arguments
        block_widths_ns : numpy.ndarray
            Block widths in nanoseconds.
        effective_samples_per_microsecond : numpy.ndarray
            Effective sample rates for each block width.
        relative_tolerance : float
            Maximum relative deviation allowed from the terminal estimate.
        minimum_plateau_points : int
            Minimum number of values required in the plateau.

    Returns
        tuple[float, float, bool]
            Plateau effective samples per microsecond, plateau start in
            nanoseconds, and whether a stable plateau was detected.
    """

    n_values = len(effective_samples_per_microsecond)
    terminal_point_count = max(minimum_plateau_points, int(np.ceil(n_values / 4)))
    terminal_start_index = n_values - terminal_point_count
    terminal_values = effective_samples_per_microsecond[terminal_start_index:]
    terminal_effective_samples = np.median(terminal_values)
    relative_deviations = np.abs(
        effective_samples_per_microsecond - terminal_effective_samples
    ) / terminal_effective_samples
    plateau_start_index = None

    for start_index in range(n_values - minimum_plateau_points + 1):
        if np.all(relative_deviations[start_index:] <= relative_tolerance):
            plateau_start_index = start_index
            break

    if plateau_start_index is None:
        plateau_start_index = terminal_start_index
        plateau_detected = False
    else:
        plateau_detected = True

    plateau_values = effective_samples_per_microsecond[plateau_start_index:]
    plateau_effective_samples = np.median(plateau_values)
    plateau_start_ns = block_widths_ns[plateau_start_index]

    return float(plateau_effective_samples), float(plateau_start_ns), plateau_detected


def plot_statistical_inefficiency(block_widths_ns, statistical_inefficiencies, output_file):
    """
    Plots statistical inefficiency versus block width.

    Arguments
        block_widths_ns : numpy.ndarray
            Block widths in nanoseconds.
        statistical_inefficiencies : numpy.ndarray
            Values of 1 + 2 * tau_ac for each block width.
        output_file : str
            Path to the output PDF.

    Returns
        None
    """

    figure, axis = plt.subplots(figsize=(8, 6))
    axis.plot(block_widths_ns, statistical_inefficiencies, "o-")
    axis.set_xlabel("Block width (ns)", fontsize=18)
    axis.set_ylabel(r"$1 + 2\tau_{ac}$ (dimensionless)", fontsize=18)
    axis.tick_params(axis="both", labelsize=14)
    axis.set_xlim(block_widths_ns[0], block_widths_ns[-1])
    figure.tight_layout()
    figure.savefig(output_file)
    plt.close(figure)


def plot_effective_samples_per_microsecond(block_widths_ns, effective_samples_per_microsecond, plateau_effective_samples, output_file):
    """
    Plots effective independent samples per microsecond versus block width.

    The median effective sample rate from the automatically detected plateau
    is shown as a red dotted horizontal line.

    Arguments
        block_widths_ns : numpy.ndarray
            Block widths in nanoseconds.
        effective_samples_per_microsecond : numpy.ndarray
            Effective sample rates for each block width.
        plateau_effective_samples : float
            Effective sample rate estimated from the plateau.
        output_file : str
            Path to the output PDF.

    Returns
        None
    """

    figure, axis = plt.subplots(figsize=(8, 6))
    axis.plot(block_widths_ns, effective_samples_per_microsecond, "o-")
    axis.axhline(
        plateau_effective_samples,
        color="red",
        linestyle=":",
        label=f"Plateau = {plateau_effective_samples:.1f} samples/μs",
    )
    axis.set_xlabel("Block width (ns)", fontsize=18)
    axis.set_ylabel("Effective samples per μs", fontsize=18)
    axis.tick_params(axis="both", labelsize=14)
    axis.set_xlim(block_widths_ns[0], block_widths_ns[-1])
    axis.legend()
    figure.tight_layout()
    figure.savefig(output_file)
    plt.close(figure)


def plot_all_sequences_statistical_inefficiency(sequence_results, output_file):
    """
    Plots statistical inefficiency for all sequences on one figure.

    Arguments
        sequence_results : dict
            Mapping of sequence names to their calculated block widths and
            statistical inefficiencies.
        output_file : str
            Path to the output PDF.

    Returns
        None
    """

    figure, axis = plt.subplots(figsize=(8, 6))

    for sequence, results in sequence_results.items():
        sequence_color = PROTEIN_COLOR[sequence[0]]
        axis.plot(
            results["block_widths_ns"],
            results["statistical_inefficiencies"],
            "o-",
            color=sequence_color,
            label=sequence,
        )

    axis.set_xlabel("Block width (ns)", fontsize=18)
    axis.set_ylabel(r"$1 + 2\tau_{ac}$ (dimensionless)", fontsize=18)
    axis.tick_params(axis="both", labelsize=14)
    axis.legend()
    figure.tight_layout()
    figure.savefig(output_file)
    plt.close(figure)


def plot_all_sequences_effective_samples_per_microsecond(sequence_results, output_file):
    """
    Plots effective samples per microsecond for all sequences on one figure.

    Arguments
        sequence_results : dict
            Mapping of sequence names to their calculated block widths and
            effective sample rates.
        output_file : str
            Path to the output PDF.

    Returns
        None
    """

    figure, axis = plt.subplots(figsize=(8, 6))

    for sequence, results in sequence_results.items():
        sequence_color = PROTEIN_COLOR[sequence[0]]
        axis.plot(
            results["block_widths_ns"],
            results["effective_samples_per_microsecond"],
            "o-",
            color=sequence_color,
            label=f"{sequence} (n = {results['plateau_effective_samples']:.0f})",
        )
        axis.axhline(
            results["plateau_effective_samples"],
            color=sequence_color,
            linestyle=":",
            linewidth=1,
        )

    axis.set_xlabel("Block width (ns)", fontsize=18)
    axis.set_ylabel("Effective samples per μs", fontsize=18)
    axis.tick_params(axis="both", labelsize=14)
    axis.legend()
    figure.tight_layout()
    figure.savefig(output_file)
    plt.close(figure)


def calculate_sequence_results(args, sequence):
    """
    Calculates autocorrelation-based statistics for one sequence.

    Arguments
        args : argparse.Namespace
            Command-line arguments controlling the analysis.
        sequence : str
            Sequence name.

    Returns
        dict
            Block widths, statistical inefficiencies, effective sample rates,
            plateau estimate, and plateau status for the sequence.
    """

    rg_values = read_data(args.input_path, sequence)
    autocorrelation = calculate_autocorrelation(rg_values)
    block_widths_ns, block_sizes = get_block_widths(
        args.minimum_block_width_ns,
        args.maximum_block_width_ns,
        args.block_width_step_ns,
        args.nanoseconds_per_frame,
    )
    statistical_inefficiencies = np.array([
        calculate_statistical_inefficiency(autocorrelation, block_size)
        for block_size in block_sizes
    ])
    effective_samples_per_microsecond = calculate_effective_samples_per_microsecond(
        statistical_inefficiencies,
        args.microseconds_per_frame,
    )
    plateau_effective_samples, plateau_start_ns, plateau_detected = find_plateau_effective_samples(
        block_widths_ns,
        effective_samples_per_microsecond,
        args.plateau_relative_tolerance,
        args.minimum_plateau_points,
    )

    return {
        "block_widths_ns": block_widths_ns,
        "block_sizes": block_sizes,
        "statistical_inefficiencies": statistical_inefficiencies,
        "effective_samples_per_microsecond": effective_samples_per_microsecond,
        "plateau_effective_samples": plateau_effective_samples,
        "plateau_start_ns": plateau_start_ns,
        "plateau_detected": plateau_detected,
        "trajectory_time_microseconds": len(rg_values) * args.microseconds_per_frame,
    }


def main(args):
    """
    Calculates statistical inefficiency and generates both output plots.

    Arguments
        args : argparse.Namespace
            Command-line arguments controlling the input, analysis, and output.

    Returns
        None

    Output files
        statistical_inefficiency_<sequence>.pdf
            Statistical inefficiency versus block width.
        effective_samples_per_microsecond_<sequence>.pdf
            Effective samples per microsecond versus block width, including
            the automatically estimated plateau as a red dotted line.
    """

    sequence_results = {
        sequence: calculate_sequence_results(args, sequence)
        for sequence in args.sequences
    }
    results = sequence_results[args.sequence]
    block_widths_ns = results["block_widths_ns"]
    block_sizes = results["block_sizes"]
    statistical_inefficiencies = results["statistical_inefficiencies"]
    effective_samples_per_microsecond = results["effective_samples_per_microsecond"]
    plateau_effective_samples = results["plateau_effective_samples"]
    plateau_start_ns = results["plateau_start_ns"]
    plateau_detected = results["plateau_detected"]
    trajectory_time_microseconds = results["trajectory_time_microseconds"]
    plateau_effective_sample_size = plateau_effective_samples * trajectory_time_microseconds
    statistical_inefficiency_output_file = f"{args.output_path}/statistical_inefficiency_{args.sequence}.pdf"
    effective_samples_output_file = f"{args.output_path}/effective_samples_per_microsecond_{args.sequence}.pdf"

    plot_statistical_inefficiency(
        block_widths_ns,
        statistical_inefficiencies,
        statistical_inefficiency_output_file,
    )
    plot_effective_samples_per_microsecond(
        block_widths_ns,
        effective_samples_per_microsecond,
        plateau_effective_samples,
        effective_samples_output_file,
    )

    print(f"Block widths (ns): {block_widths_ns.tolist()}")
    print(f"Block sizes (frames): {block_sizes.tolist()}")
    print(f"Statistical inefficiencies: {statistical_inefficiencies.tolist()}")
    print(f"Effective samples per μs: {effective_samples_per_microsecond.tolist()}")
    print(f"Automatically selected plateau start (ns): {plateau_start_ns:.6g}")
    print(f"Stable plateau detected: {plateau_detected}")
    print(f"Plateau relative tolerance: {args.plateau_relative_tolerance:.6g}")
    print(f"Plateau effective samples per μs: {plateau_effective_samples:.0f}")
    print(f"Trajectory duration (μs): {trajectory_time_microseconds:.6g}")
    print(f"Plateau effective sample size: {plateau_effective_sample_size:.0f}")

    if len(sequence_results) > 1:
        all_statistical_inefficiency_output_file = f"{args.output_path}/statistical_inefficiency_all_sequences.pdf"
        all_effective_samples_output_file = f"{args.output_path}/effective_samples_per_microsecond_all_sequences.pdf"
        plot_all_sequences_statistical_inefficiency(
            sequence_results,
            all_statistical_inefficiency_output_file,
        )
        plot_all_sequences_effective_samples_per_microsecond(
            sequence_results,
            all_effective_samples_output_file,
        )


if __name__ == "__main__":

    NANOSECONDS_PER_FRAME = 0.1
    MICROSECONDS_PER_FRAME = NANOSECONDS_PER_FRAME / 1000
    MINIMUM_BLOCK_WIDTH_NS = 1
    MAXIMUM_BLOCK_WIDTH_NS = 20
    BLOCK_WIDTH_STEP_NS = 1
    PLATEAU_RELATIVE_TOLERANCE = 0.1
    MINIMUM_PLATEAU_POINTS = 3
    PROTEIN_COLOR = {
        "F": "gray",
        "M": "deepskyblue",
        "V": "tomato",
        "L": "darkcyan",
        "A": "orchid",
        "Y": "limegreen",
        "I": "darkviolet",
    }

    parser = argparse.ArgumentParser(description="Calculate effective sample rates from Rg autocorrelation.")
    parser.add_argument("--input_path", required=True, help="Path to the input data.")
    parser.add_argument("--output_path", required=True, help="Path to the output data.")
    parser.add_argument("--sequence", default="M66", help="Sequence name.")
    parser.add_argument("--sequences", required=False, default="V66 M66 F66 A66 Y66 L66 I66", help="Sequence names.")
    parser.add_argument("--minimum_block_width_ns", type=float, default=MINIMUM_BLOCK_WIDTH_NS, help="Minimum block width in nanoseconds.")
    parser.add_argument("--maximum_block_width_ns", type=float, default=MAXIMUM_BLOCK_WIDTH_NS, help="Maximum block width in nanoseconds.")
    parser.add_argument("--block_width_step_ns", type=float, default=BLOCK_WIDTH_STEP_NS, help="Spacing between block widths in nanoseconds.")
    parser.add_argument("--nanoseconds_per_frame", type=float, default=NANOSECONDS_PER_FRAME, help="Simulation time per frame in nanoseconds.")
    parser.add_argument("--microseconds_per_frame", type=float, default=MICROSECONDS_PER_FRAME, help="Simulation time per frame in microseconds.")
    parser.add_argument("--plateau_relative_tolerance", type=float, default=PLATEAU_RELATIVE_TOLERANCE, help="Maximum relative deviation allowed within the automatically detected plateau.")
    parser.add_argument("--minimum_plateau_points", type=int, default=MINIMUM_PLATEAU_POINTS, help="Minimum number of points required in the automatically detected plateau.")
    args = parser.parse_args()
    args.sequences = args.sequences.split()
    if args.sequence not in args.sequences:
        args.sequences.insert(0, args.sequence)
    main(args)
