import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


"""
Plots cumulative Rg while preserving the times of frames removed for PIFs.

The Rg files contain only retained values.  The corresponding
``*_mindist_frames_to_delete.txt`` files identify the original frame numbers
that were removed before equilibration was discarded.
"""


SEQ_COLORS = {
    "F66": "gray",
    "M66": "deepskyblue",
    "V66": "tomato",
    "L66": "darkcyan",
    "A66": "orchid",
    "Y66": "limegreen",
    "I66": "darkviolet",
}


def read_rg(sequence, input_path):
    """
    Reads the retained, post-equilibration Rg values.

    Arguments:
        sequence : str
            The sequence name.
        input_path : str
            Path to the input data.

    Returns:
        rg_data : numpy.ndarray
            Array of retained Rg values.
    """
    filename = f"{input_path}/rg_{sequence}_wrapped_centered_PIF_equil_del.txt"
    return np.atleast_1d(np.loadtxt(filename, dtype=float))


def read_deleted_frames(sequence, deleted_frames_path):
    """
    Reads original frame numbers removed because of periodic images.

    Arguments:
        sequence : str
            The sequence name.
        deleted_frames_path : str
            Path to the deleted-frame data.

    Returns:
        deleted_frames : numpy.ndarray
            Array of deleted original frame numbers.
    """
    filename = f"{deleted_frames_path}/{sequence}_mindist_frames_to_delete.txt"
    try:
        deleted = np.atleast_1d(np.loadtxt(filename, dtype=int))
    except OSError as error:
        raise FileNotFoundError(f"Could not read deleted-frame file: {filename}") from error

    return np.unique(deleted)


def get_post_equilibration_frame_numbers(sequence, n_retained, deleted_frames_path, equilibration_frames):
    """
    Reconstructs original frame numbers for the retained Rg values.

    The processed Rg file has had both PIF frames and the initial equilibration
    frames removed.  Equilibration is applied after PIF removal, so the first
    ``equilibration_frames`` retained frames must be discarded.  It is not
    correct to discard all original frames with numbers below
    ``equilibration_frames`` because some of those original frames may already
    have been deleted as PIF frames.

    Arguments:
        sequence : str
            The sequence name.
        n_retained : int
            Number of Rg values in the processed data.
        deleted_frames_path : str
            Path to the deleted-frame data.
        equilibration_frames : int
            Number of retained frames removed for equilibration.

    Returns:
        retained_frames : numpy.ndarray
            Original frame numbers corresponding to the retained Rg values.
    """
    deleted_frames = read_deleted_frames(sequence, deleted_frames_path)

    # Reconstruct the full original trajectory, remove PIF frames, and then
    # remove equilibration frames from the resulting retained trajectory.
    n_original_frames = n_retained + len(deleted_frames) + equilibration_frames
    original_frames = np.arange(n_original_frames)
    retained_mask = ~np.isin(original_frames, deleted_frames)
    all_retained_frames = original_frames[retained_mask]
    retained_frames = all_retained_frames[equilibration_frames:]

    print(f"{sequence}:")
    print(f"  Rg values: {n_retained}")
    print(f"  Deleted frames: {len(deleted_frames)}")
    print(f"  Equilibration frames removed: {equilibration_frames}")
    print(f"  Reconstructed frames: {len(retained_frames)}")
    print(f"  First reconstructed frame: {retained_frames[0]}")
    print(f"  Last reconstructed frame: {retained_frames[-1]}")

    if len(retained_frames) != n_retained:
        raise ValueError(f"{sequence}: reconstructed {len(retained_frames)} retained frames, but the Rg file contains {n_retained} values. Check the frame spacing, equilibration time, and deleted-frame list.")

    return retained_frames


def cumulative_rg_plot_data(sequence, input_path, deleted_frames_path, ps_per_frame, equilibration_frames):
    """
    Calculates cumulative Rg values and their corrected times.

    Arguments:
        sequence : str
            The sequence name.
        input_path : str
            Path to the input data.
        deleted_frames_path : str
            Path to the deleted-frame data.
        ps_per_frame : float
            Time between original frames in picoseconds.
        equilibration_frames : int
            Number of retained frames removed for equilibration.

    Returns:
        time_in_ns : numpy.ndarray
            Corrected time values in nanoseconds.
        cumulative_rg : numpy.ndarray
            Array of cumulative average Rg values.
    """
    rg_data = read_rg(sequence, input_path)

    frame_numbers = get_post_equilibration_frame_numbers(sequence, len(rg_data), deleted_frames_path, equilibration_frames)

    cumulative_rg = np.cumsum(rg_data) / np.arange(1, len(rg_data) + 1)

    # Plot time since the beginning of the retained/equilibrated trajectory.
    deleted_frames = read_deleted_frames(sequence, deleted_frames_path)
    n_original_frames = len(rg_data) + len(deleted_frames) + equilibration_frames
    plot_frames = np.arange(equilibration_frames, n_original_frames)
    plot_cumulative_rg = np.full(len(plot_frames), np.nan)
    plot_cumulative_rg[frame_numbers - equilibration_frames] = cumulative_rg
    first_valid = np.flatnonzero(~np.isnan(plot_cumulative_rg))[0]
    for index in range(first_valid + 1, len(plot_cumulative_rg)):
        if np.isnan(plot_cumulative_rg[index]):
            plot_cumulative_rg[index] = plot_cumulative_rg[index - 1]
    time_in_ns = (plot_frames - equilibration_frames) * ps_per_frame / 1000.0

    return time_in_ns, plot_cumulative_rg


def make_plot(sequences, input_path, deleted_frames_path, output_path, ps_per_frame, equilibration_frames):
    """
    Generates a step plot of cumulative Rg over corrected time values.

    Arguments:
        sequences : list
            Sequence names to plot.
        input_path : str
            Path to the input data.
        deleted_frames_path : str
            Path to the deleted-frame data.
        output_path : str
            Path for the output PDF.
        ps_per_frame : float
            Time between original frames in picoseconds.
        equilibration_frames : int
            Number of retained frames removed for equilibration.

    Returns:
        None
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    for sequence in sequences:
        time_in_ns, cumulative_rg = cumulative_rg_plot_data(sequence, input_path, deleted_frames_path, ps_per_frame, equilibration_frames)

        # A step plot keeps the cumulative average constant during a missing
        # frame and updates it when the next retained frame is reached.
        ax.step(time_in_ns, cumulative_rg, where="post", color=SEQ_COLORS.get(sequence, "black"), label=sequence, linewidth=2)

    ax.set_ylim(14.5, 22)
    # ax.set_xlim(-20, 1200)
    ax.set_xlabel("t (ns)", fontsize=20)
    ax.set_ylabel(r"$R_g$ ($\AA$)", fontsize=20)
    ax.tick_params(axis="both", labelsize=18)
    ax.legend(fontsize=15, title_fontsize=15, loc="upper right", ncol=2)

    output_filename = f"{output_path}/cumulative_rg_all_sequences.pdf"
    plt.savefig(output_filename, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Plot cumulative Rg using original frame times after PIF removal.")
    
    parser.add_argument("--input_path", required=True, help="Directory containing Rg files.")
    parser.add_argument("--deleted_frames_path", required=True, help="Directory containing *_mindist_frames_to_delete.txt files.")
    parser.add_argument("--output_path", required=True, help="Directory for the PDF output.")
    parser.add_argument("--sequences", default="F66 M66 V66 L66 A66 Y66 I66", help="Space-separated sequence names.")
    parser.add_argument("--ps_per_frame", type=float, default=100.0, help="Time between original frames in ps (default: 100).")
    parser.add_argument("--equilibration_frames", type=int, default=8000, help="Number of retained frames removed for equilibration (default: 8000).")

    args = parser.parse_args()

    # 

    make_plot(args.sequences.split(), args.input_path, args.deleted_frames_path, args.output_path, args.ps_per_frame, args.equilibration_frames)
