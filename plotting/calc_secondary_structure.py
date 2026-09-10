import subprocess
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from matplotlib.gridspec import GridSpec
import argparse


"""
    This script plots helical and beta sheet secondary structure at the residue level for all sequences using GROMACS dssp
"""
def calc_ss(input_trajectory, topology, output_dssp, stats, group_number=1):
    """
    Detects specific patterns of hydrogen bonds between amino acid residues to determine the secondary structure of a protein.

    Arguments:
        input_trajectory : Path
            Input trajectory file.

        topology : Path
            GROMACS .tpr file.

        output_dssp : Path
            Output dssp file.

        stats : Path
            Statistics file.

        group_number : int, optional
            GROMACS index number corresponding to the Protein group.
            Default is 1.

    Returns:
        None
    """

    command = [
        "gmx",
        "dssp",
        "-f", str(input_trajectory),
        "-s", str(topology),
        "-o", str(output_dssp),
        "-num", str(stats),
    ]

    selections = f"{group_number}\n{group_number}\n"

    subprocess.run(command, input=selections, text=True, check=True)

def plot_helix_beta():
    """
    Calculates and plots the percentage of residues assigned to helices
    and beta sheets for each sequence.

    Arguments:
        None

    Returns:
        None
    """

    structure_labels = ["Helices", "Beta sheets"]

    fig, axs = plt.subplots(nrows=2,figsize=(20, 15),sharex=True)

    for seq in sequences:
        input_trajectory = os.path.join(data_path,f"{seq}_{state}.gro")
        topology = os.path.join(raw_traj_path, f"{seq}66", f"{seq}66.tpr")
        output_dssp = os.path.join(data_path, f"dssp_{seq}_{state}.dat")
        stats = os.path.join(data_path, f"stats_{seq}_{state}.xvg" )

        calc_ss(input_trajectory, topology, output_dssp, stats)

        with open(output_dssp, "r") as file:
            lines = [
                line.strip()
                for line in file
                if line.strip()
            ]

        data = [list(line) for line in lines]
        df = pd.DataFrame.from_records(data)
        df_reduced = df.replace(ss_dict)
        total_frames = len(df_reduced)
        position_variation = (df_reduced.apply(lambda column: column.value_counts()).fillna(0) / total_frames * 100)

        # Remove the first and last DSSP columns.
        position_variation = position_variation.iloc[:, 1:-1]
        residue_positions = list(position_variation.columns)
        residue_labels = [int(position) + start_residue - 1 for position in residue_positions]
        tick_indices = list(range(0, len(residue_positions), tick_step))

        for ax, structure in zip(axs, structure_labels):
            ax.plot(
                position_variation.loc[structure, :],
                label=f"{seq}66",
                color=protein_color[seq],
                linewidth=4
            )

            if structure == "Helices":
                ax.set_ylim(0, ymax_helix)
            else:
                ax.set_ylim(0, ymax_beta)

            for boundary in blob_boundaries[1:-1]:
                if boundary < len(residue_positions):
                    ax.axvline(
                        x=boundary,
                        color="grey",
                        linewidth=0.5
                    )

            ax.set_ylabel(f"{structure} (%)", fontsize=30)
            ax.set_xticks([residue_positions[i] for i in tick_indices])
            ax.set_xticklabels([residue_labels[i] for i in tick_indices],fontsize=25)

    axs[0].legend( loc="upper right", fontsize=20, title_fontsize=20)
    axs[1].legend(loc="upper right", fontsize=20, title_fontsize=20)
    axs[-1].set_xlabel("Residue", fontsize=30)

    for ax in axs:
        ax.tick_params(axis="y", labelsize=25)

    fig.tight_layout()
    output_file = os.path.join(output_path, fig_name)
    fig.savefig(output_file, bbox_inches="tight")
    plt.close(fig)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Plots helices and beta sheets at residue level.")

    parser.add_argument("--raw_traj_path", required=True, help="Path to the raw trajectory.")
    parser.add_argument("--input_path", required=True, help="Path to the data.")
    parser.add_argument("--output_path", required=True, help="Path to the figures.")
    parser.add_argument("--state", required=True, help="Contact state name of the protein.")
    parser.add_argument("--sequences", required=False, default="F M V L A Y I", help="Sequence names.")

    args = parser.parse_args()

    raw_traj_path = args.raw_traj_path
    data_path = args.input_path
    output_path = args.output_path
    state = args.state
    sequences = args.sequences.split()

    fig_name=(f"ss_{state}.pdf")

    ymax_helix = 75
    ymax_beta = 75
    start_residue = 23

    ss_dict = {
        'H': 'Helices', 'G': 'Helices', 'I': 'Helices', 'P': 'Helices',
        'B': 'Beta sheets', 'E': 'Beta sheets',
        'T': 'Coils and turns', 'S': 'Coils and turns', '~': 'Coils and turns'}

    protein_color = {'F': 'gray', 
                     'M': 'deepskyblue', 
                     'V': 'tomato', 
                     'L': 'darkcyan', 
                     'A': 'orchid', 
                     'Y': 'limegreen', 
                     'I': 'darkviolet'}


    blob_boundaries = [
        0, 1, 8, 16, 19, 25, 32, 41, 43, 50,
        65, 69, 70, 75, 76, 81, 82, 89, 91
    ]

    tick_step = 4

    plot_helix_beta()


