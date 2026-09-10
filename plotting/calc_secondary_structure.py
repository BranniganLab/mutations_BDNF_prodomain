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
def calc_ss(group_number=1):
    """
    Detects specific patterns of hydrogen bonds between amino acid residues to determine the secondary structure of a protein.

    Arguments:
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

    # calc_ss()

    dat_files = [f for f in os.listdir(data_path) if f.endswith(".dat")]

    # Create a plot with 2 subplots (one for each secondary structure)
    fig, axs = plt.subplots(nrows=2, figsize=(20, 15), sharex=True)
    structure_labels = ['Helices', 'Beta sheets']

    for file_name in dat_files:
        file_path = os.path.join(data_path, file_name)

        # Read and process the DSSP file
        with open(file_path, "r") as file:
            lines = file.readlines()

        data = [[char for char in line.strip()] for line in lines]
        df = pd.DataFrame.from_records(data)

        # Replace DSSP labels with the categories
        df_reduced = df.replace(ss_dict)

        # Calculate position variation and normalize by the total frames
        total_frames = len(df_reduced)
        position_variation = (df_reduced.apply(lambda x: x.value_counts()).fillna(0) / total_frames * 100)

        # Drop first and last residues (columns)
        position_variation = position_variation.iloc[:, 1:-1]

        # Plot data for each secondary structure label
        for ax, structure in zip(axs, structure_labels):
            if structure in position_variation.index:
                mutation = f"{file_name[5]}"
                print(file_name)
                print(f"{file_name[5]}")
                ax.plot(position_variation.loc[structure, :], label=file_name[5]+"66", color=protein_color[mutation], linewidth=4)
                blob_boundaries = [0, 1, 8, 16, 19, 25, 32, 41, 43, 50, 65, 69, 70, 75, 76, 81, 82, 89, 91, len(position_variation.columns)]

                # Set specific y-axis limits and blob boundary vline height
                if structure == "Helices":
                    ax.set_ylim(0, ymax_helix)
                    for x in blob_boundaries[1:-1]: 
                        ax.vlines(x=x, ymin=0, ymax=ymax_helix, color="grey", lw=0.5)
                elif structure == "Beta sheets":
                    ax.set_ylim(0, ymax_beta)
                    for x in blob_boundaries[1:-1]: 
                        ax.vlines(x=x, ymin=0, ymax=ymax_beta, color="grey", lw=0.5)

                ax.set_ylabel(f"{structure} (%)", fontsize=30)
                ax.legend(loc='upper right', title="Variant", fontsize=20, title_fontsize=20)

                # Set x-axis ticks and labels only for the kept residues
                residue_positions = list(position_variation.columns)  # these are column indices after slicing
                residue_labels = [f"{int(i)+22}" for i in residue_positions]  # adjust to real residue numbers

                # Optional: select every 4th residue for cleaner labeling
                tick_indices = list(range(0, len(residue_positions), 4))
                ax.set_xticks([residue_positions[i] for i in tick_indices])
                ax.set_xticklabels([residue_labels[i] for i in tick_indices], fontsize=25)

    ax.set_xlabel("Residue", fontsize=30)  # Set the x-axis label
    
    for ax in axs:
        ax.tick_params(axis='y', labelsize=25)

    plt.tight_layout()

    # Save the plot
    output_file = os.path.join(output_path, fig_name)
    plt.savefig(output_file)
    print(f"Plot saved to {output_file}")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Plots helices and beta sheets at residue level.")

    parser.add_argument("--raw_traj_path", required=True, help="Path to the raw trajectory.")
    parser.add_argument("--input_path", required=True, help="Path to the data.")
    parser.add_argument("--output_path", required=True, help="Path to the figures.")
    parser.add_argument("--state", required=True, help="Contact state name of the protein.")
    parser.add_argument("--sequences", required=False, default="F M V L A Y I", help="Sequence names.")

    args = parser.parse_args()

    raw_traj_path = args.raw_traj_path
    # processed_traj_path = args.processed_traj_path
    data_path = args.input_path
    output_path = args.output_path
    state = args.state
    sequences = args.sequences.split()

    fig_name=(f"{output_path}/ss_{state}.pdf")

    ymax_helix = 75
    ymax_beta = 75

    ss_dict = {
        'H': 'Helices', 'G': 'Helices', 'I': 'Helices', 'P': 'Helices',
        'B': 'Beta sheets', 'E': 'Beta sheets',
        'T': 'Coils and turns', 'S': 'Coils and turns', '~': 'Coils and turns'}

    protein_color = {'F': 'gray', 'M': 'deepskyblue', 'V': 'tomato', 'L': 'darkcyan', 'A': 'orchid', 'Y': 'limegreen', 'I': 'darkviolet'}


    for seq in sequences:
        input_trajectory = (f"{data_path}/{seq}_{state}.gro")
        topology = (f"{raw_traj_path}/{seq}66/{seq}66.tpr")
        output_dssp = (f"{data_path}/dssp_{seq}_{state}.dat")
        stats = (f"{data_path}/stats_{seq}_{state}.xvg")


    plot_helix_beta()


