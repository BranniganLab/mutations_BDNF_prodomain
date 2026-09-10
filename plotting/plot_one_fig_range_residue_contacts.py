import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse

"""
This script plots the contact frequency between a selection of residues in a protein. 
The data used for the contact frequency was generated from "calc_sidechain_contacts.tcl"
"""

def read_file(sequence):
    """
    Reads residue-pair side chain contact frequencies from a text file.
    
    Arguments:
        sequence : str
            The sequence name.
    
    Returns:
        contact_freq_data : pandas.DataFrame
            Array of contact frequencies.
    """
    contact_freq_data = pd.read_csv(f'{input_path}/{res_type}_contactfreqs_6A_{sequence}_{state}.txt', sep=r'\s+', header=None, names=["resid1", "resid2", "freq"])

    return contact_freq_data

def keep_only_desired_residue_range(sequence):
    """
    Filters the contact frequency data to keep only the rows where the first column is within [x_axis_beg_range, x_axis_end_range] and the second column is within [y_axis_beg_range, y_axis_end_range].
    
    Arguments:
        sequence : str
            The sequence name.
    
    Returns:
        contact_range : pandas.DataFrame
            DataFrame containing only rows where the first column is within [x_axis_beg_range, x_axis_end_range] and the second column is within [y_axis_beg_range, y_axis_end_range].
    """
    contact_freq_data = read_file(sequence)

    if x_axis_beg_range < y_axis_beg_range:
        contact_range = contact_freq_data[
            contact_freq_data["resid1"].between(x_axis_beg_range, x_axis_end_range) & 
            contact_freq_data["resid2"].between(y_axis_beg_range, y_axis_end_range)].copy()
    else:
        contact_range = contact_freq_data[
            contact_freq_data["resid1"].between(y_axis_beg_range, y_axis_end_range) & 
            contact_freq_data["resid2"].between(x_axis_beg_range, x_axis_end_range)].copy()

        contact_range[["resid1", "resid2"]] = contact_range[["resid2", "resid1"]]

    return contact_range.reset_index(drop=True)    

def create_heatmap_matrix(sequence):
    """
    Converts residue pair contact frequencies into a matrix for heatmap plotting.
    
    Arguments:
        sequence : str
            The sequence name.
    
    Returns:
        contact_freq_matrix : pandas.DataFrame
            DataFrame containing the contact frequencies in a matrix format, with rows corresponding to y-axis residues and columns corresponding to x-axis residues.
    """

    contact_freq_data = keep_only_desired_residue_range(sequence)

    # Define residue ranges
    residue_range_x = list(range(x_axis_beg_range, x_axis_end_range + 1))
    residue_range_y = list(range(y_axis_beg_range, y_axis_end_range + 1))

    # Map residues to matrix indices
    x_label_to_index = {label: i for i, label in enumerate(residue_range_x)}
    y_label_to_index = {label: i for i, label in enumerate(residue_range_y)}

    # Initialize matrix (rows = y-axis, cols = x-axis)
    contact_freq_for_plot = np.zeros((len(residue_range_y), len(residue_range_x)))

    # Fill matrix
    for _, row in contact_freq_data.iterrows():
        if row["resid1"] in x_label_to_index and row["resid2"] in y_label_to_index:
            x_idx = x_label_to_index[row["resid1"]]
            y_idx = y_label_to_index[row["resid2"]]
            contact_freq_for_plot[y_idx, x_idx] = row["freq"]

    # Convert to DataFrame for easier plotting / saving
    contact_freq_matrix = pd.DataFrame(contact_freq_for_plot, index=residue_range_y, columns=residue_range_x)
    # contact_freq_matrix.to_csv("matrix.csv")

    return contact_freq_matrix

def plot_heatmap(ax, sequence, fig): 
    """
    Generates a plot of the contact frequencies on a heatmap.
    
    Arguments:
        sequence : str
            The sequence name.
        ax : matplotlib.axes.Axes
            The axes on which to plot the heatmap.
        fig : matplotlib.figure.Figure
            The figure to which the axes belong.
    
    Returns:
        None
    """
    
    contact_freq_matrix = create_heatmap_matrix(sequence)
    
    # Generate the heatmap
    imgp = ax.imshow(contact_freq_matrix, origin='lower', aspect='equal', cmap='Purples', vmin=0, vmax=12)
    
    # Title for the plot
    ax.set_title(f"V66", loc='center', fontsize='30')

    # Protein sequence with mutation character preserved
    protein_seq_w_mutation = [sequence if aa == '*' else aa for aa in protein_seq]

    # Residue numbering
    residue_numbers = np.arange(start_residue, start_residue + len(protein_seq_w_mutation))

    # X-axis labels
    x_range = protein_seq_w_mutation[x_axis_beg_range - start_residue:x_axis_end_range - start_residue + 1]
    x_numbers = residue_numbers[x_axis_beg_range - start_residue:x_axis_end_range - start_residue + 1]

    x_labels = [f"{aa}{num}" for aa, num in zip(x_range, x_numbers)]
    ax.set_xticks(np.arange(len(x_range)))
    ax.set_xticklabels( x_labels, fontsize=12, ha='center', va='top', rotation=90)

    # Color the x tick labels
    for i, label in enumerate(ax.get_xticklabels()):
        aa = x_range[i]
        if aa in acidic:
            label.set_color(acidic[aa])
        elif aa in basic:
            label.set_color(basic[aa])
        elif aa in special:
            label.set_color(special[aa])

    # Y-axis labels
    y_range = protein_seq_w_mutation[y_axis_beg_range - start_residue:y_axis_end_range - start_residue + 1]
    y_numbers = residue_numbers[y_axis_beg_range - start_residue:y_axis_end_range - start_residue + 1]

    y_labels = [f"{aa}{num}" for aa, num in zip(y_range, y_numbers)]
    ax.set_yticks(np.arange(len(y_range)))
    ax.set_yticklabels(y_labels, fontsize=12, ha='right', va='center')

    # Color the y tick labels
    for i, label in enumerate(ax.get_yticklabels()):
        aa = y_range[i]
        if aa in acidic:
            label.set_color(acidic[aa])
        elif aa in basic:
            label.set_color(basic[aa])
        elif aa in special:
            label.set_color(special[aa])
    
    # Add axis labels and title
    ax.set_title(f"{sequence}66", loc='center', fontsize='30')

    fig.tight_layout()
    
    # Adjust the right side of the subplots to make space for the colorbar
    fig.subplots_adjust(right=0.91)  

    # Draw the colorbar right next to the subplot
    if sequence == "L" or sequence == "I":
        box = ax.get_position()  # [x0, y0, width, height]
        cbar_width = 0.02
        cbar_pad = 0.01
        cbar_x = box.x1 + cbar_pad

        if sequence == 'L_hid_11' or sequence == "L":
            cbar_ax = fig.add_axes([cbar_x, box.y0, cbar_width, box.height])
        elif sequence == 'I_hid_11' or sequence == "I":
            cbar_ax = fig.add_axes([cbar_x, box.y0, cbar_width, box.height])

        cbar = fig.colorbar(imgp, cax=cbar_ax)
        cbar.ax.tick_params(labelsize=10)
        cbar.set_label("Contact frequency (%)", fontsize=15, rotation=270, labelpad=20)
    
def plot_combined():
    """
    Generates a combined figure of the contact frequencies on heat maps.
    
    Arguments:
        None
    
    Returns:
        None
    """
    fig, axes = plt.subplots(2, 4, figsize=(10, 8))

    # TOP ROW
    plot_heatmap(axes[0, 0], "F", fig)
    plot_heatmap(axes[0, 1], "M", fig)
    plot_heatmap(axes[0, 2], "V", fig)
    plot_heatmap(axes[0, 3], "L", fig)
    # axes[0, 3].axis('off')

    # BOTTOM ROW (Beta)
    plot_heatmap(axes[1, 0], "A", fig)
    plot_heatmap(axes[1, 1], "Y", fig)
    plot_heatmap(axes[1, 2], "I", fig)
    axes[1, 3].axis('off')

    plt.savefig(f"{output_path}/{output_filename}")
    plt.close()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Plot side-chain contact frequency maps.")
    
    parser.add_argument("--state", required=True, help="Contact state name used in the contact-frequency filename.")
    parser.add_argument("--x-range", nargs=2, type=int, required=True, metavar=("START", "END"), help="Residue range for the x-axis.")
    parser.add_argument("--y-range", nargs=2, type=int, required=True, metavar=("START", "END"), help="Residue range for the y-axis.")
    parser.add_argument("--res_type", required=True, help="Residue type (e.g., 'bb' for backbone or 'sc' for side chain) used in the contact-frequency filename.")
    parser.add_argument("--input_path", required=True, help="Path to the input data.") 
    parser.add_argument("--output_path", required=True, help="Path to the output data.")
    parser.add_argument("--start_residue", required=False, default="23", help="First residue of the protein.")
    
    args = parser.parse_args()

    state = args.state
    x_axis_beg_range, x_axis_end_range = args.x_range
    y_axis_beg_range, y_axis_end_range = args.y_range
    res_type = args.res_type
    input_path = args.input_path
    output_path = args.output_path
    start_residue = int(args.start_residue)

    acidic = {
        'D': 'red',
        'E': 'red'
    }

    basic = {
        'R': 'blue',
        'H': 'blue',
        'K': 'blue'
    }

    special = {
        'M': 'gold',
        'F': 'gold',
        'Y': 'gold',
        'W': 'gold'
    }

    protein_seq = (
        'EANIRGQGGLAYPGVRTHGTLESVNGPKAGSRGLTSLADTFEH*IEELLDEDQKVRPNEENNKDADLYTSRVMLSSQVPLEPPLLFLLEEY'
    )

    residue_range_str_x = f'{x_axis_beg_range}-{x_axis_end_range}'
    residue_range_str_y = f'{y_axis_beg_range}-{y_axis_end_range}'
    residue_range_str = f'{residue_range_str_x}_vs_{residue_range_str_y}'

    output_filename = (f'{res_type}_{residue_range_str}_contactfreqs_6A_{state}_all_seqs.pdf')

    plot_combined()