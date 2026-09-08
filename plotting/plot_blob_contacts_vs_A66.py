import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
import linecache
import argparse

"""
This script plots the difference of blob contact frequency between the sequence versus the A66 reference sequence.
The data used for the contact frequency was generated from "calc_midpoints_of_all_blobs.tcl" and "calc_rg_of_all_blobs.tcl"
"""

def read_center_file(sequence):
    """
    Reads the XYZ coordinates associated with the center of a blob. All coordinates of all blobs are stored in this file.
    Each row in the text file corresponds to a frame number and each column represents x, y, or z coordinates. 
    Columns read as: column 1:x data, column 2:y data, column 3:z data (of the first blob of the sequence),
                     column 4:x data, column 5:y data, column 6:z data (of the second blob of the sequence)
                     column 7:x data, column 8:y data, column 9:z data (of the third blob of the sequence and so on)
    Arguments:
        sequence : str
            The sequence name.
    
    Returns:
        blob_center_coords : numpy.ndarray
            Array of blob center coordinates.
    """
    center_coords_filename = f"midpoint_{sequence}_wrapped_centered_PIF_equil_del.txt" 
    blob_center_coords = np.loadtxt(f'{input_path}/{center_coords_filename}')

    return blob_center_coords

def calc_excess_distance_and_contacts_for_specific_blob_groups(sequence):
    """
    Calculates the excess distance between two blobs and if the distance between the two blobs is less than the cutoff distance, a counter will update the contact_probability array.
    
    Arguments:
        sequence : str
            The sequence name.    
    
    Returns:
        contact_probability : numpy.ndarray
            Array of contact probabilities.
    """      
    rgyr_filename = f"{input_path}/blob_Rg_{sequence}_wrapped_centered_PIF_equil_del.txt"
    blob_center_coords = read_center_file(sequence)
    
    number_of_frames = blob_center_coords.shape[0]
    contact_probability = np.zeros((number_of_frames, len(GROUP_A), len(GROUP_B)))

    for i, blob_i in enumerate(GROUP_A):
        for j, blob_j in enumerate(GROUP_B):

            # ---- RGYR (use REAL blob IDs) ----
            rgyr_i = np.array([float(v) for v in linecache.getline(f"{rgyr_filename}", blob_i + 1).split()])
            rgyr_j = np.array([float(v) for v in linecache.getline(f"{rgyr_filename}", blob_j + 1).split()])

            # ---- coordinate indices (use REAL blob IDs) ----
            i_x, i_y, i_z = blob_i*3, blob_i*3+1, blob_i*3+2
            j_x, j_y, j_z = blob_j*3, blob_j*3+1, blob_j*3+2

            dist = np.sqrt(
                (blob_center_coords[:, i_x] - blob_center_coords[:, j_x])**2 +
                (blob_center_coords[:, i_y] - blob_center_coords[:, j_y])**2 +
                (blob_center_coords[:, i_z] - blob_center_coords[:, j_z])**2)

            excess_distance = dist - (rgyr_i + rgyr_j)

            contact_probability[:, i, j] = (excess_distance < cutoff_distance).astype(int)

    return contact_probability

def calc_contact_frequency(sequence):  
    """
    Calculates the average of the contact probability array and expresses it as a percentage by mulitplying by 100.

    Arguments:
        sequence : str
            The sequence name.

    Returns:
        contact_frequency : numpy.ndarray
            Array of contact frequencies.        
    """ 
    contact_probability = calc_excess_distance_and_contacts_for_specific_blob_groups(sequence)
    contact_frequency = np.mean(contact_probability, axis=0) * 100 # Averages the values across each time step for each pair of blobs and expresses it as a percentage

    return contact_frequency

def calc_contact_difference_between_variants(sequence1, sequence2):
    """
    Calculates the difference between the contact frequencies of two sequences.

    Arguments:
        sequence1 : str
            The sequence name.
        sequence2 : str
            The sequence name of the reference.

    Returns:
        contact_difference : numpy.ndarray
            Array of contact differences.
    """
    contact_frequency1 = calc_contact_frequency(sequence1)
    contact_frequency2 = calc_contact_frequency(sequence2)

    contact_difference = contact_frequency1 - contact_frequency2

    return contact_difference

def fig_gen_difference_between_variants(ax, sequence1, sequence2, fig):
    """
    Generates a plot of the contact frequencies on a heat map.

    Arguments:
        ax : matplotlib.axes._subplots.AxesSubplot
            The subplot.
        sequence1 : str
            The sequence name.
        sequence2 : str
            The sequence name of the reference.
        fig : matplotlib.figure.Figure
            The figure.

    Returns:
        None
    """
    probability_mean = calc_contact_difference_between_variants(sequence1, sequence2)
    imgp = ax.imshow(probability_mean, origin='lower', aspect='equal', cmap='PRGn', vmin=-40, vmax=40)
    ax.set_title(fr'{sequence1}-A66', loc='center', fontsize='20')

    # Draw the colorbar right next to the subplot
    if sequence1 in ["V66", "I66"]:
        box = ax.get_position()
        cbar_width = 0.02
        cbar_pad = 0.01
        cbar_x = box.x1 + cbar_pad
        cbar_ax = fig.add_axes([cbar_x, box.y0, cbar_width, box.height])

        cbar = fig.colorbar(imgp, cax=cbar_ax)
        cbar.ax.tick_params(labelsize=10)
        cbar.set_label("Contact frequency (%)", fontsize=15, rotation=270, labelpad=20)

    ax.xaxis.set_major_locator(plticker.MultipleLocator(1))
    ax.yaxis.set_major_locator(plticker.MultipleLocator(1))
    ax.set_xticklabels([])
    ax.set_yticklabels([])

def combine_blob_contact_plots_difference():
    """
    Generates a plot that contains blob contact frequencies for all sequences

    Arguments:
        None

    Returns:
        None
    """
    fig, axes = plt.subplots(2, 3, figsize=(10, 10))

    # TOP ROW
    fig_gen_difference_between_variants(axes[0, 0], "F66", "A66", fig)
    fig_gen_difference_between_variants(axes[0, 1], "M66", "A66", fig)
    fig_gen_difference_between_variants(axes[0, 2], "V66", "A66", fig)

    # BOTTOM ROW (Beta)
    fig_gen_difference_between_variants(axes[1, 0], "L66", "A66", fig)
    fig_gen_difference_between_variants(axes[1, 1], "Y66", "A66", fig)
    fig_gen_difference_between_variants(axes[1, 2], "I66", "A66", fig)

    fig.subplots_adjust(wspace=0)
    plt.savefig(f"{output_path}/blob_contacts_A66_vs_sequences.pdf", bbox_inches='tight')
    plt.close()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Plot the differences in contact frequencies between a given sequence and the reference.")

    parser.add_argument("--input_path", required=True, help="Path to the input data.") 
    parser.add_argument("--output_path", required=True, help="Path to the output data.")
    parser.add_argument("--n_independent_samples", required=False, default=1088, help="Number of samples (1088).")
    parser.add_argument("--cutoff_distance", required=False, default=0.55, help="Cutoff distance.")
    parser.add_argument("--GROUP_A", nargs="+", required=False, default=[0, 1, 2, 3, 4, 5], help="Blob indices for the y axis.")
    parser.add_argument("--GROUP_B", nargs="+", required=False, default=[7, 8, 9, 10], help="Blob indices for the x axis.")

    args = parser.parse_args()

    input_path = args.input_path
    output_path = args.output_path
    n_independent_samples = int(args.n_independent_samples)
    cutoff_distance = float(args.cutoff_distance)
    GROUP_A = args.GROUP_A
    GROUP_B = args.GROUP_B

    combine_blob_contact_plots_difference()