import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
import linecache
import argparse

"""
This script plots the contact frequency between all blobs of a protein.
The data used for the contact frequency was generated from "calc_midpoints_of_all_blobs.tcl" or "calc_geometric_centers_of_all_blobs.tcl" and "calc_rg_of_all_blobs.tcl"
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

def get_num_of_blobs(sequence):
    """
    Counts the number of lines in the Rg file which is equivalent to the number of blobs in the protein.

    Arguments:
        sequence : str
            The sequence name.
    
    Returns:
        num_of_blobs : int
            The number of blobs in the protein.
        rgyr_filename : str
            The filename of the Rg file.
    """
    rgyr_filename = f"blob_Rg_{sequence}_wrapped_centered_PIF_equil_del.txt"
    with open(f"{input_path}/{rgyr_filename}", 'r') as fp:
        num_of_blobs = len(fp.readlines()) # The number of blobs is equal to the number of lines in the rgyr file

    return num_of_blobs, rgyr_filename

def calc_excess_distance_and_contacts(sequence):
    """
    Calculates the excess distance between two blobs and if the distance between the two blobs is less than the cutoff distance, a counter will update the contact_probability array.

    Arguments:
        sequence : str
            The sequence name.
    
    Returns:
        contact_probability : numpy.ndarray
            Array of contact probabilities.
    """    
    blob_center_coords = read_center_file(sequence) 
    num_of_blobs, rgyr_filename = get_num_of_blobs(sequence)

    number_of_frames = blob_center_coords.shape[0]
    contact_probability = np.zeros((number_of_frames, num_of_blobs, num_of_blobs)) # Constructs an array for the contact probability that is the length of the time frame x length of order_of_blob x length order_of_blob

    for blob_i in range(num_of_blobs):
        blob_i_x = blob_i * 3 # column position for x coords of blob i 
        blob_i_y = blob_i * 3 + 1 # column position for y coords of blob i
        blob_i_z = blob_i * 3 + 2 # column position for z coords of blob i

        rgyr_values_blob_i = linecache.getline(f"{input_path}/{rgyr_filename}", blob_i + 1).strip()
        list_of_rgyr_i = np.array([float(value) for value in rgyr_values_blob_i.split()])  # Convert to NumPy array

        for blob_j in range(num_of_blobs):
            blob_j_x = blob_j * 3 # column position for x coords of blob j
            blob_j_y = blob_j * 3 + 1 # column position for y coords of blob j
            blob_j_z = blob_j * 3 + 2 # column position for z coords of blob j

            rgyr_values_blob_j = linecache.getline(f"{input_path}/{rgyr_filename}", blob_j + 1).strip()
            list_of_rgyr_j = np.array([float(value) for value in rgyr_values_blob_j.split()])  # Convert to NumPy array

            excess_distance = ((blob_center_coords[:, blob_i_x] - blob_center_coords[:, blob_j_x])**2 +
                                (blob_center_coords[:, blob_i_y] - blob_center_coords[:, blob_j_y])**2 +
                                (blob_center_coords[:, blob_i_z] - blob_center_coords[:, blob_j_z])**2)**0.5 - (list_of_rgyr_i + list_of_rgyr_j)  

            frame = 0
            for excess_distance_value in excess_distance: # Check if there are contacts by the excess distance value
                if excess_distance_value < cutoff_distance: # A blob-blob contact is made if the excess distance is less than the user defined cufoff distance.
                    counter =  1
                else:
                    counter = 0

                contact_probability[frame, blob_i, blob_j] = counter # Contact probability array will be updated with a 1 if there is a contact or 0 if there is no contact at the time (frame) and blob pair position in the array
                frame = frame + 1

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
    contact_probability = calc_excess_distance_and_contacts(sequence)

    contact_frequency = np.mean(contact_probability, axis=0) * 100 # Averages the values across each time step for each pair of blobs and expresses it as a percentage

    return contact_frequency

def fig_gen(ax, sequence, fig):
    """
    Generates a figure of the contact frequencies on a heat map.
    
    Arguments: 
        ax : matplotlib.axes._subplots.AxesSubplot
            The subplot.
        sequence : str
            The sequence name. 
        fig : matplotlib.figure.Figure
            The figure.
    
    Returns:
        None
    """ 
    contact_frequency = calc_contact_frequency(sequence)
    imgp = ax.imshow(contact_frequency, origin='lower', aspect='equal', cmap='Purples', vmin=0, vmax=50)
    ax.set_title(f'{sequence}', loc='center', fontsize='20')

    # Colorbar properties
    if sequence in ["L66", "I66"]:
        box = ax.get_position()
        cbar_width = 0.008
        cbar_pad = 0.01
        cbar_x = box.x1 + cbar_pad
        cbar_ax = fig.add_axes([cbar_x, box.y0, cbar_width, box.height])

        cbar = fig.colorbar(imgp, cax=cbar_ax)
        cbar.ax.tick_params(labelsize=10)
        cbar.set_label("Contact frequency (%)", fontsize=15, rotation=270, labelpad=20)

    # Section the plot
    ax.vlines(x=6-0.5, ymin=0-0.5, ymax=6-0.5, color="black", lw=4)
    ax.vlines(x=7-0.5, ymin=7-0.5, ymax=11-0.5, color="black", lw=4)
    ax.hlines(y=6-0.5, xmin=0-0.5, xmax=6-0.5, color="black", lw=4)
    ax.hlines(y=7-0.5, xmin=7-0.5, xmax=11-0.5, color="black", lw=4)
    ax.vlines(x=6-0.5, ymin=7-0.5, ymax=11-0.5, color="black", lw=4)
    ax.vlines(x=7-0.5, ymin=0-0.5, ymax=6-0.5, color="black", lw=4)
    ax.hlines(y=6-0.5, xmin=7-0.5, xmax=11-0.5, color="black", lw=4)
    ax.hlines(y=7-0.5, xmin=0-0.5, xmax=6-0.5, color="black", lw=4)

    ax.xaxis.set_major_locator(plticker.MultipleLocator(1))
    ax.yaxis.set_major_locator(plticker.MultipleLocator(1))
    ax.set_xticklabels([])
    ax.set_yticklabels([])

def combine_blob_contact_plots():
    """
    Generates one figure that displays blob contact frequencies for all sequences.
    
    Arguments:
        None
    
    Returns:
        None
    """     
    fig, axes = plt.subplots(2, 4, figsize=(20, 10))

    # TOP ROW
    fig_gen(axes[0, 0], "F66", fig)
    fig_gen(axes[0, 1], "M66", fig)
    fig_gen(axes[0, 2], "V66", fig)
    fig_gen(axes[0, 3], "L66", fig)

    # BOTTOM ROW
    fig_gen(axes[1, 0], "A66", fig)
    fig_gen(axes[1, 1], "Y66", fig)
    fig_gen(axes[1, 2], "I66", fig)
    axes[1, 3].axis('off')

    plt.savefig(f"{output_path}/all_blob_contacts_total_ensemble.pdf", bbox_inches='tight')
    plt.close()
  
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Plot the blob-blob contact frequencies.")

    parser.add_argument("--input_path", required=True, help="Path to the input data.") 
    parser.add_argument("--output_path", required=True, help="Path to the output data.")
    parser.add_argument("--cutoff_distance", required=False, default=0.55, help="Cutoff distance.")

    args = parser.parse_args()

    input_path = args.input_path
    output_path = args.output_path
    cutoff_distance = float(args.cutoff_distance)
    
    combine_blob_contact_plots()