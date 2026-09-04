import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
import linecache
import argparse

"""
This script plots the difference in contact frequencies between a given sequence and the SAHP sequence on a heat map. 
The SAHP simulation was generated from "SAHP_polymer.py" and the data to calculate the contact frequencies
was generated from "calc_midpoints_of_all_blobs.tcl" and "calc_rg_of_all_blobs.tcl"
"""

def read_reference_file():
    """
    Reads the XYZ coordinates associated with the center of a blob. All coordinates of all blobs are stored in this file.
    Each row in the text file corresponds to a frame number and each column represents x, y, or z coordinates. 
    Columns read as: column 1:x data, column 2:y data, column 3:z data (of the first blob of the sequence),
                     column 4:x data, column 5:y data, column 6:z data (of the second blob of the sequence)
                     column 7:x data, column 8:y data, column 9:z data (of the third blob of the sequence and so on)
    Arguments:
        None
    
    Returns:
        blob_center_coords : numpy.ndarray
            Array of blob center coordinates.
    """
    center_coords_filename = f"midpoint_{ref_seq}_wrapped_centered_PIF_equil_del.txt" 
    blob_center_coords = np.loadtxt(f'{input_path}/{center_coords_filename}')

    return blob_center_coords

def get_num_of_blobs():
    """
    Counts the number of lines in the Rg file which is equivalent to the number of 
    blobs in the protein and reads the name of the Rg file.

    Arguments:
        None
    
    Returns:
        num_of_blobs : int
            The number of blobs in the protein.
        rgyr_filename : str
            The filename of the Rg file.
    """
    rgyr_filename = f"blob_Rg_{ref_seq}_wrapped_centered_PIF_equil_del.txt"
    
    with open(f"{input_path}/{rgyr_filename}", 'r') as fp:
        num_of_blobs = len(fp.readlines()) # The number of blobs is equal to the number of lines in the rgyr file

    return num_of_blobs, rgyr_filename

def calc_excess_distance_and_contacts():
    """
    Calculates the excess distance between two blobs and if the distance between the two blobs is less than the cutoff distance, a counter will update the contact_probability array.

    Arguments:
        None
    
    Returns:
        contact_probability : numpy.ndarray
            Array of contact probabilities.
    """ 
    blob_center_coords = read_reference_file() 
    num_of_blobs, rgyr_filename = get_num_of_blobs()

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

def calc_contact_frequency():  
    """
    Calculates the average of the contact probability array and expresses it as a percentage by mulitplying by 100.

    Arguments:
        None  
    
    Returns:
        contact_frequency : numpy.ndarray
            Array of contact frequencies.
    """ 
    contact_probability = calc_excess_distance_and_contacts()
    
    contact_frequency = np.mean(contact_probability, axis=0) * 100 # Averages the values across each time step for each pair of blobs and expresses it as a percentage
    
    return contact_frequency

def read_SAHP_file():
    """
    Reads the contact frequency of blobs within the SAHP.
    
    Arguments:
        None
    
    Returns:
        SAHP_data : numpy.ndarray
            Array of contact frequencies of the blobs within the SAHP
    """

    SAHP_path = f"{input_path}/A66_parameritization_SAHP.txt"
    SAHP_data = np.loadtxt(SAHP_path)

    return SAHP_data

def calc_variant_SAHP_difference():
    """
    Calculates the difference in contact frequencies between a given sequence and the SAHP sequence.
    
    Arguments:
        None
    
    Returns:
        difference : numpy.ndarray
            Array of contact frequency differences.
    """
    contact_frequency_variant = calc_contact_frequency()
    contact_frequency_SAHP = read_SAHP_file()
    difference = contact_frequency_variant - contact_frequency_SAHP

    return difference

def fig_gen():
    """
    Generates a figure of the contact frequencies on a heat map.
    
    Arguments:
        None
    
    Returns:
        None
    """
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 4))
    probability_mean_difference = calc_variant_SAHP_difference()
    imgp = ax.imshow(probability_mean_difference, origin='lower', aspect='equal', cmap='PRGn', vmin=-100, vmax=100)
    ax.set_title(rf'{ref_seq}-SAHP', loc='center', fontsize='20')

    # Draw the colorbar right next to the subplot
    box = ax.get_position()  # [x0, y0, width, height]
    cbar_width = 0.02
    cbar_pad = 0.01
    cbar_x = box.x1 + cbar_pad
    cbar_ax = fig.add_axes([cbar_x, box.y0, cbar_width, box.height])

    cbar = fig.colorbar(imgp, cax=cbar_ax)
    cbar.ax.tick_params(labelsize=10)
    cbar.set_label("Contact frequency (%)", fontsize=15, rotation=270, labelpad=14)

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

    output_filename = f"{ref_seq}_vs_SAHP.pdf"
    plt.savefig(f"{output_path}/{output_filename}", bbox_inches='tight')
    plt.close()
    
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Plot the differences in contact frequencies between a given sequence and the SAHP.")

    parser.add_argument("--input_path", required=True, help="Path to the input data.") 
    parser.add_argument("--output_path", required=True, help="Path to the output data.")
    parser.add_argument("--ref_seq", required=False, default="A66", help="Reference sequence.")
    parser.add_argument("--cutoff_distance", required=False, default=0.55, help="Cutoff distance.")

    args = parser.parse_args()

    input_path = args.input_path
    output_path = args.output_path
    ref_seq = args.ref_seq
    cutoff_distance = float(args.cutoff_distance)

    fig_gen()