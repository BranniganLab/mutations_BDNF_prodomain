import numpy as np
import matplotlib.pyplot as plt
import linecache
import argparse

"""
This script plots the cumulative average contact frequency over time for two specified blobs.
The data used for the contact frequency was generated from "calc_midpoints_of_all_blobs.tcl" and "calc_rg_of_all_blobs.tcl"
"""

def read_center_file():
    """
    Reads the coordinates of the blob centers from a file.
    
    Arguments:
        None
    
    Returns:
        blob_center_coords : numpy.ndarray
            Array of blob center coordinates.
    """

    center_coords_filename = f"midpoint_{sequence}_wrapped_centered_PIF_equil_del.txt"
    blob_center_coords = np.loadtxt(f"{input_path}/{center_coords_filename}")

    return blob_center_coords

def calc_contacts():
    """
    Reads the coordinates of the blob centers from a file.
    
    Arguments:
        None

    Returns:
        contacts : numpy.ndarray
            Array of contacts.
    """

    blob_center_coords = read_center_file()

    # Columns corresponding to each blob
    i_cols = [blob_i * 3, blob_i * 3 + 1, blob_i * 3 + 2]
    j_cols = [blob_j * 3, blob_j * 3 + 1, blob_j * 3 + 2]

    # Get XYZ coordinates
    coords_i = blob_center_coords[:, i_cols]
    coords_j = blob_center_coords[:, j_cols]

    # Read Rg
    rgyr_filename = f"blob_Rg_{sequence}_wrapped_centered_PIF_equil_del.txt"
    rgyr_values_blob_i = linecache.getline(f"{input_path}/{rgyr_filename}",blob_i + 1).strip()
    rgyr_values_blob_j = linecache.getline(f"{input_path}/{rgyr_filename}",blob_j + 1).strip()

    list_of_rgyr_i = np.array([float(value) for value in rgyr_values_blob_i.split()])
    list_of_rgyr_j = np.array([float(value) for value in rgyr_values_blob_j.split()])

    # Distance between blob centers 
    center_distance = np.linalg.norm(coords_i - coords_j,axis=1) 

    # Excess distance
    excess_distance = (center_distance - list_of_rgyr_i - list_of_rgyr_j)

    # Convert to contact or no contact
    contacts = (excess_distance < cutoff_distance).astype(int)

    return contacts

def calc_cumulative_contact_frequency():
    """
    Calculates the cumulative average contact frequency at every frame.
    
    Arguments:
        None

    Returns:
        cumulative_contact_frequency : numpy.ndarray
            Array of cumulative contact frequencies.
    """

    contacts = calc_contacts()

    cumulative_contact_frequency = (np.cumsum(contacts) / np.arange(1, len(contacts) + 1))

    return cumulative_contact_frequency

def fig_gen():
    """
    Generates a scatter plot of the cumulative average contact frequency over time for two specified blobs.
    
    Arguments:
        None

    Returns:
        None
    """

    fig, ax = plt.subplots(figsize=(10, 6))

    cumulative_frequency = calc_cumulative_contact_frequency()

    time_in_ps = np.arange(len(cumulative_frequency)) * ps_per_frame # Convert frames to time in picoseconds
    time_in_ns = time_in_ps / 1000  # Convert time to nanoseconds

    ax.plot(time_in_ns, cumulative_frequency, color="black", linewidth=2)
    ax.set_ylim(0, 1.1)

    ax.set_xlabel(r"t (ns)", fontsize=20)
    ax.set_ylabel(r"$f_{v\mathrm{-}m}(t)$", fontsize=20)
    ax.tick_params(axis='both', labelsize=18)

    blob_name_i = id_to_name[blob_i]
    blob_name_j = id_to_name[blob_j]
    output_filename = (f"cumulative_avg_{blob_name_i}_{blob_name_j}_seq_{sequence}.pdf")

    plt.savefig(f"{output_path}/{output_filename}",bbox_inches="tight")
    plt.close()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Plot cumulative frequency contacts between the Variant and Mediator blobs.")
    
    parser.add_argument("--input_path", required=True, help="Path to the input data.") 
    parser.add_argument("--output_path", required=True, help="Path to the output data.")
    parser.add_argument("--blob_i", required=False, default=5, help="Index of the first blob.")
    parser.add_argument("--blob_j", required=False, default=8, help="Index of the second blob.")
    parser.add_argument("--sequence", required=False, default="M66", help="Sequence name.")
    parser.add_argument("--cutoff_distance", required=False, default=0.55, help="Cutoff distance.")
    parser.add_argument("--ps_per_frame", required=False, default=100, help="Time step in picoseconds.")
    
    args = parser.parse_args()

    input_path = args.input_path
    output_path = args.output_path
    blob_i = int(args.blob_i)
    blob_j = int(args.blob_j)
    sequence = args.sequence
    cutoff_distance = float(args.cutoff_distance)
    ps_per_frame = float(args.ps_per_frame)

    id_to_name = {
        0: "p1",
        1: "h1a",
        2: "h1b",
        3: "p2",
        4: "h2a",
        5: "h2b",
        6: "p3",
        7: "h3a",
        8: "h3b",
        9: "h3c",
        10: "h3d",
    }
    
    fig_gen()


