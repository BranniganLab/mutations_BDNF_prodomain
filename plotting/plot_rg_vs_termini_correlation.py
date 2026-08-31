import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy import stats
from matplotlib.ticker import MultipleLocator
import argparse

"""
This script plots the radius of gyration vs a termini association (e.g. end-to-end distance and terminal blob contacts)
The data used for the Rg was generated from "calc_rg_of_protein.tcl"
The data used for the Ree was generated from "calc_ree_of_protein.tcl" which calculates the distance between the of the Oxygen and Nitrogen of the termini of the blob 
The data used for the terminal blob distances was generated from "calc_one_pair_blob_contact.tcl"
"""

def read_data(sequence, data_type):
    """
    Reads a data file that the user specifies ('rg' or 'ree' or 'single_pair_blob_contact_2_16').
    
    Arguments
        sequence : str
            The sequence name.
    
    Returns
        data : pandas.DataFrame
            A DataFrame containing the data.
    """

    data = pd.read_csv(f'{input_path}/{data_type}_{sequence}_resid_23-113-capped.txt', sep=r'\s', header=None, names=['data'], engine='python')

    return data

def average_data(data):
    """
    Computes the average of the given data.

    Arguments
        data : pandas.DataFrame
            Contains the data for each frame of the given blob.

    Returns
        float
            The average value for the given blob.
    """
    
    return data['data'].mean()

def do_stats(data):
    """
    Compute the standard error of the mean (SEM) for the given data.
    
    Arguments
        data : pandas.DataFrame
            A DataFrame containing the data.

    Returns
        sem : float
            The standard error of the mean.
    """
    data = data['data'].to_numpy()
    std = np.std(data, ddof=1)
    sem = std / np.sqrt(n_independent_samples)
    
    return float(sem)

def plot_rg_vs_ree(ax):
    """
    Plot average Rg versus average Ree.
        
    Arguments
        ax : matplotlib.axes._subplots.AxesSubplot
            The axes on which to plot the data.
    
    Returns
        None
    """
    rg_values = []
    ree_values = []

    for sequence in sequences:
        rg = read_data(sequence, 'rg') 
        sem_rg = do_stats(rg)
        avg_rg = average_data(read_data(sequence, 'rg'))
        avg_ree = average_data(read_data(sequence, 'ree'))

        ax.scatter(avg_ree, avg_rg, s=50, color='black')
        ax.errorbar(avg_ree, avg_rg, yerr=sem_rg, capsize=3, elinewidth=1, color='black')
        ax.annotate(sequence, (avg_ree, avg_rg), textcoords="offset points", xytext=(15, 30), ha="right", fontsize=17)
        
        rg_values.append(avg_rg)
        ree_values.append(avg_ree)

    r, p = stats.pearsonr(ree_values, rg_values)

    ax.set_xlim(30, 38)
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.set_xlabel(r'$\langle R_{ee} \rangle$ (Å)', fontsize=20)

def plot_rg_vs_n_c_contacts(ax):
    """
    Plot average Rg versus frequency of states with termini contacts.
        
    Arguments
        ax : matplotlib.axes._subplots.AxesSubplot
            The axes on which to plot the data.

    Returns
        None
    """
    rg_values = []
    termini_values = []

    for sequence in sequences:
        rg = read_data(sequence, 'rg')
        sem_rg = do_stats(rg)
        avg_rg = average_data(read_data(sequence, 'rg'))
        contact_frequency = (average_data(read_data(sequence, 'single_pair_blob_contact_2_16'))) * 100
        
        ax.scatter(contact_frequency, avg_rg, s=50, color='black')
        ax.errorbar(contact_frequency, avg_rg, yerr=sem_rg, capsize=3, elinewidth=1, color='black')
        ax.annotate(sequence, (contact_frequency, avg_rg), textcoords="offset points", xytext=(0, 28), ha="center", fontsize=17)

        rg_values.append(avg_rg)
        termini_values.append(contact_frequency)

    r, p = stats.pearsonr(termini_values, rg_values)

    ax.set_xlabel(r'$n$-$c$ contact frequency (%)',fontsize=20)
    ax.set_xlim(26, 52)
    ax.xaxis.set_major_locator(MultipleLocator(3))

def plot_ree_vs_n_c_contacts(ax):
    """
    Plot average Ree versus frequency of states that display termini contacts.
        
    Arguments
        ax : matplotlib.axes._subplots.AxesSubplot
            The axes on which to plot the data.

    Returns
        None
    """
    ree_values = []
    termini_values = []

    for sequence in sequences:
        ree = read_data(sequence, 'ree')
        sem_ree = do_stats(ree)
        avg_ree = average_data(read_data(sequence, 'ree'))
        n_c_contact_frequency = (average_data(read_data(sequence, 'single_pair_blob_contact_2_16'))) * 100
        
        ax.scatter(n_c_contact_frequency, avg_ree, s=50, color='black')
        ax.errorbar(n_c_contact_frequency, avg_ree, yerr=sem_ree, capsize=3, elinewidth=1, color='black')
        ax.annotate(sequence, (n_c_contact_frequency, avg_ree), textcoords="offset points", xytext=(0, 28), ha="center", fontsize=17)
        ree_values.append(avg_ree)
        termini_values.append(n_c_contact_frequency)

    r, p = stats.pearsonr(termini_values, ree_values)

    ax.set_xlabel(r'$n$-$c$ contact frequency (%)',fontsize=20)
    ax.set_xlim(26, 52)
    ax.xaxis.set_major_locator(MultipleLocator(3))

def plot_multipanel():
    """
    Generates a figure with:
        Top left:  Rg vs Ree
        Top right: Rg vs Termini-contact frequency
        Bottom right: Ree vs Termini-contact frequency
        
    Arguments
        None

    Returns
        None
    """
    fig, axes = plt.subplots(2,2,figsize=(16, 12),sharey='row')

    # Left panel
    plot_rg_vs_ree(axes[0, 0])

    # Right panel
    plot_rg_vs_n_c_contacts(axes[0, 1])

    # Bottom right panel
    plot_ree_vs_n_c_contacts(axes[1, 0])

    axes[0, 0].set_ylabel(r'$\langle R_{g} \rangle$ (Å)', fontsize=20)
    axes[0, 0].set_ylim(17.7, 19.4)
    axes[0, 0].tick_params(axis='both', labelsize=18)
    axes[0, 0].yaxis.set_major_locator(MultipleLocator(0.2))

    axes[0, 1].tick_params(axis='both', labelsize=18)

    axes[1, 0].set_ylabel(r'$\langle R_{ee} \rangle$ (Å)', fontsize=20)
    axes[1, 0].set_ylim(30, 42)
    axes[1, 0].tick_params(axis='both', labelsize=18)

    axes[1, 1].axis('off')  # Hide the bottom right panel

    fig.subplots_adjust(wspace=0.08)
    plt.savefig(f"{output_path}/rg_ree_n_c_contact_correlations.pdf",bbox_inches='tight')
    plt.close()

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Plots rg vs ree, rg vs termini contact frequency and ree vs termini contact frequency.")

    parser.add_argument("--input_path", required=True, help="Path to the input data.") 
    parser.add_argument("--output_path", required=True, help="Path to the output data.")
    parser.add_argument("--n_samples", required=False, default=1088, help="Number of samples (1088).")
    parser.add_argument("--sequences", required=False, default="V66 M66 F66 A66 Y66 L66 I66", help="Sequence names.")

    args = parser.parse_args()

    input_path = args.input_path
    output_path = args.output_path
    n_independent_samples = int(args.n_samples)
    sequences = args.sequences.split()

    plot_multipanel()       


