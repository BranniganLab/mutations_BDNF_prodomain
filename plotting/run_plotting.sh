#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

ROOT_DIR="$(dirname "$SCRIPT_DIR")"

REPLICA="0"
TRAJ_DIR="$ROOT_DIR/trajectories/replica_$REPLICA/raw"
INPUT_DIR="$ROOT_DIR/data/replica_$REPLICA/Lmin4_Hstar0.37_original_params/cutoff_5.5A_original_params"
FIGURE_DIR="$ROOT_DIR/figures/replica_$REPLICA/Lmin4_Hstar0.37_original_params/cutoff_5.5A_original_params"

# These must match the parameters used by run_analysis.sh.
L_MIN=4 #4 5
H_STAR=0.37 #0.33, 0.40, 0.37

BLOB_INDICES_FILE="$INPUT_DIR/blob_indices_Lmin${L_MIN}_Hstar${H_STAR}.txt"

if [[ ! -f "$BLOB_INDICES_FILE" ]]; then
	echo "ERROR: Blob-index file not found: $BLOB_INDICES_FILE"
	echo "Run analysis first for the selected L_MIN/H_STAR parameters."
	exit 1
fi

#Load the file that contains index data of blobs
source "$BLOB_INDICES_FILE"

VARIANT_BLOB_INDEX="$variant_blob_index"
MEDIATOR_BLOB_INDEX="$mediator_blob_index"
N_BLOB_INDEX="$N_blob_index"
C_BLOB_INDEX="$C_blob_index"
VARIANT_BLOB_INDEX_NO_SBLOBS="$variant_blob_index_no_s_blobs"
MEDIATOR_BLOB_INDEX_NO_SBLOBS="$mediator_blob_index_no_s_blobs"

mkdir -p "$FIGURE_DIR"

# # Plots the radius of gyration over time for each sequence  
# python3 plot_rg_over_time.py \
# 	--input_path "$INPUT_DIR" \
# 	--output_path "$FIGURE_DIR" 	

# # Plots the cumulative average radius of gyration over time for all sequences
# python3 plot_cumulative_average_rg.py \
# 	--input_path "$INPUT_DIR" \
# 	--output_path "$FIGURE_DIR"  	

# Plots the SEM of Rg vs the block size and estimates the blocking effective sample size
python3 get_effective_sample_size.py \
	--input_path "$INPUT_DIR" \
	--output_path "$FIGURE_DIR" \
	--plateau_block_size 500

# Compares odds ratios using the original and blocking effective sample sizes
python3 plot_odds_ratios_blocking_comparison.py \
	--input_path "$INPUT_DIR" \
	--output_path "$FIGURE_DIR" \
	--effective_sample_size_file "$FIGURE_DIR/effective_sample_size_M66.txt"

# # Plots the cumulatiave average of Variant-Mediator contacts for the M66 sequence  
# python3 plot_cumulative_average_blob_contact.py \
# 	--input_path "$INPUT_DIR" \
# 	--output_path "$FIGURE_DIR" \
#  	--variant_blob_index_no_s_blobs "$VARIANT_BLOB_INDEX_NO_SBLOBS" \
#  	--mediator_blob_index_no_s_blobs "$MEDIATOR_BLOB_INDEX_NO_SBLOBS" 

# # Plots the average radius of gyration for each sequence on a bar plot  
# python3 plot_avg_rg.py \
# 	--input_path "$INPUT_DIR" \
# 	--output_path "$FIGURE_DIR" 							

# # Simulates the SAHP and plots the blob-blob contacts
# python3 SAHP_polymer.py \
# 	--input_path "$INPUT_DIR" \
# 	--output_path "$FIGURE_DIR"

# # Plots the differences of blob-blob contacts for A66 versus the SAHP  
# python3 protein_vs_SAHP.py \
# 	--input_path "$INPUT_DIR" \
# 	--output_path "$FIGURE_DIR" 						

# # Plots the differences of blob-blob contacts for all sequences versus A66  
# python3 plot_blob_contacts_vs_A66.py \
# 	--input_path "$INPUT_DIR" \
# 	--output_path "$FIGURE_DIR" 			

# # Plots radius of gyration vs Ree and N-C contacts  
# python3 plot_rg_vs_termini_correlation.py \
# 	--input_path "$INPUT_DIR" \
# 	--output_path "$FIGURE_DIR" \
# 	--N_blob_index "$N_BLOB_INDEX" \
# 	--C_blob_index "$C_BLOB_INDEX" 		

# # Plot the frequencies of contact states: V-N -> N-C
# python3 plot_state_freq.py \
# 	--association "$VARIANT_BLOB_INDEX"_"$N_BLOB_INDEX"_to_"$N_BLOB_INDEX"_"$C_BLOB_INDEX" \
# 	--bar_plot_color darkblue \
# 	--input_path "$INPUT_DIR" \
# 	--output_path "$FIGURE_DIR" \
# 	--variant_blob_index "$VARIANT_BLOB_INDEX" \
# 	--mediator_blob_index "$MEDIATOR_BLOB_INDEX" \
# 	--N_blob_index "$N_BLOB_INDEX" \
# 	--C_blob_index "$C_BLOB_INDEX" 

# # Plot the frequencies of contact states: V-C -> N-C
# python3 plot_state_freq.py \
# 	--association "$VARIANT_BLOB_INDEX"_"$C_BLOB_INDEX"_to_"$N_BLOB_INDEX"_"$C_BLOB_INDEX" \
# 	--bar_plot_color maroon \
# 	--input_path "$INPUT_DIR" \
# 	--output_path "$FIGURE_DIR" \
# 	--variant_blob_index "$VARIANT_BLOB_INDEX" \
# 	--mediator_blob_index "$MEDIATOR_BLOB_INDEX" \
# 	--N_blob_index "$N_BLOB_INDEX" \
# 	--C_blob_index "$C_BLOB_INDEX" 

# # Plot the frequencies of contact states: V-M -> M-N
# python3 plot_state_freq.py \
# 	--association "$VARIANT_BLOB_INDEX"_"$MEDIATOR_BLOB_INDEX"_to_"$MEDIATOR_BLOB_INDEX"_"$N_BLOB_INDEX" \
# 	--bar_plot_color darkgreen \
# 	--input_path "$INPUT_DIR" \
# 	--output_path "$FIGURE_DIR" \
# 	--variant_blob_index "$VARIANT_BLOB_INDEX" \
# 	--mediator_blob_index "$MEDIATOR_BLOB_INDEX" \
# 	--N_blob_index "$N_BLOB_INDEX" \
# 	--C_blob_index "$C_BLOB_INDEX" 

# # Plot the frequencies of contact states: M-N -> N-C
# python3 plot_state_freq.py \
# 	--association "$MEDIATOR_BLOB_INDEX"_"$N_BLOB_INDEX"_to_"$N_BLOB_INDEX"_"$C_BLOB_INDEX" \
# 	--bar_plot_color indigo \
# 	--input_path "$INPUT_DIR" \
# 	--output_path "$FIGURE_DIR" \
# 	--variant_blob_index "$VARIANT_BLOB_INDEX" \
# 	--mediator_blob_index "$MEDIATOR_BLOB_INDEX" \
# 	--N_blob_index "$N_BLOB_INDEX" \
# 	--C_blob_index "$C_BLOB_INDEX" 

# # Plot the frequencies of contact states: V-M -> M-C
# python3 plot_state_freq.py \
# 	--association "$VARIANT_BLOB_INDEX"_"$MEDIATOR_BLOB_INDEX"_to_"$MEDIATOR_BLOB_INDEX"_"$C_BLOB_INDEX" \
# 	--bar_plot_color palevioletred \
# 	--input_path "$INPUT_DIR" \
# 	--output_path "$FIGURE_DIR" \
# 	--variant_blob_index "$VARIANT_BLOB_INDEX" \
# 	--mediator_blob_index "$MEDIATOR_BLOB_INDEX" \
# 	--N_blob_index "$N_BLOB_INDEX" \
# 	--C_blob_index "$C_BLOB_INDEX" 

# # Plot odds ratios: V-N -> N-C
# python3 plot_odds_ratios.py \
# 	--association "$VARIANT_BLOB_INDEX"_"$N_BLOB_INDEX"_to_"$N_BLOB_INDEX"_"$C_BLOB_INDEX" \
# 	--input_path "$INPUT_DIR" \
# 	--output_path "$FIGURE_DIR" \
# 	--variant_blob_index "$VARIANT_BLOB_INDEX" \
# 	--mediator_blob_index "$MEDIATOR_BLOB_INDEX" \
# 	--N_blob_index "$N_BLOB_INDEX" \
# 	--C_blob_index "$C_BLOB_INDEX" 	

# # Plot odds ratios: V-C -> N-C
# python3 plot_odds_ratios.py \
# 	--association "$VARIANT_BLOB_INDEX"_"$C_BLOB_INDEX"_to_"$N_BLOB_INDEX"_"$C_BLOB_INDEX" \
# 	--input_path "$INPUT_DIR" \
# 	--output_path "$FIGURE_DIR" \
# 	--variant_blob_index "$VARIANT_BLOB_INDEX" \
# 	--mediator_blob_index "$MEDIATOR_BLOB_INDEX" \
# 	--N_blob_index "$N_BLOB_INDEX" \
# 	--C_blob_index "$C_BLOB_INDEX" 	

# # Plot odds ratios: V-M -> M-N
# python3 plot_odds_ratios.py \
# 	--association "$VARIANT_BLOB_INDEX"_"$MEDIATOR_BLOB_INDEX"_to_"$MEDIATOR_BLOB_INDEX"_"$N_BLOB_INDEX" \
# 	--input_path "$INPUT_DIR" \
# 	--output_path "$FIGURE_DIR" \
# 	--variant_blob_index "$VARIANT_BLOB_INDEX" \
# 	--mediator_blob_index "$MEDIATOR_BLOB_INDEX" \
# 	--N_blob_index "$N_BLOB_INDEX" \
# 	--C_blob_index "$C_BLOB_INDEX" 	

# # Plot odds ratios: M-N -> N-C
# python3 plot_odds_ratios.py \
# 	--association "$MEDIATOR_BLOB_INDEX"_"$N_BLOB_INDEX"_to_"$N_BLOB_INDEX"_"$C_BLOB_INDEX" \
# 	--input_path "$INPUT_DIR" \
# 	--output_path "$FIGURE_DIR" \
# 	--variant_blob_index "$VARIANT_BLOB_INDEX" \
# 	--mediator_blob_index "$MEDIATOR_BLOB_INDEX" \
# 	--N_blob_index "$N_BLOB_INDEX" \
# 	--C_blob_index "$C_BLOB_INDEX" 	

# # Plot odds ratios: V-M -> M-C
# python3 plot_odds_ratios.py \
# 	--association "$VARIANT_BLOB_INDEX"_"$MEDIATOR_BLOB_INDEX"_to_"$MEDIATOR_BLOB_INDEX"_"$C_BLOB_INDEX" \
# 	--input_path "$INPUT_DIR" \
# 	--output_path "$FIGURE_DIR" \
# 	--variant_blob_index "$VARIANT_BLOB_INDEX" \
# 	--mediator_blob_index "$MEDIATOR_BLOB_INDEX" \
# 	--N_blob_index "$N_BLOB_INDEX" \
# 	--C_blob_index "$C_BLOB_INDEX" 	

# # Plot side-chain contacts for V-M M-N contact state
# python3 plot_one_fig_range_residue_contacts.py \
# 	--state VM_MN \
# 	--x-range 93 97 \
# 	--y-range 31 38 \
# 	--res_type "sc" \
# 	--input_path "$INPUT_DIR" \
# 	--output_path "$FIGURE_DIR"

# # Plot side-chain contacts for M-N N-C contact state
# python3 plot_one_fig_range_residue_contacts.py \
# 	--state MN_NC \
# 	--x-range 105 111 \
# 	--y-range 31 38 \
# 	--res_type "sc" \
# 	--input_path "$INPUT_DIR" \
# 	--output_path "$FIGURE_DIR"

# # Plot blob-blob contacts  
# python3 plot_blob_contacts_one_fig.py \
# 	--input_path "$INPUT_DIR" \
# 	--output_path "$FIGURE_DIR" 	

# # Plot secondary structure for V-M M-N contact state, raw_traj_path is need for .tpr file
# python3 calc_secondary_structure.py \
# 	--raw_traj_path "$TRAJ_DIR" \
# 	--input_path "$INPUT_DIR" \
# 	--output_path "$FIGURE_DIR" \
# 	--state VM_NC

# # Plot backbone contacts for V-M M-N contact state
# python3 plot_one_fig_range_residue_contacts.py \
# 	--state VM_NC \
# 	--x-range 93 97 \
# 	--y-range 65 72 \
# 	--res_type "bb" \
# 	--input_path "$INPUT_DIR" \
# 	--output_path "$FIGURE_DIR"		

# # Plot backbone contacts for V-M M-N contact state
# python3 plot_one_fig_range_residue_contacts.py \
# 	--state VM_NC \
# 	--x-range 105 111 \
# 	--y-range 31 38 \
# 	--res_type "bb" \
# 	--input_path "$INPUT_DIR" \
# 	--output_path "$FIGURE_DIR"	

# python3 plot_all_residue_contacts.py \
# 	--state "total_ensemble" \
# 	--res_type "ca" \
# 	--input_path "$INPUT_DIR" \
# 	--output_path "$FIGURE_DIR"	

# python3 plot_all_residue_contacts.py \
# 	--state "total_ensemble" \
# 	--res_type "bb" \
# 	--input_path "$INPUT_DIR" \
# 	--output_path "$FIGURE_DIR"	

# python3 plot_all_residue_contacts.py \
# 	--state "total_ensemble" \
# 	--res_type "sc" \
# 	--input_path "$INPUT_DIR" \
# 	--output_path "$FIGURE_DIR"	
