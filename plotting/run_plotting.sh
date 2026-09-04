#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

ROOT_DIR="$(dirname "$SCRIPT_DIR")"

INPUT_DIR="$ROOT_DIR/data11"
FIGURE_DIR="$ROOT_DIR/figures11"

mkdir -p "$FIGURE_DIR"

# # Plots the radius of gyration over time for each sequence  
# python3 plot_rg_over_time.py --input_path "$INPUT_DIR" --output_path "$FIGURE_DIR" 					

# # Plots the cumulatiave average of Variant-Mediator contacts for the M66 sequence  
# python3 plot_cumulative_average_blob_contact.py --input_path "$INPUT_DIR" --output_path "$FIGURE_DIR" 

# # Plots the average radius of gyration for each sequence on a bar plot  
# python3 plot_avg_rg.py --input_path "$INPUT_DIR" --output_path "$FIGURE_DIR" 							

# # Simulates the SAHP and plots the blob-blob contacts ---FIX---
# python3 SAHP_polymer.py --input_path "$INPUT_DIR" --output_path "$FIGURE_DIR"

# # Plots the differences of blob-blob contacts for A66 versus the SAHP  
# python3 protein_vs_SAHP.py --input_path "$INPUT_DIR" --output_path "$FIGURE_DIR" 						

# Plots the differences of blob-blob contacts for all sequences versus A66  
# python3 plot_blob_contacts_vs_A66.py --input_path "$INPUT_DIR" --output_path "$FIGURE_DIR" 			

# # Plots radius of gyration vs Ree and N-C contacts  
# python3 plot_rg_vs_termini_correlation.py --input_path "$INPUT_DIR" --output_path "$FIGURE_DIR" 		

# # Plot the frequencies of contact states
# python3 plot_state_freq.py --association 8_2_to_2_16   --bar_plot_color darkblue --input_path "$INPUT_DIR" --output_path "$FIGURE_DIR"
# python3 plot_state_freq.py --association 8_16_to_2_16  --bar_plot_color maroon --input_path "$INPUT_DIR" --output_path "$FIGURE_DIR"
# python3 plot_state_freq.py --association 8_12_to_12_2  --bar_plot_color darkgreen --input_path "$INPUT_DIR" --output_path "$FIGURE_DIR"
# python3 plot_state_freq.py --association 12_2_to_2_16  --bar_plot_color indigo --input_path "$INPUT_DIR" --output_path "$FIGURE_DIR"
# python3 plot_state_freq.py --association 8_12_to_12_16 --bar_plot_color palevioletred --input_path "$INPUT_DIR" --output_path "$FIGURE_DIR"

# # # Plot odds ratios
# python3 plot_odds_ratios.py --association 8_2_to_2_16   --input_path "$INPUT_DIR" --output_path "$FIGURE_DIR"
# python3 plot_odds_ratios.py --association 8_16_to_2_16	--input_path "$INPUT_DIR" --output_path "$FIGURE_DIR"
# python3 plot_odds_ratios.py --association 8_12_to_12_2   --input_path "$INPUT_DIR" --output_path "$FIGURE_DIR"
# python3 plot_odds_ratios.py --association 12_2_to_2_16   --input_path "$INPUT_DIR" --output_path "$FIGURE_DIR"
# python3 plot_odds_ratios.py --association 8_12_to_12_16   --input_path "$INPUT_DIR" --output_path "$FIGURE_DIR"

# # Plot side-chain contacts  
# python3 plot_one_fig_range_residue_contacts.py --state VM_MN --x-range 93 97 --y-range 31 38 --input_path "$INPUT_DIR" --output_path "$FIGURE_DIR"
# python3 plot_one_fig_range_residue_contacts.py --state MN_NC --x-range 105 111 --y-range 31 38 --input_path "$INPUT_DIR" --output_path "$FIGURE_DIR"

# # Plot blob-blob contacts  
python3 plot_blob_contacts_one_fig.py --input_path "$INPUT_DIR" --output_path "$FIGURE_DIR" 			