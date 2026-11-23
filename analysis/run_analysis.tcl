# Set arguments
set structure [lindex $argv 0]
set full_traj [lindex $argv 1]
set traj_without_equilibration [lindex $argv 2]
set output_dir [lindex $argv 3]
set blobulator_path [lindex $argv 4]
set Lmin [lindex $argv 5]
set H_star [lindex $argv 6]
set hydrophobicity_scale [lindex $argv 7]
set N_terminal_index [lindex $argv 8]
set C_terminal_index [lindex $argv 9]
set blob_distance_cutoff [lindex $argv 10]
set variant_index [lindex $argv 11]
set mediator_index [lindex $argv 12]
set state1 [lindex $argv 13]
set state2 [lindex $argv 14]
set state3 [lindex $argv 15]
set residue_contact_distance_cutoff [lindex $argv 16]
set seq [lindex $argv 17]
set ref_seq [lindex $argv 18]

set full_traj_base [file rootname [file tail $full_traj]]
set traj_without_equilibration_base [file rootname [file tail $traj_without_equilibration]]

#----------Run procedures----------
# Load sequence into VMD and remove the first frame since it is the gro structure
mol new $structure type gro waitfor all
mol addfile $full_traj type xtc waitfor all
animate delete beg 0 end 0 

# Measures the radius of gyration over time of the protein that includes equilibration time
source calc_rg_of_protein.tcl
calcRgofProtein $output_dir $full_traj_base

# Load the trajectory that doesnt include equilibration time into VMD
mol new $structure type gro waitfor all
mol addfile $traj_without_equilibration type xtc waitfor all
animate delete beg 0 end 0

# Measures the radius of gyration over time of the protein that does not include equilibration time
source calc_rg_of_protein.tcl
calcRgofProtein	$output_dir $traj_without_equilibration_base 

# Measure distance between alpha carbons, sidechains and backbones
source calc_alpha_carbon_contacts.tcl 
measureContactsCA $residue_contact_distance_cutoff $output_dir "total_ensemble"
source calc_sidechain_contacts.tcl
measureContactsSC $residue_contact_distance_cutoff $output_dir "total_ensemble"
source calc_backbone_contacts.tcl
measureContactsBB $residue_contact_distance_cutoff $output_dir "total_ensemble"

if {$seq eq $ref_seq} {
	# Blobulate protein  
	source blobulate_protein.tcl
	blobulate_protein $Lmin $H_star $hydrophobicity_scale $blobulator_path

	# Get and export blob assignments of V66 to apply to other variants
	source export_blob_assignments.tcl
	export_blob_assignments $output_dir $Lmin $H_star

} else {
	source apply_ref_seq_blob_assignments.tcl
	apply_ref_seq_blob_assignments $output_dir $Lmin $H_star
}

# Scripts for SAHP parametrization  
source calc_avg_rg_of_blobs.tcl
calcRgofBlobs $output_dir $traj_without_equilibration_base
source calc_avg_ree_of_blobs.tcl
calcReeofBlobs $output_dir $traj_without_equilibration_base		 		

# Scripts for blob-blob contacts  
source calc_rg_of_all_blobs.tcl
calcAllBlobRg $output_dir $traj_without_equilibration_base
source calc_midpoints_of_all_blobs.tcl
calcAllBlobMidpoints $output_dir $traj_without_equilibration_base								

# Measures the end-to-end distance of a protein over time  
source calc_ree_of_protein.tcl
calcReeofProtein $output_dir $traj_without_equilibration_base	

# Measures contacts for one blob pair  
source calc_one_pair_blob_contact.tcl
calc_one_pair_blob_contact $N_terminal_index $C_terminal_index $blob_distance_cutoff $output_dir $traj_without_equilibration_base

# Measures odds ratios  
source calc_odds_ratios.tcl 
calc_odds_ratios $variant_index $N_terminal_index $N_terminal_index $C_terminal_index $blob_distance_cutoff $output_dir $traj_without_equilibration_base
calc_odds_ratios $variant_index $C_terminal_index $N_terminal_index $C_terminal_index $blob_distance_cutoff $output_dir $traj_without_equilibration_base
calc_odds_ratios $variant_index $mediator_index $mediator_index $N_terminal_index $blob_distance_cutoff $output_dir $traj_without_equilibration_base
calc_odds_ratios $mediator_index $N_terminal_index $N_terminal_index $C_terminal_index $blob_distance_cutoff $output_dir $traj_without_equilibration_base
calc_odds_ratios $variant_index $mediator_index $mediator_index $C_terminal_index $blob_distance_cutoff $output_dir $traj_without_equilibration_base

# Classfies conformational ensemble into contact states and measures side-chain contacts
source classify_contact_state.tcl
classify_contact_state $variant_index $mediator_index $N_terminal_index $C_terminal_index $blob_distance_cutoff $state1 $output_dir
source calc_sidechain_contacts.tcl
measureContactsSC $residue_contact_distance_cutoff $output_dir $state1
mol delete top

# Load sequence into VMD
mol new $structure type gro waitfor all
mol addfile $traj_without_equilibration type xtc waitfor all
animate delete beg 0 end 0

# Blobulate protein by applying the blob assignments of the reference sequence
source apply_ref_seq_blob_assignments.tcl
apply_ref_seq_blob_assignments $output_dir $Lmin $H_star

# Classfies conformational ensemble into contact states and measures side-chain contacts
classify_contact_state $variant_index $mediator_index $N_terminal_index $C_terminal_index $blob_distance_cutoff $state2 $output_dir
measureContactsSC $residue_contact_distance_cutoff $output_dir $state2
mol delete top

# Load sequence into VMD
mol new $structure type gro waitfor all
mol addfile $traj_without_equilibration type xtc waitfor all
animate delete beg 0 end 0

# Blobulate protein by applying the blob assignments of the reference sequence
source apply_ref_seq_blob_assignments.tcl
apply_ref_seq_blob_assignments $output_dir $Lmin $H_star

# Classfies conformational ensemble into contact states and measures backbone contacts
classify_contact_state $variant_index $mediator_index $N_terminal_index $C_terminal_index $blob_distance_cutoff $state3 $output_dir
source calc_backbone_contacts.tcl
measureContactsBB $residue_contact_distance_cutoff $output_dir $state3
			
quit