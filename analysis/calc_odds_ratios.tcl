proc calc_odds_ratios {contactA_blobi contactA_blobj contactB_blobi contactB_blobj cutoff output_dir traj_without_equilibration_base} {

    #   Writes the number of frames that display contact-states and calculates the odds ratio (OR)
    #   with these values for all blob pairs involving the Variant (V) and a termini blob to a file.
    # 	The 'X' in the variables below refers to a blob of the protein.
    #
    #   Arguments:
    # 		contactA_blobi (int): Index of blobi in contact A of the association, defined by blobulation
    #       contactA_blobj (int): Index of blobj in contact A of the association, defined by blobulation
    # 		contactB_blobi (int): Index of blobi in contact B of the association, defined by blobulation
    # 		contactB_blobj (int): Index of blobj in contact B of the association, defined by blobulation
    # 		lMin (int): Minimum length of a blob (4)
    #       H (float): Hydropathy threshold (0.37)
    #		dictInput (string): Hydrophobicity scale ("Kyte-Doolittle")
    # 		cutoff (float): Cutoff distance between two blobs for a contact in Angstroms (5.5)
    #
    #   Returns:
    #       None

    # Calculate excess distance between midpoints of blobs
    source calc_excess_dist_of_blob_sel.tcl

	set filename "${output_dir}/${traj_without_equilibration_base}_odds_ratios_${contactA_blobi}_${contactA_blobj}_to_${contactB_blobi}_${contactB_blobj}.txt"
    set fp [open $filename "w"]

	set Aij_Bij_states   [get_Aij_Bij_state   $contactA_blobi $contactA_blobj $contactB_blobi $contactB_blobj $cutoff]
	set 0Aij_Bij_states  [get_0Aij_Bij_state  $contactA_blobi $contactA_blobj $contactB_blobi $contactB_blobj $cutoff]
	set Aij_0Bij_states  [get_Aij_0Bij_state  $contactA_blobi $contactA_blobj $contactB_blobi $contactB_blobj $cutoff] 
	set 0Aij_0Bij_states [get_0Aij_0Bij_state $contactA_blobi $contactA_blobj $contactB_blobi $contactB_blobj $cutoff] 

	set a $Aij_Bij_states
	set b $0Aij_Bij_states
	set c $Aij_0Bij_states
	set d $0Aij_0Bij_states

	if {$b==0 || $c ==0} {
		set OR "NA"
	} else {
		set OR [expr {1.0 * ($a * $d) / ($b * $c)}]
	}

	set Ai [expr {int($contactA_blobi)}]
	set Aj [expr {int($contactA_blobj)}]
	set Bi [expr {int($contactB_blobi)}]
	set Bj [expr {int($contactB_blobj)}]
	
	puts $fp "${Ai}_${Aj}->${Bi}_${Bj} ${a} ${b} ${c} ${d} ${OR}" 

	close $fp
}

proc get_Aij_Bij_state {contactA_blobi contactA_blobj contactB_blobi contactB_blobj cutoff} {

    #   Writes the number of frames that display the Aij_Bij contact state.
    #
    #   Arguments:
    # 		contactA_blobi (int): Index of blobi in contact A of the association, defined by blobulation
    #       contactA_blobj (int): Index of blobj in contact A of the association, defined by blobulation
    # 		contactB_blobi (int): Index of blobi in contact B of the association, defined by blobulation
    # 		contactB_blobj (int): Index of blobj in contact B of the association, defined by blobulation
    # 		cutoff (float): Cutoff distance between two blobs for a contact in Angstroms (5.5)
    #
    #   Returns:
    #       num_contact_states (int): The number of contact states (frames) within an ensemble
	
    set numFrames [molinfo top get numframes]
	set num_contact_states 0

	for {set i 0} {$i < $numFrames} {incr i} {
		set distance_between_Aij [calcExcessDistMid $contactA_blobi $contactA_blobj $i]
		set distance_between_Bij [calcExcessDistMid $contactB_blobi $contactB_blobj $i]

		if {$distance_between_Aij < $cutoff  && $distance_between_Bij < $cutoff} {
			incr num_contact_states
		}
	}

	return $num_contact_states
}

proc get_0Aij_Bij_state {contactA_blobi contactA_blobj contactB_blobi contactB_blobj cutoff} {

    #   Writes the number of frames that display the 0Aij_Bij contact state.
    #
    #   Arguments:
    # 		contactA_blobi (int): Index of blobi in contact A of the association, defined by blobulation
    #       contactA_blobj (int): Index of blobj in contact A of the association, defined by blobulation
    # 		contactB_blobi (int): Index of blobi in contact B of the association, defined by blobulation
    # 		contactB_blobj (int): Index of blobj in contact B of the association, defined by blobulation
    # 		cutoff (float): Cutoff distance between two blobs for a contact in Angstroms (5.5)	
    #
    #   Returns:
    #       num_contact_states (int): The number of contact states (frames) within an ensemble
	
    set numFrames [molinfo top get numframes]
	set num_contact_states 0

	for {set i 0} {$i < $numFrames} {incr i} {
		set distance_between_Aij [calcExcessDistMid $contactA_blobi $contactA_blobj $i]
		set distance_between_Bij [calcExcessDistMid $contactB_blobi $contactB_blobj $i]

		if {$distance_between_Aij > $cutoff  && $distance_between_Bij < $cutoff} {
			incr num_contact_states
		}
	}

	return $num_contact_states
}

proc get_Aij_0Bij_state {contactA_blobi contactA_blobj contactB_blobi contactB_blobj cutoff} {
    
    #   Writes the number of frames that display the Aij_0Bij contact state.
    #
    #   Arguments:
    # 		contactA_blobi (int): Index of blobi in contact A of the association, defined by blobulation
    #       contactA_blobj (int): Index of blobj in contact A of the association, defined by blobulation
    # 		contactB_blobi (int): Index of blobi in contact B of the association, defined by blobulation
    # 		contactB_blobj (int): Index of blobj in contact B of the association, defined by blobulation
    # 		cutoff (float): Cutoff distance between two blobs for a contact in Angstroms (5.5)
    #
    #   Returns:
    #       num_contact_states (int): The number of contact states (frames) within an ensemble	
    
    set numFrames [molinfo top get numframes]
	set num_contact_states 0

	for {set i 0} {$i < $numFrames} {incr i} {
		set distance_between_Aij [calcExcessDistMid $contactA_blobi $contactA_blobj $i]
		set distance_between_Bij [calcExcessDistMid $contactB_blobi $contactB_blobj $i]

		if {$distance_between_Aij < $cutoff  && $distance_between_Bij > $cutoff} {
			incr num_contact_states
		}
	}

	return $num_contact_states
}

proc get_0Aij_0Bij_state {contactA_blobi contactA_blobj contactB_blobi contactB_blobj cutoff} {

    #   Writes the number of frames that display the 0Aij_0Bij contact state.
    #
    #   Arguments:
    # 		contactA_blobi (int): Index of blobi in contact A of the association, defined by blobulation
    #       contactA_blobj (int): Index of blobj in contact A of the association, defined by blobulation
    # 		contactB_blobi (int): Index of blobi in contact B of the association, defined by blobulation
    # 		contactB_blobj (int): Index of blobj in contact B of the association, defined by blobulation
    # 		cutoff (float): Cutoff distance between two blobs for a contact in Angstroms (5.5)
    #
    #   Returns:
    #       num_contact_states (int): The number of contact states (frames) within an ensemble	
	
    set numFrames [molinfo top get numframes]
	set num_contact_states 0

	for {set i 0} {$i < $numFrames} {incr i} {
		set distance_between_Aij [calcExcessDistMid $contactA_blobi $contactA_blobj $i]
		set distance_between_Bij [calcExcessDistMid $contactB_blobi $contactB_blobj $i]

		if {$distance_between_Aij > $cutoff  && $distance_between_Bij > $cutoff} {
			incr num_contact_states
		}
	}

	return $num_contact_states
}