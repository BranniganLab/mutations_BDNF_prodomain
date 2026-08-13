proc calc_odds_ratios_all_blob_pairs_VX_XT {variant_index termini_index lMin H dictInput cutoff} {

    #   Writes the number of frames that display contact-states and calculates the odds ratio (OR)
    #   with these values for all blob pairs involving the Variant (V) and a termini blob to a file.
    # 	The 'X' in the variables below refers to a blob of the protein.
    #
    #   Arguments:
    # 		variant_index (int): Index of Mutant blob, defined by blobulation
    #       termini_index (int): Index of a terminal blob, defined by blobulation
    # 		lMin (int): Minimum length of a blob (4)
    #       H (float): Hydropathy threshold (0.37)
    #		dictInput (string): Hydrophobicity scale ("Kyte-Doolittle")
    # 		cutoff (float): Cutoff distance between two blobs for a contact in Angstroms (5.5)
    #
    #   Returns:
    #       None

    # Calculate excess distance between midpoints of blobs
    source calcExcessdistMid.tcl

    # Blobulate protein
    source blobulate_all.tcl

    # Blobulate protein
    blobulate_protein $lMin $H $dictInput

    set structureName [molinfo top get name]
    set baseName [file rootname [file tail $structureName]]
	set filename "${baseName}_odds_ratios_${variant_index}_X_to_X_${termini_index}.txt"
    set fp [open $filename "w"]

    # Get all blobs
    set prot [atomselect top "protein"]
    set blob_list [lsort -real -unique [$prot get user2]]
    $prot delete

	foreach blob_x_index $blob_list {
		set VX_XT_states [get_VX_XT_state $variant_index $termini_index $blob_x_index]
		set 0VX_XT_states [get_0VX_XT_state $variant_index $termini_index $blob_x_index]
		set VX_0XT_states [get_VX_0XT_state $variant_index $termini_index $blob_x_index] 
		set 0VX_0XT_states [get_0VX_0XT_state $variant_index $termini_index $blob_x_index] 

		set a $VX_XT_states
		set b $0VX_XT_states
		set c $VX_0XT_states
		set d $0VX_0XT_states

		if {$b==0 || $c ==0} {
			set OR "NA"
		} else {
			set OR [expr {1.0 * ($a * $d) / ($b * $c)}]
		}

		set m [expr {int($variant_index)}]
		set x [expr {int($blob_x_index)}]
		set t [expr {int($termini_index)}]
		
		puts $fp "${m}_${x}->${x}_${t} ${a} ${b} ${c} ${d} ${OR}" 

	}
	close $fp
}

proc get_VX_XT_state {variant_index termini_index blob_x_index cutoff} {

    #   Writes the number of frames that display the VX_XT contact state, where X is a blob of the sequence.
    #
    #   Arguments:
    # 		variant_index (int): Index of Mutant blob, defined by blobulation
    #       termini_index (int): Index of a terminal blob, defined by blobulation
    # 		blob_x_index (int): Index of a blob, defined by blobulation
    # 		cutoff (float): Cutoff distance between two blobs for a contact in Angstroms (5.5)
    #
    #   Returns:
    #       num_contact_states (int): The number of contact states (frames) within an ensemble
	
    set numFrames [molinfo top get numframes]
	set num_contact_states 0

	for {set i 0} {$i < $numFrames} {incr i} {
		set distance_between_VX [calcExcessDistMid $variant_index $blob_x_index $i]
		set distance_between_XT [calcExcessDistMid $blob_x_index $termini_index $i]

		if {$distance_between_VX < $cutoff  && $distance_between_XT < $cutoff} {
			incr num_contact_states
		}
	}

	return $num_contact_states
}

proc get_0VX_XT_state {variant_index termini_index blob_x_index cutoff} {

    #   Writes the number of frames that display the 0VX_XT contact state, where X is a blob of the sequence.
    #
    #   Arguments:
    # 		variant_index (int): Index of Mutant blob, defined by blobulation
    #       termini_index (int): Index of a terminal blob, defined by blobulation
    # 		blob_x_index (int): Index of a blob, defined by blobulation
    # 		cutoff (float): Cutoff distance between two blobs for a contact in Angstroms (5.5)	
    #
    #   Returns:
    #       num_contact_states (int): The number of contact states (frames) within an ensemble
	
    set numFrames [molinfo top get numframes]
	set num_contact_states 0

	for {set i 0} {$i < $numFrames} {incr i} {
		set distance_between_VX [calcExcessDistMid $variant_index $blob_x_index $i]
		set distance_between_XT [calcExcessDistMid $blob_x_index $termini_index $i]

		if {$distance_between_VX > $cutoff  && $distance_between_XT < $cutoff} {
			incr num_contact_states
		}
	}

	return $num_contact_states
}

proc get_VX_0XT_state {variant_index termini_index blob_x_index cutoff} {
    
    #   Writes the number of frames that display the VX_0XT contact state, where X is a blob of the sequence.
    #
    #   Arguments:
    # 		variant_index (int): Index of Mutant blob, defined by blobulation
    #       termini_index (int): Index of a terminal blob, defined by blobulation
    # 		blob_x_index (int): Index of a blob, defined by blobulation
    # 		cutoff (float): Cutoff distance between two blobs for a contact in Angstroms (5.5)
    #
    #   Returns:
    #       num_contact_states (int): The number of contact states (frames) within an ensemble	
    
    set numFrames [molinfo top get numframes]
	set num_contact_states 0

	for {set i 0} {$i < $numFrames} {incr i} {
		set distance_between_VX [calcExcessDistMid $variant_index $blob_x_index $i]
		set distance_between_XT [calcExcessDistMid $blob_x_index $termini_index $i]

		if {$distance_between_VX < $cutoff  && $distance_between_XT > $cutoff} {
			incr num_contact_states
		}
	}

	return $num_contact_states
}

proc get_0VX_0XT_state {variant_index termini_index blob_x_index cutoff} {

    #   Writes the number of frames that display the 0VX_0XT contact state, where X is a blob of the sequence.
    #
    #   Arguments:
    # 		variant_index (int): Index of Mutant blob, defined by blobulation
    #       termini_index (int): Index of a terminal blob, defined by blobulation
    # 		blob_x_index (int): Index of a blob, defined by blobulation
    # 		cutoff (float): Cutoff distance between two blobs for a contact in Angstroms (5.5)
    #
    #   Returns:
    #       num_contact_states (int): The number of contact states (frames) within an ensemble	
	
    set numFrames [molinfo top get numframes]
	set num_contact_states 0

	for {set i 0} {$i < $numFrames} {incr i} {
		set distance_between_VX [calcExcessDistMid $variant_index $blob_x_index $i]
		set distance_between_XT [calcExcessDistMid $blob_x_index $termini_index $i]

		if {$distance_between_VX > $cutoff  && $distance_between_XT > $cutoff} {
			incr num_contact_states
		}
	}

	return $num_contact_states
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

proc calc_odds_ratios_all_blob_pairs_NX_NC {termini_N_index termini_C_index lMin H dictInput cutoff} {

    #   Writes the number of frames that display contact-states and calculates the odds ratio (OR)
    #   with these values for all blob pairs involving the Variant (V) and a termini blob to a file.
    # 	The 'X' in the variables below refers to a blob of the protein.
    #
    #   Arguments:
    # 		termini_N_index (int): Index of the N-terminal blob, defined by blobulation
    #       termini_C_index (int): Index of the C-terminal blob, defined by blobulation
    # 		lMin (int): Minimum length of a blob (4)
    #       H (float): Hydropathy threshold (0.37)
    #		dictInput (string): Hydrophobicity scale ("Kyte-Doolittle")
    # 		cutoff (float): Cutoff distance between two blobs for a contact in Angstroms (5.5)
    #
    #   Returns:
    #       None


    # Calculate excess distance between midpoints of blobs
    source calcExcessdistMid.tcl

    # Blobulate protein
    source blobulate_all.tcl

    # Blobulate protein
    blobulate_protein $lMin $H $dictInput

    set structureName [molinfo top get name]
    set baseName [file rootname [file tail $structureName]]
	set filename "${baseName}_odds_ratios_${termini_C_index}_X_to_${termini_N_index}_${termini_C_index}.txt"
    set fp [open $filename "w"]

    # Get all blobs
    set prot [atomselect top "protein"]
    set blob_list [lsort -real -unique [$prot get user2]]
    $prot delete

	foreach blob_x_index $blob_list {
		set NX_NC_states [get_NX_NC_state $termini_N_index $termini_C_index $blob_x_index]
		set 0NX_NC_states [get_0NX_NC_state $termini_N_index $termini_C_index $blob_x_index]
		set NX_0NC_states [get_NX_0NC_state $termini_N_index $termini_C_index $blob_x_index] 
		set 0NX_0NC_states [get_0NX_0NC_state $termini_N_index $termini_C_index $blob_x_index] 

		set a $NX_NC_states
		set b $0NX_NC_states
		set c $NX_0NC_states
		set d $0NX_0NC_states

		if {$b==0 || $c ==0} {
			set OR "NA"
		} else {
			set OR [expr {1.0 * ($a * $d) / ($b * $c)}]
		}
		
		set nt [expr {int($termini_N_index)}]
		set x [expr {int($blob_x_index)}]
		set ct [expr {int($termini_C_index)}]
		
		puts $fp "${ct}_${x}->${nt}_${ct} ${a} ${b} ${c} ${d} ${OR}" 
	}
	close $fp
}

proc get_NX_NC_state {termini_N_index termini_C_index blob_x_index cutoff} {

    #   Writes the number of frames that display the NX_NC contact state, where X is a blob of the sequence.
    #
    #   Arguments:
    # 		variant_index (int): Index of Mutant blob, defined by blobulation
    #       termini_index (int): Index of a terminal blob, defined by blobulation
    # 		blob_x_index (int): Index of a blob, defined by blobulation
    # 		cutoff (float): Cutoff distance between two blobs for a contact in Angstroms (5.5)
    #
    #   Returns:
    #       num_contact_states (int): The number of contact states (frames) within an ensemble		
	
    set numFrames [molinfo top get numframes]
	set num_contact_states 0

	for {set i 0} {$i < $numFrames} {incr i} {
		set distance_between_NX [calcExcessDistMid $termini_C_index $blob_x_index $i]
		set distance_between_NC [calcExcessDistMid $termini_N_index $termini_C_index $i]

		if {$distance_between_NX < cutoff  && $distance_between_NC < cutoff} {
			incr num_contact_states
		}
	}

	return $num_contact_states
}

proc get_0NX_NC_state {termini_N_index termini_C_index blob_x_index cutoff} {

    #   Writes the number of frames that display the 0NX_NC contact state, where X is a blob of the sequence.
    #
    #   Arguments:
    # 		termini_N_index (int): Index of the N-terminal blob, defined by blobulation
    #       termini_C_index (int): Index of the C-terminal blob, defined by blobulation
    # 		blob_x_index (int): Index of a blob, defined by blobulation
    # 		cutoff (float): Cutoff distance between two blobs for a contact in Angstroms (5.5)
	# 
    #   Returns:
    #       num_contact_states (int): The number of contact states (frames) within an ensemble	
	
    set numFrames [molinfo top get numframes]
	set num_contact_states 0

	for {set i 0} {$i < $numFrames} {incr i} {
		set distance_between_NX [calcExcessDistMid $termini_C_index $blob_x_index $i]
		set distance_between_NC [calcExcessDistMid $termini_N_index $termini_C_index $i]

		if {$distance_between_NX > cutoff  && $distance_between_NC < cutoff} {
			incr num_contact_states
		}
	}

	return $num_contact_states
}

proc get_NX_0NC_state {termini_N_index termini_C_index blob_x_index cutoff} {

    #   Writes the number of frames that display the NX_0NC contact state, where X is a blob of the sequence.
    #
    #   Arguments:
    # 		termini_N_index (int): Index of the N-terminal blob, defined by blobulation
    #       termini_C_index (int): Index of the C-terminal blob, defined by blobulation
    # 		blob_x_index (int): Index of a blob, defined by blobulation
    # 		cutoff (float): Cutoff distance between two blobs for a contact in Angstroms (5.5)
    #
    #   Returns:
    #       num_contact_states (int): The number of contact states (frames) within an ensemble	
	
    set numFrames [molinfo top get numframes]
	set num_contact_states 0

	for {set i 0} {$i < $numFrames} {incr i} {
		set distance_between_NX [calcExcessDistMid $termini_C_index $blob_x_index $i]
		set distance_between_NC [calcExcessDistMid $termini_N_index $termini_C_index $i]

		if {$distance_between_NX < cutoff  && $distance_between_NC > cutoff} {
			incr num_contact_states
		}
	}

	return $num_contact_states
}

proc get_0NX_0NC_state {termini_N_index termini_C_index blob_x_index cutoff} {

    #   Writes the number of frames that display the 0NX_0NC contact state, where X is a blob of the sequence.
    #
    #   Arguments:
    # 		termini_N_index (int): Index of the N-terminal blob, defined by blobulation
    #       termini_C_index (int): Index of the C-terminal blob, defined by blobulation
    # 		blob_x_index (int): Index of a blob, defined by blobulation
    # 		cutoff (float): Cutoff distance between two blobs for a contact in Angstroms (5.5)
    #   
    # 	Returns:
    #       num_contact_states (int): The number of contact states (frames) within an ensemble	
	
    set numFrames [molinfo top get numframes]
	set num_contact_states 0

	for {set i 0} {$i < $numFrames} {incr i} {
		set distance_between_NX [calcExcessDistMid $termini_C_index $blob_x_index $i]
		set distance_between_NC [calcExcessDistMid $termini_N_index $termini_C_index $i]

		if {$distance_between_NX > cutoff  && $distance_between_NC > cutoff} {
			incr num_contact_states
		}
	}

	return $num_contact_states
}