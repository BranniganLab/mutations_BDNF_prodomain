proc calc_one_pair_blob_contact {blobIndexi blobIndexj cutoff output_dir} {
	
	# 	Calculates the contacts between two blobs in a trajectory
    #       
    #   Arguments:
	# 		blobIndexi (int): index of a blob in a sequence
	# 		blobIndexj (int): index of a blob in a sequence
	# 		cutoff (float): distance cut off between two blobs in angstroms (5.5)
    #
    #   Returns:
    #       None
    
    source calc_excess_dist_of_blob_sel.tcl

    set structureName [molinfo top get name]
    set baseName [file rootname [file tail $structureName]]
    set outfileName "${output_dir}/single_pair_blob_contact_${blobIndexi}_${blobIndexj}_${baseName}.txt"
    set f [open $outfileName "w"]

	set nframes [molinfo top get numframes]

	for {set i 0} {$i < $nframes} {incr i} {

        set excess12 [calcExcessDistMid $blobIndexi $blobIndexj $i]

		if {$excess12 < $cutoff} {
			puts $f 1
		} else {
			puts $f 0
		}
	}

	close $f

	
}