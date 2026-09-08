proc calc_one_pair_blob_contact {blobIndexi blobIndexj cutoff output_dir traj_without_equilibration_base} {
	
	# 	Calculates the contacts between two blobs in a trajectory
    #       
    #   Arguments:
	# 		blobIndexi (int): index of a blob in a sequence
	# 		blobIndexj (int): index of a blob in a sequence
	# 		cutoff (float): distance cut off between two blobs in angstroms (5.5)
    #       output_dir (str): Directory for output 
    #       traj_without_equilibration_base (str): Name of the trajectory that does not contain equilibration frames    
    # 
    #   Returns:
    #       None
    
    source calc_excess_dist_of_blob_sel.tcl

    set outfileName "${output_dir}/single_pair_blob_contact_${blobIndexi}_${blobIndexj}_${traj_without_equilibration_base}.txt"
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