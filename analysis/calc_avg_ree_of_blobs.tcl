proc calcReeofBlobs {output_dir traj_without_equilibration_base} {

    # Calculates the average end to end distance of each blob in a protein sequence and outputs it into a file.
    # 
    #   Arguments:
    #       output_dir (str): Directory for output 
    #       traj_without_equilibration_base (str): Name of the trajectory that does not contain equilibration frames
    #
    #   Returns:
    #       None

    set outputFile "${output_dir}/ree_per_blob_${traj_without_equilibration_base}.txt"
    set fh [open $outputFile "w"]

    # Get the blobs in the sequence
    set sel [atomselect top "all"]
    set user2_vals [$sel get user2]
    set blobs [lsort -unique -real $user2_vals]
    $sel delete

    set numFrames [molinfo top get numframes]  

    foreach blob $blobs {
        set sum_dist 0

        # Select each blob in the protein sequence
        set sel [atomselect top "user2 == $blob"]
        set resids [lsort -integer [$sel get resid]]
        $sel delete

        set firstResid [lindex $resids 0]
        set lastResid  [lindex $resids end]

        for {set frame 0} {$frame < $numFrames} {incr frame} {

            # Get N and O atoms of the first and last resid of the protein sequence respectively
    		set selN [atomselect top "resid $firstResid and name N" frame $frame]
    		set selC [atomselect top "resid $lastResid and name O" frame $frame]

    	    set coordN [lindex [$selN get {x y z}] 0]
    	    set coordC [lindex [$selC get {x y z}] 0]
    	    
    	    $selN delete
    	    $selC delete

            # Calculate the end-to-end distance
        	set dist [vecdist $coordN $coordC]

            # Add distance to the sums (to calculate the average in the next step)
            set sum_dist [expr {$sum_dist + $dist}]
        }

        # Average the end-to-end distance
        set avg_ree [expr {$sum_dist / double($numFrames)}]
        
        puts $fh $avg_ree
    }
    
    close $fh
}