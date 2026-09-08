proc calcReeofProtein {output_dir traj_without_equilibration_base} {

    # Calculates the end to end distance of the entire protein sequence for each frame and outputs it into a file.
    # 
    #   Arguments:
    #       output_dir (str): Directory for output 
    #       traj_without_equilibration_base (str): Name of the trajectory that does not contain equilibration frames 
    # 
    #   Returns:
    #       None
    
    set numFrames [molinfo top get numframes]

    # Set first and last resid of the protein
    set sel [atomselect top "protein"]
    set firstResid [lindex [$sel get resid] 0]
    set lastResid  [lindex [$sel get resid] end]
    $sel delete

    set outfileName "${output_dir}/ree_${traj_without_equilibration_base}.txt"
    set f [open $outfileName "w"]

    set ree {}
    for {set frame 0} {$frame < $numFrames} {incr frame} {
		set selN [atomselect top "resid $firstResid and name N" frame $frame]
		set selO [atomselect top "resid $lastResid and name O" frame $frame]
	    set coordN [lindex [$selN get {x y z}] 0]
	    set coordO [lindex [$selO get {x y z}] 0]
	    $selN delete
	    $selO delete

        # Calculate end-to-end distance
    	set ree [vecdist $coordN $coordO]
    	puts $f $ree

	}

    close $f
}