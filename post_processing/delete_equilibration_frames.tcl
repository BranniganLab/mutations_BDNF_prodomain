proc delete_equilibration_frames {equil_time processed_dir} {

	#	Deletes the frames that correspond to the equilibration run in a trajectory
	#	This proc should be performed AFTER removing periodic image (mindist) frames 
	#	
	#	Arguments:
	#		equil_time (int): Time of equilibration converted to frame (ie 800ns = 8,000 frames)
    #       processed_dir (str): Directory that contains the processed trajectories
    # 
    #   Returns:
    #       None

	set runTime [expr $equil_time - 1]; #Subtract 1 frame by the time we want to delete up to (so that we include this frame in the trajectory)

    for {set i $runTime} {$i >= 0} {incr i -1} {
        set frame_to_be_deleted $i
        animate delete beg $frame_to_be_deleted end $frame_to_be_deleted
    }

	# Write output
    set structureName [molinfo top get name] 
    set seq [string index $structureName 0] 
    animate write gro "${processed_dir}/${seq}66_wrapped_centered_PIF_equil_del.gro" beg 0 end -1 top
 
}
