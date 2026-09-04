proc calcRgofProtein {output_dir traj_name} {

    #   Calculates the radius of gyration (Rg) of a protein sequence in each frame and outputs it into a file.
    # 
    #   Arguments:
    #       None
    #
    #   Returns:
    #       None

    set numFrames [molinfo top get numframes]

    # Output file named based on basename
    set outputFile "${output_dir}/rg_${traj_name}.txt"
    set f [open $outputFile "w"]

    for {set i 0} {$i < $numFrames} {incr i} {
        set sel [atomselect top "all" frame $i]
        set rg_sel [measure rgyr $sel weight mass]
        puts $f "$rg_sel"
        $sel delete
    }

    close $f
}