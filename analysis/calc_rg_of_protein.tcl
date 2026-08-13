proc calcRgofProtein {} {

    #   Calculates the radius of gyration (Rg) of a protein sequence in each frame and outputs it into a file.
    # 
    #   Arguments:
    #       None
    #
    #   Returns:
    #       None

    set numFrames [molinfo top get numframes]

    # Get the basename of the loaded molecule file (without path or extension)
    set structureName [molinfo top get name]
    set baseName [file rootname [file tail $structureName]]

    # Output file named based on basename
    set outputFile "rg_${baseName}.txt"
    set f [open $outputFile "w"]

    for {set i 0} {$i < $numFrames} {incr i} {
        set sel [atomselect top "protein" frame $i]
        set rg_sel [measure rgyr $sel weight mass]
        puts $f "$rg_sel"
        $sel delete
    }

    close $f
}