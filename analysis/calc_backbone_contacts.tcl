proc measureContactsBB {cutoff output_dir state} {
    
    #   Calculates the contact frequencies between the backbone of residues.
    #   It reports the frequencies as percentages of the entire analyzed trajectory
    #
    #   Arguments:
    #       cutoff (float/int): Distance between two residues that defines a contact (Angstroms)
    #       output_dir (str): Directory for output 
    #       state (str): Contact state
    #
    #   Returns:
    #       None

    set structureName [molinfo top get name]
    set baseName [file rootname [file tail $structureName]]
    set seq [string index $baseName 0]

    set outfileName "${output_dir}/bb_contactfreqs_${cutoff}A_${seq}_${state}.txt"
    set f [open $outfileName "w"]

    # Get list of residue IDs via CA atoms
    set resids [[atomselect top "protein and name CA"] get resid]
    set numRes [llength $resids]

    # Number of unique residue pairs (excluding i==j, and avoiding double counting)
    set numPairs [expr ($numRes * $numRes - ($numRes * ($numRes + 1) / 2))]
    set contactCounts [lrepeat $numPairs 0]
    set numFrames [molinfo top get numframes]

    # Loop over frames
    for {set frame 0} {$frame < $numFrames} {incr frame} {
        set coords {}
        # Collect center-of-mass coordinates for each residue backbone
        foreach resid $resids {
            set sel [atomselect top "backbone and resid $resid" frame $frame]
            lappend coords [measure center $sel weight mass]
            $sel delete
        }

        # Measure pairwise distances
        set index 0
        for {set i 0} {$i < [llength $coords] - 1} {incr i} {
            set p1 [lindex $coords $i]
            for {set j [expr $i + 1]} {$j < [llength $coords]} {incr j} {
                set p2 [lindex $coords $j]
                set dist [vecdist $p1 $p2]
                if {$dist <= $cutoff} {
                    set currentCount [lindex $contactCounts $index]
                    set contactCounts [lreplace $contactCounts $index $index [expr {$currentCount + 1}]]
                }
                incr index
            }
        }
    }

    # Write contact frequencies
    set index 0
    for {set i 0} {$i < $numRes - 1} {incr i} {
        set resid1 [lindex $resids $i]
        for {set j [expr $i + 1]} {$j < $numRes} {incr j} {
            set resid2 [lindex $resids $j]
            set contactCount [lindex $contactCounts $index]
            set freq [expr {($contactCount * 100.0) / $numFrames}]
            puts $f "$resid1 $resid2 $freq"
            incr index
        }
    }

    close $f

    # Delete atomselections
    foreach sel $scSelectors {
        if {$sel ne ""} {
            $sel delete
        }
    }
}
