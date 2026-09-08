proc measureContactsSC {cutoff output_dir state} {

    #   Calculates the contact frequencies between sidechains (excluding GLY) of residues.
    #
    #   Arguments:
    #       cutoff (float/int): Distance between two residues that defines a contact
    #       output_dir (str): Directory for output 
    #       state (str): Contact state
    #
    #   Returns:
    #       None

    set structureName [molinfo top get name]
    set baseName [file rootname [file tail $structureName]]
    set seq [string index $baseName 0]

    set outfileName "${output_dir}/sc_contactfreqs_${cutoff}A_${seq}_${state}.txt"
    set f [open $outfileName "w"]

    # # Get list of residue IDs via CA atoms
    set resids [[atomselect top "protein and name CA"] get resid]
    set numRes [llength $resids]

    # Number of unique residue pairs (excluding i==j, and avoiding double counting)
    set numPairs [expr {$numRes * $numRes - ($numRes * ($numRes + 1) / 2)}]
    set contactCounts [lrepeat $numPairs 0]
    set numFrames [molinfo top get numframes]

    # Precompute sidechain selectors and residue names
    set scSelectors {}
    set isGlycine {}
    foreach resid $resids {
        set resSel [atomselect top "resid $resid"]
        set resname [[lindex $resSel 0] get resname]
        $resSel delete

        if {[lindex $resname 0] == "GLY"} {
            lappend isGlycine 1
            lappend scSelectors ""
        } else {
            set sel [atomselect top "protein and resid $resid and noh and not name N CA C O"]
            lappend scSelectors $sel
            lappend isGlycine 0
        }
    }

    # Loop over frames
    for {set frame 0} {$frame < $numFrames} {incr frame} {
        set coords {}

        # Update selections and collect centers of mass
        for {set i 0} {$i < $numRes} {incr i} {
            if {[lindex $isGlycine $i]} {
                lappend coords "GLY"
            } else {
                set sel [lindex $scSelectors $i]
                $sel frame $frame
                $sel update
                lappend coords [measure center $sel weight mass]
            }
        }

        # Measure pairwise distances
        set index 0
        for {set i 0} {$i < $numRes - 1} {incr i} {
            set p1 [lindex $coords $i]
            for {set j [expr {$i + 1}]} {$j < $numRes} {incr j} {
                set p2 [lindex $coords $j]

                if {![string match "GLY*" $p1] && ![string match "GLY*" $p2]} {
                    set dist [vecdist $p1 $p2]
                    if {$dist <= $cutoff} {
                        set currentCount [lindex $contactCounts $index]
                        set contactCounts [lreplace $contactCounts $index $index [expr {$currentCount + 1}]]
                    }
                }
                incr index
            }
        }
    }

    # Write contact frequencies
    set index 0
    for {set i 0} {$i < $numRes - 1} {incr i} {
        set resid1 [lindex $resids $i]
        for {set j [expr {$i + 1}]} {$j < $numRes} {incr j} {
            set resid2 [lindex $resids $j]
            set count [lindex $contactCounts $index]
            set freq [expr {($count * 100.0) / $numFrames}]
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
