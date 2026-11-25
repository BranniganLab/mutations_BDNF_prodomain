proc apply_ref_seq_blob_assignments {blobulated_ref_seq output_dir Lmin H_star} {

    #   Applies the reference sequence Blobulator assignments to the top VMD molecule
    #       
    #   Arguments:
    #       blobulated_ref_seq (str): The reference blobulated sequence
    #       output_dir (str): Directory for output 
    #       L_Min (int): Minimum length of a blob
    #       H_star (float): Hydropathy threshold
    #
    #   Returns:
    #       None

    set infile "${output_dir}/${blobulated_ref_seq}_blob_assignments_Lmin${Lmin}_Hstar${H_star}.txt"
    set in [open $infile r]

    # Skip header
    gets $in header

    set numFrames [molinfo top get numframes]

    while {[gets $in line] >= 0} {

        set fields [split $line "\t"]

        set resid   [lindex $fields 0]
        set resname [lindex $fields 1]
        set user    [lindex $fields 2]
        set user2   [lindex $fields 3]

        # Select residue independently of frame
        set sel [atomselect top "protein and resid $resid"]

        set current_resname [lindex [lsort -unique [$sel get resname]] 0]

        if {[$sel num] > 0} {

            # Apply assignments to every trajectory frame
            for {set i 0} {$i < $numFrames} {incr i} {

                $sel frame $i
                $sel set user $user
                $sel set user2 $user2
            }

            if {$current_resname ne $resname} {
                puts "NOTE: Residue $resid is $current_resname, reference is $resname. Applied reference blob assignment."
            }
        }

        $sel delete
    }

    close $in

    puts "Applied V66 Blobulator assignments to all $numFrames frames"
}
