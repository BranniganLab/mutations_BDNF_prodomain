proc apply_ref_seq_blob_assignments {output_dir Lmin H_star} {
    
    #   Applies the reference sequence Blobulator assignments to the top VMD molecule
    #       
    #   Arguments:
    #       output_dir (str): Directory for output 
    #       L_Min (int): Minimum length of a blob (4)
    #       H_star (float): Hydropathy threshold (0.37)
    #
    #   Returns:
    #       None

    set infile "${output_dir}/V66_blob_assignments_Lmin${Lmin}_Hstar${H_star}.txt"
    set in [open $infile r]

    # Skip header
    gets $in header

    while {[gets $in line] >= 0} {
        set fields [split $line "\t"]

        set resid  [lindex $fields 0]
        set resname [lindex $fields 1]
        set user   [lindex $fields 2]
        set user2  [lindex $fields 3]

        set sel [atomselect top "protein and resid $resid"]

        # Check that the residue matches the reference sequence
        set current_resname [lindex [lsort -unique [$sel get resname]] 0]

        if {[$sel num] > 0} {
            $sel set user $user
            $sel set user2 $user2

            if {$current_resname ne $resname} {
                puts "NOTE: Residue $resid is $current_resname, reference is $resname. Applied reference blob assignment."
            }
        }

        $sel delete
    }

    close $in
    puts "Applied V66 Blobulator assignments"
}