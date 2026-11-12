proc export_blob_assignments {output_dir Lmin H_star} {
    
    # Export the blob assignments for a reference sequence (V66)
    # 
    #   Arguments:
    #       output_dir (str): Directory for output 
    #       L_Min (int): Minimum length of a blob (4)
    #       H_star (float): Hydropathy threshold (0.37)
    #
    #   Returns:
    #       None


    set outfile "${output_dir}/V66_blob_assignments_Lmin${Lmin}_Hstar${H_star}.txt"

    # Select protein residues
    set protein [atomselect top "protein"]
    set resids [lsort -integer -unique [$protein get resid]]

    set out [open $outfile w]

    # Write header
    puts $out "resid\tresname\tuser\tuser2"

    foreach resid $resids {

        set sel [atomselect top "protein and resid $resid"]

        # Retrieve residue information
        set resname [lindex [lsort -unique [$sel get resname]] 0]
        set user_values [lsort -unique [$sel get user]]
        set user2_values [lsort -unique [$sel get user2]]

        # Ensure every atom in the residue has the same assignments
        if {[llength $user_values] != 1 || [llength $user2_values] != 1} {
            puts "ERROR: Inconsistent assignments for residue $resid"
            $sel delete
            close $out
            $protein delete
            error "Assignment export stopped."
        }

        set user [lindex $user_values 0]
        set user2 [lindex $user2_values 0]

        # Write one line per residue
        puts $out "$resid\t$resname\t$user\t$user2"

        $sel delete
    }

    close $out
    $protein delete
}