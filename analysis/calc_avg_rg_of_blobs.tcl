proc calcRgofBlobs {output_dir traj_without_equilibration_base} {

    # Calculates the average radius of gyration (Rg) of each blob in a protein sequence and outputs it into a file.
    # 
    #   Arguments:
    #       lMin (int): Minimum length of a blob (4)
    #       H (float): Hydropathy threshold (0.37)
    #       dictInput (string): Hydrophobicity scale ("Kyte-Doolittle")
    #
    #   Returns:
    #       None
    
    # Get the basename of the loaded molecule file (without path or extension)
    set outputFile "${output_dir}/rg_per_blob_${traj_without_equilibration_base}.txt"
    set fh [open $outputFile "w"]

    # Get the blobs in the sequence
    set sel [atomselect top "protein"]
    set user2_vals [$sel get user2]
    set blobs [lsort -unique -real $user2_vals]
    $sel delete

    set numFrames [molinfo top get numframes]  

    foreach blob $blobs {
        set sum_rg 0
        for {set i 0} {$i < $numFrames} {incr i} {
            set sel [atomselect top "protein and user2 == $blob" frame $i]
            set rgyrinFrame [measure rgyr $sel weight mass]
            $sel delete

            # Add Rg values
            set sum_rg [expr {$sum_rg + $rgyrinFrame}]
        }

        # Calculate average Rg
        set avg_rg [expr {$sum_rg / double($numFrames)}]
        
        puts $fh $avg_rg
    }

    close $fh
}