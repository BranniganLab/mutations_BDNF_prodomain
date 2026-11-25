proc calcAllBlobRg {output_dir traj_without_equilibration_base Lmin} {

    #   This script outputs one file that contains Rg values of each frame for each blob. 
    #   The number of blobs corresponds to the number of rows in the file.
    #   Each row contains Rg values of one blob for all frames.
    #   Note: Blobs containing fewer than Lmin residues are omitted from the
    #         output. The output row order is the sequence order of the
    #         remaining blobs.
    #       
    #   Arguments:
    #       output_dir (str): Directory for output 
    #       traj_without_equilibration_base (str): Name of the trajectory that does not contain equilibration frames
    #       Lmin (int): Minimum length of a blob to include in the output.
    #
    #   Returns:
    #       None
    
    set numOfFrames [molinfo top get numframes]

    # Get list of blobs
    set sel [atomselect top "all and alpha"]
    set list_of_blobs [lsort -unique [$sel get user2]]
    $sel delete

    # Initialize a list for blob residue ranges
    set list_of_blob_resid_ranges {}

    foreach blob $list_of_blobs {
        set blob_sel [atomselect top "user2 $blob"]
        set resids_in_blob [$blob_sel get resid]
        $blob_sel delete
        set blobRange [list [lindex $resids_in_blob 0] [lindex $resids_in_blob end]]
        lappend list_of_blob_resid_ranges $blobRange
    }

    # Sort list_of_blob_resid_ranges by the first element of each range
    set list_of_blob_resid_ranges [lsort -index 0 -integer $list_of_blob_resid_ranges]
    set num_of_blobs [llength $list_of_blob_resid_ranges]

    # Open file for writing blob data
    set fp [open "${output_dir}/blob_Rg_${traj_without_equilibration_base}.txt" w]

    # Iterate over each blob residue range and calculate radius of gyration
    set count 1
    foreach blob_resid_range $list_of_blob_resid_ranges {

        # Get the start and end residue IDs for the current blob
        set firstResid [lindex $blob_resid_range 0]
        set lastResid [lindex $blob_resid_range 1]

        # Skip blobs shorter than Lmin residues.
        if {[expr {$lastResid - $firstResid + 1}] < $Lmin} {
            continue
        }

        set Rg {}
        set rangeString "$firstResid to $lastResid"
        set sel [atomselect top "resid $rangeString"]

        for {set i 0} {$i < $numOfFrames} {incr i} {
            $sel frame $i
            set rg_in_A [measure rgyr $sel weight mass]
            set rg_in_nm [expr $rg_in_A / 10.0] ;# Convert to nm
            lappend Rg $rg_in_nm

        }

        set blob_Rg($count) $Rg
        puts $fp $Rg
        $sel delete
        incr count
    }

    close $fp
    set currentDir [pwd]
}

