proc calcAllBlobMidpoints {output_dir traj_without_equilibration_base} {

    #   This script outputs one file that contains xyz coordinates of the midpoint for each blob
    #   Each row corresponds to a frame. The first 3 values are xyz coordinates of the first blob 
    #   in the protein sequence (e.g. p1 of the BDNF prodomain) and so on. The values are calculated in nm.
    #   Note: This script does not consider s-blobs!
    #       
    #   Arguments:
    #       output_dir (str): Directory for output 
    #       traj_without_equilibration_base (str): Name of the trajectory that does not contain equilibration frames
    #
    #   Returns:
    #       None

    set numOfFrames [molinfo top get numframes]

    # Get list of blobs
    set sel [atomselect top "protein"]
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

    # Open file to store midpoints
    set fp [open "${output_dir}/midpoint_${traj_without_equilibration_base}.txt" w]

    # Loop over each frame
    for {set i 0} {$i < $numOfFrames} {incr i} {
        set line ""

        # Loop over each blob
        foreach blob_resid_range $list_of_blob_resid_ranges {

            # Get the start and end residue IDs for the current blob
            set firstResid [lindex $blob_resid_range 0]
            set lastResid [lindex $blob_resid_range 1]

            # Skip if the blob_resid_range is smaller than lMin residues (we don't want to calculate midpoint for s blobs)
            if {[expr {$lastResid - $firstResid + 1}] < 4} {
                continue
            }

            # Select N atom from start residue and O atom from end residue
            set selN [atomselect top "resid $firstResid and name N" frame $i]
            set selO [atomselect top "resid $lastResid and name O" frame $i]

            set coordN [lindex [$selN get {x y z}] 0]
            set coordO [lindex [$selO get {x y z}] 0]
            set midpoint [vecscale 0.05 [vecadd $coordN $coordO]]

            append line "$midpoint " 
            
            $selN delete
            $selO delete
        }


        puts $fp $line
    }

    close $fp

    set currentDir [pwd]
}