proc calcAllBlobRg {output_dir} {

    #   This script outputs one file that contains Rg values of each frame for each blob. 
    #   The number of blobs corresponds to the number of rows in the file. 
    #   Note: This script does not consider s-blobs!
    #       
    #   Arguments:
    #       None
    #
    #   Returns:
    #       None
    
    set numOfFrames [molinfo top get numframes]

    # Get list of blobs
    set sel [atomselect top "protein and alpha"]
    set list_of_blobs [lsort -unique [$sel get user2]]
    $sel delete

    # Initialize a list for blob residue ranges
    set list_of_blob_resid_ranges {}

    foreach blob $list_of_blobs {
        set blob_sel [atomselect top "protein and user2 $blob"]
        set resids_in_blob [$blob_sel get resid]
        $blob_sel delete
        set blobRange [list [lindex $resids_in_blob 0] [lindex $resids_in_blob end]]
        lappend list_of_blob_resid_ranges $blobRange
    }

    # Sort list_of_blob_resid_ranges by the first element of each range
    set list_of_blob_resid_ranges [lsort -index 0 -integer $list_of_blob_resid_ranges]
    set num_of_blobs [llength $list_of_blob_resid_ranges]

    # Open file for writing blob data
    set structureFullName [molinfo top get name]
    set structureBaseName [file rootname [file tail $structureFullName]]
    set fp [open "${output_dir}/blob_Rg_${structureBaseName}.txt" w]

    # Iterate over each blob residue range and calculate radius of gyration
    set count 1
    foreach blob_resid_range $list_of_blob_resid_ranges {

        # Get the start and end residue IDs for the current blob
        set firstResid [lindex $blob_resid_range 0]
        set lastResid [lindex $blob_resid_range 1]

        # Skip if the blob range is smaller than lMin residues (we don't want to calculate the Rg for s blobs)
        if {[expr {$lastResid - $firstResid + 1}] < 4} {
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

