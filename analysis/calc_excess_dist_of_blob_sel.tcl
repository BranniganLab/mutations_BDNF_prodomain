proc calcExcessDistMid {blobIndexi blobIndexj frame} {

    #   Calculates the excess distance between two blobs.
    #
    #   Arguments:
    #       blobIndexi (int): The index of a blob
    #       blobIndexj (int): The index of a blob
    #       frame (int): A frame of the trajectory
    #
    #   Returns:
    #       excess (float): Excess distance between two blobs i and j
    
    source calc_midpoint_of_blob_sel.tcl
    source calc_rg_of_blob_sel.tcl

    # Blob i
    set blobSel [atomselect top "protein and user2 $blobIndexi"]
    set residsi [$blobSel get resid]
    $blobSel delete
    set firstResidi [lindex $residsi 0]
    set lastResidi  [lindex $residsi end]

    # Midpoint for blob i of pair x
    set selN1 [atomselect top "resid $firstResidi and name N" frame $frame]
    set selO1 [atomselect top "resid $lastResidi and name O" frame $frame]
    set mid1 [calcBlobMidpoint $selN1 $selO1]

    $selN1 delete
    $selO1 delete

    # Radius of gyration of blob i
    set blobSel1 [atomselect top "resid $firstResidi to $lastResidi" frame $frame]
    set rg1 [calcBlobRg $blobSel1]
    $blobSel1 delete

    # Blob j
    set blobSel [atomselect top "protein and user2 $blobIndexj"]
    set residsj [$blobSel get resid]
    $blobSel delete
    set firstResidj [lindex $residsj 0]
    set lastResidj  [lindex $residsj end]

    # Midpoint for blob j of pair x
    set selN2 [atomselect top "resid $firstResidj and name N" frame $frame]
    set selO2 [atomselect top "resid $lastResidj and name O" frame $frame]
    set mid2 [calcBlobMidpoint $selN2 $selO2]
    $selN2 delete
    $selO2 delete

    # Radius of gyration of blob j
    set blobSel2 [atomselect top "resid $firstResidj to $lastResidj" frame $frame]
    set rg2 [calcBlobRg $blobSel2]
    $blobSel2 delete

    # Excess distance between blob i and j
    set dist [vecdist $mid1 $mid2]
    set excess [expr {$dist - ($rg1 + $rg2)}]
    
    return $excess
}

