proc calcExcessDist {blobIndex1 blobIndex2 frame} {

    #   Calculates the excess distance between two blobs.
    #
    #   Arguments:
    #       blobIndex1 (int): A blob index, defined by blobulation
    #       blobIndex2 (int): A blob index, defined by blobulation
    #       frame (int): A frame of a trajectory
    #
    #   Returns:
    #       excess (float): Excess distance
    
    source calcMid.tcl
    source calcRg.tcl

    # Get blobs from user2 field
    set sel [atomselect top "protein and alpha"]
    set user2list [lsort -unique [$sel get user2]]
    $sel delete

    # Get corresponding user2 values based on the blob indices specified by user
    set blobi $blobIndexi
    set blobj $blobIndexj

    # Print and extract ranges for just blobs i and j
    foreach {blob residRange} [list $blobi residRangei $blobj residRangej] {
        set blobSel [atomselect top "protein and name CA and user2 $blob"]
        set resids [$blobSel get resid]
        $blobSel delete
        set first [lindex $resids 0]
        set last  [lindex $resids end]
        set $residRange [list $first $last]
    }

    # Extract residue ranges for blobIndex1 and blobIndex2
    set firstResidi [lindex $residRangei 0]
    set lastResidi  [lindex $residRangei 1]
    set firstResidj [lindex $residRangej 0]
    set lastResidj  [lindex $residRangej 1]
    
    # Midpoint for blob 1 of pair x
    set selN1 [atomselect top "resid $firstResidi and name N" frame $frame]
    set selO1 [atomselect top "resid $lastResidi and name O" frame $frame]
    set mid1 [calcBlobMidpoint $selN1 $selO1]
    $selN1 delete
    $selO1 delete

    # Midpoint for blob 2 of pair x
    set selN2 [atomselect top "resid $firstResidj and name N" frame $frame]
    set selO2 [atomselect top "resid $lastResidj and name O" frame $frame] 
    set mid2 [calcBlobMidpoint $selN2 $selO2]
    $selN2 delete
    $selO2 delete

    # Radius of gyration
    set blobSel1 [atomselect top "resid $firstResidi to $lastResidi" frame $frame]
    set rg1 [calcBlobRg $blobSel1]
    $blobSel1 delete

    set blobSel2 [atomselect top "resid $firstResidj to $lastResidj" frame $frame]
    set rg2 [calcBlobRg $blobSel2]
    $blobSel2 delete

    # Excess distance
    set dist [vecdist $mid1 $mid2]
    set excess [expr {$dist - ($rg1 + $rg2)}]
    
    return $excess
}

