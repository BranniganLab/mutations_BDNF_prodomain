proc calcBlobMidpoint {selN selO} {

    #   Calculates the midpoint between two selections (i.e. two atoms of one blob)
    #
    #   Arguments:
    #       selN (): An atom selection of the N atom of a blob at a particular frame
    #       selO (): An atom selection of the O atom of a blob at a particular frame
    #
    #   Returns:
    #       midpoint (list): Three coordinate xyz midpoint of a selection

    set coordN [lindex [$selN get {x y z}] 0]
    set coordO [lindex [$selO get {x y z}] 0]
    set midpoint [vecscale 0.5 [vecadd $coordN $coordO]]
    
    return $midpoint
}
