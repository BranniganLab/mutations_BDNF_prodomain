proc blobulate_protein {lMin H dictInput path_to_blobulator} {

    #   Blobulates a protein
    #
    #   Arguments:
    #       lMin (int): Minimum length of a blob (4)
    #       H (float): Hydropathy threshold (0.37)
    #		dictInput (string): Hydrophobicity scale ("Kyte-Doolittle")
    #
    #   Returns:
    #       None
    

    set orig_dir [pwd]
	cd $path_to_blobulator
    source blobulation.tcl

	::blobulator::blobulate top $lMin $H "all" $dictInput

    cd $orig_dir

}
