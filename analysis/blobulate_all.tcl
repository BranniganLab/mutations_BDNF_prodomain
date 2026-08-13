proc blobulate_protein {lMin H dictInput} {

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
    # Choose the directory that contains "blobulation.tcl, for example:"
    cd /home/lmr294/Brannigan/GitHub/blobulator/VMD_scripts 
	source blobulation.tcl

	::blobulator::blobulate top $lMin $H "all" $dictInput

    cd $orig_dir

}
