proc blobulate_protein {L_Min H_star hydrophobicity_scale path_to_blobulator} {

    #   Blobulates a protein
    #
    #   Arguments:
    #       L_Min (int): Minimum length of a blob (4)
    #       H_star (float): Hydropathy threshold (0.37)
    #		hydrophobicity_scale (string): Hydrophobicity scale ("Kyte-Doolittle")
    #       path_to_blobulator: Path to the VMD scripts in the blobulator repo
    #
    #   Returns:
    #       None
    

    set orig_dir [pwd]
	cd $path_to_blobulator
    source blobulation.tcl

	::blobulator::blobulate top $L_Min $H_star "all" $hydrophobicity_scale

    cd $orig_dir

}
