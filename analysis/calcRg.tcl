proc calcBlobRg {blob_sel} {
    
    #	Measures the Rg of a given blob
    #
    #   Arguments:
    # 		blob_sel (): An atom selection of a blob at a particular frame
    #
    #   Returns:
    #       rg_blob (float): The radius of gyration of a blob selection

	set rg_blob [measure rgyr $blob_sel weight mass]
	
	return $rg_blob
}