proc classify_contact_state {V_index M_index N_index C_index cutoff cluster output_dir} {
    
    #   Classifies an ensemble (trajectory) based on the presence or absence of contacts.
    #   This script deletes the frames of the original trajectory where the condition for a contact state is not satisfied,
    #   then rewrites the retained frames as a new trajectory (Note: The script does not modify the original trajectory)
    #
    #   Arguments:
    #       V_index (int): The index of Mutant blob (8.0)
    #       M_index (int): The index of Partner blob (12.0)
    #       N_index (int): The index of N-terminal blob (2.0)
    #       C_index (int): The index of C-terminal blob 2 (16.0)
    #       cutoff (float or int): The cutoff value for calculating a contact (Angstroms) (5.5)
    #       cluster (str): The desired cluster to calculate
    #       output_dir (str): Directory for output 
    #
    #   Returns:
    #       None

    # Calculate excess distance between midpoints of blobs
    source calc_excess_dist_of_blob_sel.tcl
    
    # List to store frames that should be deleted
    set framesToDelete {}

    set numFrames [molinfo top get numframes]

    for {set i 0} {$i < $numFrames} {incr i} {
        set VN [calcExcessDistMid $V_index $N_index $i]
        set VC [calcExcessDistMid $V_index $C_index $i]
        set VM [calcExcessDistMid $V_index $M_index $i]
        set MN [calcExcessDistMid $M_index $N_index $i]
        set MC [calcExcessDistMid $M_index $C_index $i]
        set NC [calcExcessDistMid $N_index $C_index $i]

        set cluster_condition [check_conditional $VN $VC $VM $MN $MC $NC $cluster $cutoff]
        
        if {$cluster_condition == 1} { 
            lappend framesToDelete $i
        }
    }

    # Delete frames directly from the list 
    set sorted_unwanted_frames_list [lsort -real $framesToDelete]
    for {set i [expr {[llength $sorted_unwanted_frames_list] - 1}]} {$i >= 0} {incr i -1} {
        set frame_to_be_deleted [lindex $sorted_unwanted_frames_list $i] 
        animate delete beg $frame_to_be_deleted end $frame_to_be_deleted top
    }

    # Rename and write output
    set structureName [molinfo top get name] 
    set seq [string index $structureName 0] 
    
    mol rename top "${seq}_${cluster}.gro" 
    animate write gro "${output_dir}/${seq}_${cluster}.gro" beg 0 end -1 top
}

proc check_conditional {VN VC VM MN MC NC cluster cutoff} {
    
    #   This script checks the conditions/rules for a cluster.
    #
    #   Arguments:
    #       VN: Distance between Variant-NTB
    #       VC: Distance between Variant-CTB
    #       VM: Distance between Variant-Mediator
    #       MN: Distance between Mediator-NTB
    #       MC: Distance between Mediator-CTB
    #       NC: Distance between NTB-CTB
    #
    #   Returns:
    #       1 if the condition is true. 0 is the condition is false

    # Define cluster rules
    set cluster_rules {

        VM_VC_0NC      {(($VM > $cutoff) || ($VC > $cutoff) || ($NC < $cutoff))}
        VM_0VC_0NC     {(($VM > $cutoff) || ($VC < $cutoff) || ($NC < $cutoff))}
        VM_VC_NC       {(($VM > $cutoff) || ($VC > $cutoff) || ($NC > $cutoff))}
        VM_0VC_NC      {(($VM > $cutoff) || ($VC < $cutoff) || ($NC > $cutoff))}
        
        0VM_VC_0NC     {(($VM < $cutoff) || ($VC > $cutoff) || ($NC < $cutoff))}
        0VM_0VC_0NC    {(($VM < $cutoff) || ($VC < $cutoff) || ($NC < $cutoff))}
        0VM_VC_NC      {(($VM < $cutoff) || ($VC > $cutoff) || ($NC > $cutoff))}
        0VM_0VC_NC     {(($VM < $cutoff) || ($VC < $cutoff) || ($NC > $cutoff))}

        VM_MN          {(($VM > $cutoff) || ($MN > $cutoff))}
        0VM_MN         {(($VM < $cutoff) || ($MN > $cutoff))}
        VM_0MN         {(($VM > $cutoff) || ($MN < $cutoff))}
        0VM_0MN        {(($VM < $cutoff) || ($MN < $cutoff))}

        VM_MN_VN       {(($VM > $cutoff) || ($MN > $cutoff) || ($VN > $cutoff) || ($MC < $cutoff))}
        VM_MN_MC       {(($VM > $cutoff) || ($MN > $cutoff) || ($MC > $cutoff) || ($VN < $cutoff))}
        VM_MN_VN_MC    {(($VM > $cutoff) || ($MN > $cutoff) || ($VN > $cutoff) || ($MC > $cutoff))}
        VM_MN_0VN_0MC  {(($VM > $cutoff) || ($MN > $cutoff) || ($VN < $cutoff) || ($MC < $cutoff))}

        0VM_MN_VN      {(($VM < $cutoff) || ($MN > $cutoff) || ($VN > $cutoff) || ($MC < $cutoff))}
        0VM_MN_MC      {(($VM < $cutoff) || ($MN > $cutoff) || ($MC > $cutoff) || ($VN < $cutoff))}
        0VM_MN_VN_MC   {(($VM < $cutoff) || ($MN > $cutoff) || ($VN > $cutoff) || ($MC > $cutoff))}
        0VM_MN_0VN_0MC {(($VM < $cutoff) || ($MN > $cutoff) || ($VN < $cutoff) || ($MC < $cutoff))}

        MN_NC          {(($NC > $cutoff) || ($MN > $cutoff))}
        MN_0NC         {(($NC < $cutoff) || ($MN > $cutoff))}
        0MN_NC         {(($NC > $cutoff) || ($MN < $cutoff))}
        0MN_0NC        {(($NC < $cutoff) || ($MN < $cutoff))}

        MN_VN_0MC_NC   {(($NC > $cutoff) || ($MN > $cutoff) || ($VN > $cutoff) || ($MC < $cutoff))}
        MN_MC_0VN_NC   {(($NC > $cutoff) || ($MN > $cutoff) || ($MC > $cutoff) || ($VN < $cutoff))}
        MN_VN_MC_NC    {(($NC > $cutoff) || ($MN > $cutoff) || ($VN > $cutoff) || ($MC > $cutoff))}
        MN_0VN_0MC_NC  {(($NC > $cutoff) || ($MN > $cutoff) || ($VN < $cutoff) || ($MC < $cutoff))}
    
        0NC_MN_VN      {(($NC < $cutoff) || ($MN > $cutoff) || ($VN > $cutoff) || ($MC < $cutoff))}
        0NC_MN_MC      {(($NC < $cutoff) || ($MN > $cutoff) || ($MC > $cutoff) || ($VN < $cutoff))}
        0NC_MN_VN_MC   {(($NC < $cutoff) || ($MN > $cutoff) || ($VN > $cutoff) || ($MC > $cutoff))}
        0NC_MN_0VN_0MC {(($NC < $cutoff) || ($MN > $cutoff) || ($VN < $cutoff) || ($MC < $cutoff))}
   
        0NC_VM_MN      {(($VM > $cutoff) || ($NC < $cutoff) || ($MN > $cutoff))}
        0NC_0VM_MN     {(($VM < $cutoff) || ($NC < $cutoff) || ($MN > $cutoff))}
        NC_VM_MN       {(($VM > $cutoff) || ($NC > $cutoff) || ($MN > $cutoff))}
        NC_0VM_MN      {(($VM < $cutoff) || ($NC > $cutoff) || ($MN > $cutoff))}

        0NC_VN_MN      {(($VN > $cutoff) || ($NC < $cutoff) || ($MN > $cutoff))}
        0NC_0VN_MN     {(($VN < $cutoff) || ($NC < $cutoff) || ($MN > $cutoff))}
        NC_VN_MN       {(($VN > $cutoff) || ($NC > $cutoff) || ($MN > $cutoff))}
        NC_0VN_MN      {(($VN < $cutoff) || ($NC > $cutoff) || ($MN > $cutoff))}

        VM             {$VM > $cutoff}
        NC             {$NC > $cutoff}

        VC_0NC         {(($VC > $cutoff) || ($NC < $cutoff))}
        0VC_0NC        {(($VC < $cutoff) || ($NC < $cutoff))}
        VC_NC          {(($VC > $cutoff) || ($NC > $cutoff))}
        0VC_NC         {(($VC < $cutoff) || ($NC > $cutoff))}

        VM_NC          {(($VM > $cutoff) || ($NC > $cutoff))}

    }

    # Evaluate rule
    set rule [dict get $cluster_rules $cluster]

    if {[expr $rule]} {
        return 1
    } else {
        return 0
    }
}