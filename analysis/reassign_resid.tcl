proc reassign_resid {resid newUser newUser2} {

    #   Reassigns a residue to a desired blob. Useful for ensuring that all residues within 
    #   the BDNF prodomain sequence have the same User and User2 values. In some cases, residue H65
    #   is considered to be in the s blob but it should be in the h2b, mutation-containing blob
    #   
    #   Arguments:
    #       resid (int): the resid of the sequence that will be changed (65)
    #       newUser (int): new User group for the resid (1)
    #       newUser2 (int): new User2 group for the resid (8.0)
    #
    #   Returns:
    #       None

    # Select the entire residue
    set sel [atomselect top "protein and resid $resid"]

    # Get current user assignments
    set oldUser  [lsort -unique [$sel get user]]
    set oldUser2 [lsort -unique [$sel get user2]]
    set resname  [lsort -unique [$sel get resname]]

    puts "Reassigning $resname$resid"
    puts "    user:  $oldUser  -> $newUser"
    puts "    user2: $oldUser2 -> $newUser2"

    # Change blob type and blob identity
    $sel set user  $newUser
    $sel set user2 $newUser2

    $sel delete
}