proc resolve_blob_indices {assignment_file output_file Lmin mutation_resid mediator_resid} {

    #   Retrieves the blob indices of the reference blobulated BDNF prodomain 
    #
    #   Arguments:
    #       assignment_file (str): Path to the V66 blob assignment file.
    #       output_file (str): Path for the resolved blob-index output file.
    #       Lmin (int): Minimum blob length used when generating blob data.
    #       mutation_resid (int): Residue containing the sequence mutation (66).
    #       mediator_resid (int): Residue defining the mediator blob (95).
    #
    #   Returns:
    #       list: The 1-based blob indices in the order:
    #           N-terminal, C-terminal, variant, mediator.
    #
    #   Notes:
    #       The assignment file must contain the columns:
    #           resid, resname, user, user2
    #       A user value of 1.0 identifies an H-blob, and user2 is the
    #       1-based blob index. The output file also contains zero-based
    #       indices after removing blobs shorter than Lmin residues. These
    #       reduced indices match the blob-data files used by the plotting
    #       scripts.

    if {![file exists $assignment_file]} {
        error "Blob assignment file not found: $assignment_file"
    }

    set in [open $assignment_file r]
 
    # Skip the header.
    if {[gets $in header] < 0} {
        close $in
        error "Blob assignment file is empty: $assignment_file"
    }

    set first_h_blob ""
    set last_h_blob ""
    set variant_blob ""
    set mediator_blob ""
    set variant_matches 0
    set mediator_matches 0

    # Keep the blobs in assignment-file order and record their residue counts.
    # The blob-data producers omit every blob shorter than Lmin residues, so
    # the reduced indices must use the same rule.
    set blob_order {}
    array set blob_lengths {}

    while {[gets $in line] >= 0} {
        set line [string trim $line]
        if {$line eq ""} {
            continue
        }

        set fields [split $line "\t"]
        if {[llength $fields] < 4} {
            close $in
            error "Malformed assignment line in $assignment_file: $line"
        }

        set resid [lindex $fields 0]
        set user [lindex $fields 2]
        set blob [lindex $fields 3]

        if {![string is integer -strict $resid] ||
            ![string is double -strict $user] ||
            ![string is double -strict $blob]} {
            close $in
            error "Non-numeric assignment values in $assignment_file: $line"
        }

        set blob_index [expr {int($blob)}]
        if {$blob != $blob_index} {
            close $in
            error "Blob index is not an integer in $assignment_file: $line"
        }

        # H-blobs have user == 1.0. The assignment file is ordered by
        # residue, so the first and last H-blobs are the terminal H-blobs.
        if {$user == 1.0} {
            if {$first_h_blob eq ""} {
                set first_h_blob $blob_index
            }
            set last_h_blob $blob_index

            if {$resid == $mutation_resid} {
                set variant_blob $blob_index
                incr variant_matches
            }

            if {$resid == $mediator_resid} {
                set mediator_blob $blob_index
                incr mediator_matches
            }
        }

        if {[lsearch -exact $blob_order $blob_index] < 0} {
            lappend blob_order $blob_index
        }
        if {[info exists blob_lengths($blob_index)]} {
            incr blob_lengths($blob_index)
        } else {
            set blob_lengths($blob_index) 1
        }
    }

    close $in

    if {$first_h_blob eq "" || $last_h_blob eq ""} {
        error "No h-blobs were found in: $assignment_file"
    }
    if {$variant_matches != 1} {
        error "Expected residue $mutation_resid to occur in exactly one h-blob; found $variant_matches occurrences"
    }
    if {$mediator_matches != 1} {
        error "Expected residue $mediator_resid to occur in exactly one h-blob; found $mediator_matches occurrences"
    }

    # Build a mapping from full 1-based user2 indices to compact zero-based
    # indices after removing the short blobs omitted by the data producers.
    set compact_index 0
    array set no_s_blob_index {}
    foreach blob_index $blob_order {
        if {$blob_lengths($blob_index) < $Lmin} {
            continue
        }
        set no_s_blob_index($blob_index) $compact_index
        incr compact_index
    }

    foreach required_blob [list $variant_blob $mediator_blob] {
        if {![info exists no_s_blob_index($required_blob)]} {
            error "Variant or mediator h-blob $required_blob was omitted because it has fewer than $Lmin residues"
        }
    }

    if {[info exists no_s_blob_index($first_h_blob)]} {
        set N_terminal_index_no_s_blobs $no_s_blob_index($first_h_blob)
    } else {
        set N_terminal_index_no_s_blobs ""
    }
    if {[info exists no_s_blob_index($last_h_blob)]} {
        set C_terminal_index_no_s_blobs $no_s_blob_index($last_h_blob)
    } else {
        set C_terminal_index_no_s_blobs ""
    }
    set variant_index_no_s_blobs $no_s_blob_index($variant_blob)
    set mediator_index_no_s_blobs $no_s_blob_index($mediator_blob)

    set indices [list \
        $first_h_blob \
        $last_h_blob \
        $variant_blob \
        $mediator_blob]

    set out [open $output_file w]
    puts $out "# Resolved from: [file tail $assignment_file]"
    puts $out "# Indices are 1-based user2 blob indices."
    puts $out "N_blob_index=$first_h_blob"
    puts $out "C_blob_index=$last_h_blob"
    puts $out "variant_blob_index=$variant_blob"
    puts $out "mediator_blob_index=$mediator_blob"
    puts $out "N_blob_index_no_s_blobs=$N_terminal_index_no_s_blobs"
    puts $out "C_blob_index_no_s_blobs=$C_terminal_index_no_s_blobs"
    puts $out "variant_blob_index_no_s_blobs=$variant_index_no_s_blobs"
    puts $out "mediator_blob_index_no_s_blobs=$mediator_index_no_s_blobs"
    close $out

    return $indices
}
