# Set arguments
set equil_time [lindex $argv 0]
set gro_name [lindex $argv 1]
set xtc_name [lindex $argv 2]
set processed_dir [lindex $argv 3]

# Load sequence into VMD and remove the first frame since it is the gro structure
mol new $gro_name type gro waitfor all
mol addfile $xtc_name type xtc waitfor all
animate delete beg 0 end 0 

# Deletes periodic image frames of the protein
source delete_equilibration_frames.tcl
delete_equilibration_frames $equil_time $processed_dir

quit

