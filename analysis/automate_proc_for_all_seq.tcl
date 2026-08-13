# 	Allows user to use a proc for multiple sequences by running autpmate_script.sh
# 
#   Arguments:
#       None
#
#   Returns:
#       None

# Source a script
source calc_avg_ree_of_blobs.tcl

# Set trajectory name
set traj_name [lindex $argv 0]

# Load sequence into VMD
mol new $traj_name type gro waitfor all

# Run calculation
calcReeofBlobs
quit