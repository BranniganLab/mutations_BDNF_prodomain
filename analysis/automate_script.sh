#   Automates a proc on multiple or one sequence. Run the following command in a terminal:
#   ./automate_script.sh 
#
#   Arguments:
#       Varies by proc
#
#   Returns:
#       Varies by proc


# ---------- Uncomment to automate a script ---------- 
for seq in F M V L A Y I; do
  vmd -dispdev none \
      -e automate_proc_for_all_seq.tcl \
      -eofexit \
      -args /trajectories/${seq}_hid_11_wrapped_PIF_equil_del.gro
done


# ---------- Uncomment to automate clustering ---------- 
# for seq in F M V L A Y I; do
#   for cluster in 1MC_0T 0MC_0T 1MC_1T 0MC_1T; do
#     vmd -dispdev none \
#         -e clustering.tcl \
#         -eofexit \
#         -args /home/lmr294/Brannigan/Data/trajectories/${seq}66/${seq}_hid_11_wrapped_PIF_equil_del.gro "$cluster" \
#         &> "vmd_${seq}_${cluster}.log"
#     gmx trjconv -f ${seq}_${cluster}.gro -o ${seq}_${cluster}.xtc
#   done
# done
